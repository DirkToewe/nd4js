'use strict';

/* This file is part of ND4JS.
 *
 * ND4JS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ND4JS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ND4JS. If not, see <http://www.gnu.org/licenses/>.
 */

import {  array,
        NDArray} from '../nd_array'

import {_giv_rot_qr,
        _giv_rot_rows} from '../la/_giv_rot'
import {FrobeniusNorm} from '../la/norm'
import {_rrqr_rank,
        _rrqr_decomp_inplace,
        _norm_update,
        _norm} from '../la/rrqr'
import {_triu_solve,
        _triu_t_solve} from '../la/tri'
import {_urv_decomp_full} from '../la/urv'

import {OptimizationNoProgressError} from './optimization_error';


const REPORT_STATE_READY    = 1,
      REPORT_STATE_NA       = 2,
      REPORT_STATE_CONSIDER = 3;


/** 
 */
export class TrustRegionSolverLSQ
{
  constructor( fJ, x0 )
  {
    // EVAL fJ(x0)
    if( !(fJ instanceof Function) )
      throw new Error('new TrustRegionSolverLSQ(fJ, x0): fJ must be Function.');

    this.fJ = fJ;

    this.report_x = x0 = array('float64', x0);
    if( x0.ndim !== 1 )
      throw new Error('new TrustRegionSolverLSQ(fJ, x0): x0.ndim must be 1.');

    let [f,J] = fJ( new NDArray(x0.shape, x0.data.slice()) );
    this.report_f = f = array('float64', f);
    this.report_J = J = array('float64', J);
    if( f.ndim !== 1 ) throw new Error('new TrustRegionSolverLSQ(fJ, x0): fJ(x) returned [f,J] with f.ndim !== 1.');
    if( J.ndim !== 2 ) throw new Error('new TrustRegionSolverLSQ(fJ, x0): fJ(x) returned [f,J] with J.ndim !== 2.');
  
    const [M,N] = J.shape;

    if( f.shape[0] !== M ) throw new Error('new TrustRegionSolverLSQ(fJ, x0): fJ(x) returned [f,J] where J.shape[0] did not match f.shape[0].');
    if(x0.shape[0] !== N ) throw new Error('new TrustRegionSolverLSQ(fJ, x0): fJ(x) returned [f,J] where J.shape[1] did not match x0.shape[0].');

    this.       loss      =
    this.report_loss      = NaN;
    this.report_loss_grad = null;

    // ALLOCATE MEMORY
    this.M = M;
    this.N = N;
    const L = Math.min(M,N),
          K = Math.max(L+N, M);

    this.D = new Float64Array(N), // <- trust region scaling diagonal

    // STORED FUNCTION INPUT AND OUTPUT AT CURRENT TRUST REGION CENTER
    this.X0 = new Float64Array(N*1); // <- center of the current trust radius (point where fJ was last evaluated)
    this.F0 = new Float64Array(M*1);
    this.J0 = new Float64Array(M*N);
    this.G0 = new Float64Array(N);

    // USED BY computeNewton() and computeNewtonRegularized
    this.tmp = new Float64Array(2*N); // <- used to compute the column norms of J for the diagonal scaling
    this.TMP = new Float64Array(K*N); // <- stores R of regularized RRQR

    // USED IN computeNewton()
    this.newton_P  = new   Int32Array(N); // <- column permutation by RRQR decomp.
    this.newton_Q  = new   Int32Array(N); // <- column permutation by URV  decomp.
    this.newton_R  = new Float64Array(L*L);
    this.newton_V  = new Float64Array(N*N); // <- stores complete orthogonal decomp V TODO: In unlikely case (M < N), economic SRRQR could be used instead of full
    this.newton_UF = new Float64Array(M*1);
    this.newton_dX = new Float64Array(N*1); // <- stores gauss-newton point relative to X0

    // USED IN computeRegularized(λ)
    this.regularized_R0 = new Float64Array(L*N); // <- stores R of RRQR (backed up from Newton)
    this.regularized_QF = new Float64Array(K*1);
    this.regularized_dX = new Float64Array(N*1);

    this.rank = -1;

    this._report_state = REPORT_STATE_CONSIDER;

    Object.seal(this);

    // INIT DATA
    this._considerMove_computeMSE();
    this.makeConsideredMove();
  }


  wiggle()
  {
    throw new OptimizationNoProgressError('Too many unsuccessfull iterations.');
  }


  _considerMove_computeMSE()
  {
    const { M,N, report_f: {data: F},
                 report_J: {data: J} } = this;

    let mse = 0.0;
    for( let i=M; i-- > 0; ) {
      const  s = F[i];
      mse += s*s / M;
    }

    const mse_grad = new Float64Array(N);
    for( let i=M; i-- > 0; )
    for( let j=N; j-- > 0; )
      mse_grad[j] += J[N*i+j] * F[i] / M * 2;

    this.report_loss      = mse;
    this.report_loss_grad = new NDArray(Int32Array.of(N), mse_grad); // <- TODO: computation of loss_grad could probably be delayed until makeConsideredMove() is called
  }


  scaledNorm( X )
  {
    const {N,D} = this;

    if( X.length !== N ) throw new Error('Assertion failed.');

    const norm = new FrobeniusNorm();

    for( let i=N; i-- > 0; )
      norm.include(D[i]*X[i]);

    return norm.result;
  }


  cauchyTravel()
  {
    const {M,N, J0: J, G0: G} = this;

    // polynomial along the gradient const1 + 2ax + bx² = const2 + (a/b + x)² * b
    let a=0,
        b=0;

    for( let i=N; i-- > 0; )
      a += G[i]*G[i];

    for( let i=M; i-- > 0; )
    {
      let Jg = 0;
      for( let j=N; j-- > 0; )
        Jg += J[N*i+j] * G[j];
      b += Jg * Jg;
    }

    return 0===a ? 0 :
           0===b ? -Infinity : -a/b;
  }


  report()
  {
    if( this._report_state !== REPORT_STATE_READY )
      throw new Error('TrustRegionSolverLSQ::report: can only be called once after each makeConsideredMove() but not directly after considerMove(dX).');
    this._report_state = REPORT_STATE_NA;

    const result = [
      this.report_x,
      this.report_loss,
      this.report_loss_grad,
      this.report_f,
      this.report_J
    ];

    this.report_x        =
    this.report_f        =
    this.report_J        =
    this.report_loss_grad = null;
    this.report_loss      = NaN;

    return result;
  }


  considerMove( dX )
  {
    const {M,N, X0, F0, J0} = this;

    this._report_state = REPORT_STATE_CONSIDER;

    if( N !== dX.length )
      throw new Error('Assertion failed.');

    const X_shape = Int32Array.of(N),
          X = new Float64Array(N);

    for( let i=N; i-- > 0; )
      X[i] = X0[i] + dX[i];

    let [f,J] = this.fJ(
      new NDArray( X_shape, X.slice() )
    );

    f = array('float64', f);
    J = array('float64', J);

    if( f.ndim !== 1 ) throw new Error('Assertion failed.');
    if( J.ndim !== 2 ) throw new Error('Assertion failed.');

    if( f.shape[0] !== M ) throw new Error('Assertion failed.');
    if( J.shape[0] !== M ) throw new Error('Assertion failed.');
    if( J.shape[1] !== N ) throw new Error('Assertion failed.');

    this.report_x = new NDArray(X_shape, X);
    this.report_f = f;
    this.report_J = J;

    this._considerMove_computeMSE();

    let predict_loss = 0.0;
    for( let i=M; i-- > 0; )
    {
      let s = 0;
      for( let j=N; j-- > 0; )
        s += J0[N*i+j] * (X[j] - X0[j]);
      s += F0[i];
      predict_loss += s*s / M;
    }

    return [ predict_loss,
         this.report_loss ];
  }


  makeConsideredMove()
  {
    if( this._report_state !== REPORT_STATE_CONSIDER )
      throw new Error('Assertion failed.');

    this._report_state = REPORT_STATE_READY;
    this.rank = -1;

    this.loss = this.report_loss;

    const {M,N, D, tmp: norm, X0, F0, J0, G0} = this;
  
    ;{
      const f = this.report_f.data,
            J = this.report_J.data,
           x0 = this.report_x.data;

      // COPY (x,f,J)
      for( let i=  N; i-- > 0; ) X0[i] = x0[i];
      for( let i=M  ; i-- > 0; ) F0[i] =  f[i];
      for( let i=M*N; i-- > 0; ) J0[i] =  J[i];
    }

    // COMPUTE GRADIENT OF (HALF) SQUARED ERROR (MIND THE HALF :P)
    G0.fill(0.0);
    for( let i=M; i-- > 0; )
    for( let j=N; j-- > 0; )
      G0[j] += J0[N*i+j] * F0[i];

    // UPDATE SCALING
    // compute column norms of J
    if( norm.length !== 2*N )
      throw new Error('Assertion failed.');
    norm.fill(0);
    for( let i=0; i < M; i++ )
      _norm_update(norm, J0, N*i, 0);

    // use column norms of J to update D
    for( let j=0; j < N; j++ )
      D[j] = Math.max( D[j], _norm(norm,j) );
  }


  computeNewton()
  {
    if( 0 <= this.rank )
      return;

    const {
      M,N, D, J0, F0,

      tmp: Y,      
      TMP: T,

      newton_P : P,
      newton_Q : Q,
      newton_V : V,
      newton_dX: X,
      newton_UF: UF
    } = this;

    const L = Math.min(M,N);

    // INIT
    for( let i=  N; i-- > 0; )  P[i] =    i;
    for( let i=M*N; i-- > 0; )  T[i] = J0[i];
    for( let i=M  ; i-- > 0; ) UF[i] = F0[i];

    // SOLVE GLOBAL MINIMUM
    _rrqr_decomp_inplace(M,N,1, T,0, UF,0, P,0, Y);

    
    // backup the rrqr decomposition in R (for finding regularized minima)
    ;{
      const R0 = this.regularized_R0;
      for( let i=L*N; i-- > 0; )
        R0[i] = T[i];
    }

    const rank =
     this.rank = _rrqr_rank(M,N, T,0, Y);

    if( rank < N )
    { // RANK DEFICIENT CASE -> USE COMPLETE ORTHOGONAL DECOMPOSITION

      // factor in scaling
      for( let j=N; j-- > 0; ) {
        const  d = D[P[j]];
        if(0!==d)
          for( let i=rank; i-- > 0; )
            T[N*i+j] /= d;
      }

      // complete orthogonal decompositon
      for( let i=N; i-- > 0; )
        Q[i] = P[i];
      V.fill(0.0);
      _urv_decomp_full( rank,N, T,0, V,0, Q,0 );


      // compute least squares solution
      for( let i=rank; i-- > 0; )
        Y[i] = -UF[i];
      ;{
        const R = this.newton_R;
        // backup complete orthogonal T to R
        for( let i=rank; i-- > 0; )
        for( let j=rank; j-- > 0; )
          R[L*i+j] = T[N*i+j];

        _triu_solve(rank,L,1, R,0, Y,0); // Y = S \ Y
      }

      X.fill(0.0);
      for( let k=rank; k-- > 0; )
      for( let i=N   ; i-- > 0; )
        X[i] += V[N*k+i] * Y[k]; // X = V.T @ Y

      // factor out scaling
      for( let i=N; i-- > 0; ) {
        const     d = D[i];
        if( 0 !== d )
          X[i] /= d;
      }
    }
    else
    { // FULL RANK CASE -> SCALING NOT RELEVANT TO GLOBAL SOLUTION
      for( let i=N; i-- > 0; )
        Y[i] = -UF[i];

      _triu_solve(N,N,1, T,0, Y,0);

      for( let i=N; i-- > 0; )
        X[P[i]] = Y[i];
    }
  }

  // Regularized least squares (RLS) solver. Solves the follwing equation system:
  //
  //   (JᵀJ + λDᵀD)x = Jᵀf
  //
  // The equation above can be solved as the following least squares problem:
  //
  //   min ║Ax - z║₂
  //    x
  //
  // where
  //       ┏              ┓        ┏   ┓
  //       ┃              ┃        ┃   ┃
  //       ┃              ┃        ┃   ┃
  //       ┃       J      ┃        ┃-f ┃
  //       ┃              ┃        ┃   ┃
  //       ┃              ┃        ┃   ┃
  //   A = ┃┄┄┄┄┄┄┄┄┄┄┄┄┄┄┃    z = ┃┄┄┄┃
  //       ┃√λ D[0]       ┃        ┃ 0 ┃
  //       ┃   ⋱         ┃        ┃ ⋮ ┃
  //       ┃    √λ D[N-1] ┃        ┃ 0 ┃
  //       ┗              ┛        ┗   ┛
  //
  // On top of that we can utilize the fact that J has already
  // been RRQR decomposed in computeNewton().
  computeNewtonRegularized( λ )
  {
    this.computeNewton();
    // if( ! (0 <= this.rank) )
    //   throw new Error('TrustRegionSolverLSQ.prototype.computeMinRegularized(λ): Can only be called after computeNewton().');

    const {
      M,N, D, rank: rnk,

      tmp: Y,
      TMP: T,

      newton_P: P,
      newton_UF: QF0,
      newton_R,
      newton_V: V,
      newton_dX,

      regularized_R0: R0,
      regularized_dX: X,
      regularized_QF: QF
    } = this;

    const L = Math.min(M,N);      

    if( 0===λ )
    {
      for( let i=N; i-- > 0; )
        X[i] = newton_dX[i];

      if( rnk < N )
      {
        // UNDER-DETERMINED/RANK-DEFICIENT CASE -> USE COMPLETE ORTHOGONAL DECOMP. FROM computeNewton()

        // COMPUTE DISTANCE
        const  r = this.scaledNorm(X); // <- the scaled length
        if(0===r)
          return [0,0];

        // COMPUTE DISTANCE DERIVATIVE dr/dλ, see:
        //  "THE LEVENBERG-MARQUARDT ALGORITHM: IMPLEMENTATION AND THEORY"
        //   by Jorge J. Moré
        //   Chapter 5 "The Levenberg-Marquardt Parameter", Equation (5.8)
        Y.fill(0.0, 0,rnk);
        for( let i=rnk; i-- > 0; )
        for( let j= N ; j-- > 0; )
          Y[i] += X[j]*D[j]*V[N*i+j];

        _triu_t_solve(rnk,L,1, newton_R,0, Y,0);

        let dr = 0;
        for( let i=rnk; i-- > 0; ) {
          const d = Y[i];
          dr += d*d;
        } dr /= -r;

        return [r,dr];
      }
    }

    if( !(0 <= λ) )
      throw new Error('TrustRegionSolverLSQ.prototype.computeMinRegularized(λ): λ must be positive.');

    // START AT RRQR FROM computeNewton() (by copying R0 -> T)
    for( let i=rnk*N; i-- > 0; )
      T[i] = R0[i];

    for( let i=rnk; i-- > 0; )
      QF[i] = -QF0[i];

    // TODO scale T to norm(T,'fro')=1 to avoid underflow ?

    if( 0 !== λ )
    {
      T .fill(0.0, rnk*N,(rnk+N)*N);
      QF.fill(0.0, rnk,   rnk+N   );

      const λSqrt = Math.sqrt(λ);

      for( let i=N; i-- > 0; )
      { let    Dλ = D[P[i]];
        if(0===Dλ)
               Dλ = 1;
        else   Dλ *= λSqrt;
        if( !( Dλ > 0 ) ) throw new Error('TrustRegionSolverLSQ.prototype.computeMinRegularized(λ): λ too small (caused underflow).');

        const       j = (i-rnk+N) % N;
        T[rnk*N + N*j+i] = Dλ;

        // TODO: if we were to eliminate each regularization row right away, TMP could be smaller
      }

      // COMPLETE QR DECOMPOSITION
      for( let i=0; i < N; i++ ) { const J = Math.min(i+1,rnk)+N;
      for( let j=N; j < J; j++ )
      {
        // USE GIVENS ROTATION TO ELIMINATE ELEMENT R_ji
        const ji = N*j+i, T_ji = T[ji]; if(0 === T_ji) continue;
        const ii = N*i+i, T_ii = T[ii],
         [c,s,norm] = _giv_rot_qr(T_ii,T_ji);
            T[ji]= 0; if(0 === s) continue;
            T[ii]= norm;
        _giv_rot_rows(T, N-1-i, ii+1,
                                ji+1, c,s);
        _giv_rot_rows(QF, 1, i,j, c,s);
      }}
    }
    else if( rnk !== N )
      throw new Error('Assertion failed.');

    for( let i=N; i-- > 0; )
      Y[i] = QF[i];

    if( 0 !== λ )
    {
      _triu_solve(N,N,1, T,0, Y,0);

      for( let i=N; i-- > 0; )
        X[P[i]] = Y[i];
    }

    // COMPUTE DISTANCE
    const  r = this.scaledNorm(X); // <- the scaled length
    if(0===r)
      return [0,0];

    // COMPUTE DISTANCE DERIVATIVE dr/dλ, see:
    //  "THE LEVENBERG-MARQUARDT ALGORITHM: IMPLEMENTATION AND THEORY"
    //   by Jorge J. Moré
    //   Chapter 5 "The Levenberg-Marquardt Parameter", Equation (5.8)
    for( let i=N; i-- > 0; ) {
      const j = P[i];
      Y[i] = X[j]*D[j]*D[j];
    }

    _triu_t_solve(N,N,1, T,0, Y,0);

    let dr = 0;
    for( let i=N; i-- > 0; ) {
      const d = Y[i];
      dr += d*d;
    } dr /= -r;

    return [r,dr];
  }
}
