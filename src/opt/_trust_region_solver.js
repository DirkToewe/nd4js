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

import {NDArray} from '../nd_array'

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


/** 
 */
export class TrustRegionSolverLSQ
{
  constructor( M, N )
  {
    if( M%1 !== 0 ) throw new Error('Assertion failed.');
    if( N%1 !== 0 ) throw new Error('Assertion failed.');
    if( !(0 <  M) ) throw new Error('Assertion failed.');
    if( !(0 <  N) ) throw new Error('Assertion failed.');
    this.M = M|0;
    this.N = N|0;

    const L = Math.min(M,N),
          K = Math.max(L+N, M);

    this.P = new Int32Array(N); // <- column permutation by RRQR
    this.Q = new Int32Array(N); // <- column permutation by RRQR

    this.D = new Float64Array(N), // <- trust region scaling diagonal

    this.R = new Float64Array( L*N ); // <- stores result of RRQR
    this.S = new Float64Array( L*L ); // <- stores complete orthogonal decomp R (in rank-deficient case)
    this.T = new Float64Array( K*N ); // <- working memory for decompositions (RRQR/URV)

    this.V = new Float64Array(N*N); // <- stores complete orthogonal decomp V TODO: In unlikely case (M < N), economic SRRQR could be used instead of full

    this.F = new Float64Array( M*1 );
    this.Y = new Float64Array( K*1 ); // <- (L+N)*1 should be enough ... right?

    this.X = new Float64Array(N*1); // <- stores solution of current min. (global or regularized)
    this.Z = new Float64Array(N*1); // <- backup for the global min.

    this.norm = new Float64Array(2*N); // <- couldn't Y be used instead?

    this.G = new Float64Array(N);

    this.rank = -1;

    this._state = 0;
    Object.seal(this);
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

  cauchyPointTravel()
  {
    if( this._state !== 1 ) throw new Error('TrustRegionSolverLSQ.prototype.cauchyPointTravel(): can only be called directly after update(f,J) call.');

    const {M,N, X, F,T: J,G} = this;

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

// METHODS NEED TO BE CALLED IN THE FOLLOWING ORDER:
//   1) update(f,J)
//   2) computeMinGlobal()
//   3) computeMinRegularized() (can be called multiple times)
// computeMinGlobal() and computeMinRegularized() store their results
// in X overwriting the previous result.

  update( f, J )
  {
    this._state = 1;

    const {M,N, D, P,Q,T,F, G, norm} = this;

    // ASSERTIONS
    if( ! (f instanceof NDArray) ) throw new Error('Assertion failed.');
    if( ! (J instanceof NDArray) ) throw new Error('Assertion failed.');

    if( f.ndim !== 1 ) throw new Error('Assertion failed.');
    if( J.ndim !== 2 ) throw new Error('Assertion failed.');

    if( f.shape[0] !== M ) throw new Error('Assertion failed.');
    if( J.shape[0] !== M ) throw new Error('Assertion failed.');
    if( J.shape[1] !== N ) throw new Error('Assertion failed.');

    f = f.data;
    J = J.data;

    // COPY f -> F
    for( let i=M*1; i-- > 0; ) F[i] = f[i];

    // COPY J -> T
    for( let i=M*N; i-- > 0; ) T[i] = J[i];

    // INIT P
    for( let i=  N; i-- > 0; ) P[i] =   i;

    // UPDATE SCALING
    // compute column norms of J
    norm.fill(0);
    for( let i=0; i < M; i++ )
      _norm_update(norm, J, N*i, 0);

    // use column norms of J to update D
    for( let j=0; j < N; j++ )
      D[j] = Math.max( D[j], _norm(norm,j) );

    // COMPUTE GRADIENT OF THE (HALF) SQUARED ERROR (MIND THE HALF :P)
    G.fill(0.0);
    for( let i=M; i-- > 0; )
    for( let j=N; j-- > 0; )
      G[j] += J[N*i+j] * F[i];
  }

  computeMinGlobal()
  {
    if( this._state !== 1 &&
        this._state !== 2 &&
        this._state !== 3 ) throw new Error('TrustRegionSolverLSQ.prototype.computeMinGlobal(): can only be called once after every update(f,J) call.');

    const {M,N, P,Q, D, R,S,T, V, X,Z, F,Y, norm} = this;

    if( this._state === 1 )
    {
      const L = Math.min(M,N);

      // SOLVE GLOBAL MINIMUM
      _rrqr_decomp_inplace(M,N,1, T,0, F,0, P,0, norm);

      // backup the rrqr decomposition in R (for finding regularized minima)
      for( let i=L*N; i-- > 0; )
        R[i] = T[i];

      const rank =
      this.rank = _rrqr_rank(M,N, T,0, /*tmp=*/norm);

      if( rank < N )
      { // RANK DEFICIENT CASE -> USE COMPLETE ORTHOGONAL DECOMPOSITION

        // factor in scaling
        for( let i=L; i-- > 0; )
        for( let j=N; j-- > 0; ) {
          const         d = D[P[j]];
          if(     0 !== d )
            T[N*i+j] /= d;
        }

        // complete orthogonal decompositon
        for( let i=N; i-- > 0; )
          Q[i] = P[i];
        V.fill(0.0);
        _urv_decomp_full( rank,N, T,0, V,0, Q,0 );

        // backup complete orthogonal T to S
        for( let i=rank; i-- > 0; )
        for( let j=rank; j-- > 0; )
          S[L*i+j] = T[N*i+j];

        // compute least squares solution
        for( let i=rank; i-- > 0; )
          Y[i] = -F[i];

        _triu_solve(rank,L,1, S,0, Y,0); // Y = S \ Y

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
          Y[i] = -F[i];

        _triu_solve(N,N,1, T,0, Y,0);

        for( let i=N; i-- > 0; )
          X[P[i]] = Y[i];
      }

      for( let i=N; i-- > 0; )
        Z[i] = X[i]; // <- backup global min. to Z
    }
    else if( this._state === 3 ) {
      for( let i=N; i-- > 0; )
        X[i] = Z[i]; // <- read backup of gloabl min. from Z
    }

    this._state = 2;
  }

  // Regularized least squares (RLS) solver. Solves the follwing equation system:
  //
  //   (JᵀJ + λI)x = Jᵀf
  //
  // The equation above can be solved as the following least squares problem:
  //
  //   min ║Ax - z║₂
  //    x
  //
  // where
  //       ┏       ┓        ┏   ┓
  //       ┃       ┃        ┃   ┃
  //       ┃       ┃        ┃   ┃
  //       ┃   J   ┃        ┃-f ┃
  //       ┃       ┃        ┃   ┃
  //       ┃       ┃        ┃   ┃
  //   A = ┃┄┄┄┄┄┄┄┃    z = ┃┄┄┄┃
  //       ┃√λ     ┃        ┃ 0 ┃
  //       ┃   ⋱   ┃        ┃ ⋮ ┃
  //       ┃    √λ ┃        ┃ 0 ┃
  //       ┗       ┛        ┗   ┛
  //
  // On top of that we can utilize the fact that J has already
  // been RRQR decomposed in computeMinGlobal().
  computeMinRegularized( λ )
  {
    if(  this._state !== 2
      && this._state !== 3 ) throw new Error('TrustRegionSolverLSQ.prototype.computeMinRegularized(λ): can only be (repeatedly) called after computeMinGlobal().');

    if( 0 !== λ )
      this._state = 3;

    const {M,N, P, D, R,S,T, V, X, F,Y, rank: rnk} = this;

    const L = Math.min(M,N);

    if( 0===λ && rnk < N )
    {
      // UNDER-DETERMINED/RANK-DEFICIENT CASE -> USE COMPLETE ORTHOGONAL DECOMPOSTION
      this.computeMinGlobal();

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

      _triu_t_solve(rnk,L,1, S,0, Y,0);

      let dr = 0;
      for( let i=rnk; i-- > 0; ) {
        const d = Y[i];
        dr += d*d;
      } dr /= -r;

      return [r,dr];
    }

    if( !(0 <= λ) )
      throw new Error('TrustRegionSolverLSQ.prototype.computeMinRegularized(λ): λ must be positive.');

    for( let i=rnk; i-- > 0; )
      Y[i] = -F[i];
    Y.fill(0.0, rnk,rnk+N);

    // COPY R->T TO REUSE THE RRQR OF J
    for( let i=rnk*N; i-- > 0; )
      T[i] = R[i];
    T.fill(0.0, rnk*N,(rnk+N)*N);

    // TODO scale T to norm(T,'fro')=1 to avoid underflow ?

    if( 0 !== λ )
    {
      const λSqrt = Math.sqrt(λ);

      for( let i=N; i-- > 0; ) {
        let Dλ = D[P[i]];
        if( Dλ===0 )
            Dλ = 1;
        else
            Dλ *= λSqrt;

        if( !(0 < Dλ) ) throw new Error('TrustRegionSolverLSQ.prototype.computeMinRegularized(λ): λ too small (caused underflow).');

        const       j = (i-rnk+N) % N;
        T[rnk*N + N*j+i] = Dλ;
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
        _giv_rot_rows(Y, 1, i,j, c,s);
      }}
    }
    else if( rnk !== N )
      throw new Error('Assertion failed.');

    if( 0===λ )
      this.computeMinGlobal();
    else {
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
