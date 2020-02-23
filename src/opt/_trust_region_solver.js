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

import {_giv_rot_rows} from '../la/_giv_rot'
import {FrobeniusNorm} from '../la/norm'
import {_rrqr_rank,
        _rrqr_decomp_inplace,
        _norm_update,
        _norm} from '../la/rrqr'
import {_triu_solve} from '../la/tri'
import {_urv_decomp_full} from '../la/urv'


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

    const L = Math.min(M,N);

    this.P = new Int32Array(N); // <- column permutation by RRQR
    this.Q = new Int32Array(N); // <- column permutation by RRQR

    this.D = new Float64Array(N), // <- trust region scaling diagonal
    this.D.fill(Number.MIN_VALUE);

    this.R = new Float64Array(     L*N ); // <- stores result of RRQR
    this.T = new Float64Array( (M+N)*N ); // <- used to compute URV decomp (for rank-deficient least-squares to find global min) and regularized QR decomp (used to find regularized min.).

//    this.V = this.T.subarray(-N*N); // <- TODO: In unlikely case that M < N, use economic SRRQR instead of full SRRQR
    this.V = new Float64Array(N*N); // <- TODO: In unlikely case that M < N, use economic SRRQR instead of full SRRQR

    this.F = new Float64Array(     M*1 );
    this.Y = new Float64Array( (M+N)*1 );

    this.X = new Float64Array(N*1);

    this.norm = new Float64Array(2*N);

    this.G = new Float64Array(N);

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

    if( !(0 <= a) ) throw new Error(`Assertion failed.`);
    if( !(0 <= b) ) throw new Error(`Assertion failed.`);

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
      G[j] += T[N*i+j] * F[i];
  }

  computeMinGlobal()
  {
    if( this._state !== 1 ) throw new Error('TrustRegionSolverLSQ.prototype.computeMinGlobal(): can only be called once after every update(f,J) call.');
    this._state = 2;

    const {M,N, P,Q, D, R,T,V,X, F,Y, norm} = this;

    const L = Math.min(M,N);

    // SOLVE GLOBAL MINIMUM
    _rrqr_decomp_inplace(M,N,1, T,0, F,0, P,0, norm);

    // backup the rrqr decomposition in R (for finding regularized minima)
    for( let i=L*N; i-- > 0; )
      R[i] = T[i];

    const rank = _rrqr_rank(M,N, T,0, /*tmp=*/norm);

    if( rank < N )
    { // RANK DEFICIENT CASE -> USE COMPLETE ORTHOGONAL DECOMPOSITION

      // factor in scaling
      for( let i=L; i-- > 0; )
      for( let j=N; j-- > 0; )
        T[N*i+j] /= D[P[j]];

      // complete orthogonal decompositon
      for( let i=N; i-- > 0; )
        Q[i] = P[i];
      V.fill(0.0);
      _urv_decomp_full( rank,N, T,0, V,0, Q,0 );

      // compute least squares solution
      for( let i=N; i-- > 0; )
        Y[i] = -F[i];

      _triu_solve(rank,N,1, T,0, Y,0); // Y = T \ Y

      X.fill(0.0);
      for( let k=rank; k-- > 0; )
      for( let i=N   ; i-- > 0; )
        X[i] += V[N*k+i] * Y[k]; // X = V.T @ Y

      // factor out scaling
      for( let i=N; i-- > 0; )
        X[i] /= D[i];
    }
    else
    { // FULL RANK CASE -> SCALING NOT RELEVANT TO GLOBAL SOLUTION
      for( let i=N; i-- > 0; )
        norm[i] = -F[i];

      _triu_solve(N,N,1, T,0, /*tmp=*/norm,0);

      for( let i=N; i-- > 0; )
        X[P[i]] = /*tmp[i]=*/norm[i];
    }
  }

  computeMinRegularized( λ )
  {
    if(  this._state !== 2
      && this._state !== 3 ) throw new Error('TrustRegionSolverLSQ.prototype.computeMinRegularized(λSqrt): can only be (repeatedly) called after computeMinGlobal().');
    this._state = 3;

    if( !(0 < λ) )
      throw new Error('TrustRegionSolverLSQ.prototype.computeMinRegularized(λSqrt): λSqrt must be positive.');

    const {M,N, P, D, R,T,X, F,Y, norm} = this;

    const L = Math.min(M,N);

    for( let i=0; i < M; i++ )
      Y[i] = -F[i];
    Y.fill(0.0, M);

    for( let i=L*N; i-- > 0; )
      T[i] = R[i];
    T.fill(0.0, L*N);

    const λSqrt = Math.sqrt(λ);

    for( let i=N; i-- > 0; ) {
      const Dλ = D[P[i]] * λSqrt;

      if( !(0 < Dλ) ) throw new Error('TrustRegionSolverLSQ.prototype.computeMinRegularized(λSqrt): λSqrt too small (caused underflow).');

      T[M*N + N*i+i] = Dλ;
    }

    // COMPLETE QR DECOMPOSITION
    for( let i=0              ; i <  N  ; i++ )
    for( let j=Math.max(M,i+1); j <= M+i; j++ )
    {
      // USE GIVENS ROTATION TO ELIMINATE ELEMENT R_ji
      const ji = N*j+i, T_ji = T[ji]; if(0 === T_ji) continue;
      const ii = N*i+i, T_ii = T[ii],
                 norm = Math.hypot(T_ji,T_ii),
      c = T_ii / norm,
      s = T_ji / norm;
          T[ji]= 0; if(0 === s) continue;
          T[ii]= norm;
      _giv_rot_rows(T, N-1-i, ii+1,
                              ji+1, c,s);
      _giv_rot_rows(Y, 1, i,j, c,s);
    }

    _triu_solve(N,N,1, T,0, Y,0);

    for( let i=N; i-- > 0; )
      X[P[i]] = Y[i];

    // COMPUTE DISTANCE
    const r = this.scaledNorm(X); // <- the scaled length

    // COMPUTE DISTANCE DERIVATIVE, see:
    //  "THE LEVENBERG-MARQUARDT ALGORITHM: IMPLEMENTATION AND THEORY"
    //   by Jorge J. Moré
    //   Chapter 5 "The Levenberg-Marquardt Parameter", Eq. (5.8)
    norm.fill(0.0, 0,N);
    for( let i=N; i-- > 0; ) {
      const j = P[i];
      norm[i] = X[j]*D[j]*D[j];
    }

    // FORWARD SUBSTITUTION
    for( let i=0; i < N; i++ )
    {                          norm[i] /= T[N*i+i];
      for( let j=i; ++j < N; ) norm[j] -= T[N*i+j] * norm[i];
    }

    let dr = 0;
    for( let i=N; i-- > 0; ) {
      const d = norm[i];
      dr += d*d;
    } dr /= -r;

    return [r,dr];
  }
}
