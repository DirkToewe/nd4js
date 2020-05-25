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

import {ARRAY_TYPES} from '../dt'
import {asarray, NDArray} from '../nd_array'
import {srrqr_decomp_full} from './srrqr'
import {_giv_rot_qr,
        _giv_rot_rows} from './_giv_rot'
import {_triu_solve} from './tri'


// TODO: add economic URV decomposition


export function _urv_decomp_full( M,N, R,R_off, V,V_off, P,P_off )
{
  if( 0 !== M%1 ) throw new Error('Assertion failed.');
  if( 0 !== N%1 ) throw new Error('Assertion failed.');
  if( 0 !== R_off%1 ) throw new Error('Assertion failed.');
  if( 0 !== V_off%1 ) throw new Error('Assertion failed.');

/*DEBUG*/  for( let i=N; i-- > 0; )
/*DEBUG*/  for( let j=N; j-- > 0; )
/*DEBUG*/    if( 0 !== V[V_off + N*i+j] ) throw new Error('Assertion failed.');

  M |= 0;
  N |= 0;
  R_off |= 0;
  V_off |= 0;

  if(!(M <=N) ) throw new Error('Assertion failed.');
  if(  M===N  )
  {
    for( let i=N; i-- > 0; ) {
      const         j = P[P_off+i];
      V[V_off + N*i+j] = 1;
    }
  }
  else
  {
    // INIT V TO IDENTITY
    for( let i=N; i-- > 0; )
      V[V_off + N*i+i] = 1;

    for( let i=M; i-- > 0;     )
    for( let j=M; j   < N; j++ )
    { const ij = R_off + N*i+j, R_ij = R[ij]; if(0 === R_ij) continue;
      const ii = R_off + N*i+i, R_ii = R[ii];
      const [c,s,norm] = _giv_rot_qr(R_ii,R_ij);
      if( s !== 0 )
      { // ROT COLUMNS IN R
        for( let k=i; k-- > 0; ) // <- TODO can this be made cache friendlier?
        { const  ki = R_off + N*k+i, R_ki = R[ki],
                 kj = R_off + N*k+j, R_kj = R[kj];
            R[ki] = R_kj*s + R_ki*c;
            R[kj] = R_kj*c - R_ki*s;
        }
        // ROT ROWS IN V
        _giv_rot_rows(V, 1+j-i, V_off + N*i+i,
                                V_off + N*j+i, c,s);
        R[ii] = norm;
      } R[ij] = 0;
    }

    // apply colum permutations
    // nested loop but actually only O(N) swaps
    for( let i=0; i < N; ++i )
    for(;;) // <- start swap cycle
    { const j = P[P_off+i];
                P[P_off+i] = P[P_off+j];
                             P[P_off+j] = j;
      if( j <= i )
        break;
      // SWAP COLUMNS IN V
      for( let k=N; k-- > 0; )
      { const tmp = V[V_off + N*k+i];
                    V[V_off + N*k+i] = V[V_off + N*k+j];
                                       V[V_off + N*k+j] = tmp;
      }
    }
  }
}


export function urv_decomp_full( A ) // <- TODO add tolerance parameters to pass on to SRRQR
{
  let [U, RR, P, rnks] = srrqr_decomp_full(A);
                                           A = undefined;
  const [M,N] = RR.shape.slice(-2),
      V_shape = RR.shape.slice();
      V_shape[V_shape.length-2] = N;

        P =    P.data;
  const r = rnks.data,
        R =   RR.data,
        V = new ARRAY_TYPES[RR.dtype](R.length/M*N);

  for( let
    R_off=0,
    V_off=0,
    P_off=0,
    r_off=0;     r_off < r.length;
    R_off += M*N,
    V_off += N*N,
    P_off +=   N,
    r_off +=   1
  )
  {
    const rnk = r[r_off];
    if( rnk < N )
    { // set close-to-zero lower right region to zero
      for( let i=rnk; i < M; i++ )
      for( let j=rnk; j < N; j++ )
        R[R_off + N*i+j] = 0; // <- TODO consider doing this in SRRQR already
    }
    _urv_decomp_full( rnk,N, R,R_off, V,V_off, P,P_off ); // <- TODO if we init V to identity then we can do faster rotations of V. P would then be applied after computing V.
  }

  return [U, RR, new NDArray(V_shape, V), rnks];
}


export function _urv_lstsq( rnk, I,J,K,L,M, U,U_off, R,R_off, V,V_off, X,X_off, Y,Y_off, tmp )
{
  //   ┏        ┓ ┏        ┓ ┏        ┓ ┏        ┓   ┏        ┓
  //   ┃        ┃ ┃        ┃ ┃        ┃ ┃        ┃   ┃        ┃
  //   ┃ U[I,J] ┃ ┃ R[J,K] ┃ ┃ V[K,L] ┃ ┃ X[L,M] ┃ = ┃ Y[I,M] ┃
  //   ┃        ┃ ┃        ┃ ┃        ┃ ┃        ┃   ┃        ┃
  //   ┗        ┛ ┗        ┛ ┗        ┛ ┗        ┛   ┗        ┛

  if(      rnk%1 !== 0 ) throw new Error('Assertion failed.');
  if(        I%1 !== 0 ) throw new Error('Assertion failed.');
  if(        J%1 !== 0 ) throw new Error('Assertion failed.');
  if(        K%1 !== 0 ) throw new Error('Assertion failed.');
  if(        L%1 !== 0 ) throw new Error('Assertion failed.');
  if(        M%1 !== 0 ) throw new Error('Assertion failed.');
  if(    U_off%1 !== 0 ) throw new Error('Assertion failed.');
  if(    R_off%1 !== 0 ) throw new Error('Assertion failed.');
  if(    V_off%1 !== 0 ) throw new Error('Assertion failed.');
  if(    X_off%1 !== 0 ) throw new Error('Assertion failed.');
  if(    Y_off%1 !== 0 ) throw new Error('Assertion failed.');

  if( ! (0 <= I) ) throw new Error('Assertion failed.');
  if( ! (0 <= J) ) throw new Error('Assertion failed.');
  if( ! (0 <= K) ) throw new Error('Assertion failed.');
  if( ! (0 <= L) ) throw new Error('Assertion failed.');
  if( ! (0 <= M) ) throw new Error('Assertion failed.');

  if( ! (0 <= U_off) ) throw new Error('Assertion failed.');
  if( ! (0 <= R_off) ) throw new Error('Assertion failed.');
  if( ! (0 <= V_off) ) throw new Error('Assertion failed.');
  if( ! (0 <= X_off) ) throw new Error('Assertion failed.');
  if( ! (0 <= Y_off) ) throw new Error('Assertion failed.');

  if( ! (U_off <= U.length-I*J) ) throw new Error('Assertion failed.');
  if( ! (R_off <= R.length-J*K) ) throw new Error('Assertion failed.');
  if( ! (V_off <= V.length-K*L) ) throw new Error('Assertion failed.');
  if( ! (X_off <= X.length-L*M) ) throw new Error('Assertion failed.');
  if( ! (Y_off <= Y.length-I*M) ) throw new Error('Assertion failed.');

  if( tmp.length < I*M ) throw new Error('Assertion failed.');
  tmp.fill(0.0, 0,I*M);

  // TMP = U.T @ Y
  for( let k=0; k <  I ; k++ )
  for( let i=0; i < rnk; i++ )
  for( let j=0; j <  M ; j++ )
    tmp[M*i+j] += U[U_off + J*k+i] * Y[Y_off + M*k+j];

  // TMP = R \ TMP
  _triu_solve(rnk,K,M, R,R_off, tmp,0);

  // TMP = V.T @ TMP
  for( let k=0; k < rnk; k++ )
  for( let i=0; i <  L ; i++ )
  for( let j=0; j <  M ; j++ )
    X[X_off + M*i+j] += V[V_off + L*k+i] * tmp[M*k+j];
}


export function urv_lstsq( U,R,V, ranks, Y )
{
  if( Y == null )
  {
    if(  null != ranks
      || null != V )
      throw new Error('urv_lstsq( U,R,V,ranks, Y ): Either 2 ([U,R,V,ranks], Y) or 5 arguments (U,R,V,ranks, Y) expected.');
    Y = R;
    ([U,R,V,ranks] = U);
  }

  U = asarray(U); if( ! (U.ndim >= 2) ) throw new Error('urv_lstsq(U,R,V, Y): U.ndim must be at least 2.')
  R = asarray(R); if( ! (R.ndim >= 2) ) throw new Error('urv_lstsq(U,R,V, Y): R.ndim must be at least 2.')
  V = asarray(V); if( ! (V.ndim >= 2) ) throw new Error('urv_lstsq(U,R,V, Y): V.ndim must be at least 2.')
  Y = asarray(Y); if( ! (Y.ndim >= 2) ) throw new Error('urv_lstsq(U,R,V, Y): Y.ndim must be at least 2.')
  ranks = asarray(ranks);

  const U_shape =     U.shape,
        R_shape =     R.shape,
        V_shape =     V.shape,
        Y_shape =     Y.shape,
    ranks_shape = ranks.shape;

  const dtype = [U,R,V,Y].every(x => x.dtype==='float32') ? 'float32' : 'float64';

  const DTypeArray = ARRAY_TYPES[dtype];

  const ndim = [U,R,V,Y].reduce(
    (max,{ndim}) => Math.max(max,ndim),
    ranks.ndim+2
  );

  const X_shape = new Int32Array(ndim);
        X_shape.fill(1);

  // FIND COMMON (BROADCASTED) SHAPE
  for( let arr of [U,R,V,Y] )
    for( let i=ndim-2, j=arr.ndim-2; i-- > 0 && j-- > 0; )
      if( 1 === X_shape[i] )
        X_shape[i] = arr.shape[j];
      else if( X_shape[i] != arr.shape[j] && arr.shape[j] != 1 )
        throw new Error('urv_lstsq( U,R,V,ranks, Y ): U,R,V,ranks, Y not broadcast-compatible.');

  for( let i=ndim-2, j=ranks.ndim; i-- > 0 && j-- > 0; )
    if( 1 === X_shape[i] )
      X_shape[i] = ranks.shape[j];
    else if( X_shape[i] != ranks.shape[j] && ranks.shape[j] != 1 )
      throw new Error('urv_lstsq( U,R,V,ranks, Y ): U,R,V,ranks, Y not broadcast-compatible.');

  // CHECK MATRIX SHAPES
  //   ┏        ┓ ┏        ┓ ┏        ┓ ┏        ┓   ┏        ┓
  //   ┃        ┃ ┃        ┃ ┃        ┃ ┃        ┃   ┃        ┃
  //   ┃ U[I,J] ┃ ┃ R[J,K] ┃ ┃ V[K,L] ┃ ┃ X[L,M] ┃ = ┃ Y[I,M] ┃
  //   ┃        ┃ ┃        ┃ ┃        ┃ ┃        ┃   ┃        ┃
  //   ┗        ┛ ┗        ┛ ┗        ┛ ┗        ┛   ┗        ┛
  const [I,J]= U.shape.slice(-2),
        [K,L]= V.shape.slice(-2),
           M = Y.shape[Y.ndim-1];

  if( R.shape[R.ndim-2] !== J ) throw new Error('urv_lstsq( U,R,V,ranks, Y ): Matrix dimensions incompatible.');
  if( R.shape[R.ndim-1] !== K ) throw new Error('urv_lstsq( U,R,V,ranks, Y ): Matrix dimensions incompatible.');
  if( Y.shape[Y.ndim-2] !== I ) throw new Error('urv_lstsq( U,R,V,ranks, Y ): Matrix dimensions incompatible.');

  if( J !== K ) {
    if( I < L ) if( I !== J ) throw new Error('Assertion failed.');
    else        if( K !== L ) throw new Error('Assertion failed.');
  }
  if( !(I >= J) ) throw new Error('Assertion failed.');
  if( !(K <= L) ) throw new Error('Assertion failed.');

  X_shape[ndim-2] = L;
  X_shape[ndim-1] = M;

  ranks = ranks.data;
  U = U.data;
  R = R.data;
  V = V.data;
  Y = Y.data;
  const X = new DTypeArray( X_shape.reduce((m,n) => m*n) ),
      tmp = new DTypeArray( I*M );

  let   U_off =     U.length,     U_stride = 1,
        R_off =     R.length,     R_stride = 1,
        V_off =     V.length,     V_stride = 1,
        Y_off =     Y.length,     Y_stride = 1,
    ranks_off = ranks.length, ranks_stride = 1,
        X_off =     X.length;

  function solv(d) {
    if( d === ndim-2 ) {
          U_stride = I*J;
          R_stride = J*K;
          V_stride = K*L;
          Y_stride = I*M;
      ranks_stride = 1;

          U_off -=     U_stride;
          R_off -=     R_stride;
          V_off -=     V_stride;
          Y_off -=     Y_stride;
      ranks_off -= ranks_stride;
          X_off -= L*M;

      const rnk = ranks[ranks_off];

      _urv_lstsq( rnk, I,J,K,L,M, U,U_off, R,R_off, V,V_off, X,X_off, Y,Y_off, tmp )
    }
    else {
      for( let l=X_shape[d]; ; l-- ) {
        solv(d+1);
        if( l == 1 ) break;
        if( ! (    U_shape[ d - ndim   +     U_shape.length ] > 1) )     U_off +=     U_stride;
        if( ! (    R_shape[ d - ndim   +     R_shape.length ] > 1) )     R_off +=     R_stride;
        if( ! (    V_shape[ d - ndim   +     V_shape.length ] > 1) )     V_off +=     V_stride;
        if( ! (    Y_shape[ d - ndim   +     Y_shape.length ] > 1) )     Y_off +=     Y_stride;
        if( ! (ranks_shape[ d - ndim+2 + ranks_shape.length ] > 1) ) ranks_off += ranks_stride;
      }
          U_stride *=     U_shape[ d - ndim   +     U_shape.length ] || 1;
          R_stride *=     R_shape[ d - ndim   +     R_shape.length ] || 1;
          V_stride *=     V_shape[ d - ndim   +     V_shape.length ] || 1;
          Y_stride *=     Y_shape[ d - ndim   +     Y_shape.length ] || 1;
      ranks_stride *= ranks_shape[ d - ndim+2 + ranks_shape.length ] || 1;
    }
  }
  solv(0);

  return new NDArray(X_shape, X);
}
