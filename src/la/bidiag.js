'use strict';

/* This file is part of ND.JS.
 *
 * ND.JS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ND.JS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ND.JS. If not, see <http://www.gnu.org/licenses/>.
 */

import {asarray, NDArray} from '../nd_array'
import {super_dtype, ARRAY_TYPES} from '../dt'
import {_giv_rot_rows} from './_giv_rot'
import {FrobeniusNorm} from './norm'
import {_transpose_inplace} from './transpose_inplace'


function _bidiag_decomp_vert( M,N, U,U_off, B,V,BV_off )
{
       M |= 0;
       N |= 0;
   U_off |= 0;
  BV_off |= 0;

  const NORM = new FrobeniusNorm();

  if( M < N ) throw new Error('Assertion failed.');

  // init V to identity
  for( let i=N; i-- > 0; ) V[BV_off + N*i+i] = 1;

  for( let i=0; i < N; i++ )
  { // ELIMINATE ELEMENTS BELOW (i,i) USING GIVENS
                               const ii = U_off + N*i+i;
    for( let j=i; ++j < M; ) { const ji = U_off + N*j+i;
      const   B_ji = U[ji];
      if( 0===B_ji ) continue;
      const   B_ii = U[ii];
      let            norm = Math.hypot(B_ii,B_ji),
          c = B_ii / norm,
          s = B_ji / norm;
      if( s !== 0 ) { // <- might happen in case of underflow
        if( c  <  0 ) {
            c *= -1;
            s *= -1;
         norm *= -1;
        }
        _giv_rot_rows(U, N-1-i, ii+1,
                                ji+1, c,s);
        U[ii] = norm;
      } U[ji] = s;
    }
    if( i < N-2 )
    { // ELIMINATE ELEMENTS RIGHT OF (i,i+1) USING HOUSEHOLDER
      NORM.reset();
      for( let j=N; --j > i+1; )
          NORM.include(U[U_off + N*i+j]);
      if( NORM.max === 0 ) continue;
      const norm = NORM.resultIncl(U[ii+1]) * (U[ii+1] > 0 ? -1 : +1),
             div = NORM.resultIncl(U[ii+1] -= norm);
      for( let j=i; ++j < N; )
        U[U_off + N*i+j] /= div;
      // apply householder to right of V
      for( let j=0; j < N; j++ ) {
        let sum = 0; for(let k=i; ++k < N;) sum += V[BV_off + N*j+k] * U[U_off + N*i+k];
            sum*= 2; for(let k=i; ++k < N;)        V[BV_off + N*j+k]-= U[U_off + N*i+k] * sum;
      }
      // apply householder to right of B
      for( let j=i; ++j < M; ) {
        let sum = 0; for(let k=i; ++k < N;) sum += U[U_off + N*j+k] * U[U_off + N*i+k];
            sum*= 2; for(let k=i; ++k < N;)        U[U_off + N*j+k]-= U[U_off + N*i+k] * sum;
      }
                U[ii+1] = norm;
      U.fill(0.0, ii+2, U_off+N*(i+1));
    }
  }
  _transpose_inplace(N, V,BV_off);
  // COPY U -> B
  for( let i=N; --i > 0; ) {
    B[BV_off + N* i   +i] = U[U_off + N* i   +i];
    B[BV_off + N*(i-1)+i] = U[U_off + N*(i-1)+i];
  } B[BV_off            ] = U[U_off];
  // COMPUTE U
  // init upper right triangle of U
  for( let i=0; i < N; i++ )
  for( let j=i; j < N; j++ ) U[U_off + N*i+j] = +(i===j);

  for( let i=N; i-- > 0; )
  for( let j=M; --j > i; ) {
    const s = U[U_off + N*j+i]; if(0 === s) continue;
              U[U_off + N*j+i]  =  0;
    const c = Math.sqrt( (1+s)*(1-s) );
    _giv_rot_rows(U, N-i, U_off + N*j+i, 
                          U_off + N*i+i, c,s);
  }
}


function _bidiag_decomp_square( N, U,B,V, off )
{
  N   |= 0;
  off |= 0;
  const NORM = new FrobeniusNorm();
  // INIT U & V TO IDENTITY
  for( let i=N; i-- > 0; ) U[off + N*i+i] = 1;
  for( let i=N; i-- > 0; ) V[off + N*i+i] = 1;

  for( let i=0; i < N-1; i++ )
  { // ELIMINATE ELEMENTS BELOW (i,i) USING GIVENS
                               const ii = off + N*i+i;
    for( let j=i; ++j < N; ) { const ji = off + N*j+i;
      const  B_ji = B[ji];
      if(0===B_ji) continue;
      const  B_ii = B[ii], norm = Math.hypot(B_ii,B_ji),
                c = B_ii / norm,
                s = B_ji / norm; if(0===s) continue; // <- can happen due to underflow
                    B[ii]= norm;
                    B[ji]= 0;
      _giv_rot_rows(B, N-1-i, ii+1,
                              ji+1,      c,s);
      _giv_rot_rows(U,   1+j, off + N*i,
                              off + N*j, c,s);
    }
    // ELIMINATE ELEMENTS RIGHT OF (i,i+1) USING HOUSEHOLDER
    NORM.reset();
    for( let j=N; --j > i+1; )
        NORM.include(B[off + N*i+j]);
    if( NORM.max === 0 ) continue;
    const norm = NORM.resultIncl(B[ii+1]) * (B[ii+1] > 0 ? -1 : +1),
           div = NORM.resultIncl(B[ii+1]-= norm);
    for( let j=i; ++j < N; )
      B[off + N*i+j] /= div;
    // apply householder to right of V
    for( let j=0; j < N; j++ ) {
      let sum = 0; for(let k=i; ++k < N;) sum += V[off + N*j+k] * B[off + N*i+k];
          sum*= 2; for(let k=i; ++k < N;)        V[off + N*j+k]-= B[off + N*i+k] * sum;
    }
    // apply householder to right of B
    for( let j=i; ++j < N; ) {
      let sum = 0; for(let k=i; ++k < N;) sum += B[off + N*j+k] * B[off + N*i+k];
          sum*= 2; for(let k=i; ++k < N;)        B[off + N*j+k]-= B[off + N*i+k] * sum;
    }
              B[ii+1] = norm;
    B.fill(0.0, ii+2, off+N*(i+1));
  }
  _transpose_inplace(N, V,off);
  _transpose_inplace(N, U,off);
}


export function _bidiag_decomp_horiz( M,N, U,U_off, B,B_off, V,V_off )
{
  M              |=0;
  N              |=0;
  U_off          |=0;
  B_off          |=0;
  V_off = V_off+N| 0;

  if( M > N ) throw Error('Assertion failed.');

  const NORM = new FrobeniusNorm();
  // INIT U TO IDENTITY
  for( let i=M; i-- > 0; ) U[U_off + M*i+i] = 1;

  for( let i=0; i < M && i < N-1; i++ )
  { // ELIMINATE ELEMENTS BELOW (i,i) USING GIVENS
                               const ii = V_off + N*i+i;
    for( let j=i; ++j < M; ) { const ji = V_off + N*j+i;
      const  B_ji = V[ji];
      if(0===B_ji) continue
      const  B_ii = V[ii], norm = Math.hypot(B_ii,B_ji),
                c = B_ii / norm,
                s = B_ji / norm;
                    V[ji]= 0; if(0===s) continue;
                    V[ii]= norm;
      _giv_rot_rows(V, N-1-i, ii+1,
                              ji+1,        c,s);
      _giv_rot_rows(U,   1+j, U_off + M*i,
                              U_off + M*j, c,s);
    }
    // ELIMINATE ELEMENTS RIGHT OF (i,i+1) USING HOUSEHOLDER

    // compute householder vector
    const above = V_off + N*(i-1); // <- write householder to row (i-1)
    NORM.reset();
    for( let j=N-1;; ) {
      const        V_j = V[above+j] = V[V_off + N*i+j]; if(--j <= i ) break;
      NORM.include(V_j);
    }
    if(    0 === NORM.max  )   { V[above+(i+1)] = 0; continue; }
    const norm = NORM.resultIncl(V[above+(i+1)]) * (V[above+(i+1)] > 0 ? -1 : +1),
         denom = NORM.resultIncl(V[above+(i+1)] -= norm);
    for( let j=i; ++j < N; )     V[above+ j   ] /= denom;

    if(0===NORM.max) continue;

    // apply householder to right of V
    for( let j=i; ++j < M; ){
    let sum = 0; for( let k=i; ++k < N; ) sum += V[V_off + N*j+k] * V[above+k];
        sum*= 2; for( let k=i; ++k < N; )        V[V_off + N*j+k]-= V[above+k] * sum;
    }
    V[ii+1] = norm;
  }

  // TRANSPOSE U
  _transpose_inplace(M, U,U_off);

  // COPY V -> B
  for( let i=M; i-- > 0; ) {
    B[B_off + (M+1)*i+(i+1)] = i < N-1 ? V[V_off + N*i+(i+1)] : 0;
    B[B_off + (M+1)*i+ i   ] =           V[V_off + N*i+ i   ];
  }

  // COMPUTE V
  for( let i=Math.min(M,N-1);; ) {
         --i;
    V.fill(0.0, V_off+N*i, V_off+N*(i+1));
              V[V_off+N*i+(i+1)] = 1;

    if(0 > i) break;

    const above = V_off + N*(i-1); // <- read householder from row (i-1)
    // apply householder to right of V
    for( let j=i; j < M; j++ ) {
      let sum = 0; for(let k=i; ++k < N;) sum += V[V_off + N*j+k] * V[above+k];
          sum*= 2; for(let k=i; ++k < N;)        V[V_off + N*j+k]-= V[above+k] * sum;
    }
  }
}


// TODO consider Householder reflections for eliminating right side of a row
export function bidiag_decomp(A)
{
  A = asarray(A);
  if( A.ndim < 2                    ) throw new Error('bidiag_decomp(A): A must be at least 2D.');
  if( A.dtype.startsWith('complex') ) throw new Error('bidiag_decomp(A): complex A not yet supported.');
  const
    DType = A.dtype === 'float32' ? 'float32' : 'float64',
    DTypeArray = ARRAY_TYPES[DType], // <- ensure at least double precision
    U_shape = Int32Array.from(A.shape),
    B_shape = Int32Array.from(U_shape),
    V_shape = Int32Array.from(U_shape),
    [N,M] = A.shape.slice(-2), // <- TODO should be [M,N]
    len   = A.data.length / (M*N)| 0,
    I     = Math.min(N,M)        | 0,// #rows    of bidiag. matrix B
    J     = (N >= M ? I : I+1)   | 0;// #columns of bidiag. matrix B
  U_shape[U_shape.length-1] = I;
  B_shape[B_shape.length-2] = I;
  B_shape[B_shape.length-1] = J;
  V_shape[V_shape.length-2] = J;

  const [U,B,V] = function(){
    if( M === N ) {
      const B =     DTypeArray.from(A.data),
            U = new DTypeArray(B.length),
            V = new DTypeArray(B.length);

      for( let off=0; off < B.length; off += N*N )
        _bidiag_decomp_square(N, U,B,V, off);

      return [U,B,V];
    }
    else if( N > M ) {
      const U =     DTypeArray.from(A.data); A = undefined;
      const B = new DTypeArray(len*M*M),
            V = new DTypeArray(len*M*M);
      for(
        let U_off=0,
           BV_off=0; BV_off < B.length;
            U_off += N*M,
           BV_off +=   M*M
      ) _bidiag_decomp_vert(N,M, U,U_off, B,V,BV_off);

      return [U,B,V];
    }
    else {/*N < M*/
      const V = new DTypeArray(len*(N+1)*M);

      A = A.data;
      for( let j=0,i=M; i < V.length; i += M )
      for( let I=i+N*M; i < I;   j++, i++    ) V[i] = A[j];
      A = undefined;

      const U = new DTypeArray(len * N * N   ),
            B = new DTypeArray(len * N *(N+1));
      for(
        let U_off=0,
            B_off=0,
            V_off=0; B_off < B.length;
            U_off += N*N,
            B_off +=   N*(N+1),
            V_off +=     (N+1)*M
      ) _bidiag_decomp_horiz(N,M, U,U_off, B,B_off, V,V_off);

      return [U,B,V];
    }
  }();

  return [
    new NDArray(U_shape, U),
    new NDArray(B_shape, B),
    new NDArray(V_shape, V)
  ];
}
