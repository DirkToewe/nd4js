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
import {ARRAY_TYPES, eps} from '../dt'
import {matmul2} from './matmul'
import {qr_decomp} from './qr'
import {transpose_inplace} from './transpose_inplace'
import {_svd_jac_rot_rows,
        _svd_jac_rot_cols,
        _svd_jac_angles,
        _svd_jac_post_skip1 } from './_svd_jac_utils'


const B =  4, // <- block size should match cache line size
     BB = B*B;

// The memory of `S` is organized in a tiled order with a tile size
// of [B,B]. If for example BR=2 and BC=3, a matrix of shape [4,6]
// would have the following memory order.
// ┌                       ┐
// │  0  1    4  5    8  9 │
// │  2  3    6  7   10 11 │
// │                       │
// │ 12 13   16 17   20 21 │
// │ 14 15   18 19   22 23 │
// └                       ┘
// The purpose of this memory order is to make column (Givens) rotations
// faster.
//
// Let's say the cache line size is L. If we had a fortran-style,
// column major memory order that would correspond to a tile size of
// [L,1]. Say we have an [L,L] block and we want to access one row and
// one column of it. Accessing the column requires L cache lines to be
// read/written. Accessing the row only takes one cache line for a
// total of (L+1) cache lines.
//
// If we, on the other hand, used a tile size of [√(L),√(L)], both
// row and column access would require √(L) cache lines to be
// read/written each for a total of only 2*√(L) cache lines.
//
// ('B' as in Block)
//
// A power of two square block size is chosen such that it barely fits
// at least on cache line (64 bytes).
//
// Some quick and dirty tests have shown that is this method is faster
// for matrix sizes 800 and above.

export function svd_jac_2sided_blocked(A)
{
  A = asarray(A);
  if( A.dtype.startsWith('complex') )
    throw new Error('svd_jac_1sided(A): A.dtype must be float.');
  const shape = A.shape,
    N = shape[shape.length-2] | 0;

  // QR DECOMPOSE RECTANGULAR MATRICES AND WORK ON (QUADRATIC) R
  {
    const M = shape[shape.length-1];
    // if A is not square use QR Decomposition
    if( N > M ) {
      const [Q,R] = qr_decomp(A),
         [U,sv,V] = svd_jac_2sided_blocked(R)
      return [matmul2(Q,U), sv, V]
    }
    if( N < M ) {
      const [Q,R] = qr_decomp(A.T),
         [U,sv,V] = svd_jac_2sided_blocked(R)
      transpose_inplace(V)
      return [V, sv, matmul2(Q,U).T]
    }
  }
  // ALLOCATE RESULT DATA
  const DType = A.dtype==='float32' ? 'float32' : 'float64',
        DTypeArray = ARRAY_TYPES[DType],
        TOL = (N * eps(DType))**2,
        n = Math.ceil(N/B) | 0,
        U = DTypeArray.from(A.data); A = undefined; // <- potentially allow GC
  const S = new DTypeArray(n*n*B*B), // <- tempory storage for decomposition
        V = new DTypeArray(U.length),
       sv = new DTypeArray(U.length/N),
      ord = Int32Array.from({length: N}, (_,i) => i);

  if( 1 >  N ) throw new Error('Assertion failed.');
  if( 1 == N ) {
    for( let i=U.length; i-- > 0; )
      if( U[i]  < +0.0 ) {
          U[i] *= -1.0;
         sv[i]  = -1.0;
      }
      else sv[i] = +1.0;
    return [
      new NDArray(shape,            sv ),
      new NDArray(shape.slice(0,-1), U ),
      new NDArray(shape,      V.fill(1))
    ];
  }

  for( let UV_off=0,
           sv_off=0; sv_off < sv.length; UV_off += N*N,
                                         sv_off += N )
  {
    // MOVE FROM U TO S
    for( let i=0; i < N; i++ )
    for( let j=0; j < N; j++ ) {
      const k = BB*(n*(i/B | 0) + (j/B | 0)) + B*(i%B) + (j%B) | 0;
      S[k] = U[UV_off + N*i+j];
             U[UV_off + N*i+j] = +(i===j);
    };
    // INIT V TO IDENTITY
    for( let i=0; i < N; i++ )
    for( let j=0; j < N; j++ ) V[UV_off + N*i+j] = +(i===j);

     //
    // JACOBI SVD ITERATIONS
   //
    for( let finished = false;
           ! finished; )
    {        finished = true;
      for( let Q = n             ; Q-- >  0;    )
      for( let P = Q             ; P   >= 0; P--)
      for( let q = B             ; q-- >  0;    )
      for( let p =(P===Q) ? q : B; p-- >  0;    )
      {
        const S_pp = S[BB*(n*P+P) + B*p+p],
              S_pq = S[BB*(n*P+Q) + B*p+q],
              S_qp = S[BB*(n*Q+P) + B*q+p],
              S_qq = S[BB*(n*Q+Q) + B*q+q];
        // stopping criterion inspiredy by:
        //  "Jacobi's Method is More Accurate than QR"
        //   by James Demmel
        //   SIAM J. Matrix Anal. Appl, vol. 13, pp. 1204-1245, 1992
        if( ! ( S_pq*S_pq + S_qp*S_qp > Math.abs(S_pp*S_qq) * TOL ) )
          continue;

        finished = false;
 
        const [cα,sα, cβ,sβ] = _svd_jac_angles(S_pp,S_pq,
                                               S_qp,S_qq);

        // ROTATE ROWS IN S
        for( let i = P*n*BB + p*B,
                 j = Q*n*BB + q*B,
                 I = i + n*BB; i < I; i += BB,
                                      j += BB )
        {
          for( let k=0; k < B; k++ ) {
            const S_ki   =  cα*S[k+i] + sα*S[k+j];
                  S[k+j] = -sα*S[k+i] + cα*S[k+j];
                  S[k+i] = S_ki;
          }
        }
        // ROTATE COLUMNS IN S
        for( let i = P*BB + p + n*n*B*B,
                 j = Q*BB + q + n*n*B*B;
                (i -= n*BB) >= 0; )
        {        j -= n*BB;
          for( let k=BB; (k -= B) >= 0; ) {
            const S_ik   = S[i+k]*cβ + S[j+k]*-sβ;
                  S[j+k] = S[i+k]*sβ + S[j+k]* cβ;
                  S[i+k] = S_ik;
          }
        }

        // ENTRIES (k,l) AND (l,k) ARE REMAINDERS (CANCELLATION ERROR) FROM ELIMINATION => SHOULD BE SAFELY ZEROABLE
        S[BB*(n*P+Q) + B*p+q] =
        S[BB*(n*Q+P) + B*q+p] = 0.0;

        const rowP = UV_off + N*(B*P + p),
              rowQ = UV_off + N*(B*Q + q);
        // ROTATE U & V
        _svd_jac_rot_rows(U, N, rowP,rowQ, cα, sα);
        _svd_jac_rot_rows(V, N, rowP,rowQ, cβ,-sβ);
      }
    }

    for( let i=0; i < N; i++ ) {
      const k = BB*(n+1)*(i/B | 0) + (B+1)*(i%B) | 0;
      sv[sv_off + i] = S[k];
    }

    _svd_jac_post_skip1( N, U,V, UV_off, sv, sv_off, ord );
  }

  return [
    new NDArray(shape,             U),
    new NDArray(shape.slice(0,-1),sv),
    new NDArray(shape,             V)
  ];
}
