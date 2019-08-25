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

import {ARRAY_TYPES, eps} from '../dt'
import {asarray, NDArray} from '../nd_array'
import {matmul2} from './matmul'
import {qr_decomp} from './qr'
import {transpose_inplace} from './transpose_inplace'
import {_svd_jac_rot_rows,
        _svd_jac_rot_cols,
        _svd_jac_angles,
        _svd_jac_post } from './_svd_jac_utils'


export function svd_jac_2sided(A)
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
         [U,sv,V] = svd_jac_2sided(R)
      return [matmul2(Q,U), sv, V]
    }
    if( N < M ) {
      const [Q,R] = qr_decomp(A.T),
         [U,sv,V] = svd_jac_2sided(R)
      transpose_inplace(V)
      return [V, sv, matmul2(Q,U).T]
    }
  }
  // ALLOCATE RESULT DATA
  const DType = A.dtype==='float32' ? 'float32' : 'float64',
        DTypeArray = ARRAY_TYPES[DType],
        TOL = (N * eps(DType))**2,
        B = DType === 'float32' ? 64/4 : 64/8, // <- block size should match cache line size
        U = DTypeArray.from(A.data); A = undefined; // <- potentially allow GC
  const S = new DTypeArray(N*N), // <- tempory storage for decomposition
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
    for( let j=0; j < N; j++ ) { S[N*i+j] = U[UV_off + N*i+j];
                                            U[UV_off + N*i+j] = (i===j)*1 };
    // INIT V TO IDENTITY
    for( let i=0; i < N; i++ )
    for( let j=0; j < N; j++ ) V[UV_off + N*i+j] = (i===j)*1;

     //
    // JACOBI SVD ITERATIONS
   //
    for( let finished = false;
           ! finished; )
    {        finished = true;
      // SWEEP
      for( let Q=0; Q <  N; Q += B )
      for( let P=0; P <= Q; P += B )
      for( let q=Q; q < Q+B && q < N; q++ )
      for( let p=P; p < P+B && p < q; p++ )
      {
        const S_pp = S[N*p+p],
              S_pq = S[N*p+q],
              S_qp = S[N*q+p],
              S_qq = S[N*q+q];
        // stopping criterion inspiredy by:
        //  "Jacobi's Method is More Accurate than QR"
        //   by James Demmel
        //   SIAM J. Matrix Anal. Appl, vol. 13, pp. 1204-1245, 1992
//        if( ! ( Math.max( Math.abs(S_pq),
//                          Math.abs(S_qp) ) > Math.sqrt(Math.abs(S_pp*S_qq)) * TOL ) )
        if( ! ( S_pq*S_pq + S_qp*S_qp > Math.abs(S_pp*S_qq) * TOL ) )
          continue;

        finished = false;
 
        const [cα,sα, cβ,sβ] = _svd_jac_angles(S_pp,S_pq,
                                               S_qp,S_qq);

        // ROTATE S
        _svd_jac_rot_rows(S, N, N*p,N*q, cα,sα);
        _svd_jac_rot_cols(S, N,   p,  q, cβ,sβ);

        // ENTRIES (k,l) AND (l,k) ARE REMAINDERS (CANCELLATION ERROR) FROM ELIMINATION => SHOULD BE SAFELY ZEROABLE
        S[N*p+q] = 0.0;
        S[N*q+p] = 0.0;
  
        // ROTATE U & V
        _svd_jac_rot_rows(U, N, UV_off + N*p,
                                UV_off + N*q, cα, sα);
        _svd_jac_rot_rows(V, N, UV_off + N*p,
                                UV_off + N*q, cβ,-sβ);
      }
    }

    _svd_jac_post( N, U,S,V, UV_off, sv, sv_off, ord );
  }

  return [
    new NDArray(shape,             U),
    new NDArray(shape.slice(0,-1),sv),
    new NDArray(shape,             V)
  ];
}
