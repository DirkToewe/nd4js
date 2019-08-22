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


function rot_rows( A, N, i, j, c, s )
{
  i |= 0;
  j |= 0;
  c = +c;
  s = +s;
  for( let k=N | 0; k-- > 0; ) {
    const A_i = A[i],
          A_j = A[j];
    A[i] = c*A_i + s*A_j;
    A[j] = c*A_j - s*A_i;
    i = i+1 | 0;
    j = j+1 | 0;
  }
}


function rot_cols( A, N, i, j, c, s )
{
  i |= 0;
  j |= 0;
  c = +c;
  s = +s;
  for( let k=N | 0; k-- > 0; ) {
    const A_i = A[i],
          A_j = A[j];
    A[i] = c*A_i - s*A_j;
    A[j] = c*A_j + s*A_i;
    i = i+N | 0;
    j = j+N | 0;
  }
}


//function angles( S_pp, S_pq,
//                 S_qp, S_qq )
//{
//  // determine rotation angles such that
//  // ┌               ┐ ┌            ┐ ┌               ┐   ┌        ┐
//  // │ cos(α) sin(α) │ │ S_pp  S_pq │ │ cos(β) sin(β) │ ! │ s1   0 │
//  // │               │ │            │ │               │ = │        │
//  // │-sin(α) cos(α) │ │ S_qp  S_qq │ │-sin(β) cos(β) │   │  0  s2 │
//  // └               ┘ └            ┘ └               ┘   └        ┘ 
//  const d = S_qp - S_pq;
//  let  ca = 1,
//       sa = 0,
//        x = (S_qq - S_pp) / (2*S_pq);
//  if( 0 !== d )
//  { // SYMMETRIZE
//    x = (S_pp + S_qq) / d;
//    const    y = Math.sqrt(1 + x*x);
//    sa = 1 / y;
//    ca = x / y;
//    x = ( x*(S_qq - S_pp) - (S_pq + S_qp) )
//      / ( x*(S_pq + S_qp) + (S_qq - S_pp) );
//  }
//
//  // DIAGONALIZE
//  let   y = Math.abs(x) + Math.sqrt(1 + x*x);
//  const s = x < 0 ? -1 : +1,
//        z = 1 / y;
//  let  cb = 1 / Math.sqrt(1 + z*z),
//       sb = s / Math.sqrt(1 + y*y);
//
//   x = ca*cb + sa*sb;
//  sa = sa*cb - ca*sb;
//  ca = x;
//
//  x = cb*(sa*S_qp + ca*S_pp) - sb*(sa*S_qq + ca*S_pq),
//  y = sb*(ca*S_qp - sa*S_pp) + cb*(ca*S_qq - sa*S_pq);
//
//  if( Math.abs(x) < Math.abs(y) ) {
//    [ca,sa] = [-sa,ca];
//    [sb,cb] = [-cb,sb]; x = y;
//  }
//
//  if( x < 0 ) {
//    cb = -cb;
//    sb = -sb;
//  }
//
//  return [ca,sa, cb,sb];
//}


function angles( S_pp, S_pq,
                 S_qp, S_qq )
{
  // determine rotation angles such that
  // ┌               ┐ ┌            ┐ ┌               ┐   ┌        ┐
  // │ cos(α) sin(α) │ │ S_pp  S_pq │ │ cos(β) sin(β) │ ! │ s1   0 │
  // │               │ │            │ │               │ = │        │
  // │-sin(α) cos(α) │ │ S_qp  S_qq │ │-sin(β) cos(β) │   │  0  s2 │
  // └               ┘ └            ┘ └               ┘   └        ┘ 
  let x = Math.atan2( S_qp - S_pq, S_qq + S_pp ),
      y = Math.atan2( S_qp + S_pq, S_qq - S_pp );

  const a = (x-y)/2,
        b = (x+y)/2;

  let ca = Math.cos(a), sa = Math.sin(a),
      cb = Math.cos(b), sb = Math.sin(b);

  x = cb*(sa*S_qp + ca*S_pp) - sb*(sa*S_qq + ca*S_pq);
  y = sb*(ca*S_qp - sa*S_pp) + cb*(ca*S_qq - sa*S_pq);

  if( Math.abs(x) < Math.abs(y) ) {
    [sa,ca] = [ca,-sa];
    [cb,sb] = [sb,-cb]; x = y;
  }

  if( x < 0 ) {
    cb = -cb;
    sb = -sb;
  }

  return [ca,sa, cb,sb];
}


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
 
        const [cα,sα, cβ,sβ] = angles(S_pp,S_pq,
                                      S_qp,S_qq);

        // ROTATE S
        rot_rows(S, N, N*p,N*q, cα,sα);
        rot_cols(S, N,   p,  q, cβ,sβ);

        // ENTRIES (k,l) AND (l,k) ARE REMAINDERS (CANCELLATION ERROR) FROM ELIMINATION => SHOULD BE SAFELY ZEROABLE
        S[N*p+q] = 0.0;
        S[N*q+p] = 0.0;
  
        // ROTATE U & V
        rot_rows(U, N, UV_off + N*p,
                       UV_off + N*q, cα, sα);
        rot_rows(V, N, UV_off + N*p,
                       UV_off + N*q, cβ,-sβ);
      }
    }

    // MOVE S TO sv (AND MAKE POSITIVE)
    for( let i=0; i < N; i++ ) {
      const sv_i = S[N*i+i];
      // flip sign if necessary
      if(   sv_i < 0 || Object.is(sv_i,-0) )
        for( let j=0; j < N; j++ )
          U[UV_off + N*i+j] *= -1;
      sv[sv_off + i] = Math.abs(sv_i);
    }

    // SORT SINGULAR VALUES
    ord.sort( (i,j) => sv[sv_off+j] - sv[sv_off+i] );

    for( let i=0; i < N;  ++i  )
    for( let j=i;; ) // <- start swap cycle
    {
      let tmp = ord[j];
                ord[j] = j;
                         j = tmp;
      if( j <= i )
        break;
      // SWAP ROWS IN U AND V
      const rowI = UV_off + ord[j]*N,
            rowJ = UV_off +     j *N;
      for( let k=0; k < N; k++ ) { const U_ik = U[rowI+k]; U[rowI+k] = U[rowJ+k]; U[rowJ+k] = U_ik; }
      for( let k=0; k < N; k++ ) { const V_ik = V[rowI+k]; V[rowI+k] = V[rowJ+k]; V[rowJ+k] = V_ik; }
      // SWAP SV
      tmp = sv[sv_off + ord[j]];
            sv[sv_off + ord[j]] = sv[sv_off + j];
                                  sv[sv_off + j] = tmp;
    }

    // TRANSPOSE U IN-PLACE (TRANSPOSED FOR CACHE LOCALITY REASONS)
    for( let i=0;   i < N-1; i++ )
    for( let j=i; ++j < N  ; ) {
      const ij = UV_off + N*i+j,
            ji = UV_off + N*j+i, U_ij = U[ij];
                                        U[ij] = U[ji];
                                                U[ji] = U_ij;
    }
  }

  return [
    new NDArray(shape,             U),
    new NDArray(shape.slice(0,-1),sv),
    new NDArray(shape,             V)
  ]
}
