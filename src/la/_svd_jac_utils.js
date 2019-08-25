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
// THIS FILE IS FOR INTERNAL USE ONLY

import {asarray, NDArray} from '../nd_array'
import {ARRAY_TYPES, eps} from '../dt'
import {matmul2} from './matmul'
import {qr_decomp} from './qr'
import {transpose_inplace} from './transpose_inplace'


/** Applies givens roation to rows i and j.
 */
export function _svd_jac_rot_rows( W, N, i, j, c, s )
{
  i |= 0;
  j |= 0;
  c = +c;
  s = +s;
  for( let k=N | 0; k-- > 0; ) {
    const W_i = W[i],
          W_j = W[j];
    W[i] = c*W_i + s*W_j;
    W[j] = c*W_j - s*W_i;
    i = i+1 | 0;
    j = j+1 | 0;
  }
}


/** Applies givens rotation to columns i and j.
 */
export function _svd_jac_rot_cols( W, N, i, j, c, s )
{
  i |= 0;
  j |= 0;
  c = +c;
  s = +s;
  for( let k=N | 0; k-- > 0; ) {
    const W_i = W[i],
          W_j = W[j];
    W[i] = c*W_i - s*W_j;
    W[j] = c*W_j + s*W_i;
    i = i+N | 0;
    j = j+N | 0;
  }
}


//function _svd_jac_angles( S_pp, S_pq,
//                          S_qp, S_qq )
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


export function _svd_jac_angles( S_pp, S_pq,
                                 S_qp, S_qq )
{
  // determine rotation angles such that
  // ┌               ┐ ┌            ┐ ┌               ┐   ┌        ┐
  // │ cos(α) sin(α) │ │ S_pp  S_pq │ │ cos(β) sin(β) │ ! │ s1   0 │
  // │               │ │            │ │               │ = │        │
  // │-sin(α) cos(α) │ │ S_qp  S_qq │ │-sin(β) cos(β) │   │  0  s2 │
  // └               ┘ └            ┘ └               ┘   └        ┘
  //
  // such that: s1 ≥ |s2|
  //
  // FIXME: update the following documentation
  //
  // => 0 = cos(α)*{D_kl*cos(β) - D_ll*sin(β)} + sin(α)*{D_kk*cos(β) - D_lk*sin(β)} = ½⋅{(D_ll - D_kk)⋅sin(α-β) - (D_kk + D_ll)⋅sin(α+β) + (D_kl + D_lk)⋅cos(α-β) + (D_kl - D_lk)⋅cos(α+β)}
  //    0 = cos(α)*{D_kk*sin(β) + D_lk*cos(β)} - sin(α)*{D_kl*sin(β) + D_ll*cos(β)} = ½⋅{(D_ll - D_kk)⋅sin(α-β) + (D_kk + D_ll)⋅sin(α+β) + (D_kl + D_lk)⋅cos(α-β) - (D_kl - D_lk)⋅cos(α+β)}
  //   d1 = cos(α)*{D_kk*cos(β) - D_lk*sin(β)} - sin(α)*{D_kl*cos(β) - D_ll*sin(β)}
  //   d2 = cos(α)*{D_kl*sin(β) + D_ll*cos(β)} + sin(α)*{D_kk*sin(β) + D_lk*cos(β)}
  //
  // => 0 = {D_kl⋅cos(α) - D_kk⋅sin(α)}⋅cos(β) + {D_lk⋅sin(α) - D_ll⋅cos(α)}⋅sin(β)
  //    0 = {D_kl⋅sin(α) + D_kk⋅cos(α)}⋅sin(β) + {D_lk⋅cos(α) + D_ll⋅sin(α)}⋅cos(β)
  //   d1 = {D_kl⋅sin(α) + D_kk⋅cos(α)}⋅cos(β) + {D_lk⋅cos(α) - D_ll⋅sin(α)}⋅sin(β)
  //   d2 = {D_kl⋅cos(α) - D_kk⋅sin(α)}⋅sin(β) + {D_lk⋅sin(α) + D_ll⋅cos(α)}⋅cos(β)
  //
  // => 0 = (D_kk+D_ll)⋅sin(α-β) + (D_kl-D_lk)⋅cos(α-β)  
  //    0 = (D_kk-D_ll)⋅sin(α+β) + (D_kl+D_lk)⋅cos(α+β) 
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


/** Postprocessing of the jacobi iteration result:
 *    1) Copy singular values from work matrix S to result array sv
 *    2) Make singular values positive
 *    3) Sort singular values
 *    4) Transpose U (was transposed for cache friendliness reasons)
 */
export function _svd_jac_post( N, U,S,V, UV_off, sv, sv_off, ord )
{
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

  // nested loop but actually O(N) operation
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

  // TRANSPOSE U
  for( let i=0;   i < N-1; i++ )
  for( let j=i; ++j < N  ;     ) {
    const ij = UV_off + N*i+j,
          ji = UV_off + N*j+i, U_ij = U[ij];
                                      U[ij] = U[ji];
                                              U[ji] = U_ij;
  }
}
