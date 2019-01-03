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
import {ARRAY_TYPES} from '../dt'


export function hessenberg_decomp(A)
{
  A = asarray(A);
  if( A.ndim < 2 ) throw new Error('A must be at least 2D.');
  const
    DTypeArray = ARRAY_TYPES[A.dtype==='float32' ? 'float32' : 'float64'], // <- ensure at least double precision
    shape = A.shape,
    N = shape[shape.length-1],
    w   = new DTypeArray(N-1),
    tmp = new DTypeArray(N);
  if( N != shape[shape.length-2] ) throw new Error('A is not square.');

  const H = DTypeArray.from(A.data); A = undefined;
  const U = new DTypeArray(H.length);

  for( let off=0; off < H.length; off += N*N )
  { // INIT P TO IDENTITY
    for( let i=0; i < N; i++ ) U[off + N*i+i] = 1.0;

    // ELIMINATION BEGIN
    for( let i=N; --i > 1; )
    { 
      let row_norm; { // BEGIN CREATE HOUSEHOLDER VECTOR
        let
          norm_sum = 0.0,
          norm_max = 0.0;
        for( let j=i; j-- > 0; ) {
          let H_ij = H[off + N*i+j];
          w[j] = H_ij;
          H_ij = Math.abs(H_ij);
          if(   H_ij > 0 ) {
            if( H_ij > norm_max ) {
              const scale = norm_max / H_ij; norm_max = H_ij;
              norm_sum *= scale*scale;
            }
            const ratio = H_ij / norm_max;
            norm_sum += ratio*ratio;
          }
        }
        row_norm = isFinite(norm_max) ? Math.sqrt(norm_sum)*norm_max : norm_max;
        if( w[i-1] > 0 )
          row_norm *= -1 // <- avoid losing signifciance (other people do, so ... I should do it as well) TODO verify usefulness
        w[i-1] -= row_norm
        norm_sum = 0.0;
        norm_max = 0.0;
        for( let j=i; j-- > 0; ) {
          const w_j = Math.abs(w[j]);
          if(   w_j > 0 ) {
            if( w_j > norm_max ) {
              const scale = norm_max / w_j; norm_max = w_j;
              norm_sum *= scale*scale;
            }
            const ratio = w_j / norm_max;
            norm_sum += ratio*ratio;
          }
        }
        norm_sum = isFinite(norm_max) ? Math.sqrt(norm_sum)*norm_max : norm_max;
        if( norm_sum == 0 ) continue;
        for( let j=i; j-- > 0; ) w[j] /= norm_sum;
      } // END CREATE HOUSEHOLDER VECTOR

      // APPLY w TO RIGHT OF H
      for( let j=N; j-- > 0; )
        if( j != i )
        { let sum=0.0
          for( let k=i; k-- > 0; ) sum += w[k] * H[off + N*j+k];
          sum *= 2;
          for( let k=i; k-- > 0; ) H[off + N*j+k] -= w[k] * sum;
        }
      H.fill( 0.0, off + N*i+0, off + N*i+i-1 )
      H[off + N*i+i-1] = row_norm;

      // APPLY w TO LEFT OF H
      tmp.fill(0.0);
      for( let j=i; j-- > 0; )
      for( let k=N; k-- > 0; ) tmp[k] += w[j] * H[off + N*j+k];
      for( let j=i; j-- > 0; )
      for( let k=N; k-- > 0; ) H[off + N*j+k] -= 2 * w[j] * tmp[k];

      // APPLY w TO RIGHT OF U
      for( let j=N; j-- > 0; )
      { let sum=0.0
        for( let k=i; k-- > 0; ) sum += w[k] * U[off + N*j+k];
        sum *= 2;
        for( let k=i; k-- > 0; ) U[off + N*j+k] -= w[k] * sum;
      }
    }
    // ELIMINATION END
  }

  return [
    new NDArray(shape,U),
    new NDArray(shape,H)
  ];
};


//la.hessenberg_decomp = A => {
//  A = nd.asarray(A);
//  if( A.ndim < 2 ) throw new Error('A must be at least 2D.');
//  const
//    DTypeArray = nd.dtypes[ nd.super_dtype(A.dtype,'float64') ], // <- ensure at least double precision
//    shape = A.shape,
//    N = shape[shape.length-1];
//  if( N != shape[shape.length-2] ) throw new Error('A is not square.');
//
//  const H = DTypeArray.from(A.data); A = undefined;
//  const U = new DTypeArray(H.length);
//
//  for( let off=0; off < H.length; off += N*N )
//  { // INIT P TO IDENTITY
//    for( let i=0; i < N; i++ ) U[off + N*i+i] = 1.0;
//
//    // { ELIMINATION BEGIN
//    for( let i=N-1; --i > 0; )
//    for( let j=i  ; j-- > 0; )
//    { // ROTATE COLUMS/ROWS j and i
//      // determine angle
//      const H_Ij = H[off + N*(i+1)+j]; if( H_Ij == 0 ) continue;
//      const H_Ii = H[off + N*(i+1)+i],
//            norm = Math.hypot(H_Ii,H_Ij),
//            c = H_Ii / norm,
//            s = H_Ij / norm;
//
//      // rotate columns in H
//      for( let k=i+2; k-- > 0; )
//      {
//        const kj = off + N*k+j, H_kj = H[kj],
//              ki = off + N*k+i, H_ki = H[ki];
//        H[kj] = c*H_kj - s*H_ki;
//        H[ki] = s*H_kj + c*H_ki;
//      }
//
//      // rotate rows in H
//      for( let k=N; k-- > 0; )
//      {
//        const jk = off + N*j+k, H_jk = H[jk],
//              ik = off + N*i+k, H_ik = H[ik];
//        H[jk] = c*H_jk - s*H_ik;
//        H[ik] = s*H_jk + c*H_ik;
//      }
//
//      // rotate rows in U (transposed for cache locality reasons)
//      for( let k=N; k-- > 0; ) // <- FIXME are there operations to be saved here?
//      {
//        const jk = off + N*j+k, U_jk = U[jk],
//              ik = off + N*i+k, U_ik = U[ik];
//        U[jk] = c*U_jk - s*U_ik;
//        U[ik] = s*U_jk + c*U_ik;
//      }
//    }
//    // } ELIMINATION END
//
//    // TRANSPOSE U BACK (transposed for cache locality reasons)
//    for( let i=0; i < N; i++ )
//    for( let j=0; j < i; j++ ) {
//      const ij = off + N*i+j,
//            ji = off + N*j+i, U_ij = U[ij]; U[ij] = U[ji]; U[ji] = U_ij;
//    }
//  }
//
//  return [
//    new nd.Array(shape,U),
//    new nd.Array(shape,H)
//  ]
//};
