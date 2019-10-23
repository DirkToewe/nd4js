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
import {qr_decomp} from './qr'


export const det_tri = a =>
{
  a = asarray(a);
  if( a.ndim < 2 )
    throw new Error(`det_tri(a): a.shape=[${a.shape}]; a.ndim must be at least 2.`);

  const a_shape = a.shape,
          [M,N] = a_shape.slice(-2);
  if( M !== N )
    throw new Error(`det_tri(a): a must be square matrices.`);

  const DTypeArray = ARRAY_TYPES[a.dtype];

  const A = a.data; a = undefined;
  const D = new DTypeArray( A.length / (N*N) );

  for( let  D_off=0,
            A_off=0; D_off < D.length;
                     D_off++ )
  { let i = A_off;
            A_off += N*N;
    let d = 1;
    for( ; i < A_off; i += N+1 )
      d *= A[i];
    D[D_off] = d;
  }

  return new NDArray(a_shape.slice(0,-2), D);
}


export const slogdet_tri = A =>
{
  A = asarray(A);
  if( A.ndim < 2 )
    throw new Error(`det_tri(A): A.ndim must be at least 2.`);

  const A_shape = A.shape,
          [M,N] = A_shape.slice(-2);
  if( M !== N )
    throw new Error(`det_tri(A): A must be square matrices.`);

  const DTypeArray = ARRAY_TYPES[A.dtype==='float32' ? 'float32' : 'float64'];

        A = A.data;
  const S = new DTypeArray( A.length / (N*N) ),
        D = new DTypeArray( S.length );

  for( let SD_off=0,
            A_off=0; SD_off < D.length;
                     SD_off++ )
  { let i = A_off;
            A_off += N*N;
    let s = 1,
        d = 0;
    for( ; i < A_off; i += N+1 ) {
      s *=           Math.sign(A[i]);
      d += Math.log( Math. abs(A[i]) );
    }
    S[SD_off] = s;
    D[SD_off] = d;
  }

  const SD_shape = A_shape.slice(0,-2);
  return [
    new NDArray(SD_shape, S),
    new NDArray(SD_shape, D)
  ];
}


export const det = A =>
{
  const [Q,R] = qr_decomp(A);
  return det_tri(R);
}


export const slogdet = A =>
{
  const [Q,R] = qr_decomp(A);
  return slogdet_tri(R);
}
