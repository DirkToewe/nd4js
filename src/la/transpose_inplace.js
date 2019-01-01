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

import {asarray} from '../nd_array'


export function transpose_inplace(A)
{
  const [N,M] = A.shape.slice(-2);
  if( N != M ) throw new Error('In-place transposition is only supported for square matrices.');
  A = A.data;
  for( let off=0; off < A.length; off += N*N )
  for( let i=0  ; i < N-1; i++ )
  for( let j=1+i; j < N  ; j++ ) {
    const
      ij = off + N*i+j,
      ji = off + N*j+i,
      A_ij = A[ij]; A[ij] = A[ji]; A[ji] = A_ij;
  }
}
