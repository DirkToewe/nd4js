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


export function diag_mat(diag)
{
  if( diag.ndim < 1 )
    throw new Error(`diag_mat(diag): diag.ndim=${diag.ndim} is less than 1.`)
  diag = asarray(diag)
  const
    shape = diag.shape,
    N = shape[shape.length-1],
    d = diag.data,
    D = ARRAY_TYPES[diag.dtype].from({ length: shape.reduce((a,b) => a*b)*N }, ()=>0)
  diag = undefined

  if( shape.length < 1 ) throw new Error('diag_mat(diag): diag must be at least 1d.')
  if( N <= 0 ) throw new Error('Assertion Failed!');

  for( let d_off=0,
           D_off=0;; )
  {
    const d_end = d_off+N;
    while(true) {
      D[D_off] = d[d_off];
      if( ++d_off === d_end ) break;
      D_off += N+1;
    }
    if( ++D_off == D.length ) break;
  }
  return new NDArray(Int32Array.of(...shape,N), D)
}


export function diag(A, offset=0)
{
  if( A.ndim < 2 )
    throw new Error(`diag(A, offset): A.ndim=${A.ndim} is less than 2.`)
    
  A = asarray(A)
  const  [N,M] = A.shape.slice(-2),
         shape = A.shape.slice(0,-1),
    DTypeArray = ARRAY_TYPES[A.dtype]
  A = A.data

  if( ! (offset > -N
      && offset < +M) )
    throw new Error(`diag(A, offset): offset=${offset} out of A.shape=[${A.shape}].`)

  const L = Math.min(
    N, N+offset,
    M, M-offset
  )
  shape[shape.length-1] = L;

  const
    d = new DTypeArray( shape.reduce((a,b) => a*b) ),
    a_off = Math.max( 0, offset, -offset*M );

  for( let d_off=0,
           A_off=0; A_off < A.length;
           d_off += L,
           A_off += N*M
  )
  {
    for( let i=0; i < L; i++ )
      d[d_off+i] = A[A_off + a_off + M*i+i]
  }
  return new NDArray(shape, d)
}
