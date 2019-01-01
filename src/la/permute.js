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


export function permute_rows( A, P )
{
  A = asarray(A); if( A.ndim < 2 ) throw new Error('permute_rows(A,P): A.ndim must be at least 2.')
  P = asarray(P); if( P.ndim < 1 ) throw new Error('permute_rows(A,P): P.ndim must be at least 1.')
  if( P.dtype !== 'int32' )        throw new Error('permute_rows(A,P): P.dtype must be "int32".')

  const ndim = Math.max(A.ndim,P.ndim+1),
     B_shape = Int32Array.from({length: ndim}, () => 1),
       [M,N] = A.shape.slice(-2)

  if( M !== P.shape[P.ndim-1] )
    throw new Error('permute_rows(A,P): A.shape[-2] and P.shape[-1] must match.')

  for( let i=  ndim,
           j=A.ndim; i-- > 0 &&
                     j-- > 0; )
    B_shape[i] = A.shape[j];

  // FIND COMMON (BROADCASTED) SHAPE
  for( let i=  ndim-2,
           j=P.ndim-1; i-- > 0 &&
                       j-- > 0; )
    if( 1 === B_shape[i] )
      B_shape[i] = P.shape[j];
    else if( B_shape[i] != P.shape[j] && P.shape[j] != 1 )
      throw new Error('permute_rows(A,P): A and P not broadcast-compatible.');

  // GENERATE RESULT DATA
  const B = new ARRAY_TYPES[A.dtype]( B_shape.reduce((a,b) => a*b) ),
        A_shape = A.shape, A_ndim = A.ndim,
        P_shape = P.shape, P_ndim = P.ndim
  A = A.data
  P = P.data

  let A_off = 0, A_stride = -1,
      P_off = 0, P_stride = -1,
      B_off = 0

  function perm(d) {
    if( d === ndim-2 ) {
      // Q.T @ y
      for( let i=0; i < M; i++ )
      {
        const k = P[P_off + i]
        if( ! (k >= 0) ) throw new Error('permute_rows(A,P): P.contains invalid indices.')
        if( ! (k <  M) ) throw new Error('permute_rows(A,P): P.contains invalid indices.')

        for( let j=0; j < N; j++ )
          B[B_off + N*i+j] = A[A_off + N*k+j]
      }

      A_off += A_stride = M*N;
      P_off += P_stride = M;
      B_off += A_stride;

      return;
    }
    for( let l=B_shape[d]; ; l-- ) {
      perm(d+1);
      if( l === 1 ) break;
      if( ! (A_shape[ d - ndim   + A_ndim ] > 1) ) A_off -= A_stride;
      if( ! (P_shape[ d - ndim+1 + P_ndim ] > 1) ) P_off -= P_stride;
    }
    A_stride *= A_shape[ d - ndim   + A_ndim ] || 1;
    P_stride *= P_shape[ d - ndim+1 + P_ndim ] || 1;
  }
  perm(0);

  return new NDArray(B_shape,B);
}


export function permute_cols( A, P )
{
  A = asarray(A); if( A.ndim < 2 ) throw new Error('permute_cols(A,P): A.ndim must be at least 2.')
  P = asarray(P); if( P.ndim < 1 ) throw new Error('permute_cols(A,P): P.ndim must be at least 1.')
  if( P.dtype !== 'int32' )        throw new Error('permute_cols(A,P): P.dtype must be "int32".')

  const ndim = Math.max(A.ndim,P.ndim+1),
     B_shape = Int32Array.from({length: ndim}, () => 1),
       [M,N] = A.shape.slice(-2)

  if( N !== P.shape[P.ndim-1] )
    throw new Error('permute_cols(A,P): A.shape[-1] and P.shape[-1] must match.')

  for( let i=  ndim,
           j=A.ndim; i-- > 0 &&
                     j-- > 0; )
    B_shape[i] = A.shape[j];

  // FIND COMMON (BROADCASTED) SHAPE
  for( let i=  ndim-2,
           j=P.ndim-1; i-- > 0 &&
                       j-- > 0; )
    if( 1 === B_shape[i] )
      B_shape[i] = P.shape[j];
    else if( B_shape[i] != P.shape[j] && P.shape[j] != 1 )
      throw new Error('permute_cols(A,P): A and P not broadcast-compatible.');

  // GENERATE RESULT DATA
  const B = new ARRAY_TYPES[A.dtype]( B_shape.reduce((a,b) => a*b) ),
        A_shape = A.shape, A_ndim = A.ndim,
        P_shape = P.shape, P_ndim = P.ndim
  A = A.data
  P = P.data

  let A_off = 0, A_stride = -1,
      P_off = 0, P_stride = -1,
      B_off = 0

  function perm(d) {
    if( d === ndim-2 ) {
      // Q.T @ y
      for( let i=0; i < M; i++ )
      for( let j=0; j < N; j++ )
      {
        const k = P[P_off + j]
        if( ! (k >= 0) ) throw new Error('permute_cols(A,P): P.contains invalid indices.')
        if( ! (k <  N) ) throw new Error('permute_cols(A,P): P.contains invalid indices.')
        B[B_off + N*i+j] = A[A_off + N*i+k]
      }

      A_off += A_stride = M*N;
      P_off += P_stride =   N;
      B_off += A_stride;

      return;
    }
    for( let l=B_shape[d]; ; l-- ) {
      perm(d+1);
      if( l === 1 ) break;
      if( ! (A_shape[ d - ndim   + A_ndim ] > 1) ) A_off -= A_stride;
      if( ! (P_shape[ d - ndim+1 + P_ndim ] > 1) ) P_off -= P_stride;
    }
    A_stride *= A_shape[ d - ndim   + A_ndim ] || 1;
    P_stride *= P_shape[ d - ndim+1 + P_ndim ] || 1;
  }
  perm(0);

  return new NDArray(B_shape,B);
}


export function unpermute_rows( A, P )
{
  A = asarray(A); if( A.ndim < 2 ) throw new Error('unpermute_rows(A,P): A.ndim must be at least 2.')
  P = asarray(P); if( P.ndim < 1 ) throw new Error('unpermute_rows(A,P): P.ndim must be at least 1.')
  if( P.dtype !== 'int32' )        throw new Error('unpermute_rows(A,P): P.dtype must be "int32".')

  const ndim = Math.max(A.ndim,P.ndim+1),
     B_shape = Int32Array.from({length: ndim}, () => 1),
       [M,N] = A.shape.slice(-2)

  if( M !== P.shape[P.ndim-1] )
    throw new Error('unpermute_rows(A,P): A.shape[-2] and P.shape[-1] must match.')

  for( let i=  ndim,
           j=A.ndim; i-- > 0 &&
                     j-- > 0; )
    B_shape[i] = A.shape[j];

  // FIND COMMON (BROADCASTED) SHAPE
  for( let i=  ndim-2,
           j=P.ndim-1; i-- > 0 &&
                       j-- > 0; )
    if( 1 === B_shape[i] )
      B_shape[i] = P.shape[j];
    else if( B_shape[i] != P.shape[j] && P.shape[j] != 1 )
      throw new Error('unpermute_rows(A,P): A and P not broadcast-compatible.');

  // GENERATE RESULT DATA
  const B = new ARRAY_TYPES[A.dtype]( B_shape.reduce((a,b) => a*b) ),
        A_shape = A.shape, A_ndim = A.ndim,
        P_shape = P.shape, P_ndim = P.ndim
  A = A.data
  P = P.data

  let A_off = 0, A_stride = -1,
      P_off = 0, P_stride = -1,
      B_off = 0

  function perm(d) {
    if( d === ndim-2 ) {
      // Q.T @ y
      for( let i=0; i < M; i++ )
      {
        const k = P[P_off + i]
        if( ! (k >= 0) ) throw new Error('unpermute_rows(A,P): P.contains invalid indices.')
        if( ! (k <  M) ) throw new Error('unpermute_rows(A,P): P.contains invalid indices.')

        for( let j=0; j < N; j++ )
          B[B_off + N*k+j] = A[A_off + N*i+j]
      }

      A_off += A_stride = M*N;
      P_off += P_stride = M;
      B_off += A_stride;

      return;
    }
    for( let l=B_shape[d]; ; l-- ) {
      perm(d+1);
      if( l === 1 ) break;
      if( ! (A_shape[ d - ndim   + A_ndim ] > 1) ) A_off -= A_stride;
      if( ! (P_shape[ d - ndim+1 + P_ndim ] > 1) ) P_off -= P_stride;
    }
    A_stride *= A_shape[ d - ndim   + A_ndim ] || 1;
    P_stride *= P_shape[ d - ndim+1 + P_ndim ] || 1;
  }
  perm(0);

  return new NDArray(B_shape,B);
}


export function unpermute_cols( A, P )
{
  A = asarray(A); if( A.ndim < 2 ) throw new Error('unpermute_cols(A,P): A.ndim must be at least 2.')
  P = asarray(P); if( P.ndim < 1 ) throw new Error('unpermute_cols(A,P): P.ndim must be at least 1.')
  if( P.dtype !== 'int32' )        throw new Error('unpermute_cols(A,P): P.dtype must be "int32".')

  const ndim = Math.max(A.ndim,P.ndim+1),
     B_shape = Int32Array.from({length: ndim}, () => 1),
       [M,N] = A.shape.slice(-2)

  if( N !== P.shape[P.ndim-1] )
    throw new Error('unpermute_cols(A,P): A.shape[-1] and P.shape[-1] must match.')

  for( let i=  ndim,
           j=A.ndim; i-- > 0 &&
                     j-- > 0; )
    B_shape[i] = A.shape[j];

  // FIND COMMON (BROADCASTED) SHAPE
  for( let i=  ndim-2,
           j=P.ndim-1; i-- > 0 &&
                       j-- > 0; )
    if( 1 === B_shape[i] )
      B_shape[i] = P.shape[j];
    else if( B_shape[i] != P.shape[j] && P.shape[j] != 1 )
      throw new Error('unpermute_cols(A,P): A and P not broadcast-compatible.');

  // GENERATE RESULT DATA
  const B = new ARRAY_TYPES[A.dtype]( B_shape.reduce((a,b) => a*b) ),
        A_shape = A.shape, A_ndim = A.ndim,
        P_shape = P.shape, P_ndim = P.ndim
  A = A.data
  P = P.data

  let A_off = 0, A_stride = -1,
      P_off = 0, P_stride = -1,
      B_off = 0

  function perm(d) {
    if( d === ndim-2 ) {
      // Q.T @ y
      for( let i=0; i < M; i++ )
      for( let j=0; j < N; j++ )
      {
        const k = P[P_off + j]
        if( ! (k >= 0) ) throw new Error('unpermute_cols(A,P): P.contains invalid indices.')
        if( ! (k <  N) ) throw new Error('unpermute_cols(A,P): P.contains invalid indices.')
        B[B_off + N*i+k] = A[A_off + N*i+j]
      }

      A_off += A_stride = M*N;
      P_off += P_stride =   N;
      B_off += A_stride;

      return;
    }
    for( let l=B_shape[d]; ; l-- ) {
      perm(d+1);
      if( l === 1 ) break;
      if( ! (A_shape[ d - ndim   + A_ndim ] > 1) ) A_off -= A_stride;
      if( ! (P_shape[ d - ndim+1 + P_ndim ] > 1) ) P_off -= P_stride;
    }
    A_stride *= A_shape[ d - ndim   + A_ndim ] || 1;
    P_stride *= P_shape[ d - ndim+1 + P_ndim ] || 1;
  }
  perm(0);

  return new NDArray(B_shape,B);
}
