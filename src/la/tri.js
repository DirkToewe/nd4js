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

import {ARRAY_TYPES} from '../dt'
import {asarray, NDArray} from '../nd_array'


export function tril(m, k=0)
{
  m = asarray(m);
  if( m.ndim < 2 ) throw new Error('Input must be at least 2D.');
  return m.mapElems( m.dtype, (x,...indices) => {
    const [i,j] = indices.slice(-2);
    return i < j-k ? 0 : x
  });
}


export function triu(m, k=0)
{
  m = asarray(m);
  if( m.ndim < 2 ) throw new Error('Input must be at least 2D.');
  return m.mapElems( m.dtype, (x,...indices) => {
    const [i,j] = indices.slice(-2);
    return i > j-k ? 0 : x
  });
}


export function tril_solve(L,y)
{
  L = asarray(L); if( L.ndim < 2 ) throw new Error('tril_solve(L,y): L.ndim must be at least 2.');
  y = asarray(y); if( y.ndim < 2 ) throw new Error('tril_solve(L,y): y.ndim must be at least 2.');

  const
    [N,M] = L.shape.slice(-2),
    [I,J] = y.shape.slice(-2);
  if( N != M ) throw new Error('Last two dimensions of L must be quadratic.')
  if( I != M ) throw new Error("L and y don't match.");

  const
    ndim = Math.max(L.ndim, y.ndim),
    shape = Int32Array.from({ length: ndim }, () => 1 );
  shape[ndim-2] = I;
  shape[ndim-1] = J;

  // FIND COMMON (BROADCASTED) SHAPE
  for( let arr of [L,y] )
    for( let i=ndim-2, j=arr.ndim-2; i-- > 0 && j-- > 0; )
      if( 1 === shape[i] )
        shape[i] = arr.shape[j];
      else if( shape[i] != arr.shape[j] && arr.shape[j] != 1 )
        throw new Error('Shapes are not broadcast-compatible.');

  // GENERATE RESULT DATA
  const
    dtype = [L,y].every(a => a.dtype==='float32') ? 'float32' : 'float64',
    x_dat = new ARRAY_TYPES[dtype]( shape.reduce((a,b) => a*b) ),
    L_dat = L.data,
    y_dat = y.data;
  let
    L_off = 0, L_stride = 1,
    y_off = 0, y_stride = 1,
    x_off = 0;

  function solv(d) {
    if( d === ndim-2 ) {
      L_stride = N*N;
      y_stride = N*J;

      // COPYING y
      for( let i=0; i < y_stride; i++ )
        x_dat[x_off+i] = y_dat[y_off+i]

      // FORWARD SUBSTITUTION
      for( let i=0; i < I; i++ )
      for( let j=0; j < J; j++ )
      {
        for( let k=0; k < i; k++ )
          x_dat[x_off+i*J+j] -= L_dat[L_off+N*i+k] * x_dat[x_off+k*J+j]
        x_dat[x_off+i*J+j] /= L_dat[L_off+N*i+i]
      }

      L_off += L_stride;
      y_off += y_stride;
      x_off += y_stride;

      return;
    }
    for( let l=shape[d]; ; l-- ) {
      solv(d+1);
      if( l == 1 ) break;
      if( ! (L.shape[ d - ndim + L.ndim ] > 1) ) L_off -= L_stride;
      if( ! (y.shape[ d - ndim + y.ndim ] > 1) ) y_off -= y_stride;
    }
    L_stride *= L.shape[ d - ndim + L.ndim ] || 1;
    y_stride *= y.shape[ d - ndim + y.ndim ] || 1;
  }
  solv(0);

  return new NDArray(shape,x_dat);
}


export function triu_solve (U,y)
{
  U = asarray(U); if( U.ndim < 2 ) throw new Error('triu_solve(U,y): U.ndim must be at least 2.');
  y = asarray(y); if( y.ndim < 2 ) throw new Error('triu_solve(U,y): y.ndim must be at least 2.');

  const
    [K,N] = U.shape.slice(-2),
    [I,J] = y.shape.slice(-2);
  if( K != N ) throw new Error('Last two dimensions of U must be quadratic.')
  if( I != N ) throw new Error("U and y don't match.");

  const
    ndim = Math.max(U.ndim, y.ndim),
    shape = Int32Array.from({length: ndim}, ()=>1);
  shape[ndim-2] = I;
  shape[ndim-1] = J;

  // FIND COMMON (BROADCASTED) SHAPE
  for( let arr of [U,y] )
    for( let i=ndim-2, j=arr.ndim-2; i-- > 0 && j-- > 0; )
      if( 1 === shape[i] )
        shape[i] = arr.shape[j];
      else if( shape[i] != arr.shape[j] && arr.shape[j] != 1 )
        throw new Error('Shapes are not broadcast-compatible.');

  // GENERATE RESULT DATA
  const
    dtype = [U,y].every(a => a.dtype==='float32') ? 'float32' : 'float64',
    x_dat = new ARRAY_TYPES[dtype]( shape.reduce((a,b) => a*b) ),
    U_dat = U.data,
    y_dat = y.data;
  let
    U_off = 0, U_stride = 1,
    y_off = 0, y_stride = 1,
    x_off = 0;

  function solv(d) {
    if( d === ndim-2 ) {
      U_stride = N*N;
      y_stride = N*J;

      // COPYING y
      for( let i=0; i < y_stride; i++ )
        x_dat[x_off+i] = y_dat[y_off+i]

      // BACKWARD SUBSTITUTION
      for( let i=I; i-- > 0; )
      for( let j=J; j-- > 0; )
      {
        for( let k=K; --k > i; )
          x_dat[x_off+i*J+j] -= U_dat[U_off+N*i+k] * x_dat[x_off+k*J+j]
        x_dat[x_off+i*J+j] /= U_dat[U_off+N*i+i]
      }

      U_off += U_stride;
      y_off += y_stride;
      x_off += y_stride;

      return;
    }
    for( let l=shape[d]; ; l-- ) {
      solv(d+1);
      if( l == 1 ) break;
      if( ! (U.shape[ d - ndim + U.ndim ] > 1) ) U_off -= U_stride;
      if( ! (y.shape[ d - ndim + y.ndim ] > 1) ) y_off -= y_stride;
    }
    U_stride *= U.shape[ d - ndim + U.ndim ] || 1;
    y_stride *= y.shape[ d - ndim + y.ndim ] || 1;
  }
  solv(0);

  return new NDArray(shape,x_dat);
}
