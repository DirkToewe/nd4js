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
import {KahanSum} from '../kahan_sum'
import {asarray, NDArray} from '../nd_array'

import {_tril_solve,
        _tril_t_solve} from './tri'


export function _cholesky_decomp(M,N, L,L_off)
{
  if( ! (M <= N) ) throw new Error('Assertion failed.');

  const kahan = new KahanSum();

  // https://de.wikipedia.org/wiki/Cholesky-Zerlegung
  for( let i=0; i<M; i++ )
  for( let j=0; j<M; j++ )
    if( i < j ) {
      L[L_off+N*i+j] = 0;
    } else {
      kahan.set( L[L_off+N*i+j] );
      for( let k=0; k<j; k++ )
        kahan.add( - L[L_off+N*i+k] * L[L_off+N*j+k] );

      if( i > j ) L[L_off+N*i+j] =           kahan.sum / L[L_off+N*j+j];
      else {      L[L_off+N*i+i] = Math.sqrt(kahan.sum);
        if( isNaN(L[L_off+N*i+i]) )
          throw new Error('Matrix contains NaNs or is (near) singular.');
      }
    }
}


export function cholesky_decomp(S)
{
  S = asarray(S);
  const
    dtype = S.dtype === 'float32' ? 'float32' : 'float64',
    shape = S.shape,
    [N,M] =   shape.slice(-2),
    L = ARRAY_TYPES[dtype].from(S.data);
  S = undefined;

  if( N != M )
    throw new Error('Last two dimensions must be quadratic.')

  for( let L_off=0; L_off < L.length; L_off += N*N )
    _cholesky_decomp(N,N, L,L_off);

  return new NDArray(shape,L);
}


export function cholesky_solve(L,y)
{
  L = asarray(L);
  y = asarray(y);
  if( L.ndim < 2 ) throw new Error('L must be at least 2D.');
  if( y.ndim < 2 ) throw new Error('y must be at least 2D.');

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
    dtype = L.dtype === 'float32'
         && y.dtype === 'float32' ? 'float32' : 'float64',
    x_dat = new ARRAY_TYPES[dtype]( shape.reduce((a,b) => a*b, 1) ),
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

      _tril_solve  (I,I,J, L_dat,L_off, x_dat,x_off);
      _tril_t_solve(I,I,J, L_dat,L_off, x_dat,x_off);

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
