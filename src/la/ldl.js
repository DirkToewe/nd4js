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


// THE FOLLOWING IMPLEMENTATION HAS BETTER CACHE ALIGNMENT BUT REQUIRES AN EXTRA ARRAY OF LENGTH M-1
//export function _ldl_decomp(M,N, LD,LD_off, V)
//{
//  if( ! (M   <= N       ) ) throw new Error('Assertion failed.');
//  if( ! (M-1 <= V.length) ) throw new Error('Assertion failed.');
//
//  // https://arxiv.org/abs/1111.4144
//  for( let j=0; j < M; j++ )
//  {
//    for( let k=0; k < j; k++ )
//      V[k] = LD[LD_off + N*j+k] * LD[LD_off + N*k+k];
//
//    for( let i=j; i < M; i++ )
//    {
//      for( let k=0; k < j; k++ )
//        LD[LD_off + N*i+j] -= LD[LD_off + N*i+k] * V[k];
//
//      if( j < i )
//        LD[LD_off + N*i+j] /= LD[LD_off + N*j+j];
//    }
//  }
//}


export function _ldl_decomp(M,N, LD,LD_off)
{
  if( ! (M <= N) ) throw new Error('Assertion failed.');

  // https://arxiv.org/abs/1111.4144
  for( let j=0; j < M; j++ )
  {
    for( let k=0; k < j; k++ )
    {
      const V_k = LD[LD_off + N*j+k] * LD[LD_off + N*k+k];

      for( let i=j; i < M; i++ )
        LD[LD_off + N*i+j] -= LD[LD_off + N*i+k] * V_k;
    }

    for( let i=j; ++i < M; )
      LD[LD_off + N*i+j] /= LD[LD_off + N*j+j];
  }
}


export function ldl_decomp(S)
{
  S = asarray(S);
  const
    DTypeArray = ARRAY_TYPES[S.dtype === 'float32' ? 'float32' : 'float64'],
    shape = S.shape,
    [N,M] =   shape.slice(-2);
  S = S.data;

  if( N !== M )
    throw new Error('Last two dimensions must be quadratic.')

  const LD = new DTypeArray(S.length);


  for( let LD_off=0; LD_off < LD.length; LD_off += N*N )
  {
    for( let i=0; i < N; i++ )
    for( let j=0; j <=i; j++ )
      LD[LD_off + N*i+j] = S[LD_off + N*i+j];

    _ldl_decomp(N,N, LD,LD_off);
  }

  return new NDArray(shape,LD);
}


export function _ldl_solve(M,N,O, LD,LD_off, X,X_off)
{
  if( 0 !==LD.length%1 ) throw new Error('Assertion failed.');
  if( 0 !== X.length%1 ) throw new Error('Assertion failed.');
  if( 0 !==LD_off%1 ) throw new Error('Assertion failed.');
  if( 0 !== X_off%1 ) throw new Error('Assertion failed.');
  if( 0 !== M%1 ) throw new Error('Assertion failed.');
  if( 0 !== N%1 ) throw new Error('Assertion failed.');
  if( 0 !== O%1 ) throw new Error('Assertion failed.');

  if(LD.length -LD_off < M*N ) throw new Error('Assertion failed.');
  if( X.length - X_off < M*O ) throw new Error('Assertion failed.');

  if( !(0 < M) ) throw new Error('Assertion failed.');
  if( !(0 < N) ) throw new Error('Assertion failed.');
  if( !(0 < O) ) throw new Error('Assertion failed.');

  if( !(M <= N) ) throw new Error('Assertion failed.');

  // FORWARD SUBSTITUTION
  for( let i=1; i < M; i++ )
  for( let k=0; k < i; k++ )
  for( let j=0; j < O; j++ )
    X[X_off+O*i+j] -= LD[LD_off+N*i+k] * X[X_off+O*k+j];

  // SCALING
  for( let i=0; i < M; i++ )
  for( let j=0; j < O; j++ )
    X[X_off + O*i+j] /= LD[LD_off + N*i+i];

  // BACKWARD SUBSTITUTION
  for( let k=M; k-- > 1; )
  for( let i=k; i-- > 0; )
  for( let j=O; j-- > 0; )
    X[X_off + O*i+j] -= LD[LD_off + N*k+i] * X[X_off + O*k+j]
}


export function ldl_solve(LD,y)
{
  LD= asarray(LD);
  y = asarray(y );
  if(LD.ndim < 2 ) throw new Error('ldl_solve(LD,y): LD must be at least 2D.');
  if( y.ndim < 2 ) throw new Error('ldl_solve(LD,y): y must be at least 2D.');

  const
    [N,M] = LD.shape.slice(-2),
    [I,J] =  y.shape.slice(-2);
  if( N != M ) throw new Error('ldl_solve(LD,y): Last two dimensions of LD must be quadratic.')
  if( I != M ) throw new Error("ldl_solve(LD,y): LD and y don't match.");

  const
    ndim = Math.max(LD.ndim, y.ndim),
    shape = Int32Array.from({ length: ndim }, () => 1 );
  shape[ndim-2] = I;
  shape[ndim-1] = J;

  // FIND COMMON (BROADCASTED) SHAPE
  for( let arr of [LD,y] )
    for( let i=ndim-2, j=arr.ndim-2; i-- > 0 && j-- > 0; )
      if( 1 === shape[i] )
        shape[i] = arr.shape[j];
      else if( shape[i] != arr.shape[j] && arr.shape[j] != 1 )
        throw new Error('Shapes are not broadcast-compatible.');

  // GENERATE RESULT DATA
  const
    dtype = LD.dtype === 'float32'
          && y.dtype === 'float32' ? 'float32' : 'float64',
    x_dat = new ARRAY_TYPES[dtype]( shape.reduce((a,b) => a*b, 1) ),
   LD_dat = LD.data,
    y_dat =  y.data;
  let
   LD_off = 0, LD_stride = 1,
    y_off = 0,  y_stride = 1,
    x_off = 0;

  function solv(d) {
    if( d === ndim-2 ) {
      LD_stride = N*N;
       y_stride = N*J;

      // COPY y
      for( let i=0; i < y_stride; i++ )
        x_dat[x_off+i] = y_dat[y_off+i];

      _ldl_solve(N,N,J, LD_dat,LD_off, x_dat,x_off);

      LD_off += LD_stride;
       y_off += y_stride;
       x_off += y_stride;

      return;
    }
    for( let l=shape[d]; ; l-- ) {
      solv(d+1);
      if( l == 1 ) break;
      if( ! (LD.shape[ d - ndim + LD.ndim ] > 1) ) LD_off -= LD_stride;
      if( ! ( y.shape[ d - ndim +  y.ndim ] > 1) )  y_off -=  y_stride;
    }
    LD_stride *= LD.shape[ d - ndim + LD.ndim ] || 1;
     y_stride *=  y.shape[ d - ndim +  y.ndim ] || 1;
  }
  solv(0);

  return new NDArray(shape,x_dat);
}
