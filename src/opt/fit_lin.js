'use strict';

/* This file is part of ND4JS.
 *
 * ND4JS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ND4JS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ND4JS. If not, see <http://www.gnu.org/licenses/>.
 */

import {asarray, NDArray} from '../nd_array'

import {lstsq} from '../la/lstsq'


export function fit_lin(x,y, regularization, funcs)
{
  if(funcs == null) {  
     funcs = regularization;
             regularization = 0;
  }
  else if( null == regularization )
    regularization = 0;
  else if( !(0 <= regularization && regularization < Infinity) )
    throw new Error(`fit_lin(x,y, regularization, funcs): Invalid regularization: ${regularization}.`);

  funcs = [...funcs];

  x = asarray('float64', x);
  y = asarray('float64', y);

  if( y.ndim !== 1 )
    throw new Error('fit_lin(x,y, (regularization,) funcs): y.ndim must be 1.');

  const M = y.shape[0],
        N = funcs.length;

  const x_ndim = x.ndim;

  if( x.ndim !== 1 &&
      x.ndim !== 2 )
    throw new Error('fit_lin(x,y, (regularization,) funcs): x.ndim must be 1 or 2.');

  if( x.shape[0] !== y.shape[0] )
    throw new Error('fit_lin(x,y, (regularization,) funcs): x.shape[0] and y.shape[0] must equal.');

  if( x.ndim === 1 )
    x = x.data;
  else {
    const   x_shape = x.shape.slice(1),
      [L] = x_shape;

    x = x.data.slice();
    Object.freeze(x.buffer);
    x = Array.from(
      {length: M},
      (_,i) => Object.freeze(
        new NDArray(x_shape, x.subarray( L*i, L*(i+1) ))
      )
    );
  }

  const A_shape = Int32Array.of( (regularization!==0)*N + M, N )
  let   A = new Float64Array(A_shape[0]*A_shape[1]),
        z = y;

  if( regularization !== 0 )
  {
    regularization = Math.sqrt(regularization);
    for( let i=N; i-- > 0; )
      A[M*N + N*i+i] = regularization;
  }

  for( let i=M; i-- > 0; )
  for( let j=N; j-- > 0; )
    A[N*i+j] = funcs[j](x[i]);

  if( regularization !== 0 )
  {
    y = y.data;
    const z = new Float64Array(M+N);
    for( let i=M; i-- > 0; )
      z[i] = y[i];
    y = new NDArray(Int32Array.of(M+N,1), z);
  }
  else
    y = new NDArray(Int32Array.of(M,1), y.data);

  A = new NDArray(A_shape, A);

  const {data: coeffs} = lstsq(A,y);

  const param_lin_func = x =>
  {
    const {coeffs, funcs} = param_lin_func;

    if( coeffs.length !== funcs.length )
      throw new Error('Assertion failed.');

    let result = 0;
    for( let i=coeffs.length; i-- > 0; )
      result += coeffs[i] * funcs[i](x);

    return result;
  };

  param_lin_func.coeffs = coeffs;
  param_lin_func.funcs  =  funcs;

  return param_lin_func;
}
