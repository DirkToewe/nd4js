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
 * along with ND.JS. If not, see <http://www.gnu.org/licenses/>.
 */

import {asarray, NDArray} from '../../nd_array'
import {tabulate} from '../../tabulate'


export function rosenbrock( x )
{
  // https://en.wikipedia.org/wiki/Rosenbrock_function
  x = asarray(x)

  if(         x.ndim    < 1 ) throw new Error('rosenbrock(x): x.ndim must be at least 1.');
  if( x.shape[x.ndim-1] < 2 ) throw new Error('rosenbrock(x): x.shape[-1] must be at least 2.');

  const N = x.shape[x.ndim-1];

  const F_shape = x.shape.slice(0,-1),
        F       = new (x.dtype==='float32' ? Float32Array
                                           : Float64Array)(x.data.length/N);
  x = x.data;

  for( let x_off=x.length,
           F_off=F.length;
           F_off-- > 0; )
  {        x_off -= N;
    for( let i=N-1; i-- > 0; ) {
      const xj = x[x_off + i+1],
            xi = x[x_off + i],
             u = xj - xi*xi,
             v = 1-xi;
      F[F_off] += 100*u*u + v*v;
    }
  }

  return new NDArray(F_shape, F);
}


export function rosenbrock_grad( x )
{
  x = asarray(x)

  if(         x.ndim    < 1 ) throw new Error('rosenbrock(x): x.ndim must be at least 1.');
  if( x.shape[x.ndim-1] < 2 ) throw new Error('rosenbrock(x): x.shape[-1] must be at least 2.');

  const N = x.shape[x.ndim-1];

  const G_shape = x.shape,
        G       = new (x.dtype==='float32' ? Float32Array
                                           : Float64Array)(x.data.length);
  x = x.data;

  for( let off=x.length;
           off >= 0; )
  {        off -= N;
    for( let i=N; i-- > 0; )
    {
      const xh = x[off + i-1],
            xi = x[off + i  ],
            xj = x[off + i+1];

      if(     i < N-1 ) G[off + i]  = 400*(xi*xi - xj   )*xi - 2*(1-xi);
      if( 0 < i       ) G[off + i] += 200*(xi    - xh*xh);
    }
  }

  return new NDArray(G_shape, G);
}


export function rosenbrock_hess( x )
{
  x = asarray(x)

  if(         x.ndim    < 1 ) throw new Error('rosenbrock(x): x.ndim must be at least 1.');
  if( x.shape[x.ndim-1] < 2 ) throw new Error('rosenbrock(x): x.shape[-1] must be at least 2.');

  const N = x.shape[x.ndim-1];

  const H_shape = Int32Array.of(...x.shape, N),
        H = new (x.dtype==='float32' ? Float32Array
                                     : Float64Array)(x.data.length*N);
  x = x.data;

  for( let H_off=H.length,
           x_off=x.length;
           x_off >= 0; )
  {        x_off -= N;
           H_off -= N*N;
    for( let i=N; i-- > 0; ) {
      if( 0 < i ) {
        H[H_off + N* i   + i   ] = 200;
        H[H_off + N*(i-1)+ i   ] =
        H[H_off + N* i   +(i-1)] = -400*x[x_off + i-1];
      }
      if( i < N-1 ) {
        const xi = x[x_off + i];
        H[H_off + N*i+i] -= 400*x[x_off + i+1] - 1200*xi*xi  -  2;
      }
    }
  }

  return new NDArray(H_shape, H);
}


export function rosenbrock_lsq( x )
{
  x = asarray(x)

  if(         x.ndim    < 1 ) throw new Error('rosenbrock_lsq(x): x.ndim must be at least 1.');
  if( x.shape[x.ndim-1] < 2 ) throw new Error('rosenbrock_lsq(x): x.shape[-1] must be at least 2.');

  const N = x.shape[x.ndim-1],
        M = (N-1)*2;

  const F_shape = x.shape.slice(),
        F       = new (x.dtype==='float32' ? Float32Array
                                           : Float64Array)(x.data.length/N*M);
  F_shape[F_shape.length-1] = M;
  x = x.data;

  for( let x_off=x.length,
           F_off=F.length;
          (F_off -= M) >= 0; )
  {        x_off -= N;
    for( let i=N-1; i-- > 0; ) {
      const xj = x[x_off + i+1],
            xi = x[x_off + i],
             u = xj - xi*xi,
             v = 1-xi;
      F[F_off + 2*i+0] = u*10;
      F[F_off + 2*i+1] = v;
    }
  }

  return new NDArray(F_shape, F);
}


export function rosenbrock_lsq_jac( x )
{
  x = asarray(x)

  if(         x.ndim    < 1 ) throw new Error('rosenbrock(x): x.ndim must be at least 1.');
  if( x.shape[x.ndim-1] < 2 ) throw new Error('rosenbrock(x): x.shape[-1] must be at least 2.');

  const N = x.shape[x.ndim-1],
        M = (N-1)*2;

  const J_shape = Int32Array.of(...x.shape.slice(0,-1),M,N),
        J       = new (x.dtype==='float32' ? Float32Array
                                           : Float64Array)(x.data.length*M);
  x = x.data;

  for( let x_off=x.length,
           J_off=J.length;
          (J_off -= M*N) >= 0; )
  {        x_off -= N;
    for( let i=N-1; i-- > 0; ) {
      const  j = i+1,
            xj = x[x_off + j],
            xi = x[x_off + i];
      J[J_off + N*(2*i+0)+i] += -20*xi;
      J[J_off + N*(2*i+0)+j] +=  10;
      J[J_off + N*(2*i+1)+i] -=  1;
    }
  }

  return new NDArray(J_shape, J);
}
