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


// REFERENCES
// ----------
// .. [1] "Testing Unconstrained Optimization Software"
//         Jorge J. Moré, Burton S. Garbow, Kenneth E. Hillstrom
//         ACM Transactions on Mathematical Software, Vol 7, No. 1, March 1982, pp. 17-41
// .. [2] "A rapidly convergent descent method for minimization"
//         R. Fletcher, M.J.D. Powell
//         Comput. J. 6 (1963), 163-168


const atan2 = (y,x) => Math.atan2(
  y * ( (x < 0 || Object.is(x,-0)) ? -1 : +1 ),
  Math.abs(x)
);


export const helical_valley = x =>
{
  x = asarray(x)

  if(         x.ndim     <  1 ) throw new Error('helical_valley(x): x.ndim must be at least 1.');
  if( x.shape[x.ndim-1] !== 3 ) throw new Error('helical_valley(x): x.shape[-1] must be 3.');

  const F_shape = x.shape.slice(0,-1),
        F       = new (x.dtype==='float32' ? Float32Array
                                           : Float64Array)(x.data.length/3);
  x = x.data;

  for( let x_off=x.length,
           F_off=F.length;
           F_off-- > 0; )
  {        x_off -= 3;
    const x1 = x[x_off+0],
          x2 = x[x_off+1],
          x3 = x[x_off+2];

    const f1 = x3*10 - 50*( atan2(x2,x1) / Math.PI + (x1 < 0) ),
          f2 =    10 * Math.hypot(x2,x1) - 10,
          f3 = x3;

    F[F_off] = f1*f1 + f2*f2 + f3*f3;
  }

  return new NDArray(F_shape, F);
}


helical_valley.nIn =
helical_valley.nOut= 3;


helical_valley.minima       =
helical_valley.minima_global=
helical_valley.roots        = [[1,0,0]];


helical_valley.grad = x =>
{
  x = asarray(x)

  if(         x.ndim     <  1 ) throw new Error('helical_valley.grad(x): x.ndim must be at least 1.');
  if( x.shape[x.ndim-1] !== 3 ) throw new Error('helical_valley.grad(x): x.shape[-1] must be 3.');

  const G_shape = x.shape,
        G       = new (x.dtype==='float32' ? Float32Array
                                           : Float64Array)(x.data.length);
  x = x.data;

  for( let off=x.length; (off-=3) >= 0; )
  {
    const x1 = x[off+0],
          x2 = x[off+1],
          x3 = x[off+2];

    const f1 = x3*10 - 50*( atan2(x2,x1) / Math.PI + (x1 < 0) );

    G[off+0] = 100*( x1*(2  -  2/Math.hypot(x2,x1))  +  x2*f1 / (Math.PI * (x2*x2 + x1*x1)) );
    G[off+1] = 100*( x2*(2  -  2/Math.hypot(x2,x1))  -  x1*f1 / (Math.PI * (x2*x2 + x1*x1)) );
    G[off+2] = f1*20 + x3*2;
  }

  return new NDArray(G_shape, G);
}


helical_valley.hess = x =>
{
  x = asarray(x)

  if(         x.ndim     <  1 ) throw new Error('helical_valley.hess(x): x.ndim must be at least 1.');
  if( x.shape[x.ndim-1] !== 3 ) throw new Error('helical_valley.hess(x): x.shape[-1] must be 2.');

  const H_shape = Int32Array.of(...x.shape, 3),
        H = new (x.dtype==='float32' ? Float32Array
                                     : Float64Array)(x.data.length*3);
  x = x.data;

  for( let H_off=H.length,
           x_off=x.length;
          (x_off -= 3) >= 0; )
  {        H_off -= 9;
    const x1 = x[x_off+0],
          x2 = x[x_off+1],
          x3 = x[x_off+2];

    const f1 = x3*10 - 50*( atan2(x2,x1) / Math.PI + (x1 < 0) );

    const dot = x1*x1 + x2*x2;

    H[H_off+0] = 200 * (  1  -  x2*x2 / dot**1.5  +  (25*x2/Math.PI - x1*f1) * x2 / (Math.PI * dot*dot)  );
    H[H_off+1] = 100/Math.PI * f1/dot + 200*x1*x2 / dot**1.5 - (5000*x1/Math.PI + 200*x2*f1) * x2 / (Math.PI * dot*dot);
    H[H_off+2] = +1000/Math.PI * x2/dot;

    H[H_off+3] = H[H_off+1];
    H[H_off+4] = 200 * ( 1  -  x1*x1 / dot**1.5  +  (25*x1/Math.PI + x2*f1) * x1 / (Math.PI * dot*dot) );
    H[H_off+5] = -1000/Math.PI * x1/dot;

    H[H_off+6] = H[H_off+2];
    H[H_off+7] = H[H_off+5];
    H[H_off+8] = 200 + 2;
  }

  return new NDArray(H_shape, H);
}


helical_valley.lsq = x =>
{
  x = asarray(x)

  if(         x.ndim     <  1 ) throw new Error('helical_valley.lsq(x): x.ndim must be at least 1.');
  if( x.shape[x.ndim-1] !== 3 ) throw new Error('helical_valley.lsq(x): x.shape[-1] must be 3.');

  const F_shape = x.shape.slice(),
        F       = new (x.dtype==='float32' ? Float32Array
                                           : Float64Array)(x.data.length);
  x = x.data;

  for( let off = x.length;
          (off-= 3) >= 0; )
  {
    const x1 = x[off+0],
          x2 = x[off+1],
          x3 = x[off+2];

    F[off+0] = x3*10 - 50*( atan2(x2,x1) / Math.PI + (x1 < 0) );
    F[off+1] =    10 * Math.hypot(x2,x1) - 10;
    F[off+2] = x3;
  }

  return new NDArray(F_shape, F);
}


helical_valley.lsq_jac = x =>
{
  x = asarray(x)

  if(         x.ndim     <  1 ) throw new Error('helical_valley.lsq_jac(x): x.ndim must be at least 1.');
  if( x.shape[x.ndim-1] !== 3 ) throw new Error('helical_valley.lsq_jac(x): x.shape[-1] must be 3.');

  const J_shape = Int32Array.of(...x.shape, 3),
        J       = new (x.dtype==='float32' ? Float32Array
                                           : Float64Array)(x.data.length*3);
  x = x.data;

  for( let x_off = x.length,
           J_off = J.length;
          (x_off-= 3) >= 0; )
  {        J_off-= 9;
    const x1 = x[x_off+0],
          x2 = x[x_off+1];

    J[J_off+0] = +50/Math.PI * x2 / (x2*x2 + x1*x1);
    J[J_off+1] = -50/Math.PI * x1 / (x2*x2 + x1*x1);
    J[J_off+2] = 10;
    J[J_off+3] = 10 * x1 / Math.hypot(x2,x1);
    J[J_off+4] = 10 * x2 / Math.hypot(x2,x1);
    J[J_off+5] = 0;
    J[J_off+6] = 0;
    J[J_off+7] = 0;
    J[J_off+8] = 1;
  }

  return new NDArray(J_shape, J);
}


Object.freeze(helical_valley);
