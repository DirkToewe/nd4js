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

import {array, NDArray} from '../nd_array'
import {ARRAY_TYPES, eps} from '../dt'


export const num_grad = ( f, {h_rel=undefined, h_abs=undefined}={} ) => x => {
  x = array(x);

  const dtype = x.dtype==='float32' ? 'float32' : 'float64',
        DTypeArray = ARRAY_TYPES[dtype],
        shape = x.shape;

  x = x.data;

  const  g = new DTypeArray(x.length),
    epsRel = null != h_rel ? h_rel : eps(dtype)**(1/3),
    epsAbs = null != h_abs ? h_abs : eps(dtype)**(1/3);

  // https://en.wikipedia.org/wiki/Finite_difference#Forward,_backward,_and_central_differences
  // https://www.geometrictools.com/Documentation/FiniteDifferences.pdf
  for( let i=x.length; i-- > 0; )
  {
    const h = Math.max(Math.abs(x[i])*epsRel, epsAbs); // TODO maybe a check for h===0 might be sensible

    for( const [d,w] of [[+2*h,-1],
                         [+1*h,+8],
                         [-1*h,-8],
                         [-2*h,+1]] )
    {
      const     X = new NDArray( shape, DTypeArray.from(x) );
                X.data[i] += d;
      g[i] += f(X) * w;
    }
    g[i] /= 12*h;
  }

  return new NDArray(shape, g)
}


export const num_grad_forward = ( f, {h_rel=undefined, h_abs=undefined}={} ) => x => {
  x = array(x)

  const dtype = x.dtype==='float32' ? 'float32' : 'float64',
        DTypeArray = ARRAY_TYPES[dtype],
        shape = x.shape
  x = x.data;

  const  g = new DTypeArray(x.length),
    epsRel = null != h_rel ? h_rel : eps(dtype)**(1/3),
    epsAbs = null != h_abs ? h_abs : eps(dtype)**(1/3);

  const f0 = f( new NDArray(shape, x.slice()) );
  g.fill(-1.5*f0);

  // http://macs.citadel.edu/chenm/344.dir/temp.dir/lect4_1.pdf
  for( let i=x.length; i-- > 0; )
  {
    const h = Math.max(Math.abs(x[i])*epsRel, epsAbs) // TODO maybe a check for h===0 might be sensible
    for( const [d,w] of [[1*h,+2.0],
                         [2*h,-0.5]] )
    {
      const     X = new NDArray( shape, DTypeArray.from(x) );
                X.data[i] += d;
      g[i] += f(X) * w;
    }
    g[i] /= h;
  }

  return new NDArray(shape, g)
}
