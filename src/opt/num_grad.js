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


export const num_grad = (f,eps=2**-12) => x => {
  x = asarray(x)

  const DTypeArray = ARRAY_TYPES[x.dtype==='float32' ? 'float32' : 'float64'],
        shape = x.shape
  x = x.data

  const g = new DTypeArray(x.length),
      eps = Math.sqrt(Number.EPSILON)

  // https://en.wikipedia.org/wiki/Finite_difference#Forward,_backward,_and_central_differences
  // https://www.geometrictools.com/Documentation/FiniteDifferences.pdf
  for( let i=x.length; i-- > 0; )
  {
    const h = Math.max(Math.abs(x[i])*eps, eps)
    for( const [d,w] of [[+2*h,-1],
                         [+1*h,+8],
                         [-1*h,-8],
                         [-2*h,+1]] )
    {
      const X = new NDArray(shape, DTypeArray.from( x, (xj,j) => xj + d*(i===j) ) )
      g[i] += w * f(X)
    }
    g[i] /= 12*h 
  }

  return new NDArray(shape, g)
}
