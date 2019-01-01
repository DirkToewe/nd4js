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

import {asarray} from '../../nd_array'
import {tabulate} from '../../tabulate'


export function rosenbrock( x )
{
  x = asarray(x)

  if(         x.ndim    < 1 ) throw new Error('rosenbrock(x): x.ndim must be at least 1.');
  if( x.shape[x.ndim-1] < 2 ) throw new Error('rosenbrock(x): x.shape[-1] must be at least 2.');

  const N = x.shape[x.ndim-1],
    dtype = x.dtype==='float32' ? 'float32' : 'float64'

  return tabulate(x.shape.slice(0,-1), dtype, (...i) => {
    let sum = 0
    for( let j=N-1; j-- > 0; )
      sum += 100*( x(...i,j+1) - x(...i,j)**2 )**2 + ( 1 - x(...i,j) )**2
    return sum
  })
}


export function rosenbrock_grad( x )
{
  x = asarray(x)

  if(         x.ndim    < 1 ) throw new Error('rosenbrock(x): x.ndim must be at least 1.');
  if( x.shape[x.ndim-1] < 2 ) throw new Error('rosenbrock(x): x.shape[-1] must be at least 2.');

  const N = x.shape[x.ndim-1],
    dtype = x.dtype==='float32' ? 'float32' : 'float64'

  return x.mapElems(dtype, (xj,...i) => {
    const j = i.pop()

    let result = 0
    if(j < N-1) result += 400*( xj**2 - x(...i,j+1)    )*xj - 2*(1 - xj)
    if(j >  0 ) result += 200*( xj    - x(...i,j-1)**2 )
    return result
  })
}
