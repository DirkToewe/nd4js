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
import {NDArray} from '../nd_array'


export function eye(...shape)
{
  const dtype = shape[0] in ARRAY_TYPES ? shape.shift() : 'float64'

  if( shape.length <  1 ) throw new Error('Size parameter missing.')
  if( shape.length == 1 ) shape.push( shape[shape.length-1] )

  shape = Int32Array.from(shape, s => {
    if( 0 !== s%1 )
      throw new Error(`eye(): Invalid shape [${shape}].`)
    return s
  })

  const I = new ARRAY_TYPES[dtype]( shape.reduce((m,n) => m*n) ),
     [M,N]= shape.slice(-2)

  for( let off = I.length; (off -= M*N) >= 0; )
  for( let i=Math.min(M,N); i-- > 0; )
    I[off + N*i+i] = 1

  const result = new NDArray(shape, I)

  return result
}
