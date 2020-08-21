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

import {asarray, NDArray} from './nd_array'
import {ARRAY_TYPES, super_dtype} from './dt'


export function concat(axis, dtype, ndarrays)
{
  if( null == ndarrays )
  {
    if( null == dtype ){ ndarrays = axis; axis = undefined }
    else {
      ndarrays = dtype; dtype = undefined
      if( 'string' === typeof axis ) { dtype = axis; axis = undefined }
    }
  }

  ndarrays = Array.from( ndarrays, arr => asarray(dtype,arr) );

  if( null == axis ) axis = 0
  if( null == dtype) dtype = super_dtype( ...ndarrays.map( a => a.dtype ) );

  if( 0 > axis )  axis += ndarrays[0].shape.length
  if( 0 > axis || axis >= ndarrays[0].shape.length ) throw new Error('Axis out of bounds.')

  const newShape = Int32Array.from(ndarrays[0].shape)
  for( let i=ndarrays.length; --i > 0; )
  {
    const shape = ndarrays[i].shape

    if(   newShape.length != shape.length )
      throw new Error('All shapes must have the same length.')
    if( ! newShape.every( (len,d) => newShape[d] === shape[d] || axis === d ) )
      throw new Error('Shape along all axes but the concatentation axis must match.')

    newShape[axis] += shape[axis]
  }

  const
    rest = newShape.slice(axis+1).reduce((a,b) => a*b, 1),
    indices = Int32Array.from(ndarrays, ndarr => ndarr.data.length),
    newData = new ARRAY_TYPES[dtype]( newShape.reduce((a,b) => a*b, 1) )
  let newIdx = newData.length

  function fill(d)
  {
    if( d === axis )
      for( let i=ndarrays.length; i-- > 0; )
        for( let j=ndarrays[i].shape[d]*rest; j-- > 0; )
          newData[--newIdx] = ndarrays[i].data[--indices[i]]
    else for( let i=newShape[d]; i-- > 0; )
      fill(d+1)
  }
  fill(0)

  return new NDArray(newShape, newData)
}
