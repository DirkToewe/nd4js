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

import {concat} from './concat'


export function stack(axis, dtype, ndarrays)
{
  if( null == ndarrays )
  {
    if( null == dtype ){ ndarrays = axis; axis = undefined }
    else {
      ndarrays = dtype; dtype = undefined
      if( 'string' === typeof axis ) { dtype = axis; axis = undefined }
    }
  }
  if( ! ('length' in ndarrays) ) ndarrays = [...ndarrays]
  if( null == axis ) axis = 0

  if( 0 > axis )  axis += ndarrays[0].shape.length+1
  if( 0 > axis || axis >  ndarrays[0].shape.length ) throw new Error('Axis out of bounds.')

  ndarrays = ndarrays.map( arr => arr.reshape(
    ...arr.shape.slice(0,axis),
    1,
    ...arr.shape.slice(axis)
  ))
  return concat(axis,dtype,ndarrays)
}
