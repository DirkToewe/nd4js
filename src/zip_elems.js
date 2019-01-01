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

import {ARRAY_TYPES, _check_dtype} from './dt'
import {asarray, NDArray} from './nd_array'


export function zip_elems( ndarrays, dtype, zip_fn )
{
  // PREPROCESS ARGUMENTS
  if( null == zip_fn && dtype instanceof Function ) { zip_fn = dtype; dtype = undefined }

  if( ! (ndarrays instanceof Array) )
    ndarrays = [...ndarrays]

  ndarrays = ndarrays.map(asarray)

  if( null == dtype ) dtype = 'object'
  if( null == zip_fn ) {
    if( 'object' !== dtype )
      throw new Error('If zip_fn is undefined, dtype must be "object" or undefined.')
    const L = ndarrays.length
    zip_fn = (...x) => x.slice(0,L) 
  }

  _check_dtype(dtype)

  const
    ndim = ndarrays.reduce( (ndim,arr) => Math.max(ndim,arr.ndim), 0),
    shape = Int32Array.from({ length: ndim }, () => 1 )

  // FIND COMMON (BROADCASTED) SHAPE
  for( let arr of ndarrays )
    for( let i=ndim, j=arr.ndim; i-- > 0 && j-- > 0; )
      if( 1 === shape[i] )
        shape[i] = arr.shape[j]
      else if( shape[i] != arr.shape[j] && arr.shape[j] != 1 )
        throw new Error('Shapes are not broadcast-compatible.')

  // GENERATE RESULT DATA
  const
    multi_idx = new Int32Array(ndim), // <- index in result
    data = new ARRAY_TYPES[dtype]( shape.reduce((a,b) => a*b, 1) ),

    values  = new      Array(ndarrays.length), // <- cache of ndarrays[indices]
    indices = new Int32Array(ndarrays.length), // <- indices in ndarrays(s)
    strides = new Int32Array(ndarrays.length);

  let flat_idx = 0

  function write(d) {
    if( d === ndim ) {
      strides.fill(1)
      for( let i=ndarrays.length; i-- > 0; )
        values[i] = ndarrays[i].data[indices[i]++]
      data[flat_idx++] = zip_fn(...values, ...multi_idx)
      return
    }
    for( multi_idx[d] = 0;; )
    {
      write(d+1)
      if( ++multi_idx[d] >= shape[d] )
        break
      for( let i=ndarrays.length; i-- > 0; )
        if( ! (ndarrays[i].shape[ d - ndim + ndarrays[i].ndim ] > 1) ) // <- handles undefined (index out of bounds)
          indices[i] -= strides[i]
    }
    for( let i=ndarrays.length; i-- > 0; )
      strides[i] *= ( ndarrays[i].shape[ d - ndim + ndarrays[i].ndim ] || 1 )
  }
  write(0)

  return new NDArray(shape,data)
}
