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

import {NDArray} from './nd_array'
import {ARRAY_TYPES, _check_dtype} from './dt'


export function tabulate(shape, dtype, idx2val)
{
  shape = Int32Array.from(shape, s => {
    if( s%1 !== 0 )
      throw new Error(`tabulate(shape, dtype='object', idx2val): Invalid shape [${shape}].`)
    return s
  })
  if( null == idx2val){ idx2val = dtype; dtype = undefined }
  if( null == idx2val) throw new Error("tabulate(shape, dtype='object', idx2val): idx2val missing.")
  if( null == dtype ) dtype = 'object'
  _check_dtype(dtype)

  const
    multi_idx = new Int32Array(shape.length), // <- index in result
    data = new ARRAY_TYPES[dtype]( shape.reduce((a,b) => a*b, 1) )

  let flat_idx = 0

  function write(d) {
    if( d === shape.length )
      data[flat_idx++] = idx2val(...multi_idx)
    else for(
      multi_idx[d] = 0;
      multi_idx[d] < shape[d];
      multi_idx[d]++
    )
      write(d+1)
  }
  write(0)

  return new NDArray(shape,data)
}
