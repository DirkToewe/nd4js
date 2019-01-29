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

import {Complex} from './complex'
import {Complex128Array} from './complex_array'

export * from './complex'
export * from './complex_array'

export const ARRAY_TYPES = {
       'int32':      Int32Array,
     'float32':    Float32Array,
     'float64':    Float64Array,
  'complex128': Complex128Array,
     'object' :           Array
}

export function eps( dtype )
{
  _check_dtype(dtype)

  switch(dtype)
  {
    case 'complex128': return new Complex(Number.EPSILON)
    case    'float32': return 1.1920928955078125e-7
    default          : return Number.EPSILON
  }
}

export function cast_scalar(x, dtype)
{
  if( dtype ===      'int32' ) return x & 0xFFFFFFFF;
  if( dtype ===    'float32' ) return Math.fround(x);
  if( dtype === 'complex128' ) return x instanceof Complex ? x : new Complex(x);
  return x*1;
}

export function _check_dtype(dtype)
{
  if( ! ARRAY_TYPES.hasOwnProperty(dtype) )
    throw new Error("Invalid dtype '" + dtype + "'. Must be one of {'" + Object.getOwnPropertyNames(nd.dtypes).join("', '") + "'}.")
}

export function dtypeof(value)
{
  if( value % 1 === 0 )
  {
    if(    value <= ~(1 << 31)
        && value >=  (1 << 31) )
      return 'int32'
    return 'object'
  }
  if( value * 1 == value ) return 'float64'
  if( value instanceof Complex ) return 'complex128';
  return 'object'
}
  
export const super_dtype = (...dtypes) => dtypes.reduce( (dtype1,dtype2) => {
  _check_dtype(dtype1)
  _check_dtype(dtype2)
  if( dtype1 ===     'object' || dtype2 ===     'object' ) return 'object'
  if( dtype1 === 'complex128' || dtype2 === 'complex128' ) return 'complex128'
  if( dtype1 ===    'float64' || dtype2 ===    'float64' ) return 'float64'
  if( dtype1 ===    'float32' || dtype2 ===    'float32' ) return 'float32'
  return 'int32'
})
  
export function is_subdtype(sub_dtype, sup_dtype)
{
  _check_dtype(sub_dtype)
  _check_dtype(sup_dtype)
  const rank = {
        'int32': 0,
      'float32': 1,
      'float64': 2,
   'complex128': 3,
       'object': 4
  }
  return rank[sub_dtype] <= rank[sup_dtype]
}
