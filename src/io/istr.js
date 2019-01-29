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

import {ARRAY_TYPES, _check_dtype} from '../dt'
import {b64_encode,
        b64_decode} from './b64'
import {IS_LITTLE_ENDIAN} from '.'
import {asarray, NDArray} from '../nd_array'


export function istr_to_nd( b64_chars )
{  
  b64_chars = b64_chars[Symbol.iterator]()

  let dtype = []
  for(;;){
    const {value, done} = b64_chars.next()
    if(done) throw new Error('b64_parse(b64_chars): Invalid b64_chars.')
    if(value === '[') break
    dtype.push(value)
  }
  dtype = dtype.join('').trim()
  if(dtype === '') throw new Error('b64_parse(b64_chars): dtype=object not (yet) supported.')
  _check_dtype(dtype)

  let shape = []
  for(let s=[];;)
  {
    const {value, done} = b64_chars.next()
    if(done) throw new Error('b64_parse(b64_chars): Invalid b64_chars.')
    if(value===']' && s.length=== 0 ) break
    if(value===']' ||    value===',') {
      const d = s.join('')
      if( isNaN(d) )
        throw new Error(`b64_parse(b64_chars): Invalid shape entry: "${d}".`)
      shape.push(1*d)
      if(value===']') break
      s = []
      continue
    }
    s.push(value)
  }
  shape = Int32Array.from(shape)

  b64_chars = b64_decode(b64_chars)

  const DTypeArray = ARRAY_TYPES[dtype],
    data = new Uint8Array(shape.reduce((m,n) => m*n, DTypeArray.BYTES_PER_ELEMENT) ),
    word = new Uint8Array(                           DTypeArray.BYTES_PER_ELEMENT)
        
  for( let i=0; i < data.length; )
  {
    for( let j=0; j < word.length; j++ )
    {
      const {value, done} = b64_chars.next()
      if(done) throw new Error('b64_parse(b64_chars): b64_chars too short.')
      word[j] = value
    }

    if( ! IS_LITTLE_ENDIAN )
      word.reverse()

    for( let j=0; j < word.length; j++,i++ )
      data[i] = word[j]
  }

  return new NDArray(shape, new DTypeArray(data.buffer))
}


export function* nd_to_istr( ndarray )
{
  ndarray = asarray(ndarray)
  const data = ndarray.data
  if( ndarray.dtype === 'object' )
    throw new Error(`b64_format(A): A.dtype='${ndarray.dtype}' not supported.`)

  if( ! IS_LITTLE_ENDIAN )
    throw new Error('Big endianness not (yet) supported.')

  yield* `${ndarray.dtype}[${ndarray.shape.join(',')}]\n`
  yield* b64_encode(
    new Uint8Array(data.buffer, data.byteOffset, data.byteLength),
    {lineLimit: 128}
  )
}
