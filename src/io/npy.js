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

import {pyon_parse} from './pyon'
import {ARRAY_TYPES} from '../dt'
import {asarray, NDArray} from '../nd_array'
import {IS_LITTLE_ENDIAN} from '.'


const MAGIC_STRING = '\u0093NUMPY';


export function npy_serialize( A )
{
  return Uint8Array.from( npy_serialize_gen(A) )
}


export function* npy_serialize_gen( A )
{
  A = asarray(A)

  let dt = function(){
    switch(A.dtype) {
      default: throw new Error(`nd_to_npy: A.dtype=${A.dtype} not yet supported.`);
      case      'int32': return  'i4'
      case    'float32': return  'f4'
      case    'float64': return  'f8'
      case 'complex128': return 'c16'
    }
  }();
  if( IS_LITTLE_ENDIAN ) dt = '<'+dt;
  else                   dt = '>'+dt;

  const header = `{"descr": "${dt}", "fortran_order": False, "shape": (${A.shape}${A.ndim > 0 ? ',': ''})}`;

  // MAGIC STRING
  for( let i=0; i < MAGIC_STRING.length; i++ )
    yield MAGIC_STRING.codePointAt(i);

  // VERSION
  yield 1
  yield 0

  const headerLen = header.length+11+63 >>> 6 << 6;
  if( headerLen > 0xFFFF )
    throw new Error('nd_to_npy: Header too large.');

  yield headerLen-10>>>0 & 255
  yield headerLen-10>>>8 & 255

  if( headerLen%64 !== 0 )
    throw new Error('Assertion failed.');

  for( let i=0; i < header.length; i++ )
    yield header.codePointAt(i)

  for( let i=header.length+11; i < headerLen; i++ )
    yield ' '.codePointAt(0)
  yield '\n'.codePointAt(0)

  yield* new Uint8Array(
    A.data.buffer,
    A.data.byteOffset,
    A.data.byteLength
  );
}


export function npy_deserialize( npy_bytes )
{
  let nRead = 0;
  const next = function(){
    const iter = npy_bytes[Symbol.iterator]();
    return () => {
      const {value, done} = iter.next();
      ++nRead;
      if(done) throw new Error('npy_to_nd: byte sequence ended unexpectedly.');
      return value;
    };
  }()
  npy_bytes = undefined;

  let magic_str = '';
  while( magic_str.length < 6 )
    magic_str += String.fromCodePoint( next() );

  if( magic_str !== MAGIC_STRING )
    throw new Error("npy_to_nd: byte sequence does not start with '\u0093NUMPY'.");

  const version = `${next()}.${next()}`;
  switch(version) {
    default: throw new Error(`npy_bytes: npy-file version ${version} not supported.`);
    case '1.0':
    case '2.0':
  }

  const header = function(){
    let headerLen  = next() * (1<< 0)
                  +  next() * (1<< 8)
    if( version === '2.0' ) {
        headerLen += next() * (1<<16)
                  +  next() * (1<<24)
    }
    if( headerLen < 0 ) throw new Error('Assertion failed.');

    return pyon_parse(
      function*(){
        while( headerLen-- > 0 ) yield String.fromCodePoint( next() )
      }()
    );
  }();

  const LE = function(){
    switch(header.descr[0])
    {
      default : throw new Error(`npy_to_nd: dtype '${header.descr}' not yet supported.`);
      case '<': return true;
      case '>': return false;
    }
  }();

  const dtype = function(){
    switch(header.descr.slice(1))
    {
      default : throw new Error(`npy_to_nd: dtype '${header.descr}' not yet supported.`);
      case  'i4': return      'int32';
      case  'f4': return    'float32';
      case  'f8': return    'float64';
      case 'c16': return 'complex128';
    }
  }();

  const DTypeArray = ARRAY_TYPES[dtype];

  const data = Uint8Array.from(
    { length: header.shape.reduce((m,n) => m*n, DTypeArray.BYTES_PER_ELEMENT) },
    next
  );

  if( LE !== IS_LITTLE_ENDIAN )
  {
    const wordLen = DTypeArray.BYTES_PER_ELEMENT / (dtype.startsWith('complex') ? 2 : 1);
    for( let off = 0; off < data.length; off += wordLen )
    for( let i=off, j=off+wordLen; i < --j; i++ ) {
      const data_i = data[i];
                     data[i] = data[j];
                               data[j] = data_i;
    }
  }

  const shape = Int32Array.from(header.shape);
  if( header.fortran_order )
    shape.reverse()

  let result = new NDArray( shape, new DTypeArray(data.buffer) );

  if( header.fortran_order && result.ndim > 1 )
    result = result.transpose(
      ...function*(){
        for( let i=result.ndim; i-- > 0; ) yield i;
      }()
    );

  return result;
}
