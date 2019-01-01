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


const
  BYTE_TO_CHAR = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/',
  WHITESPACES  = '\f\n\r\t\v ',
  CHAR_TO_BYTE = Int32Array.from({length: 256}, () => -2)

if( BYTE_TO_CHAR.length !== 64 )
  throw new Error('Assertion failed.')

for( let i=BYTE_TO_CHAR.length; i-- > 0; ) CHAR_TO_BYTE[BYTE_TO_CHAR.charCodeAt(i)] =  i
for( let i= WHITESPACES.length; i-- > 0; ) CHAR_TO_BYTE[ WHITESPACES.charCodeAt(i)] = -1

export function* b64_decode( b64_chars )
{
  b64_chars = b64_chars[Symbol.iterator]()

  let done

  function next6() {
    while(true)
    {
      let value
      ({value,done} = b64_chars.next())
      done = done || '=' === value
      if(done) return NaN

      if('string' !== typeof value || value.length !== 1)
        throw new Error(`_b64_decode(b64_chars): Invalid token in character sequence: ${value}.`)

      value = CHAR_TO_BYTE[value.charCodeAt(0)]
      // skip whitespaces
      if(value === -1) continue
  
      if( ! (value >= 0) ) // <- handles undefined
        throw new Error(`_b64_decode(b64_chars): Invalid base64 character '${value}'.`)
  
      return value
    }
  }

  function warn() {
    console.warn(`_b64_decode(b64_chars): Base64 sequence ended unexpectedly.`)
  }

  while(true) {
    const bit0to5  = next6(); if(done) return
    const bit6to11 = next6(); if(done) return warn(); yield ( bit0to5             << 2) | (bit6to11  >> 4)
    const bit12to17= next6(); if(done) return       ; yield ((bit6to11  & 0b1111) << 4) | (bit12to17 >> 2)
    const bit18to23= next6(); if(done) return       ; yield ((bit12to17 & 0b0011) << 6) | (bit18to23 & 0b111111)
  }
}


export function* b64_encode( bytes, { pad=true, lineLimit=Infinity }={})
{
  if( ! (0 < lineLimit) )
    throw new Error(`_b64_encode(bytes,{lineLimit}): Invalid lineLimit: ${lineLimit}.`)
  const rawBytes = bytes[Symbol.iterator]()
  bytes = function*(){
    for( let byte of rawBytes ) {
      byte *= 1
      if( ! (0 <= byte && byte < 256) )
        throw new Error(`_b64_encode(bytes): Invalid token in byte sequence: ${byte}`)
      yield byte
    }
  }()

  const sixPacks = function*()
  {  
    let value,done;
    while(true) {
      let six
      ({value,done} = bytes.next()); if(done)           return; yield value >> 2;       six = (value & 0b11  ) << 4;
      ({value,done} = bytes.next()); if(done){yield six;return} yield value >> 4 | six; six = (value & 0b1111) << 2;
      ({value,done} = bytes.next()); if(done){yield six;return} yield value >> 6 | six; yield  value & 0b111111
    }
  }()

  let i=0
  for( const six of sixPacks )
  {
    if( 'number' !== typeof six  ) throw new Error('Assertion failed.')
    if( ! (0 <= six && six < 64) ) throw new Error('Assertion failed.')
    yield BYTE_TO_CHAR[six]
    if( 0 == ++i % lineLimit )
      yield '\n'
  }

  // padding
  if(pad)
    while( 0 !== i++ % 4 )
      yield '='
}
