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

import {forEachItemIn, CUSTOM_MATCHERS} from '../jasmine_utils'
import {b64_decode,
        b64_encode} from './b64'
import {tabulate} from '../tabulate'
import {WHITESPACES} from '.'


describe('b64', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  const randInt = until => Math.trunc(Math.random()*until)


  forEachItemIn(
    function*(){
      for( let length=0; length < 1337; length++ )
        yield Array.from({length}, () => randInt(256))
                   .map( c => String.fromCharCode(c) )
                   .join('')
    }()
  ).it('b64_decode works for random examples with padding', STR => {
    let B64 = Array.from( btoa(STR) )

    let nInsert = Math.max( 0, randInt(B64.length*2) - B64.length )
    while( nInsert-- > 0 )
    {
      const i = randInt(B64.length+1),
           ws = WHITESPACES[randInt(WHITESPACES.length)]
      B64.splice(i,0,ws)
    }
    
    const bytes = b64_decode(B64),
            str = Array.from(bytes, c => String.fromCharCode(c) ).join('')

    expect(str).toBe(STR)
  })


  forEachItemIn(
    function*(){
      for( let length=0; length < 1337; length++ )
        yield Array.from({length}, () => Math.trunc(Math.random()*256))
                   .map( c => String.fromCharCode(c) )
                   .join('')
    }()
  ).it('b64_decode works for random examples with padding', STR => {
    let B64 = Array.from( btoa(STR).replace('=','') )

    let nInsert = Math.max( 0, randInt(B64.length*2) - B64.length )
    while( nInsert-- > 0 )
    {
      const i = randInt(B64.length+1),
           ws = WHITESPACES[randInt(WHITESPACES.length)]
      B64.splice(i,0,ws)
    }
    
    const bytes = b64_decode(B64),
            str = Array.from(bytes, c => String.fromCharCode(c) ).join('')

    expect(str).toBe(STR)
  })


  forEachItemIn(
    function*(){
      for( let length=0; length < 1337; length++ )
        yield Array.from({length}, () => Math.trunc(Math.random()*256))
                   .map( c => String.fromCharCode(c) )
                   .join('')
    }()
  ).it('b64_encode works for random examples', STR => {
    const bytes = Array.from(STR, c => c.charCodeAt(0))
    const B64 = b64_encode(bytes),
          str = atob(B64)

    // check proper padding
    expect([...B64].filter(s => ! WHITESPACES.includes(s)).length%4).toBe(0)

    expect(str).toBe(STR)
  })
})
