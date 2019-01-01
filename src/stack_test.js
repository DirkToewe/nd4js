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

import {array} from './nd_array'
import {stack} from './stack'


describe('nd.stack', () => {
  it('works on 3 arrays, each of shape [2,3]', () => {
    const
      a = array([
        [11,12,13],
        [21,22,23]
      ]),
      b = array([
        [31,32,33],
        [41,42,43]
      ]),
      c = array([
        [51,52,53],
        [61,62,63]
      ])
  
    let
    d = stack(0,[a,b,c])
    expect(d.shape).toEqual( Int32Array.of(3,2,3) )
    expect(d.data ).toEqual( Int32Array.of(
      11, 12, 13,
      21, 22, 23,

      31, 32, 33,
      41, 42, 43,

      51, 52, 53,
      61, 62, 63
    ))
    
    d = stack(1,[a,b,c])
    expect(d.shape).toEqual( Int32Array.of(2,3,3) )
    expect(d.data ).toEqual( Int32Array.of(
      11,12,13, 31,32,33, 51,52,53,
      21,22,23, 41,42,43, 61,62,63
    ))
    
    d = stack(2,[a,b,c])
    expect(d.shape).toEqual( Int32Array.of(2,3,3) )
    expect(d.data ).toEqual( Int32Array.of(
      11,31,51,
      12,32,52,
      13,33,53,

      21,41,61,
      22,42,62,
      23,43,63
    ))
  })
})
