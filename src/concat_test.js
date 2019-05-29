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
import {array} from './nd_array'


describe('nd.concat', () => {
  it('works on arrays of shapes [2,2,3], [2,1,3] and [2,3,3] along axis 1', () => {
    const a = array(
            [[[111,112,113],
              [121,122,123]],
             [[211,212,213],
              [221,222,223]]]
          ),
          b = array(
            [[[131,132,133]],
             [[231,232,233]]]
          ),
          c = array(
            [[[141,142,143],
              [151,152,153],
              [161,162,163]],
             [[241,242,243],
              [251,252,253],
              [261,262,263]]]
          )

    for( const d of [
      concat(1,          [a,b,c]),
      concat(1,null,     [a,b,c]),
      concat(1,undefined,[a,b,c]),
      concat(1,'int32',  [a,b,c])
    ])
    {
      expect(a.dtype).toBe('int32')
      expect(b.dtype).toBe('int32')
      expect(c.dtype).toBe('int32')
      expect(d.dtype).toBe('int32')
  
      expect(d.shape).toEqual( Int32Array.of(2,6,3) )
      expect(d.data ).toEqual(
        Int32Array.of(
          111,112,113,
          121,122,123,
          131,132,133,
          141,142,143,
          151,152,153,
          161,162,163,
          211,212,213,
          221,222,223,
          231,232,233,
          241,242,243,
          251,252,253,
          261,262,263
        )
      )
    }
  })

  it('works on arrays of shapes [2,3], [1,3] and [3,3] along axis 0', () => {
    const a = array([
            [11,12,13],
            [21,22,23]
          ]),
          b = array([
            [31,32,33]
          ]),
          c = array([
            [41,42,43],
            [51,52,53],
            [61,62,63],
          ])

    for( const d of [
      concat([a,b,c]),
      
      concat(undefined,[a,b,c]),
      concat(null,[a,b,c]),
      concat(0,[a,b,c]),
      concat('int32',[a,b,c]),

      concat(undefined, undefined, [a,b,c]),
      concat(undefined, null,      [a,b,c]),
      concat(undefined, 'int32',   [a,b,c]),
      concat(null,      undefined, [a,b,c]),
      concat(null,      null,      [a,b,c]),
      concat(null,      'int32',   [a,b,c]),
      concat(0,         undefined, [a,b,c]),
      concat(0,         null,      [a,b,c]),
      concat(0,         'int32',   [a,b,c])
    ])
    {
      expect(a.dtype).toBe('int32')
      expect(b.dtype).toBe('int32')
      expect(c.dtype).toBe('int32')
      expect(d.dtype).toBe('int32')
  
      expect(d.shape).toEqual( Int32Array.of(6,3) )
      expect(d.data ).toEqual(
        Int32Array.of(
          11,12,13,
          21,22,23,
          31,32,33,
          41,42,43,
          51,52,53,
          61,62,63
        )
      )
    }
  })

  it('works on arrays of shapes [3,2], [3,1] and [3,3] along axis 1', () => {
    const a = array([
            [11,12],
            [21,22],
            [31,32]
          ]),
          b = array([
            [13],
            [23],
            [33]
          ]),
          c = array([
            [14,15,16],
            [24,25,26],
            [34,35,36],
          ])

    for( const d of [
      concat(1,          [a,b,c]),
      concat(1,null,     [a,b,c]),
      concat(1,undefined,[a,b,c]),
      concat(1,'int32',  [a,b,c])
    ])
    {
      expect(a.dtype).toBe('int32')
      expect(b.dtype).toBe('int32')
      expect(c.dtype).toBe('int32')
      expect(d.dtype).toBe('int32')
  
      expect(d.shape).toEqual( Int32Array.of(3,6) )
      expect(d.data ).toEqual(
        Int32Array.of(
          11,12,13,14,15,16,
          21,22,23,24,25,26,
          31,32,33,34,35,36
        )
      )
    }
  })
})
