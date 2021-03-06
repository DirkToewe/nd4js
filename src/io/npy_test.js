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
import {b64_decode_gen} from './b64'
import {npy_serialize_gen,
        npy_deserialize} from './npy'
import {npy_test_data} from './npy_test_data'
import {tabulate} from '../tabulate'


describe('npy', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  forEachItemIn(
    npy_test_data()
  ).it('npy_deserialize works for random examples', ([b64, reference]) => {
    const result = npy_deserialize( b64_decode_gen(b64) );

    expect(result.dtype).toBe   (reference.dtype)
    expect(result.shape).toEqual(reference.shape)

    if( result.dtype.startsWith('complex') ) expect(result.data._array).toEqual(reference.data._array)
    else                                     expect(result.data       ).toEqual(reference.data       )
  })


  for( const dtype of ['int32','float32','float64'] )
    forEachItemIn(
      function*(){
        const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

        for( let run=173; run-- > 0; )
        {
          const shape = Int32Array.from({length: randInt(0,5)}, () => randInt(1,24)),
                A = tabulate(shape, dtype, () => Math.random() < 0.1 ? 0 : Math.random()*2e3-1e3)
          Object.freeze(A.data.buffer)
          yield A
        }
      }()
    ).it(`npy_serialize works for random ${dtype.padStart(7)} examples`, (A) => {
      const B = npy_deserialize(npy_serialize_gen(A))

      expect(B.dtype).toBe   (A.dtype)
      expect(B.shape).toEqual(A.shape)

      if( A.dtype.startsWith('complex') ) expect(B.data._array).toEqual(A.data._array)
      else                                expect(B.data       ).toEqual(A.data       )
    })
})
