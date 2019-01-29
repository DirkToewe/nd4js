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
import {nd_to_istr,
        istr_to_nd} from './istr'
import {tabulate} from '../tabulate'
import {WHITESPACES} from '.'


describe('istr', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  const randInt = until => Math.trunc(Math.random()*until)


  forEachItemIn(
    function*(){
      function* shapes() {
        yield[]
        for( let i=1; i <= 6; i++ ) { yield [i]
        for( let j=1; j <= 6; j++ ) { yield [i,j]
        for( let k=1; k <= 6; k++ ) { yield [i,j,k]
        for( let l=1; l <= 6; l++ ) { yield [i,j,k,l] }}}}
      }

      for( const dtype of ['int32', 'float32', 'float64', 'complex128'] )
      for( const shape of shapes() )
      {
        const A = tabulate(shape, dtype, () => Math.random()*2048 - 1024)
        Object.freeze(A.data.buffer)
        yield A
      }
    }()
  ).it('istr_to_nd(nd_to_istr(...)) results in same array for random examples', A => {
    const a_istr = Array.from( nd_to_istr(A) ),
          a=istr_to_nd(a_istr)

    expect(a.dtype).toBe(A.dtype)
    expect(a.shape).toEqual(A.shape)
    expect(a).toBeAllCloseTo(A, {rtol:0, atol:0})
  })
})
