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

import {tabulate} from './tabulate'
import {forEachItemIn} from './jasmine_utils'


describe('tabulate', () => {
  forEachItemIn(
    function*(){
      yield []
      for( let i=1; i < 10; i++ ) { yield [i]
      for( let j=1; j < 10; j++ ) { yield [i,j]
      for( let k=1; k < 10; k++ ) { yield [i,j,k]
      for( let l=1; l < 10; l++ ) { yield [i,j,k,l] }}}}
    }()
  ).it('works for generated examples', shape => {
    shape = Int32Array.from(shape)
    const arr = tabulate(shape, (...indices) => indices.reduce((a,b) => 10*a+b, 0) )
    expect(arr.shape).toEqual(shape)

    function test(d,...indices) {
      if( d == shape.length )
        expect( arr(...indices) ).toBe( indices.reduce((a,b) => 10*a+b, 0) )
      else for( let i=0; i < shape[d]; i++ )
        test(d+1, ...indices, i)
    }

    test(0)
  })
})
