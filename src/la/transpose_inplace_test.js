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
import {transpose_inplace} from './transpose_inplace'
import {tabulate} from '../tabulate'


describe('transpose_inplace', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  forEachItemIn(
    function*(){
      function* shapes() {
        for( let i=1; i < 7; i++ ) { yield [i,i]
        for( let j=1; j < 7; j++ ) { yield [i,j,j]
        for( let k=1; k < 7; k++ ) { yield [i,j,k,k]
        for( let l=1; l < 7; l++ ) { yield [i,j,k,l,l] }}}}
      }

      for( const shape of shapes() )
      {
        const A = tabulate( shape, 'int32', (...idx) => idx.reduce((flat,s) => 10*flat + s+1, 0) ),
           [M,N]= shape.slice(-2)
        yield A
      }
    }()
  ).it('works on generated examples', A => {
    const A_T = A.T
    transpose_inplace(A)
    expect(A).toBeAllCloseTo(A_T, {rtol:0,atol:0});
  })
})
