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
import {tabulate} from '../tabulate'
import {matmul, matmul2} from './matmul'
import {bidiag_decomp} from './bidiag'
import {eye} from './eye'


describe('bidiag', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

      for( let run=1024; run-- > 0; )
      {
        const ndim = randInt(2,5),
             shape = Int32Array.from({ length: ndim }, () => randInt(1,24) )
        yield tabulate(shape, 'float64', () => Math.random()*2 - 1 )
      }
    }()
  ).it('bidiag_decomp works on random examples', A => {
    const [N,M] = A.shape.slice(-2),
        [U,B,V] = bidiag_decomp(A)
  
    const  a = matmul(U,B,V),
      absMax = (x,y) => math.max(
        math.abs(x),
        math.abs(y)
      )

    expect(B).toBeUpperBidiagonal()

    if( N >= M ) {
      const I = eye(M)
      expect( matmul2(U.T,U) ).toBeAllCloseTo(I)
      expect( matmul2(V.T,V) ).toBeAllCloseTo(I)
      expect( matmul2(V,V.T) ).toBeAllCloseTo(I)
    }
    else {
      const I = eye(N)
      expect( matmul2(U.T,U) ).toBeAllCloseTo(I)
      expect( matmul2(U,U.T) ).toBeAllCloseTo(I)
      expect( matmul2(V,V.T) ).toBeAllCloseTo( eye(N+1) )
    }

    expect(a).toBeAllCloseTo(A)
  })
})
