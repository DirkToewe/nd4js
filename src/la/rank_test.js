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
import {matmul2} from './matmul'
import {qr_decomp} from './qr'
import {rank} from './rank'


describe('rank', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

      for( let run=1024; run-- > 0; )
      {
        let M = randInt(1,32),
            N = randInt(1,32)
        if( N > M ) [M,N] = [N,M]
        const L = Math.min(M,N),
          ndim = randInt(0,7)
        let r_shape = Array.from({ length: ndim-2 }, () => randInt(1,8) ),
            A_shape = r_shape.concat([M,N])
        r_shape = Int32Array.from(r_shape)
        A_shape = Int32Array.from(A_shape)

        const A = tabulate(A_shape, 'float64', () => Math.random() < 0.1 ? 0 : Math.random()*2-1)
        Object.freeze(A.data.buffer)
        const [Q,_]= qr_decomp(A),
               r = tabulate(    r_shape, 'int32', () => Math.random()*L),
               R = tabulate([...r_shape,L,N], 'float64', (...idx) => {
                 const j = idx.pop(),
                       i = idx.pop(), R = r(...idx)
                 if(R <= i || j < i) return 0
                 if(j == i) return (0.5 + Math.random()) * (randInt(0,2)*2 - 1)
                 return Math.random()*0.5 - 0.25
               })
        yield [r,matmul2(Q,R)]
      }
    }()
  ).it('works on random examples', ([R,A]) => {
    const r = rank(A)

    expect(r.shape).toEqual(A.shape.slice(0,-2))
    expect(R.shape).toEqual(A.shape.slice(0,-2))
    expect(r.dtype).toBe('int32')
    expect(R.dtype).toBe('int32')
    expect(r).toBeAllCloseTo(R, {rtol:0, atol:0})
  })
})
