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
import {rosenbrock, rosenbrock_grad} from './test_fn/rosenbrock'
import {LineSearchNoProgressError} from './line_search/line_search_error'
import {min_lbfgs_gen} from './lbfgs'
import math from '../math'
import {NDArray} from '../nd_array'


describe('lbfgs', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  forEachItemIn(
    function*(){
      const  S = 3.14, Δ = 0.42

      for( let w = -S; w <= +S; w+=Δ )
      for( let x = -S; x <= +S; x+=Δ ) { yield [w,x]
      for( let y = -S; y <= +S; y+=Δ ) { yield [w,x,y] }}
    }()
  ).it('min_lbfgs_gen works on rosenbrock', x0 => {
    let nCalls = 0
    const fg = x => {
      expect(++nCalls).toBeLessThan(256)
      return [
        rosenbrock(x),
        rosenbrock_grad(x)
      ]
    }

    const opt = {
      negDir0: g => g.mapElems('float64',g=>g/1024)
    }

    let x,f,g, nIter = -1
    min_lbfgs_gen(fg,x0)
    try {
      for( [x,f,g] of min_lbfgs_gen(fg, x0, opt) )
      {
        expect(x).toEqual( jasmine.any(NDArray) )
        expect(f).toEqual( jasmine.any(NDArray) )
        expect(g).toEqual( jasmine.any(NDArray) )

        expect(x.ndim).toBe(1)
        expect(f.ndim).toBe(0)
        expect(g.ndim).toBe(1)

        expect(x.shape).toEqual( Int32Array.of(x0.length) )
        expect(g.shape).toEqual( Int32Array.of(x0.length) )

        const gNorm = g.reduceElems(math.hypot)
        if( gNorm <= 1e-8 )
          break
        expect(++nIter).toBeLessThan(128)
      }
    }
    catch(err) {
      if( ! (err instanceof LineSearchNoProgressError) )
        throw err
      console.log('NO_PROGRESS')
    }

//    console.log('nCalls, nIter:', nCalls, nIter)
    expect(x).toBeAllCloseTo(1)
    expect(f).toBeAllCloseTo(0)
    expect(g).toBeAllCloseTo(0)
  })
})
