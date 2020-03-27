'use strict';

/* This file is part of ND4JS.
 *
 * ND4JS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ND4JS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ND4JS. If not, see <http://www.gnu.org/licenses/>.
 */

import {forEachItemIn, CUSTOM_MATCHERS} from '../jasmine_utils'
import {NDArray} from '../nd_array'

import {min_lbfgs_gen} from './lbfgs'

import {LineSearchNoProgressError} from './line_search/line_search_error'

import {rosenbrock, rosenbrock_grad} from './test_fn/rosenbrock'


describe('lbfgs', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

//*DEBUG*/  const samples = [];

  forEachItemIn(
    function*(){                     const n = 16;
      function*       range() { for( let i=n+1; i-- > 0; ) yield Math.PI*(1-2*(i/n)); }
      for( const x of range() )
      for( const y of range() ) { yield [x,y];
      for( const z of range() ) { yield [x,y,z]; }}

//*DEBUG*/      const avg = samples.reduce((x,y) => x+y) / samples.length,
//*DEBUG*/            std = Math.hypot( ...samples.map( x => (x-avg) / Math.sqrt(samples.length) ) );
//*DEBUG*/
//*DEBUG*/      console.log('L-BFGS')
//*DEBUG*/      console.log('------')
//*DEBUG*/      console.log('MIN:', samples.reduce((x,y) => Math.min(x,y)) );
//*DEBUG*/      console.log('MAX:', samples.reduce((x,y) => Math.max(x,y)) );
//*DEBUG*/      console.log('AVG:', avg );
//*DEBUG*/      console.log('STD:', std );
    }()
  ).it('min_lbfgs_gen works on rosenbrock', x0 => {
    let nCalls = 0
    const fg = x => {
      expect(++nCalls).toBeLessThan(96)
      return [
        rosenbrock(x),
        rosenbrock_grad(x)
      ]
    }

    let   x,f,g, nIter = -1
    for( [x,f,g] of min_lbfgs_gen(fg, x0) )
    {
      expect(x).toEqual( jasmine.any(NDArray) )
      expect(f).toEqual( jasmine.any(Number ) )
      expect(g).toEqual( jasmine.any(NDArray) )

      expect(x.ndim).toBe(1)
      expect(g.ndim).toBe(1)

      expect(x.shape).toEqual( Int32Array.of(x0.length) )
      expect(g.shape).toEqual( Int32Array.of(x0.length) )

      expect(f).toBeAllCloseTo(rosenbrock     (x), {rtol:0, atol:0})
      expect(g).toBeAllCloseTo(rosenbrock_grad(x), {rtol:0, atol:0})

      const gNorm = Math.hypot(...g.data)
      if(   gNorm <= 1e-8 )
        break
      expect(++nIter).toBeLessThan(128)
    }

//*DEBUG*/    samples.push(nCalls);

    expect(x).toBeAllCloseTo(1)
    expect(f).toBeAllCloseTo(0)
    expect(g).toBeAllCloseTo(0)
  })
})
