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

import {forEachItemIn, CUSTOM_MATCHERS} from '../../jasmine_utils'
import {rosenbrock, rosenbrock_grad} from '../test_fn/rosenbrock'
import {array, NDArray} from '../../nd_array'
import {zip_elems} from '../../zip_elems'
import math from '../../math'
import {strong_wolfe} from './strong_wolfe'
import {LineSearchNoProgressError} from './line_search_error'


describe('strong_wolfe', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  forEachItemIn(
    function*(){
      const  S = 3.14, Δ = 0.42

      for( let w = -S; w <= +S; w+=Δ )
      for( let x = -S; x <= +S; x+=Δ ) { yield [w,x]
      for( let y = -S; y <= +S; y+=Δ ) { yield [w,x,y]
      for( let z = -S; z <= +S; z+=Δ ) { yield [w,x,y,z] }}}
    }()
  ).it('works on rosenbrock', x0 => {
    let nCalls = 0
    const fg = x => {
      ++nCalls
      return [
        rosenbrock(x),
        rosenbrock_grad(x)
      ]
    }

    const [f0,g0] = fg(x0), negDir = g0.mapElems('float64', x => x*0.1),
       [c1,c2,c3] = [0.4,0.8,1.6]

    let x,f,g;
    try {
      ( [x,f,g] = strong_wolfe({c1,c2,c3})(fg)(x0,f0,g0, negDir) );
    }
    catch(err) {
      if( err instanceof LineSearchNoProgressError )
        return
      throw err
    }

    expect(x).toEqual(jasmine.any(NDArray))
    expect(f).toEqual(jasmine.any(NDArray))
    expect(g).toEqual(jasmine.any(NDArray))

    expect( rosenbrock_grad(x) ).toBeAllCloseTo(g, {atol:0, gtol:0})

    expect(x.shape).toEqual( Int32Array.of(x0.length) )
    expect(g.shape).toEqual( Int32Array.of(x0.length) )

    expect(nCalls).toBeLessThan(16)

    const [p,p0] = [g,g0].map( g => {
      let pg = 0
      for( let [i]=negDir.shape; i-- > 0; )
        pg -= negDir.data[i] * g.data[i]
      return pg
    })

    expect( Math.abs(p) ).not.toBeGreaterThan( -c2*p0 )

    let α = zip_elems([x,x0], math.sub).data.reduce(math.hypot) / negDir.data.reduce(math.hypot)
    expect( f - f0 ).not.toBeGreaterThan( c1*α*p0 )
  })
})
