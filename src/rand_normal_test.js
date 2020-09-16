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

import {forEachItemIn, CUSTOM_MATCHERS} from './jasmine_utils'

import {rand_normal} from './rand_normal'


describe('nd.rand_normal', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  it('Passes Jarque–Bera test', () => {
    // https://de.wikipedia.org/wiki/Jarque-Bera-Test
    // I'm no statistician, so ... please don't laugh if this test is implemented incorrectly :P
    const length = 4*1024*1024;

    const vals = Float64Array.from({length}, () => rand_normal());

    const E = vals.reduce((E,x) => E +  x      , 0) / length,
          s = vals.reduce((s,x) => s + (x-E)**2, 0) / length,
          S = vals.reduce((S,x) => S + (x-E)**3, 0) / length / s**1.5,
          C = vals.reduce((C,x) => C + (x-E)**4, 0) / length / s**2;

    const JB = length/6 * (S*S + (C-3)*(C-3)/4);

    // expect(JB).toBeLessThan(9.2); // α ≈ 0.01
    expect(JB).toBeLessThan(10.597); // α ≈ 0.005
  });

  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from;

      for( let run=16; run-- > 0; )
        yield randInt(4,32)*1024*1024;
    }()
  ).it('Has σ=1 and E=0', len => {
    let E = 0,
        s = 0;

    for( let i=0; i < len; i++ )
    {
      const x = rand_normal();
      E += x;
      s += x*x;
    }

    E = E/len;
    s = s/len - E*E; 

    expect(E).toBeCloseTo(0, {atol: 1e-2, rtol: 0});
    expect(s).toBeCloseTo(1, {atol: 1e-2, rtol: 0});
  })
})

