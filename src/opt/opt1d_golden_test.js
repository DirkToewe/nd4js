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

import {opt1d_golden} from './opt1d_golden'
import {CUSTOM_MATCHERS, forEachItemIn} from '../jasmine_utils'


describe('opt1d_golden', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  forEachItemIn(
    function*(){
      for( let xx=0; xx <= 8; xx += 1/1024 )
        yield xx
    }()
  ).it('computes square roots correctly', xx => {
    const off = Math.random()*0.2 - 0.1;

    const f = x => off + Math.abs(x*x - xx),
          x = opt1d_golden( f, 0, Math.max(1,xx) )
    expect(x).toBeAllCloseTo( Math.sqrt(xx));
  })

  forEachItemIn(
    function*(){
      for( let xxx=0; xxx <= 8; xxx += 1/1024 )
        yield xxx
    }()
  ).it('computes cube roots correctly', xxx => {
    const off = Math.random()*0.2 - 0.1;

    const f = x => off + Math.abs(x*x*x - xxx),
          x = opt1d_golden( f, 0, Math.max(1,xxx) )
    expect(x).toBeAllCloseTo( xxx**(1/3), {atol:1e-4} );
  })
})
