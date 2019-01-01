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

import {root1d_bisect} from './root1d_bisect'
import {forEachItemIn} from '../jasmine_utils'


describe('root1d_bisect', () => {
  const abs = Math.abs.bind(Math)

  forEachItemIn(
    function*(){
      for( let xx=0; xx <= 8; xx += 1/1024 )
        yield xx
    }()
  ).it('computes square roots correctly', xx => {
    const f = x => x*x - xx,
          x = root1d_bisect( f, 0, Math.max(1,xx) )
    expect( abs(f(x)) ).not.toBeGreaterThan( abs(f(xx**0.5)) )
  })

  forEachItemIn(
    function*(){
      for( let xxx=0; xxx <= 8; xxx += 1/1024 )
        yield xxx
    }()
  ).it('computes cube roots correctly', xxx => {
    const f = x => x*x*x - xxx,
          x = root1d_bisect( f, 0, Math.max(1,xxx) )
    expect( abs(f(x)) ).not.toBeGreaterThan( abs(f(xxx**(1/3))) )
  })
})
