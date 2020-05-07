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


import {forEachItemIn, CUSTOM_MATCHERS} from '../../jasmine_utils'
import {cartesian_prod, linspace} from '../../iter'

import {generic_test_test_fn} from './_generic_test_test_fn'
import {powell_badscale} from './powell_badscale'


describe(`test_fn.${powell_badscale.name}`, () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  forEachItemIn(
    cartesian_prod(
      linspace(-24, +24, 513),
      linspace(-24, +24, 513)
    )
  ).it('works for generated examples', ([x,y]) => {
    const f = powell_badscale([x,y]),
          f1 = 1e4*x*y - 1,
          f2 = Math.exp(-x) + Math.exp(-y) - 1.0001;

    expect(f).toBeAllCloseTo(f1*f1 + f2*f2);
  })


  forEachItemIn(
    function*(){
      for( let run=0; run++ < 173*1337; )
        yield [
          Math.random()*64 - 32,
          Math.random()*64 - 32
        ];
    }()
  ).it('works for generated examples', ([x,y]) => {
    const f = powell_badscale([x,y]),
          f1 = 1e4*x*y - 1,
          f2 = Math.exp(-x) + Math.exp(-y) - 1.0001;

    expect(f).toBeAllCloseTo(f1*f1 + f2*f2);
  })
})


generic_test_test_fn(powell_badscale, [[-8,+8],
                                       [-8,+8]]);
