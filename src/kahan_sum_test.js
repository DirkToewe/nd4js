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

import {KahanSum} from './kahan_sum'


describe('KahanSum', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  forEachItemIn(Float64Array.of(
    1e6, 3e6, 6e6, 7e6, 8e6, 9e6,
    1e5, 3e5, 6e5, 7e5, 8e5, 9e5,
    1e7, 3e7, 6e7, 7e7, 8e7, 9e7
  )).it(`n times 1/n sums up to 1.`, N => {
    const kahan = new KahanSum();

    for( let repeat=3; repeat-- > 0; )
    {
      for( let i=N; i-- > 0; )
        kahan.add(1/N);

      expect(kahan.sum).toBeAllCloseTo(1, {rtol: 0, atol: Number.EPSILON});

      kahan.set(0);
    }
  })
})
