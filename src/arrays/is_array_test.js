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


import {forEachItemIn} from '../jasmine_utils'
import {_rand_int} from '../_test_data_generators'

import {Complex128Array} from "../dt/complex_array";

import {is_array} from './is_array'


describe('arrays.compare', () => {

  forEachItemIn([
                Array,
         Float32Array,
         Float64Array,
            Int8Array,
           Int16Array,
           Int32Array,
           Uint8Array,
          Uint16Array,
          Uint32Array,
    Uint8ClampedArray,
      Complex128Array
  ]).it('returns true given Arrays and typed arrays', TypedArray => {
    expect( is_array(TypedArray.of(1,2,3)) ).toBe(true);
  })

  it('returns false for other objects', () => {
    expect( is_array({}) ).toBe(false);
    expect( is_array({length: 13}) ).toBe(false);
  })

})
