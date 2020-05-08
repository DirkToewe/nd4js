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
import {_rand_int} from '../_test_data_generators'

import {compare} from './comparator'


describe('arrays.compare', () => {

  const rand_char = () => {
    const  POOL = '0123456789ABCDEFGHIJKLMNOPQRSTUVabcdefghijklmnopqrstuv';
    return POOL[_rand_int(0,POOL.length)];
  };

  forEachItemIn(
    function*(){
      for( let run=0; run++ < 1337*1337; )
      {
        const length = _rand_int(0,32);

        const a = Array.from({length}, rand_char),
              b = Array.from({length}, rand_char);

        yield Object.freeze([a,b]);
      }
    }()
  ).it('works given random character arrays', ([a,b]) => {
    const A = a.join(''),
          B = b.join('');

    expect( compare(a,b) ).toBe( (a > b) - (a < b) );
  })

})
