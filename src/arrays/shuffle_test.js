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

import {forEachItemIn} from '../jasmine_utils'
import {shuffle} from './shuffle'


describe('shuffle', () => {

  const chi_square_table = [
    10.828,
    13.816,
    16.266,
    18.467,
    20.515,
    22.458,
    24.322,
    26.124,
    27.877,
    29.588,
    31.264,
    32.909,
    34.528,
    36.123,
    37.697,
    39.252,
    40.790,
    42.312,
    43.820,
    45.315,
    46.797,
    48.268,
    49.728,
    51.179,
    52.620,
    54.052,
    55.476,
    56.892,
    58.301,
    59.703,
    61.098,
    62.487,
    63.870,
    65.247,
    66.619,
    67.985,
    69.346,
    70.703,
    72.055,
    73.402,
    74.745,
    76.084,
    77.419,
    78.750,
    80.077,
    81.400,
    82.720,
    84.037,
    85.351
  ];

  forEachItemIn(
    function*(){
      for( let N=1; N++ < 6; ) {
        const array = Int32Array.from({length: N}, (_,i) => i);
        Object.freeze(array.buffer);
        yield array;
      }
    }()
  ).it('passes chi-squared test', array => {
    const N = array.length,
      freqs = Array.from(array, () => array.map(()=>0) );

    const M = 16*1024*1024;

    for( let repeat=0; repeat++ < M; )
    {
      const arr = array.slice();
      shuffle(arr);
      for( let i=N; i-- > 0; )
        freqs[i][arr[i]] += 1;
    }

    const E = M/N;

    for( const row of freqs )
    {
      let chi = 0;
      for( const n of row )
        chi += (n-E)**2 / E;
      expect(chi).toBeLessThan(chi_square_table[N-2]);
    }
  });

});
