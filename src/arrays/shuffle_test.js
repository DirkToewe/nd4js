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

import {shuffle} from './shuffle'


describe('shuffle', () => {

  // χ²-table; α = 0.001;
  const chiq = [
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
    62.487
  ];


  for( let N=1; N++ < 8; )
    it(`passes chi-squared test shuffling an Int32Array of length ${N}`, () => {
      const array =       Int32Array.from({length: N}, (_,i) => i),
            freqs = new Float64Array(N*N);

      const M = 17*1024*1024;

      for( let repeat=0; repeat++ < M; )
      {
        const arr = array.slice();
        shuffle(arr);
        for( let i=N; i-- > 0; )
          freqs[N*i + arr[i]] += 1;
      }

      const E = M/N;

      for( let i=N; i-- > 0; )
      {
        let chi = 0;
        for( let j=N; j-- > 0; ) {
          const   d = freqs[N*i+j] - E;
          chi += d*d / E;
        }
        expect(chi).toBeLessThan(chiq[N-2]);
      }
    });

});
