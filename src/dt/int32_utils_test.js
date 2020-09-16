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

import {forEachItemIn} from "../jasmine_utils";

import {bitCount} from "./int32_utils";


describe('int32_utils', () => {

  const test_body = bits => {
    let N_BITS = 0;
    for( let i=0; i < 32; i++ )
      N_BITS += bits>>>i & 1;

    const nBits = bitCount(bits);

    if( nBits !== N_BITS )
      expect(nBits).toBe(N_BITS);
  };


  forEachItemIn(
    function*(){
      const N = 1 << 19,
          b31 = 2**31;

      for( let bits = -b31  ; bits < N-b31; bits++ ) yield bits;
      for( let bits =     -N; bits < N    ; bits++ ) yield bits;
      for( let bits =  b31-N; bits <   b31; bits++ ) yield bits;
    }()
  ).it('bitCount counts propery given generated examples', test_body);


  forEachItemIn(
    function*(){
      for( let run=0; run++ < (1<<21); )
        yield Math.floor(Math.random() * 2**32) - 2**31;
    }()
  ).it('bitCount counts propery given random examples #1', test_body);


  forEachItemIn(
    function*(){
      for( let run=0; run++ < (1<<20); )
      {
        let bits = 0;
        for( let i=32; i-- > 0; )
          bits = bits << 1 | Math.random() < 0.5;
        yield bits;
      }
    }()
  ).it('bitCount counts propery given random examples #2', test_body);

});
