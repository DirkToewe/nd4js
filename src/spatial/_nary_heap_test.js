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

import {NAryHeap} from "./_nary_heap";


describe('arrays.compare', () => {
  forEachItemIn(
    function*(){
      for( let run=0; run++ < 7*1337; )
      {
        const length = _rand_int(1,1337);
        yield Object.freeze(
          Array.from( {length}, () => _rand_int(-1337,+1337) )
        );
      }
    }()
  ).it('works correctly given random examples', arr => {
    const heap = new NAryHeap();

    const pool = arr.slice();
    while(pool.length > 0)
    {
      const [key] = pool.splice( _rand_int(0,pool.length), 1 );
      heap.add({key});

      if( Math.random() < 0.1 ) {
        const    {key} = heap.popMin();
        pool.push(key);
      }
    }

    expect(heap.size).toBe(arr.length);

    const sorted = [];
    while( 0 < heap.size )
      sorted.push( heap.popMin().key );
    Object.freeze(sorted);

    const SORTED = Object.freeze( arr.slice().sort( (x,y) => x-y ) );

    expect(sorted).toEqual(SORTED);
  })
})