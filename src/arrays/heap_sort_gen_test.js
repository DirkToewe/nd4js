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
import {heap_sort_gen} from './heap_sort_gen'


describe('heap_sort_gen', () => {

  for( const ArrayType of [Int32Array,Float64Array] )
  for( const [suffix, ...args] of [
    [''],
    [' with descending comparator', (x,y) => y-x],
    [' with ascending comparator',  (x,y) => x-y]
  ])
    forEachItemIn(
      function*(){
        for( let length=0; length++ < 137; ) {
          const         items = ArrayType.from({length}, () => Math.random()*512 - 256 );
          Object.freeze(items.buffer);
          yield         items;
        }
      }()
    ).it(`works given random ${ArrayType.name} examples ${suffix}`, unsorted => {
      const items = unsorted.slice(),
            order = unsorted.slice();

      Object.freeze(order.buffer);
      expect(items.buffer).not.toBe(order.buffer);

      order.sort( args[0] || ((x,y) => x-y) );

      const seq = heap_sort_gen(items, ...args);

      let n=0;
      for( const x of seq ) {
        expect(x).toBe(items[n++]);
        expect( items.subarray(0,n) ).toEqual( order.subarray(0,n) );
      }
      expect(n).toBe(unsorted.length);

      expect(items).toEqual(order);
    });

});
