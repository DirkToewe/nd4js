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

import {CUSTOM_MATCHERS, forEachItemIn} from '../jasmine_utils'
import {_min1d_interp_quad,
        _heap_sort } from './_opt_utils'


describe('_heap_sort', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  for( const [suffix, ...args] of [
    [''],
    [' with comparator', (x,y) => y-x]
  ])
    forEachItemIn(
      function*(){
        for( const ArrayType of [Int32Array,Float64Array] )
          for( let length=0; length++ < 512; ) {
            const         items = ArrayType.from({length}, () => Math.random()*8192 - 4096);
            Object.freeze(items.buffer);
            yield         items;
          }
      }()
    ).it('works on random examples' + suffix, unsorted => {
      const items = unsorted.slice(),
            order = unsorted.slice();

      Object.freeze(order.buffer);
      expect(items.buffer).not.toBe(order.buffer);

      if( args.length > 0 )
        order.sort(args[0]);
      else
        order.sort( (x,y) => x-y );

      const seq = _heap_sort( items, ...args.map( f => (x,y) => f(x,y) <= 0 ) );

      let n=0;
      for( const x of seq ) {
        expect(x).toEqual(items[n++]);
        expect( items.subarray(0,n) ).toEqual( order.subarray(0,n) );
      }
      expect(n).toBe(unsorted.length);

      expect(items).toEqual(order);
    });
})


describe('_min1d_interp_quad', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  forEachItemIn(
    function*(){
      const N = 4,
            M = 2;

      for( let a = -M; a <= +M; a += 1/N )
      for( let b =1/N; b <= +M; b += 1/N )
      for( let z = -M; z <= +M; z += 1/N )
      {
        const f = x => a + b*(x-z)*(x-z),
              g = x =>   2*b*(x-z);

        for( let x1 = -M; x1 <= +M; x1 += 1/N )
        for( let x2 = -M; x2 <= +M; x2 += 1/N )
          if( x1 !== x2 )
            yield [ [x1,x2, f(x1),f(x2), g(x1)], z ];
      }
    }()
  ).it('works on generated examples.', ([args, xMin]) => {
    expect( _min1d_interp_quad(...args) ).toBeAllCloseTo(xMin);
  });

  forEachItemIn(
    function*(){
      for( let i=512; i-- > 0; )
      {
        const a = Math.random()*8 - 4,
              b = Math.random()*8 + 1/4,
              z = Math.random()*8 - 4;

        const f = x => a + b*(x-z)*(x-z),
              g = x =>   2*b*(x-z);

        for( let j=1024; j-- > 0; )
        {
          const x1 =  Math.random()*8 - 4,
           x2 = x1 + (Math.random()*4 + 1e-4) * (Math.random() < 0.5 ? +1 : -1);
          yield [ [x1,x2, f(x1),f(x2), g(x1)], z ];
        }
      }
    }()
  ).it('works on random examples.', ([args, xMin]) => {
    expect( _min1d_interp_quad(...args) ).toBeAllCloseTo(xMin, {atol: 1e-4});
  });
})
