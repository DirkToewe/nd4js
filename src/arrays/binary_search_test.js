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
import {binary_rangesearch,
        binary_search} from './binary_search'
import {_rand_int,
        _shuffle} from '../_test_data_generators'


describe('binary_search', () => {

  forEachItemIn(
    function*(){
      for( let run=0; run < 128; run++ )
      for( let len=0; len < 128; len++ )
      {
        const array = [];
        let x = _rand_int(-733,+733);
        for( let i=len; i-- > 0; )
          array.push( x += _rand_int(0,4) );
        yield Object.freeze(array);
      }
    }()
  ).it('works for random examples of integer arrays', array => {
    let [min,max] = [0,0];
    if( array.toBeLessThan !== 0 )
      [min,max] = array.reduce(
        ([min,max],x) => [Math.min(min,x),
                          Math.max(max,x)],
        [Infinity,-Infinity]
      );

    for( let key=min-2; key < max+2; key++ )
    {
      let i = binary_search(array, key);
      if( i < 0 )
      {   i = ~i;
        expect(i).toBeLessThanOrEqual(array.length);
        expect(i).toBeGreaterThanOrEqual(0);
        if( i > 0            ) expect(key).toBeGreaterThan(array[i-1]);
        if( i < array.length ) expect(key).toBeLessThan   (array[i  ]);
      }
      else {
        expect(i).toBeLessThan(array.length);
        expect(i).toBeGreaterThanOrEqual(0);
        expect(array[i]).toBe(key);
      }
    }
  })


  forEachItemIn(
    function*(){
      for( let run=0; run < 128; run++ )
      for( let len=0; len < 128; len++ )
      {
        const array = [];
        let x = _rand_int(-733,+733);
        for( let i=len; i-- > 0; )
          array.push( x -= _rand_int(0,4) );
        yield Object.freeze(array);
      }
    }()
  ).it('works for random examples of integer arrays with descending order', array => {
    let [min,max] = [0,0];
    if( array.toBeLessThan !== 0 )
      [min,max] = array.reduce(
        ([min,max],x) => [Math.min(min,x),
                          Math.max(max,x)],
        [Infinity,-Infinity]
      );

    for( let key=min-2; key < max+2; key++ )
    {
      let i = binary_search( array, key, (x,y) => (x<y) - (x>y) );
      if( i < 0 )
      {   i = ~i;
        expect(i).toBeLessThanOrEqual(array.length);
        expect(i).toBeGreaterThanOrEqual(0);
        if( i > 0            ) expect(key).toBeLessThan   (array[i-1]);
        if( i < array.length ) expect(key).toBeGreaterThan(array[i  ]);
      }
      else {
        expect(i).toBeLessThan(array.length);
        expect(i).toBeGreaterThanOrEqual(0);
        expect(array[i]).toBe(key);
      }
    }
  })

});


describe('binary_rangesearch', () => {

  forEachItemIn(
    function*(){
      for( let run=0; run < 73; run++ )
      for( let len=0; len < 73; len++ )
      {
        const array = [];
        let x = _rand_int(-337,+337);
        for( let i=len; i-- > 0; )
          array.push( x += _rand_int(0,4) );
        yield Object.freeze(array);
      }
    }()
  ).it('works for random examples of integer arrays', array => {
    let [min,max] = [0,0];
    if( array.toBeLessThan !== 0 )
      [min,max] = array.reduce(
        ([min,max],x) => [Math.min(min,x),
                          Math.max(max,x)],
        [Infinity,-Infinity]
      );

    for( let key=min-2; key < max+2; key++ )
    {
      let i = binary_rangesearch( 0, array.length, i => {
        expect(i).toBeLessThan(array.length);
        expect(i).toBeGreaterThanOrEqual(0);
        return (array[i] > key) - (array[i] < key);
      });
      if( i < 0 )
      {   i = ~i;
        expect(i).toBeLessThanOrEqual(array.length);
        expect(i).toBeGreaterThanOrEqual(0);
        if( i > 0            ) expect(key).toBeGreaterThan(array[i-1]);
        if( i < array.length ) expect(key).toBeLessThan   (array[i  ]);
      }
      else {
        expect(i).toBeLessThan(array.length);
        expect(i).toBeGreaterThanOrEqual(0);
        expect(array[i]).toBe(key);
      }
    }
  })


  forEachItemIn(
    function*(){
      for( let run=0; run < 73; run++ )
      for( let len=0; len < 73; len++ )
      {
        const array = [];
        let x = _rand_int(-337,+337);
        for( let i=len; i-- > 0; )
          array.push( x -= _rand_int(0,4) );
        yield Object.freeze(array);
      }
    }()
  ).it('works for random examples of integer arrays with descending order', array => {
    let [min,max] = [0,0];
    if( array.toBeLessThan !== 0 )
      [min,max] = array.reduce(
        ([min,max],x) => [Math.min(min,x),
                          Math.max(max,x)],
        [Infinity,-Infinity]
      );

    for( let key=min-2; key < max+2; key++ )
    {
      let i = binary_rangesearch( 0, array.length, i => {
        expect(i).toBeLessThan(array.length);
        expect(i).toBeGreaterThanOrEqual(0);
        return (array[i] < key) - (array[i] > key);
      });
      if( i < 0 )
      {   i = ~i;
        expect(i).toBeLessThanOrEqual(array.length);
        expect(i).toBeGreaterThanOrEqual(0);
        if( i > 0            ) expect(key).toBeLessThan   (array[i-1]);
        if( i < array.length ) expect(key).toBeGreaterThan(array[i  ]);
      }
      else {
        expect(i).toBeLessThan(array.length);
        expect(i).toBeGreaterThanOrEqual(0);
        expect(array[i]).toBe(key);
      }
    }
  })

});
