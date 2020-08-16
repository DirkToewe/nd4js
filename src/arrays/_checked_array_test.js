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
import {_rand_int} from "../_test_data_generators";

import {checked_array,
        IndexOutOfBoundsError} from "./_checked_array";


for( const Arr of [
         Array,
    Int32Array,
  Float64Array
])
  describe(`checked_array(_: ${Arr.name})`, () => {

    forEachItemIn(
      function*(){
        for( let length=0; length < 64; length++ )
        for( let    run=0;    run < 64;    run++ )
        {
          const arr = Arr.from({length}, () => _rand_int(-1337,+1337) ),
                idx = run < 8 ? -1-run : _rand_int(Number.MIN_SAFE_INTEGER,0);
          yield [arr,idx];
        }
      }()
    ).it('should protect from negative index read access', ([arr,idx]) => {
      expect( () => checked_array(arr)[idx] ).toThrowMatching( err => err instanceof IndexOutOfBoundsError );
    })

    forEachItemIn(
      function*(){
        for( let length=0; length < 64; length++ )
        for( let    run=0;    run < 64;    run++ )
        {
          const arr = Arr.from({length}, () => _rand_int(-1337,+1337) ),
                idx = run < 8 ? arr.length+run : _rand_int(arr.length, Number.MAX_SAFE_INTEGER);
          yield [arr,idx];
        }
      }()
    ).it('should protect from too large index read access', ([arr,idx]) => {
      expect( () => checked_array(arr)[idx] ).toThrowMatching( err => err instanceof IndexOutOfBoundsError );
    })

    forEachItemIn(
      function*(){
        for( let length=1; length < 32; length++ )
        for( let    run=0;    run < 64;    run++ )
        {
          const arr = Arr.from({length}, () => _rand_int(-1337,+1337) ),
                idx = run < arr.length ? run : _rand_int(0, arr.length);
          yield [arr,idx];
        }
      }()
    ).it('should allow correct index read access', ([arr,idx]) => {
      expect( checked_array(arr)[idx] ).toBe( arr[idx] );
    })



    forEachItemIn(
      function*(){
        for( let length=0; length < 64; length++ )
        for( let    run=0;    run < 64;    run++ )
        {
          const arr = Arr.from({length}, () => _rand_int(-1337,+1337) ),
                idx = run < 8 ? -1-run : _rand_int(Number.MIN_SAFE_INTEGER,0),
                val = _rand_int(Number.MIN_SAFE_INTEGER, Number.MAX_SAFE_INTEGER);
          yield [arr,idx,val];
        }
      }()
    ).it('should protect from negative index write access', ([arr,idx,val]) => {
      const backup = Arr.from(arr);
      expect( () => checked_array(arr)[idx] = val ).toThrowMatching( err => err instanceof IndexOutOfBoundsError );
      expect(arr).toEqual(backup);
    })

    forEachItemIn(
      function*(){
        for( let length=0; length < 64; length++ )
        for( let    run=0;    run < 64;    run++ )
        {
          const arr = Arr.from({length}, () => _rand_int(-1337,+1337) ),
                idx = run < 8 ? arr.length+run : _rand_int(arr.length, Number.MAX_SAFE_INTEGER),
                val = _rand_int(Number.MIN_SAFE_INTEGER, Number.MAX_SAFE_INTEGER);
          yield [arr,idx,val];
        }
      }()
    ).it('should protect from too large index write access', ([arr,idx,val]) => {
      const backup = Arr.from(arr);
      expect( () => checked_array(arr)[idx] = val ).toThrowMatching( err => err instanceof IndexOutOfBoundsError );
      expect(arr).toEqual(backup);
    })

    forEachItemIn(
      function*(){
        for( let length=1; length < 32; length++ )
        for( let    run=0;    run < 64;    run++ )
        {
          const arr = Arr.from({length}, () => _rand_int(-1337,+1337) ),
                idx = run < arr.length ? run : _rand_int(0, arr.length),
                val = _rand_int(Number.MIN_SAFE_INTEGER, Number.MAX_SAFE_INTEGER);
          yield [arr,idx,val];
        }
      }()
    ).it('should allow correct index write access', ([arr,idx,val]) => {
      const backup = Arr.from(arr);
      expect( checked_array(arr)[idx] = val ).toBe( val );
                          backup[idx] = val;
      expect(arr).toEqual(backup);
    })

  });
