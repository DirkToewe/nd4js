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

import {argmin,
        argmax,
           min,
           max} from './min_max'
import {forEachItemIn, CUSTOM_MATCHERS} from '../jasmine_utils'
import {_rand_int} from '../_test_data_generators'


describe('nd.iter min/max', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  });




  forEachItemIn(
    function*(){
      for( let run=0; run++ < 4096; )
        yield Object.freeze(
          Array.from({length: _rand_int(1,1337)}, () => _rand_int(-1337,+1337))
        );
    }()
  ).it('argmin works given random examples', seq => {
    const                                                      i = argmin(seq);
              expect(seq           ).toBeAllGreaterOrClose(seq[i], {rtol:0, atol:0});
    if(0!==i) expect(seq.slice(0,i)).toBeAllGreater       (seq[i], {rtol:0, atol:0});
  });


  forEachItemIn(
    function*(){
      for( let run=0; run++ < 4096; )
        yield Object.freeze(
          Array.from({length: _rand_int(1,1337)}, () => _rand_int(-1337,+1337))
        );
    }()
  ).it('argmin works given random examples and custom comparator', seq => {
    const                                                   i = argmin(seq, (x,y) => y-x);
              expect(seq           ).toBeAllLessOrClose(seq[i], {rtol:0, atol:0});
    if(0!==i) expect(seq.slice(0,i)).toBeAllLess       (seq[i], {rtol:0, atol:0});
  });


  forEachItemIn(
    function*(){
      for( let run=0; run++ < 4096; )
        yield Object.freeze(
          Array.from({length: _rand_int(1,1337)}, () => _rand_int(-1337,+1337))
        );
    }()
  ).it('argmax works given random examples', seq => {
    const                                                   i = argmax(seq);
              expect(seq           ).toBeAllLessOrClose(seq[i], {rtol:0, atol:0});
    if(0!==i) expect(seq.slice(0,i)).toBeAllLess       (seq[i], {rtol:0, atol:0});
  });


  forEachItemIn(
    function*(){
      for( let run=0; run++ < 4096; )
        yield Object.freeze(
          Array.from({length: _rand_int(1,1337)}, () => _rand_int(-1337,+1337))
        );
    }()
  ).it('argmax works given random examples and custom comparator', seq => {
    const                                                      i = argmax(seq, (x,y) => y-x);
              expect(seq           ).toBeAllGreaterOrClose(seq[i], {rtol:0, atol:0});
    if(0!==i) expect(seq.slice(0,i)).toBeAllGreater       (seq[i], {rtol:0, atol:0});
  });




  forEachItemIn(
    function*(){
      for( let run=0; run++ < 4096; )
        yield Object.freeze(
          Array.from({length: _rand_int(1,1337)}, () => _rand_int(-1337,+1337))
        );
    }()
  ).it('min works given random examples', seq => {
    const                             v = min(seq);
    expect(seq).toBeAllGreaterOrClose(v, {rtol:0, atol:0});
  })


  forEachItemIn(
    function*(){
      for( let run=0; run++ < 4096; )
        yield Object.freeze(
          Array.from({length: _rand_int(1,1337)}, () => _rand_int(-1337,+1337))
        );
    }()
  ).it('min works given random examples and custom comparator', seq => {
    const                          v = min(seq, (x,y) => y-x);
    expect(seq).toBeAllLessOrClose(v, {rtol:0, atol:0});
  })


  forEachItemIn(
    function*(){
      for( let run=0; run++ < 4096; )
        yield Object.freeze(
          Array.from({length: _rand_int(1,1337)}, () => _rand_int(-1337,+1337))
        );
    }()
  ).it('max works given random examples', seq => {
    const                          v = max(seq);
    expect(seq).toBeAllLessOrClose(v, {rtol:0, atol:0});
  })


  forEachItemIn(
    function*(){
      for( let run=0; run++ < 4096; )
        yield Object.freeze(
          Array.from({length: _rand_int(1,1337)}, () => _rand_int(-1337,+1337))
        );
    }()
  ).it('max works given random examples and custom comparator', seq => {
    const                             v = max(seq, (x,y) => y-x);
    expect(seq).toBeAllGreaterOrClose(v, {rtol:0, atol:0});
  })
});
