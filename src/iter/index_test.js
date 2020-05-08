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

import {cartesian_prod,
        enumerate,
        linspace,
        range} from '.'
import {forEachItemIn, CUSTOM_MATCHERS} from '../jasmine_utils'
import {_rand_int} from '../_test_data_generators'


describe('_iter_utils', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  forEachItemIn(
    function*(){
      for( let run=0; run++ < 16*1024; )
      {
        const start = Math.random()*512 - 256,
                end = Math.random()*512 - 256;

        for( let num=1; num++ < 16; )
          yield [start,end,num];
      }
    }()
  ).it('linspace works for random examples', ([start,end,num]) => {
    const samples = [...linspace(start,end,num) ];

    expect(samples.length).toBe(num);

    expect(samples[               0]).toBe(start);
    expect(samples[samples.length-1]).toBe(end);

    for( let i=2; i < num; i++ )
      expect(samples[i]-samples[i-1]).toBeAllCloseTo(samples[i-1]-samples[i-2])
  })


  forEachItemIn(
    function*(){
      for( let start = -64; start <=  64; start++ )
      for( let stop  = -64; stop  <=  64; stop ++ )
      for( let step  =   1; step  <= 128; step ++ )
        yield [start, stop, step];
    }()
  ).it('range works for generated examples with step > 0', ([start, stop, step]) => {
    const seq = range(start, stop, step);

    expect(step).toBeGreaterThan(0);

    for( let i=start; i < stop; i+=step )
      expect(seq.next().value).toBe(i);

    expect(seq.next().done).toBe(true);
  })


  forEachItemIn(
    function*(){
      for( let start = -64; start <= 64; start++ )
      for( let stop  = -64; stop  <= 64; stop ++ )
        yield [start, stop];
    }()
  ).it('range works for generated examples with step=undefined', ([start, stop]) => {
    for( const seq of [
      range(start, stop),
      range(start, stop, 1),
      range(start, stop, undefined)
    ])
    {
      for( let i=start; i < stop; i++ )
        expect(seq.next().value).toBe(i);

      expect(seq.next().done).toBe(true);
    }
  })


  forEachItemIn(
    function*(){
      for( let start = -64; start <=  64; start++ )
      for( let stop  = -64; stop  <=  64; stop ++ )
      for( let step  =   1; step  <= 128; step ++ )
        yield [start, stop, -step];
    }()
  ).it('range works for generated examples with step < 0', ([start, stop, step]) => {
    const seq = range(start, stop, step);

    expect(step).toBeLessThan(0);

    for( let i=start; i > stop; i+=step )
      expect(seq.next().value).toBe(i);

    expect(seq.next().done).toBe(true);
  })


  it('cartesian_prod works for examples with 1, 2, 3 and 4 args', () => {
    const a = Object.freeze([...range( 16, 32) ]),
          b = Object.freeze([...range( 64, 80) ]),
          c = Object.freeze([...range( 96, 99) ]),
          d = Object.freeze([...range(-17, -3) ]);

    const prod_a = [],
          prod_ab = [],
          prod_abc = [],
          prod_abcd = [];
    for( const x of a ) { prod_a   .push([x]);
    for( const y of b ) { prod_ab  .push([x,y]);
    for( const z of c ) { prod_abc .push([x,y,z]);
    for( const u of d ) { prod_abcd.push([x,y,z,u]); }}}}
      
    expect([...cartesian_prod(a      ) ]).toEqual( Object.freeze(prod_a   ) );
    expect([...cartesian_prod(a,b    ) ]).toEqual( Object.freeze(prod_ab  ) );
    expect([...cartesian_prod(a,b,c  ) ]).toEqual( Object.freeze(prod_abc ) );
    expect([...cartesian_prod(a,b,c,d) ]).toEqual( Object.freeze(prod_abcd) );
  })


  forEachItemIn(
    function*(){
      for( let run=0; run++ < 7*1337; )
      {
        const length = _rand_int(0,1337);

        yield Object.freeze(
          Array.from({length}, () => _rand_int(-1337,+1337) )
        );
      }
    }()
  ).it('enumerate works given random examples', arr => {
    let n=0;
    for( const [i,x] of enumerate(arr) )
    {
      expect(i).toBe(n++);
      expect(x).toBe(arr[i])
    }
    expect(n).toBe(arr.length);
  })
})
