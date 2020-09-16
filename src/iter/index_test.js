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
        range,
        repeat,
        zip} from '.'
import {forEachItemIn, CUSTOM_MATCHERS} from '../jasmine_utils'
import {_rand_int} from '../_test_data_generators'


describe('nd.iter', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  forEachItemIn(
    function*(){
      for( let run=0; run++ < 2048; )
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
      for( let start = -32; start <= 32; start++ )
      for( let stop  = -32; stop  <= 32; stop ++ )
      for( let step  =   1; step  <= 64; step ++ )
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
      for( let start = -48; start <= 48; start++ )
      for( let stop  = -48; stop  <= 48; stop ++ )
      for( let step  =   1; step  <= 48; step ++ )
        yield [start, stop, -step];
    }()
  ).it('range works for generated examples with step < 0', ([start, stop, step]) => {
    const seq = range(start, stop, step);

    expect(step).toBeLessThan(0);

    for( let i=start; i > stop; i+=step )
      expect(seq.next().value).toBe(i);

    expect(seq.next().done).toBe(true);
  })


  it('cartesian_prod works for a single  example  with 0, 1, 2, 3 and 4 args', () => {
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

    expect([...cartesian_prod(       ) ]).toEqual( Object.freeze([[]]     ) );
    expect([...cartesian_prod(a      ) ]).toEqual( Object.freeze(prod_a   ) );
    expect([...cartesian_prod(a,b    ) ]).toEqual( Object.freeze(prod_ab  ) );
    expect([...cartesian_prod(a,b,c  ) ]).toEqual( Object.freeze(prod_abc ) );
    expect([...cartesian_prod(a,b,c,d) ]).toEqual( Object.freeze(prod_abcd) );
  })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from;

      for( let i=0; i < 7; i++ ) { const a = function(){ const a0 = randInt(-1337,+1337); return Array.from({length: i}, (_,h) => h+a0); }();
      for( let j=0; j < 7; j++ ) { const b = function(){ const b0 = randInt(-1337,+1337); return Array.from({length: j}, (_,h) => h+b0); }();
      for( let k=0; k < 7; k++ ) { const c = function(){ const c0 = randInt(-1337,+1337); return Array.from({length: k}, (_,h) => h+c0); }();
      for( let l=0; l < 7; l++ ) { const d = function(){ const d0 = randInt(-1337,+1337); return Array.from({length: l}, (_,h) => h+d0); }();
        yield Object.freeze([a,b,c,d]);
      }}}}
    }()
  ).it('cartesian_prod works for generated examples with    1, 2, 3 and 4 args', ([a,b,c,d]) => {
    // console.log('a,b,c,d:\n[' + a + ']\n[' + b + ']\n[' + c + ']\n[' + d + ']');

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
      for( let run=0; run++ < 1733; )
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


  forEachItemIn(
    function*(){
      for( let run=0; run++ < 3*1337; )
        yield Object.freeze([
          Array.from({length: _rand_int(0,512)}, () => _rand_int(-1337,+1337) ),
          Array.from({length: _rand_int(0,512)}, () => _rand_int(-1337,+1337) )
        ]);
    }()
  ).it('zip works given random inputs of 2 arrays', ([L,R]) => {
    const N = Math.min(L.length, R.length);

    let n=0;

    for( const lr of zip(L,R) )
    {
      expect(n).toBeLessThan(N);
      expect(lr).toEqual([ L[n], R[n] ]);
      n++;
    }
    expect(n).toBe(N);
  })


  forEachItemIn(
    function*(){
      for( let run=0; run++ < 3*1337; )
        yield Object.freeze([
          Array.from({length: _rand_int(0,512)}, () => _rand_int(-1337,+1337) ),
          Array.from({length: _rand_int(0,512)}, () => _rand_int(-1337,+1337) ),
          Array.from({length: _rand_int(0,512)}, () => _rand_int(-1337,+1337) )
        ]);
    }()
  ).it('zip works given random inputs of 3 arrays', ([L,M,R]) => {
    const N = Math.min(L.length, M.length, R.length);

    let n=0;

    for( const lmr of zip(L,M,R) )
    {
      expect(n).toBeLessThan(N);
      expect(lmr).toEqual([ L[n], M[n], R[n] ]);
      n++;
    }
    expect(n).toBe(N);
  })


  forEachItemIn(
    cartesian_prod(
      ['Cat', 'Dog', 'Golden hamster', 'Guinea pig', 'Tarantula', 'Turtle', 'Bird', 'Gold fish'],
      range(0,17)
    )
  ).it('repeat(n,seq) works given generated inputs', ([s,n]) => {
    expect( [...repeat(n,s) ].join('') ).toBe( s.repeat(n) );
  })


  forEachItemIn(
    ['Cat', 'Dog', 'Golden hamster', 'Guinea pig', 'Tarantula', 'Turtle', 'Bird', 'Gold fish', [1,2,3], [5,6], [7]]
  ).it('repeat(  seq) works given generated inputs', s => {
    const r = repeat(s)[Symbol.iterator]();

    for( let run=0; run++ < 1337; )
      for( const x of s )
        expect(r.next().value).toBe(x);
  })
})
