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


import {forEachItemIn, CUSTOM_MATCHERS} from '../jasmine_utils'
import {_rand_rankdef} from '../_test_data_generators'
import {tabulate} from '../tabulate'
import {_rand_rows0,
        _rand_cols0} from '../_test_data_generators'

import {eye} from './eye'
import {generic_test_lstsq} from "./_generic_test_lstsq";
import {generic_test_rank } from "./_generic_test_rank";
import {matmul,
        matmul2} from './matmul'
import {urv_decomp_full,
        urv_lstsq} from './urv'


for( const lstsq of Object.values({
  [`urv_decomp + urv_lstsq   `]: (A,y) => urv_lstsq(   urv_decomp_full(A), y ),
  [`urv_decomp + urv_lstsq...`]: (A,y) => urv_lstsq(...urv_decomp_full(A), y ),
}))
  generic_test_lstsq(lstsq);


for( const rank of Object.values({
  [`urv_decomp(_)[3]`]: A => urv_decomp_full(A)[3]
}))
  generic_test_rank(rank);


describe('urv_decomp_full', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  });


  const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from;


  const test_urv_decomp_full = A =>
  {
    const lead = A.shape.slice(0,-2),
         [M,N] = A.shape.slice(  -2);

    const [U,R,V, ranks] = urv_decomp_full(A);

    expect(ranks.shape).toEqual(lead)
    expect(ranks.dtype).toBe('int32')

    expect(U.dtype).toBe(A.dtype);
    expect(R.dtype).toBe(A.dtype);
    expect(V.dtype).toBe(A.dtype);

    expect( U.shape ).toEqual( Int32Array.of(...lead, M,M) )
    expect( R.shape ).toEqual( Int32Array.of(...lead, M,N) )
    expect( V.shape ).toEqual( Int32Array.of(...lead, N,N) )

    expect(R).toBeUpperTriangular();

    let         i = 0;
    for( const Ri of R.reshape(-1,M,N) )
    {
      const  rnk = ranks.data[i++];
      if(N > rnk) { 
        const  rhs = Ri.sliceElems([,,], [rnk,,]);
        expect(rhs).toBeAllCloseTo(0, {rtol:0, atol:0});
      }
      if(M > rnk ) {
        const  bot = Ri.sliceElems([rnk,,]);
        expect(bot).toBeAllCloseTo(0, {rtol:0, atol:0});
      }
    }

    let                                     I = eye(M,M);
    expect( matmul2(U,U.T) ).toBeAllCloseTo(I);
    expect( matmul2(U.T,U) ).toBeAllCloseTo(I);
                                            I = eye(N,N);
    expect( matmul2(V,V.T) ).toBeAllCloseTo(I);
    expect( matmul2(V.T,V) ).toBeAllCloseTo(I);

    expect( matmul(U,R,V) ).toBeAllCloseTo(A, {atol: 1e-5});
  };


  for( const [rng,suffix] of [
    [() =>                           Math.random()*8 - 4, ''                      ],
    [() => Math.random() < 0.1 ? 0 : Math.random()*8 - 4, ' with occasional zeros']
  ])
  forEachItemIn(
    function*(){
      for( let run=733; run-- > 0; )
      {
        const ndim = randInt(2,6),
            shape = [
              ...Array.from({ length: ndim-2 }, () => randInt(1,4) ),
              randInt(1,32),
              randInt(1,32)
            ];

        yield tabulate(shape, 'float64', rng);
      }
    }()
  ).it('correctly decomposes random examples' + suffix, test_urv_decomp_full);


  forEachItemIn(
    function*(){
      for( let run=733; run-- > 0; )
      {
        const                         sparseness = Math.random(),
          rng = () => Math.random() < sparseness ? 0 : Math.random()*2 - 1;

        const ndim = randInt(2,6),
             shape = [
               ...Array.from({ length: ndim-2 }, () => randInt(1,4) ),
               randInt(1,32),
               randInt(1,32)
             ];

        yield tabulate(shape, 'float64', rng);
      }
    }()
  ).it('correctly decomposes random sparse examples', test_urv_decomp_full);


  forEachItemIn(
    function*(){
      function* shapes() {
        for( let M=8; M > 1; M-- )
        for( let N=8; N > 1; N-- )
          yield [M,N];

        for( let run=128; run-- > 0; )
        {
          const M = randInt(1,128),
                N = randInt(1,128);
          yield [M,N];
        }
      }

      for( const [M,N] of shapes() )
        yield _rand_rows0(M,N);
    }()
  ).it('correctly decomposes random matrices with zero rows', test_urv_decomp_full);


  forEachItemIn(
    function*(){
      function* shapes() {
        for( let M=8; M > 1; M-- )
        for( let N=8; N > 1; N-- )
          yield [M,N];

        for( let run=128; run-- > 0; )
        {
          const M = randInt(1,128),
                N = randInt(1,128);
          yield [M,N];
        }
      }

      for( const [M,N] of shapes() )
        yield _rand_cols0(M,N);
    }()
  ).it('correctly decomposes random matrices with zero columns', test_urv_decomp_full);


  // FIXME make the following test pass (infinite loop due to SRRQR not using scaling)
  forEachItemIn(
    function*(){
      for( let run=733; run-- > 0; )
      {
        const                         sparseness = Math.random(),
          rng = () => Math.random() < sparseness ? 0 : Math.random()*512*Number.MIN_VALUE;

        const ndim = randInt(2,6),
             shape = [
               ...Array.from({ length: ndim-2 }, () => randInt(1,4) ),
               randInt(1,32),
               randInt(1,32)
             ];

        yield tabulate(shape, 'float64', rng);
      }
    }()
  ).it('correctly decomposes random examples close to zero', test_urv_decomp_full);


  forEachItemIn(
    function*(){
      for( let run=512; run-- > 0; )
      {
        const ndim = randInt(2,6);

        yield _rand_rankdef(
          ...Array.from({ length: ndim-2 }, () => randInt(1,4) ),
          randInt(1,32),
          randInt(1,32)
        );
      }
    }()
  ).it(`correctly decomposes random rank-deficient examples`, ([A,RANKS]) => {
    const lead = A.shape.slice(0,-2),
         [M,N] = A.shape.slice(  -2);

    const [U,R,V, ranks] = urv_decomp_full(A);

    expect(ranks.shape).toEqual(lead)
    expect(RANKS.shape).toEqual(lead)
    expect(ranks.dtype).toBe('int32')
    expect(RANKS.dtype).toBe('int32')
    expect(ranks).toBeAllCloseTo(RANKS, {rtol:0, atol:0});

    expect(U.dtype).toBe(A.dtype);
    expect(R.dtype).toBe(A.dtype);
    expect(V.dtype).toBe(A.dtype);

    expect( U.shape ).toEqual( Int32Array.of(...lead, M,M) )
    expect( R.shape ).toEqual( Int32Array.of(...lead, M,N) )
    expect( V.shape ).toEqual( Int32Array.of(...lead, N,N) )

    expect(R).toBeUpperTriangular();

    let       i = 0;
    for( let Ri of R.reshape(-1,M,N) )
    {
      const  rnk = ranks.data[i++];
      if(N > rnk) { 
        const  rhs = Ri.sliceElems([,,], [rnk,,]);
        expect(rhs).toBeAllCloseTo(0, {rtol:0, atol:0});
      }
      if(M > rnk ) {
        const  bot = Ri.sliceElems([rnk,,]);
        expect(bot).toBeAllCloseTo(0, {rtol:0, atol:0});
      }
    }

    let                                     I = eye(M,M);
    expect( matmul2(U,U.T) ).toBeAllCloseTo(I);
    expect( matmul2(U.T,U) ).toBeAllCloseTo(I);
                                            I = eye(N,N);
    expect( matmul2(V,V.T) ).toBeAllCloseTo(I);
    expect( matmul2(V.T,V) ).toBeAllCloseTo(I);

    expect( matmul(U,R,V) ).toBeAllCloseTo(A);
  });
});
