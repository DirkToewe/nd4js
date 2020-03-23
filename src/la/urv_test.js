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
import {zip_elems} from '../zip_elems'

import {eye} from './eye'
import {matmul,
        matmul2} from './matmul'
import {permute_cols} from './permute'
import {solve} from './solve'
import {urv_decomp_full,
        urv_lstsq} from './urv'


describe('urv', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from;


  for( const [rng,suffix] of [
    [() =>                           Math.random()*8 - 4, ''                      ],
    [() => Math.random() < 0.1 ? 0 : Math.random()*8 - 4, ' with occasional zeros']
  ])
  forEachItemIn(
    function*(){
      for( let run=1024; run-- > 0; )
      {
        const ndim = randInt(2,5),
             shape = [
               ...Array.from({ length: ndim-2 }, () => randInt(1,8) ),
               randInt(1,64),
               randInt(1,64)
             ];

        yield tabulate(shape, 'float64', rng);
      }
    }()
  ).it('urv_decomp_full works on random examples' + suffix, A => {
    const lead = A.shape.slice(0,-2),
         [M,N] = A.shape.slice(-2);

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

    expect( matmul(U,R,V) ).toBeAllCloseTo(A);
  })


  forEachItemIn(
    function*(){
      for( let run=1024; run-- > 0; )
      {
        const                         sparseness = Math.random(),
          rng = () => Math.random() < sparseness ? 0 : Math.random()*2 - 1;

        const ndim = randInt(2,5),
             shape = [
               ...Array.from({ length: ndim-2 }, () => randInt(1,8) ),
               randInt(1,64),
               randInt(1,64)
             ];

        yield [
          tabulate(shape, 'float64', rng),
          sparseness
        ];
      }
    }()
  ).it('urv_decomp_full works on random sparse examples', ([A,sparseness]) => {
    const lead = A.shape.slice(0,-2),
         [M,N] = A.shape.slice(-2);

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
  })


  forEachItemIn(
    function*(){
      for( let run=1024; run-- > 0; )
      {
        const ndim = randInt(2,5);

        yield _rand_rankdef(
          ...Array.from({ length: ndim-2 }, () => randInt(1,8) ),
          randInt(1,64),
          randInt(1,64)
        );
      }
    }()
  ).it(`urv_decomp_full works on random rank-deficient examples`, ([A,RANKS]) => {
    const lead = A.shape.slice(0,-2),
         [M,N] = A.shape.slice(-2);

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
  })


  forEachItemIn(
    function*(){
      function* shapes() {
        for( let M=8; M > 1; M-- )
        for( let N=8; N > 1; N-- )
          yield [M,N];

        for( let run=1024; run-- > 0; )
        {
          const M = randInt(1,128),
                N = randInt(1,128);
          yield [M,N];
        }
      }

      for( const [M,N] of shapes() )
        yield _rand_rows0(M,N);
    }()
  ).it('urv_decomp_full works on random matrices with zero rows', A => {
    const lead = A.shape.slice(0,-2),
         [M,N] = A.shape.slice(-2);

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

    expect( matmul(U,R,V) ).toBeAllCloseTo(A);
  });


  forEachItemIn(
    function*(){
      function* shapes() {
        for( let M=8; M > 1; M-- )
        for( let N=8; N > 1; N-- )
          yield [M,N];

        for( let run=1024; run-- > 0; )
        {
          const M = randInt(1,128),
                N = randInt(1,128);
          yield [M,N];
        }
      }

      for( const [M,N] of shapes() )
        yield _rand_cols0(M,N);
    }()
  ).it('urv_decomp_full works on random matrices with zero columns', A => {
    const lead = A.shape.slice(0,-2),
         [M,N] = A.shape.slice(-2);

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

    expect( matmul(U,R,V) ).toBeAllCloseTo(A);
  });


  // FIXME make the following test pass (infinite loop due to SRRQR not using scaling)
  forEachItemIn(
    function*(){
      for( let run=1024; run-- > 0; )
      {
        const                         sparseness = Math.random(),
          rng = () => Math.random() < sparseness ? 0 : Math.random()*512*Number.MIN_VALUE;

        const ndim = randInt(2,5),
             shape = [
               ...Array.from({ length: ndim-2 }, () => randInt(1,8) ),
               randInt(1,64),
               randInt(1,64)
             ];

        yield tabulate(shape, 'float64', rng);
      }
    }()
  ).it('urv_decomp_full works on random examples close to zero', A => {
    const lead = A.shape.slice(0,-2),
         [M,N] = A.shape.slice(-2);

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

    expect( matmul(U,R,V) ).toBeAllCloseTo(A, {atol: 1e-7});
  })


  forEachItemIn(
    function*(){
      for( let run=1024; run-- > 0; )
      {
        let ndim = randInt(0,4),
          shapes = [ Array.from({length: ndim}, () => randInt(1,8)) ]
        shapes.splice( randInt(0,2), 0, shapes[0].slice( randInt(0,ndim) ) )

        for( let d=ndim; d > 0; d-- )
        for( let i=randInt(0,2); i-- > 0; ) {
          const    shape = shapes[randInt(0,2)],
               j = shape.length - d
          if(0<=j) shape[j] = 1
        }

        let M = randInt(1,24),
            N = randInt(M,24)+1;
      
        const J = randInt(1,32)
        shapes[0].push(M,N)
        shapes[1].push(M,J)

        yield shapes.map(
          s => tabulate(s,'float64', () => Math.random()*2-1)
        )
      }
    }()
  ).it('urv_decomp_full+urv_lstsq computes one solution of random under-determined examples', ([A,y]) => {
    const [U,R,V,ranks]= urv_decomp_full(A),
                     x = urv_lstsq(U,R,V,ranks, y);

    // check that it's a solution
    expect( matmul2(A,x) ).toBeAllCloseTo(y);

    // check that it's the minimum norm solution
    // https://www.math.usm.edu/lambers/mat419/lecture15.pdf
    const AAT = matmul2(A,A.T);
    expect(x).toBeAllCloseTo( matmul2(A.T, solve(AAT,y)) );
  })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

      for( let run=1024; run-- > 0; )
      {
        let ndim = randInt(0,5),
          shapes = [ Array.from({length: ndim}, () => randInt(1,8)) ]
        shapes.splice( randInt(0,2), 0, shapes[0].slice( randInt(0,ndim) ) )

        for( let d=ndim; d > 0; d-- )
        for( let i=randInt(0,2); i-- > 0; ) {
          const    shape = shapes[randInt(0,2)],
               j = shape.length - d
          if(0<=j) shape[j] = 1
        }

        let M = randInt(1,32),
            N = randInt(1,32),
            L = randInt(1,32);

        shapes[0].push(M,N)
        shapes[1].push(M,L)

        const [A]= _rand_rankdef(...shapes[0]),
               y = tabulate(shapes[1], 'float64', () => Math.random()*8-4);
        Object.freeze(y.data.buffer);
        Object.freeze(y);

        yield [A,y];
      }
    }()
  ).it('urv_decomp_full+urv_lstsq solves least squares for random rank-deficient examples', ([A,y]) => {
    const [U,R,V,ranks]= urv_decomp_full(A),
                     x = urv_lstsq(U,R,V,ranks, y),
                    Ax = matmul2(A,x);
  
    expect( matmul2(A.T, zip_elems([Ax,y], 'float64', (x,y) => x-y) ) ).toBeAllCloseTo(0)
  })
})
