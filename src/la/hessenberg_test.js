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
import {tabulate} from '../tabulate'
import {_rand_rows0,
        _rand_cols0} from '../_test_data_generators'

import {hessenberg_decomp} from './hessenberg'
import {matmul, matmul2} from './matmul'
import {eye} from './eye'


describe('hessenberg', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  const test_body = A => {
    const [N]= A.shape.slice(-1),
        [U,H]= hessenberg_decomp(A);

    expect(U.shape).toEqual(A.shape);
    expect(H.shape).toEqual(A.shape);

    expect(H).toBeUpperHessenberg();

    const                                   I = eye(N);
    expect( matmul2(U,U.T) ).toBeAllCloseTo(I);
    expect( matmul2(U.T,U) ).toBeAllCloseTo(I);

    const  a = matmul(U,H,U.T);
    expect(a).toBeAllCloseTo(A);
  };


  function* ndarray_shapes(rng) {
    for( let run=0; run++ < 32; )
    for( let   N=0;   N++ <  8; )
      yield [N,N];

    for( let run=0; run++ < 512; )
    {
      const N = rng.int(1,48),
         ndim = rng.int(2,6);
      yield [
        ...Array.from({ length: ndim-2 }, () => rng.int(1,4) ),
        N,N
      ];
    }
  }


  function* matrix_shapes(rng) {
    for( let run=0; run++ < 32; )
    for( let   N=0;   N++ <  8; )
      yield [N,N];

    for( let run=256; run-- > 0; )
    {
      const  N = rng.int(1,128); // <- TODO remove after testing
      yield [N,N];
    }
  }


  for( const [rand,suffix] of [
    [rng => () => rng.uniform(-4,+4)                           , ''                      ],
    [rng => () => rng.uniform(-4,+4) * (rng.uniform(0,1) < 0.9), ' with occasional zeros']
  ])
    forEachItemIn(
      function*(rng){
        for( const shape of ndarray_shapes(rng) )
          yield tabulate(shape, 'float64', rand(rng));
      }
    ).it('works on random examples'+suffix, test_body);


  forEachItemIn(
    function*(rng){
      for( const shape of ndarray_shapes(rng) ) {
        const                                  density = rng.uniform(0,1);
        yield tabulate(shape, 'float64', () =>(density > rng.uniform(0,1)) * rng.uniform(-4,+4) );
      }
    }
  ).it('works on random sparse examples', test_body);


  forEachItemIn(
    function*(rng){
      for( const shape of ndarray_shapes(rng) )
        yield rng.rankDef(...shape)[0];
    }
  ).it('works on random rank-deficient examples', test_body);


  forEachItemIn(
    function*(rng){
      for( const [M,N] of matrix_shapes(rng) )
        yield _rand_rows0(M,N);
    }
  ).it('works on random matrices with zero rows', test_body);


  forEachItemIn(
    function*(rng){
      for( const [M,N] of matrix_shapes(rng) )
        yield _rand_cols0(M,N);
    }
  ).it('works on random matrices with zero columns', test_body);

})
