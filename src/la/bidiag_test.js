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

import {matmul, matmul2} from './matmul'
import {bidiag_decomp} from './bidiag'
import {eye} from './eye'


describe('bidiag_decomp', () => {

  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  const test_body = A =>
  {
    const [M,N] = A.shape.slice(-2),
        [U,B,V] = bidiag_decomp(A);
    Object.freeze(U.data.buffer); Object.freeze(U);
    Object.freeze(B.data.buffer); Object.freeze(B);
    Object.freeze(V.data.buffer); Object.freeze(V);

    expect(B).toBeUpperBidiagonal();

    const I = eye( Math.min(M,N) );

    if( M >= N ) {
      expect( matmul2(U.T,U) ).toBeAllCloseTo(I);
      expect( matmul2(V.T,V) ).toBeAllCloseTo(I);
      expect( matmul2(V,V.T) ).toBeAllCloseTo(I);
    }
    else {
      expect( matmul2(U.T,U) ).toBeAllCloseTo(I);
      expect( matmul2(U,U.T) ).toBeAllCloseTo(I);
      expect( matmul2(V,V.T) ).toBeAllCloseTo( eye(M+1) );
    }

    const  a = matmul(U,B,V);
    expect(a).toBeAllCloseTo(A, {atol: 1e-7});
  };


  function* ndarray_shapes(rng) {
    for( let run=0; run++ < 32; )
    for( let   M=0;   M++ <  8; )
    for( let   N=0;   N++ <  8; )
      yield [M,N];

    for( let run=0; run++ < 733; )
    {
      const ndim = rng.int(2,6);
      yield [
        ...Array.from({ length: ndim-2 }, () => rng.int(1,4) ),
        rng.int(1,48),
        rng.int(1,48)
      ];
    }
  }


  function* matrix_shapes(rng) {
    for( let run=0; run++ < 32; )
    for( let   M=0;   M++ <  8; )
    for( let   N=0;   N++ <  8; )
      yield [M,N];

    for( let run=337; run-- > 0; )
    {
      const M = rng.int(1,128),
            N = rng.int(1,128); // <- TODO remove after testing
      yield [M,N];
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
