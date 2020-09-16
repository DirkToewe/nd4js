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

import {forEachItemIn, CUSTOM_MATCHERS} from '../jasmine_utils'
import {generic_test_lstsq} from './_generic_test_lstsq';
import {generic_test_rank } from './_generic_test_rank' ;
import {generic_test_solve} from './_generic_test_solve';
import {eps} from '../dt'
import {array} from '../nd_array'
import {tabulate} from '../tabulate';

import {diag_mat} from './diag'
import {eye} from './eye'
import {matmul, matmul2} from './matmul'
import {FrobeniusNorm, norm} from "./norm";
import {svd_lstsq,
        svd_rank,
        svd_solve} from './svd'


export function generic_test_svd_decomp( svd_decomp )
{
  const [ rank,
    solve1,solve2,
    lstsq1,lstsq2
  ] = Object.values({
    [`${svd_decomp.name} + svd_rank`]: A => {
      const        [U,sv,V] = svd_decomp(A);
      return svd_rank(sv);
    },
    [`${svd_decomp.name} + svd_solve   `]: (A,y) => svd_solve(   svd_decomp(A), y),
    [`${svd_decomp.name} + svd_solve...`]: (A,y) => svd_solve(...svd_decomp(A), y),
    [`${svd_decomp.name} + svd_lstsq   `]: (A,y) => svd_lstsq(   svd_decomp(A), y),
    [`${svd_decomp.name} + svd_lstsq...`]: (A,y) => svd_lstsq(...svd_decomp(A), y)
  });
  generic_test_lstsq(lstsq1);
  generic_test_lstsq(lstsq2);
  generic_test_rank ( rank );
  generic_test_solve(solve1);
  generic_test_solve(solve2);

  describe(`${svd_decomp.name} [generic SVD tests]`, () => {

    beforeEach( () => {
      jasmine.addMatchers(CUSTOM_MATCHERS)
    });


    const test_ndarray = A =>
    {
      expect(A.ndim).toBeGreaterThanOrEqual(2);

      const [M,N] = A.shape.slice(-2), L = Math.min(M,N),
         [U,sv,V] = svd_decomp(A),     D = diag_mat(sv);
      Object.freeze( U.data.buffer); Object.freeze( U);
      Object.freeze(sv.data.buffer); Object.freeze(sv);
      Object.freeze( V.data.buffer); Object.freeze( V);
      Object.freeze( D.data.buffer); Object.freeze( D);

      expect(sv.dtype).toMatch(/^float/)
      expect(sv.dtype).toBe(U.dtype)
      expect(sv.dtype).toBe(V.dtype)

      expect( U.shape).toEqual( Int32Array.of(...A.shape.slice(0,-2),M,L  ) );
      expect(sv.shape).toEqual( Int32Array.of(...A.shape.slice(0,-2),  L  ) );
      expect( V.shape).toEqual( Int32Array.of(...A.shape.slice(0,-2),  L,N) );

      const abs = x => x < 0 ? -x : x;

      for( let i=sv.data.length; (i-=L) >= 0; )
      {
        const  sv_i = sv.data.subarray(i,i+L);
        expect(sv_i).toEqual( sv_i.map(abs).sort((x,y) => y-x) );
      }

      const a = matmul(U,D,V);

      const U_TOL = eps(A.dtype) * M*4,
            V_TOL = eps(A.dtype) * N*4;

      const I = eye(L);
      if( M >= N ) {
        expect( matmul2(U.T,U) ).toBeAllCloseTo(I, {rtol:0, atol:U_TOL});
        expect( matmul2(V.T,V) ).toBeAllCloseTo(I, {rtol:0, atol:V_TOL});
      }
      if( M <= N ) {
        expect( matmul2(U,U.T) ).toBeAllCloseTo(I, {rtol:0, atol:U_TOL});
        expect( matmul2(V,V.T) ).toBeAllCloseTo(I, {rtol:0, atol:V_TOL});
      }
      expect(D).toBeDiagonal()
      expect(a).toBeAllCloseTo(A, {atol:1e-7})
    };


    const test_matrix = A =>
    {
      expect(A.ndim).toBe(2);

      const [M,N]= A.shape,       L = Math.min(M,N),
         [U,sv,V]= svd_decomp(A), D = diag_mat(sv);
      Object.freeze( U.data.buffer); Object.freeze( U);
      Object.freeze(sv.data.buffer); Object.freeze(sv);
      Object.freeze( V.data.buffer); Object.freeze( V);
      Object.freeze( D.data.buffer); Object.freeze( D);

      // expect(sv.dtype).toMatch(/^float/)
      expect(sv.dtype).toMatch('float64')
      expect(sv.dtype).toBe(U.dtype)
      expect(sv.dtype).toBe(V.dtype)

      expect( U.shape).toEqual( Int32Array.of(...A.shape.slice(0,-2),M,L  ) );
      expect(sv.shape).toEqual( Int32Array.of(...A.shape.slice(0,-2),  L  ) );
      expect( V.shape).toEqual( Int32Array.of(...A.shape.slice(0,-2),  L,N) );

      const abs = x => x < 0 ? -x : x;

      for( let i=sv.data.length; (i-=L) >= 0; )
      {
        const  sv_i = sv.data.subarray(i,i+L);
        expect(sv_i).toEqual( sv_i.map(abs).sort((x,y) => y-x) );
      }

      const a = matmul(U,D,V);

      expect(a.shape).toEqual(A.shape);

      const A_TOL = eps(A.dtype) *48 * Math.max(M,N) * norm(A),
            U_TOL = eps(A.dtype) * 4 * M,
            V_TOL = eps(A.dtype) * 4 * N;

      const I = eye(L);
      if( M >= N ) {
        expect( matmul2(U.T,U) ).toBeAllCloseTo(I, {rtol:0, atol:U_TOL});
        expect( matmul2(V.T,V) ).toBeAllCloseTo(I, {rtol:0, atol:V_TOL});
      }
      if( M <= N ) {
        expect( matmul2(U,U.T) ).toBeAllCloseTo(I, {rtol:0, atol:U_TOL});
        expect( matmul2(V,V.T) ).toBeAllCloseTo(I, {rtol:0, atol:V_TOL});
      }
      expect(D).toBeDiagonal();

      const NORM = new FrobeniusNorm();
      for( let i=0; i < M*N; i++ )
        NORM.include(A.data[i] - a.data[i]);

      // expect(NORM.result).toBe( norm(zip_elems([A,a], (x,y) => x-y)) );

      expect(NORM.result).toBeLessThanOrEqual(A_TOL);
    };


    it(
      ` correctly decomposes hand-crafted int32 example correctly`,
      () => test_matrix(array(
        'int32',
        [[1,1],
         [1,2],
         [1,3],
         [1,4],
         [1,5]]
      ))
    );


    forEachItemIn(
      function*(rng){
        for( let run=1024; run-- > 0; )
        {
          const N = rng.int(1,32),
             ndim = rng.int(0,3),
            shape = Array.from({ length: ndim }, () => rng.int(1,4) );
          shape.push(N);

          const S = tabulate(shape, 'float64', () => rng.uniform(-4,+4) * (rng.uniform(0,1) < 0.9) ),
                A = diag_mat(S);

          for( let i=S.data.length; i-- > 0; )
            S.data[i] = Math.abs(S.data[i]);

          for( let i=S.data.length; (i -= N) >= 0; )
            S.data.subarray(i,i+N).sort( (x,y) => y-x );

          Object.freeze(S.data.buffer)
          Object.freeze(A.data.buffer)
          yield [S,A];
        }
      }
    ).it(` correctly decomposes random batches of diagonal matrices`, ([SV,A]) => {
      const  [N] = A.shape.slice(-1),
        [U,sv,V] = svd_decomp(A);

      const I = eye(N),
          tol = svd_decomp.name.startsWith('svd_jac') ? {rtol:0, atol:0} : {};

      expect(sv).toBeAllCloseTo(SV, tol);
      expect( matmul2(U,U.T) ).toBeAllCloseTo(I, tol);
      expect( matmul2(U.T,U) ).toBeAllCloseTo(I, tol);
      expect( matmul2(V,V.T) ).toBeAllCloseTo(I, tol);
      expect( matmul2(V.T,V) ).toBeAllCloseTo(I, tol);
      expect( matmul(U,diag_mat(sv),V) ).toBeAllCloseTo(A, tol);
    });


    for( const [tab_fn,suffix] of [
      [rng => () => rng.uniform(-4,+4)                           , ''                      ],
      [rng => () => rng.uniform(-4,+4) * (rng.uniform(0,1) < 0.9), ' with occasional zeros']
    ])
      forEachItemIn(
        function*(rng){
          for( let run=1024; run-- > 0; )
          {
            const ndim = rng.int(0,3),
                 shape = Array.from({ length: ndim }, () => rng.int(1,4) );
            shape.push(
              rng.int(1,32),
              rng.int(1,32)
            );
            const A = tabulate(shape, 'float64', tab_fn(rng)); Object.freeze(A.data.buffer);
            yield A;
          }
        }
      ).it(` correctly decomposes random examples${suffix}`, test_ndarray);


    forEachItemIn(
      function*(rng){
        for( let run=1024; run-- > 0; )
        {
          const ndim = rng.int(0,3),
               shape = Array.from({ length: ndim }, () => rng.int(1,4) );
          shape.push(
            rng.int(1,32),
            rng.int(1,32)
          );
          const [A] = rng.rankDef(...shape);
          yield  A;
        }
      }
    ).it(` correctly decomposes random rank-deficient examples`, test_ndarray);
  

    forEachItemIn(
      function*(rng){
        for( let run=733; run-- > 0; )
        {
          const ndim = rng.int(0,3),
               shape = Array.from({ length: ndim }, () => rng.int(1,4) );
          shape.push(
            rng.int(1,32),
            rng.int(1,32)
          );
          const                                  sparseness = rng.uniform(0,1),
            A = tabulate(shape, 'float64', () => sparseness > rng.uniform(0,1) ? 0 : rng.uniform(-4,+4) );
          Object.freeze(A.data.buffer);
          yield A;
        }
      }
    ).it(` correctly decomposes random sparse examples`, test_ndarray);


    for( const [tab_fn,suffix] of [
      [rng => () => rng.uniform(-4,+4)                           , ''                      ],
      [rng => () => rng.uniform(-4,+4) * (rng.uniform(0,1) < 0.9), ' with occasional zeros']
    ])
      forEachItemIn(
        function*(rng){
          function* shapes() {
            for( let M=0; M++ < 16; )
            for( let N=0; N++ < 16; )
              yield [M,N];

            for( let run=0; run++ < 337; )
            {
              const M = rng.int(1,64),
                    N = rng.int(1,64);
              yield [M,N];
            }
          }

          for( const [M,N] of shapes() )
          {
            const         A = tabulate([M,N], 'float64', tab_fn(rng));
            Object.freeze(A.data.buffer);
            yield         A;
          }
        }
      ).it(
        'accurately decomposes random matrices'+suffix,
        test_matrix
      );


    for( const dr of [0,1] )
    for( const dc of [0,1] )
    if( 1 !== dr*dc )
    forEachItemIn(
      function*(rng){
        function* sizes() {
          const steps_per_binade = 3;

          for( let N=0; N++ < 16; )
            yield N;

          for( let run=4*steps_per_binade; run < 8*steps_per_binade; run++ )
            yield Math.round(2**(run/steps_per_binade))
        }

        for( const L of sizes() )
        {
          const M = L+dr,
                N = L+dc;
          const        [A] = rng.rankDef(M,N);
          Object.freeze(A);
          Object.freeze(A.data.buffer);
          yield         A;
        }
      }
    ).it(
      `accurately decomposes random rank-deficient matrices of shape [N+${dr},N+${dc}]`,
      test_matrix
    );


    forEachItemIn(
      function*(rng){
        function* shapes() {
          for( let M=0; M++ < 16; )
          for( let N=0; N++ < 16; )
            yield [M,N];

          for( let run=0; run++ < 256; )
          {
            const M = rng.int(1,64),
                  N = rng.int(1,64);
            yield [M,N];
          }
        }
  
        for( const [M,N] of shapes() )
        {
          const                                  sparseness = rng.uniform(0,1),
            A = tabulate([M,N], 'float64', () => sparseness > rng.uniform(0,1) ? 0 : rng.uniform(-4,+4) );
          Object.freeze(A.data.buffer);
          yield A;
        }
      }
    ).it(
      'accurately decomposes random sparse matrices',
      test_matrix
    );


// FIXME make svd_dc pass the following test:
//  forEachItemIn(
//    function*(){
//      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from;
//
//      function* shapes() {
//        for( let M=8; M > 1; M-- )
//        for( let N=8; N > 1; N-- )
//          yield [M,N];
//
//        for( let run=1024; run-- > 0; )
//        {
//          const M = randInt(1,64),
//                N = randInt(1,64);
//          yield [M,N];
//        }
//      }
//
//      for( const [M,N] of shapes() )
//        yield _rand_rows0(M,N);
//    }()
//  ).it(
//    'svd_dc accurately decomposes random matrices with zero rows',
//    test_accuracy
//  );


// FIXME make svd_dc pass the following test:
//  forEachItemIn(
//    function*(){
//      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from;
//
//      function* shapes() {
//        for( let M=8; M > 1; M-- )
//        for( let N=8; N > 1; N-- )
//          yield [M,N];
//
//        for( let run=1024; run-- > 0; )
//        {
//          const M = randInt(1,64),
//                N = randInt(1,64);
//          yield [M,N];
//        }
//      }
//
//      for( const [M,N] of shapes() )
//        yield _rand_cols0(M,N);
//    }()
//  ).it(
//    'svd_dc accurately decomposes random matrices with zero columns',
//    test_accuracy
//  );
  });
}
