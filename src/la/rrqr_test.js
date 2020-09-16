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
import {generic_test_lstsq} from "./_generic_test_lstsq";
import {generic_test_rank } from "./_generic_test_rank";
import {generic_test_solve} from "./_generic_test_solve";
import math from '../math'
import {tabulate} from '../tabulate'

import {eye} from './eye'
import {matmul2} from './matmul'
import {permute_cols} from './permute'
import { rrqr_decomp,
         rrqr_decomp_full,
        _rrqr_decomp_inplace,
         rrqr_rank,
         rrqr_solve,
         rrqr_lstsq}       from  './rrqr'
import {srrqr_decomp_full} from './srrqr'


for( const lstsq of Object.values({
  ['srrqr_decomp_full + rrqr_lstsq   ']: (A,y) => rrqr_lstsq(   srrqr_decomp_full(A)           , y ),
  ['srrqr_decomp_full + rrqr_lstsq...']: (A,y) => rrqr_lstsq(...srrqr_decomp_full(A).slice(0,3), y ),
  [' rrqr_decomp      + rrqr_lstsq   ']: (A,y) => rrqr_lstsq(    rrqr_decomp     (A)           , y ),
  [' rrqr_decomp      + rrqr_lstsq...']: (A,y) => rrqr_lstsq(... rrqr_decomp     (A)           , y ),
  [' rrqr_decomp_full + rrqr_lstsq   ']: (A,y) => rrqr_lstsq(    rrqr_decomp_full(A)           , y ),
  [' rrqr_decomp_full + rrqr_lstsq...']: (A,y) => rrqr_lstsq(... rrqr_decomp_full(A)           , y )
}))
  generic_test_lstsq(lstsq, 'some solution');


for( const rank of Object.values({
  ['srrqr_decomp_full(_)[3]'      ]: A =>            srrqr_decomp_full(A)[3],
  ['srrqr_decomp_full + rrqr_rank']: A => rrqr_rank( srrqr_decomp_full(A)[1] ),
  [' rrqr_decomp      + rrqr_rank']: A => rrqr_rank(  rrqr_decomp     (A)[1] ),
  [' rrqr_decomp_full + rrqr_rank']: A => rrqr_rank(  rrqr_decomp_full(A)[1] )
}))
  generic_test_rank(rank);


for( const solve of Object.values({
  ['srrqr_decomp_full + rrqr_solve   ']: (A,y) => rrqr_solve(   srrqr_decomp_full(A)           , y ),
  ['srrqr_decomp_full + rrqr_solve...']: (A,y) => rrqr_solve(...srrqr_decomp_full(A).slice(0,3), y ),
  [' rrqr_decomp      + rrqr_solve   ']: (A,y) => rrqr_solve(    rrqr_decomp     (A)           , y ),
  [' rrqr_decomp      + rrqr_solve...']: (A,y) => rrqr_solve(... rrqr_decomp     (A)           , y ),
  [' rrqr_decomp_full + rrqr_solve   ']: (A,y) => rrqr_solve(    rrqr_decomp_full(A)           , y ),
  [' rrqr_decomp_full + rrqr_solve...']: (A,y) => rrqr_solve(... rrqr_decomp_full(A)           , y )
}))
  generic_test_solve(solve);


describe('rrqr', () => {

  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  });

  // TODO: test rrqr_solve, rrqr_lstsq and rrqr_rank without using rrqr_decomp

})


describe('_rrqr_decomp_inplace', () => {

  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  });

  forEachItemIn(
    function*(rng){
      for( let run=733; run-- > 0; )
      {
        const M  = rng.int(1,64),
              N  = rng.int(1,64),
              L  = rng.int(1,64),
             [A] = rng.rankDef(M,N),
              Y  =   tabulate([M,L], 'float64', () => rng.uniform(-2,+2) *  (rng.uniform(0,1) < 0.95) );
        Object.freeze(Y.data.buffer);
        Object.freeze(Y);
        yield [A,Y];
      }
    }
  ).it(`_rrqr_decomp_inplace works on random rank-deficient examples`, ([A,Y]) => {
    const [M,N] = A.shape,
             L  = Y.shape[1];

    const [Q,R,P] = rrqr_decomp_full(A);

    const a = A.mapElems(),
          y = Y.mapElems(),
          p =   Int32Array.from({length: N  }, (_,i) => i),
       norm = Float64Array.from({length: N*2}, (_,i) => NaN);

    _rrqr_decomp_inplace(M,N,L, a.data,0, y.data,0, p,0, norm); // <- TODO test with offset

    expect(p).toBeAllCloseTo(P, {rtol: 0, atol: 0});
    expect(a).toBeAllCloseTo(R);
    expect(y).toBeAllCloseTo( matmul2(Q.T,Y) );
  })
});


for( const rrqr_deco of [
  srrqr_decomp_full,
  rrqr_decomp,
  rrqr_decomp_full
])
  describe(`${rrqr_deco.name} [generic RRQR tests]`, () => {
    beforeEach( () => {
      jasmine.addMatchers(CUSTOM_MATCHERS)
    })


    const IS_STRONG = rrqr_deco.name.startsWith('srrqr_');


    const test_body = A => {
      const [M,N] = A.shape.slice(-2),
          [Q,R,P] = rrqr_deco(A),
               L  = rrqr_deco.name.endsWith('_full') ? M : Math.min(M,N),
               a  = matmul2(Q,R)
      Object.freeze(Q.data.buffer); expect(Q.dtype).toBe('float64')
      Object.freeze(R.data.buffer); expect(R.dtype).toBe('float64')
      Object.freeze(P.data.buffer); expect(P.dtype).toBe(  'int32')
      Object.freeze(a.data.buffer)

      expect( Q.shape ).toEqual( Int32Array.of(...A.shape.slice(0,-2), M, L) )
      expect( R.shape ).toEqual( Int32Array.of(...A.shape.slice(0,-2), L, N) )
      expect( P.shape ).toEqual( Int32Array.of(...A.shape.slice(0,-2),    N) )

      A = permute_cols(A,P);

      expect(R).toBeUpperTriangular(IS_STRONG ? {tol: 1e-8} : {});
      expect(a).toBeAllCloseTo(A);

      const I = eye(L),
          Q_T = Q.T;
      Object.freeze(I.data.buffer)
                  expect( matmul2(Q_T,Q) ).toBeAllCloseTo(I)
      if( L===M ) expect( matmul2(Q,Q_T) ).toBeAllCloseTo(I)

      if( ! IS_STRONG )
      {
        const RR = R.data;
        for( let off=RR.length; (off -= L*N) >= 0; )
          for( let i=Math.min(M,N); --i > 0; )
          {
            const  R_II = math.abs(RR[off + N* i   + i   ]),
                   R_ii = math.abs(RR[off + N*(i-1)+(i-1)])
            expect(R_ii).toBeGreaterThanOrEqual(R_II)
          }
      }
    };


    forEachItemIn(
      function*(rng){
        for( let run=256; run-- > 0; ) { const                        ndim      =  rng.int(2, 5),
                                     shape = Int32Array.from({length: ndim}, () => rng.int(1,24) )
          const         A = tabulate(shape, 'float64', () => rng.uniform(-2,+2) * (rng.uniform(0,1) < 0.9) );
          Object.freeze(A.data.buffer);
          yield         A;
        }
      }
    ).it(`properly decomposes random examples`, test_body);


    forEachItemIn(
      function*(rng){
        for( let run=173; run-- > 0; ) { const                               ndim      =  rng.int(2, 5),
                                            shape = Int32Array.from({length: ndim}, () => rng.int(1,24) );
          const        [A] = rng.rankDef(...shape);
          Object.freeze(A.data.buffer);
          yield         A;
        }
      }
    ).it(`properly decomposes random rank-deficient examples`, test_body);

  })
