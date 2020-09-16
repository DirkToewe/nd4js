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

import {diag,
        diag_mat} from './diag'
import {eye} from './eye';
import {matmul, matmul2} from './matmul'
import {_svd_jac_post,
        _svd_jac_post_skip1} from './_svd_jac_utils'
import {_transpose_inplace} from './transpose_inplace'


describe('_svd_jac_utils', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  const test_svd_jac_post = ([U,sv,V, ord]) => {
    const     N = sv.length,
      sv_sorted = sv.map(s => Math.abs(s)).sort((x,y) => y-x);

    let                S = diag_mat(sv);
    const A = matmul(U,S,V);

    _transpose_inplace(N, U.data,0);
    _svd_jac_post( N, U.data,S.data,V.data,0, sv,0, ord);

    S = diag_mat(sv);

    expect( ord.slice().sort((x,y) => x-y) ).toEqual( Int32Array.from({length: N}, (_,i) => i) );

    expect(sv).toEqual(sv_sorted);

    const  a = matmul(U,S,V);
    expect(a).toBeAllCloseTo(A);

    const U_T = U.T,
          V_T = V.T,                        I = eye(N);
    expect( matmul2(U,U_T) ).toBeAllCloseTo(I);
    expect( matmul2(U_T,U) ).toBeAllCloseTo(I);
    expect( matmul2(V,V_T) ).toBeAllCloseTo(I);
    expect( matmul2(V_T,V) ).toBeAllCloseTo(I);
  };


  forEachItemIn(
    function*(rng){
      for( let run=173; run-- > 0; )
      {
        const             N = rng.int(1,96),
            U = rng.ortho(N),
            V = rng.ortho(N),
            sv = Float64Array.from({length: N},    () => rng.uniform(-4,+4)),
            ord =  Int32Array.from({length: N}, (_,i) => i);
        rng.shuffle(ord);
        Object.freeze(U.data.buffer);
        Object.freeze(V.data.buffer);
        Object.freeze(    sv.buffer);
        Object.freeze(   ord.buffer);
        yield [U,sv,V, ord];
      }
    }
  ).it('_svd_jac_post sign-flips and sorts singular values given random examples.', test_svd_jac_post);


  forEachItemIn(
    function*(rng){
      for( let run=256; run-- > 0; )
      {
        const             N = rng.int(1,96),
            U = rng.ortho(N),
            V = rng.ortho(N),
            sv = Float64Array.from({length: N},    () => rng.uniform(-4,+4)*(rng.uniform(0,1) < 0.9) ),
            ord =  Int32Array.from({length: N}, (_,i) => i);
        rng.shuffle(ord);
        Object.freeze(U.data.buffer);
        Object.freeze(V.data.buffer);
        Object.freeze(    sv.buffer);
        Object.freeze(   ord.buffer);
        yield [U,sv,V, ord];
      }
    }
  ).it('_svd_jac_post sign-flips and sorts singular values given random sparse examples.', test_svd_jac_post);
});
