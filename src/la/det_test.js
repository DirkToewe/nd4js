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
import {zip_elems} from '../zip_elems'
import {det, det_tri, slogdet, slogdet_tri} from './det'
import {matmul2} from './matmul'


describe('det', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  forEachItemIn(
    function*(){
      function* shapes() {
        for( let k=1; k <= 8; k++ ) { yield     [k,k]
        for( let j=1; j <= 8; j++ ) { yield   [j,k,k]
        for( let i=1; i <= 8; i++ ) { yield [i,j,k,k] }}}
      }

      for( const shape of shapes() )
      for( const dtype of ['int32', 'float32', 'float64'] )
      {
        const L = tabulate(shape, dtype, (...idx) => {
                const [i,j] = idx.slice(-2);
                return i < j ? 0 : Math.random()*200 - 100;
              }),
              U = L.T;
        Object.freeze(L.data.buffer); yield L;
        Object.freeze(U.data.buffer); yield U;
      }
    }()
  ).it('det_tri works for random matrices', A => {
    const N = A.shape[A.ndim-1];

    const diag_A = A.reshape(...A.shape.slice(0,-2), N*N)
                    .sliceElems('...', [0, N*N, N+1]),
           det_A = det_tri(A),
           DET_A = diag_A.reduceElems(-1, A.dtype, (x,y) => x*y);

    expect(det_A.dtype).toBe(A.dtype);
    expect(det_A).toBeAllCloseTo(DET_A);
  })


  forEachItemIn(
    function*(){
      function* shapes() {
        for( let k=1; k <= 8; k++ ) { yield     [k,k]
        for( let j=1; j <= 8; j++ ) { yield   [j,k,k]
        for( let i=1; i <= 8; i++ ) { yield [i,j,k,k] }}}
      }

      for( const shape of shapes() )
      for( const dtype of ['float32', 'float64'] )
      {
        const L = tabulate(shape, dtype, (...idx) => {
                const [i,j] = idx.slice(-2);
                return i < j ? 0 : Math.random()*200 - 100;
              }),
              U = L.T;
        Object.freeze(L.data.buffer); yield L;
        Object.freeze(U.data.buffer); yield U;
      }
    }()
  ).it('slogdet_tri works for random matrices', A => {
    const N = A.shape[A.ndim-1];

    const diag_A = A.reshape(...A.shape.slice(0,-2), N*N)
                    .sliceElems('...', [0, N*N, N+1]),
           [S,L] = slogdet_tri(A),
               s = diag_A.   mapElems(       A.dtype,  x    => Math.sign(x) )
                         .reduceElems( -1,   A.dtype, (x,y) => x*y ),
               l = diag_A.   mapElems(      'float64', x    => Math.log(Math.abs(x)) )
                         .reduceElems( -1,  'float64',(x,y) => x+y ),
           det_A =          zip_elems([S,L],'float64',(s,l) => s*Math.exp(l)),
           DET_A = diag_A.reduceElems( -1,   A.dtype, (x,y) => x*y);

    expect(S).toBeAllCloseTo(s);
    expect(L).toBeAllCloseTo(l);
    expect(det_A).toBeAllCloseTo(DET_A);
  })


  forEachItemIn(
    function*(){
      function* shapes() {
        for( let k=1; k <= 8; k++ ) { yield     [k,k]
        for( let j=1; j <= 8; j++ ) { yield   [j,k,k]
        for( let i=1; i <= 8; i++ ) { yield [i,j,k,k] }}}
      }

      for( const shape of shapes() )
      for( const dtype of ['int32', 'float32', 'float64'] )
      {
        const L = tabulate(shape, dtype, (...idx) => {
                const [i,j] = idx.slice(-2);
                return i < j ? 0 : Math.random()*2 - 1;
              }),
              U = tabulate(shape, dtype, (...idx) => {
                const [i,j] = idx.slice(-2);
                return i > j ? 0 : Math.random()*2 - 1;
              }),
              LU = matmul2(L,U);
        Object.freeze(L .data.buffer);
        Object.freeze( U.data.buffer);
        Object.freeze(LU.data.buffer);
        yield [L,U, LU];
      }
    }()
  ).it('det works for random matrices', ([L,U,A]) => {
    const tol = A.dtype==='float32' ? {rtol: 1e-3, atol: 1e-4} : {},
        det_A = zip_elems([det_tri(L),
                           det_tri(U)], (x,y) => x*y);
    expect( det(A) ).toBeAllCloseTo(det_A, tol);
  })


  forEachItemIn(
    function*(rng){
      function* shapes() {
        for( let k=1; k <= 8; k++ ) { yield     [k,k]
        for( let j=1; j <= 8; j++ ) { yield   [j,k,k]
        for( let i=1; i <= 8; i++ ) { yield [i,j,k,k] }}}
      }

      for( const shape of shapes() )
      for( const dtype of ['float64'] )
      {
        const L = tabulate(shape, dtype, (...idx) => {
                const [i,j] = idx.slice(-2);
                return i < j ? 0 : rng.uniform(-2,+2);
              }),
              U = tabulate(shape, dtype, (...idx) => {
                const [i,j] = idx.slice(-2);
                return i > j ? 0 : rng.uniform(-2,+2);
              }),
              LU = matmul2(
                L.mapElems('float64'),
                U.mapElems('float64')
              ).mapElems(dtype);
        Object.freeze(L .data.buffer);
        Object.freeze( U.data.buffer);
        Object.freeze(LU.data.buffer);
        yield [L,U, LU];
      }
    }
  ).it('slogdet works for random matrices', ([L,U,A]) => {
    const [s1,l1] = slogdet_tri(L),
          [s2,l2] = slogdet_tri(U),
                s = zip_elems([s1,s2], (x,y) => x*y),
                l = zip_elems([l1,l2], (x,y) => x+y),
            [S,D] = slogdet(A);

    expect(S).toBeAllCloseTo(s);
    expect(D).toBeAllCloseTo(l);
  })
})
