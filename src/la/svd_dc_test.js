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
import {NDArray} from '../nd_array'
import {diag_mat} from './diag'
import {tabulate} from '../tabulate'
import {zip_elems} from '../zip_elems'
import {matmul, matmul2} from './matmul'
import {eye} from './eye'
import {eps} from '../dt'
import {norm} from './norm'

import {svd_rank,
        svd_decomp,
        svd_solve,
        svd_lstsq} from './svd'
import {svd_jac_2sided        } from './svd_jac_2sided'
import {svd_jac_2sided_blocked} from './svd_jac_2sided_blocked'
import {svd_jac_classic       } from './svd_jac_classic'
import {rand_ortho} from './rand_ortho'
import {_svd_dc_1x2,
        _svd_dc_2x3,
        _svd_dc_neves} from './svd_dc'


describe('svd', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  forEachItemIn(
    function*(){
      for( let run=1024; run-- > 0; )
      {
        const rng = () => Math.random() < 0.1 ? 0 : Math.random()*2-1
        const B = new NDArray(
            Int32Array.of(1,2),
          Float64Array.of( rng(), rng() )
        );
        Object.freeze(B.data.buffer);
        yield B;
      }
    }()
  ).it('_svd_dc_1x2 works random square examples', B => {
    const U_arr = new Float64Array(1),
          B_arr =     Float64Array.of( B(0,0), B(0,1) ),
          V_arr = new Float64Array(4);

    _svd_dc_1x2(2, U_arr,0, B_arr,0, V_arr,0);

    const U = new NDArray(Int32Array.of(1,1), U_arr),
          V = new NDArray(Int32Array.of(2,2), V_arr),
          S = new NDArray(Int32Array.of(1,2), Float64Array.of(B_arr[0], 0) );

    expect(B_arr[0]).not.toBeLessThan(0);

    const I = eye(2);

    expect( matmul2(U,U.T) ).toBeAllCloseTo(1)
    expect( matmul2(U.T,U) ).toBeAllCloseTo(1)

    expect( matmul2(V,V.T) ).toBeAllCloseTo(I)
    expect( matmul2(V.T,V) ).toBeAllCloseTo(I)
  
    expect( matmul(U,S,V) ).toBeAllCloseTo(B)
  })


  forEachItemIn(
    function*(){
      for( let run=1024; run-- > 0; )
      {
        const rng = () => Math.random() < 0.1 ? 0 : Math.random()*2-1
        const B = new NDArray(
            Int32Array.of(2,3),
          Float64Array.of(
               rng(), rng(), 0,
            0, rng(), rng()
          )
        );
        Object.freeze(B.data.buffer);
        yield B;
      }
    }()
  ).it('_svd_dc_2x3 works random square examples', B => {
    const U_arr = new Float64Array(4),
          B_arr =     Float64Array.of(B(0,0), B(0,1), B(1,1), B(1,2)),
          V_arr = new Float64Array(9);

    _svd_dc_2x3(3, U_arr,0, B_arr,0, V_arr,0);

    const U = new NDArray(Int32Array.of(2,2), U_arr),
          V = new NDArray(Int32Array.of(3,3), V_arr),
          S = new NDArray(Int32Array.of(2,3), Float64Array.of(B_arr[0],       0,  0,
                                                                    0,  B_arr[2], 0) );

    expect(B_arr[2]).not.toBeLessThan(0);
    expect(B_arr[0]).not.toBeLessThan(B_arr[2]);

    const I2 = eye(2),
          I3 = eye(3);

    expect( matmul2(U,U.T) ).toBeAllCloseTo(I2)
    expect( matmul2(U.T,U) ).toBeAllCloseTo(I2)

    expect( matmul2(V,V.T) ).toBeAllCloseTo(I3)
    expect( matmul2(V.T,V) ).toBeAllCloseTo(I3)
  
    expect( matmul(U,S,V) ).toBeAllCloseTo(B)
  })


  forEachItemIn(
    function*(){
      for( let run=0; run < 64; run++ )
      for( let N=3; N < 64; N++ )
//      for( let N=3; N < 8; N++ )
      {
        const M = N-1;

        const rng = () => (Math.random()*2-1) * (Math.random() < 0.875),
             diag = Float64Array.from({length: M}, () => Math.random());
        diag[0] = 0;
        if( Math.random() >= 0.875 )
        for( let k=Math.random()*(M-1) | 0; k-- > 0; )
        {
          let i = (Math.random()*(M-1) | 0) + 1,
              j =  Math.random()* M    | 0;
          diag[i] = diag[j];
        }
        if( 0 !== diag[0] ) throw new Error('Assertion failed.');
        diag.sort((x,y) => y-x);

        const B = tabulate([M,N], () => 0);
        for( let i=0; i < M; i++ ) {
          B.set([ i,i], diag[i]);
          B.set([-1,i], rng());
        }
//        B.set([-1,-2], Math.random()*2-1)

        Object.freeze(B.data.buffer);
        yield B;
      }
    }()
  ).it('_svd_dc_neves works random square examples', B => {
    const [M,N] = B.shape.slice(-2);

    const TOL = eps(B.dtype) * 512*1024,
       B_norm = norm(B);

    const B_arr = new Float64Array(2*M);
    for( let i=0; i < M; i++ ) {
                  B_arr[2*i+1] = B(-1,i)
      if(i < M-1) B_arr[2*i  ] = B( i,i)
    }

//    console.log('\n\n');
//    console.log('B:\n' + B.toString())

    const U = eye(M),
          V = eye(N);
    U.data.fill(0);
    V.data.fill(0);

    const [Q,sv,W] = svd_decomp(B);
    Object.freeze( Q.data.buffer);
    Object.freeze(sv.data.buffer);
    Object.freeze( W.data.buffer);

//    console.log('sv:\n' + sv.mapElems(x => x/*.toFixed(6)*/) );

    const ord = Int32Array.from({length: M}, (_,i) => i);

    _svd_dc_neves(N, N, U.data,0, B_arr,0, V.data,0, ord)

//    console.log('{{{{')
//    console.log( 'Q:\n' +  Q.mapElems(x => x.toFixed(6)) );
//    console.log( 'U:\n' +  U.mapElems(x => x.toFixed(6))  + '\n');
//    console.log('sv:\n' + sv.mapElems(x => x/*.toFixed(6)*/) );
//    console.log( 'B:\n' + B_arr.filter((_,i)=>0===i%2).map(x => x/*.toFixed(6)*/) + '\n');
//    console.log( 'W:\n' +  W.mapElems(x => x.toFixed(6)) );
//    console.log( 'V:\n' +  V.mapElems(x => x.toFixed(6)) );
//    console.log('}}}}\n\n\n\n')


    // check orthogonality
    { const I = eye(M);
      expect(matmul2(U,U.T)).toBeAllCloseTo(I);
      expect(matmul2(U.T,U)).toBeAllCloseTo(I); }
    { const I = eye(N);
      expect(matmul2(V,V.T)).toBeAllCloseTo(I);
      expect(matmul2(V.T,V)).toBeAllCloseTo(I); }

    // check decomposition
    const S  = tabulate([M,N], (i,j) => B_arr[2*i]*(i===j)),
         USV = matmul(U,S,V);

    expect(USV).toBeAllCloseTo(B);

//    const  USV_err = norm( zip_elems([USV, B], B.dtype, (x,y) => x-y) );
//    expect(USV_err).not.toBeGreaterThan(B_norm*TOL);

    expect( U.data.every(x => isFinite(x)) ).toBe(true);

    // check singular values
    const  SV = B_arr.filter((_,i) => 0===i%2);
    expect(SV).toBeAllCloseTo(sv);

//    console.log('B:\n' + B.toString());
  })
})
