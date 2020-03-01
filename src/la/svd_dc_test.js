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
import {_rand_rows0,
        _rand_cols0} from '../_test_data_generators'
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
import {svd_dc,
       _svd_dc_1x2,
       _svd_dc_2x3,
       _svd_dc_neves,
       _svd_dc_bidiag} from './svd_dc'


describe('svd_dc', () => {
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
         BV_arr =     Float64Array.of( B(0,0), B(0,1), 0,0,0,0 );

    _svd_dc_1x2(2, U_arr,0, BV_arr,0,2);

    const U = new NDArray(Int32Array.of(1,1), U_arr),
          V = new NDArray(Int32Array.of(2,2), BV_arr.slice(2)),
          S = new NDArray(Int32Array.of(1,2), Float64Array.of(BV_arr[0], 0) );

    expect(BV_arr[0]).not.toBeLessThan(0);

    const I = eye(2);

    expect( matmul2(U,U.T) ).toBeAllCloseTo(1)
    expect( matmul2(U.T,U) ).toBeAllCloseTo(1)

    expect( matmul2(V,V.T) ).toBeAllCloseTo(I)
    expect( matmul2(V.T,V) ).toBeAllCloseTo(I)
  
    expect( matmul(U,S,V.T) ).toBeAllCloseTo(B)
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
         BV_arr =     Float64Array.of(B(0,0), B(0,1), B(1,1), B(1,2), ...new Float64Array(9));

    _svd_dc_2x3(3, U_arr,0, BV_arr,0,4);

    const U = new NDArray(Int32Array.of(2,2), U_arr),
          V = new NDArray(Int32Array.of(3,3), BV_arr.slice(4)),
          S = new NDArray(Int32Array.of(2,3), Float64Array.of(BV_arr[0],        0,  0,
                                                                     0,  BV_arr[2], 0) );

    expect(BV_arr[2]).not.toBeLessThan(0);
    expect(BV_arr[0]).not.toBeLessThan(BV_arr[2]);

    const I2 = eye(2),
          I3 = eye(3);

    expect( matmul2(U,U.T) ).toBeAllCloseTo(I2)
    expect( matmul2(U.T,U) ).toBeAllCloseTo(I2)

    expect( matmul2(V,V.T) ).toBeAllCloseTo(I3)
    expect( matmul2(V.T,V) ).toBeAllCloseTo(I3)
  
    expect( matmul(U,S,V.T) ).toBeAllCloseTo(B)
  })


  forEachItemIn(
    function*(){
      for( let run=0; run < 256; run++ )
      for( let N=2; N < 32; N++ )
      {
        const M = N-1,
            mid = M >>> 1;

        const rng = () => (Math.random()*2-1) * (Math.random() < 0.875),
             diag = Float64Array.from({length: M}, () => Math.random());
        diag[0] = 0;

        if( Math.random() >= 0.9375 )
          // create duplicates on diagonal
          for( let k=Math.random()*(M-1) | 0; k-- > 0; )
          {
            let i = (Math.random()*(M-1) | 0) + 1,
                j =  Math.random()* M    | 0;
            diag[i] = diag[j];
          }

        if( Math.random() >= 0.9375 )
          // create near-duplicates on diagonal
          for( let k=Math.random()*(M-1) | 0; k-- > 0; )
          {
            let i = (Math.random()*(M-1) | 0) + 1,
                j =  Math.random()* M    | 0;

            const factor = 1 + Number.EPSILON*Math.random()*(200 - 100);

            diag[i] = diag[j] * factor;

            if( 0 === diag[i] ) diag[i] = Number.MIN_VALUE*Math.random()*(200 - 100);
          }

        if(      0 !== diag[0] ) throw new Error('Assertion failed.');
        const diag_0 = diag[0];
                       diag[0] = diag[mid];
                                 diag[mid] = diag_0;

        const A = tabulate([M,N], 'float64', () => 0);
        for( let i=0; i < M; i++ ) {
          A.set([  i,i], diag[i]);
          A.set([mid,i], rng());
        }

        const L = rand_ortho('float64', M),
              R = rand_ortho('float64', N);

        Object.freeze(L.data.buffer);
        Object.freeze(A.data.buffer);
        Object.freeze(R.data.buffer);
        yield [L,A,R];
      }
    }()
  ).it('_svd_dc_neves works for random square examples', ([L,A,R]) => {
    const [M,N] = A.shape.slice(-2),
           mid  = M >>> 1;

    expect(L.dtype).toBe('float64')
    expect(A.dtype).toBe('float64')
    expect(R.dtype).toBe('float64')

    // check orthogonality
    { const I = eye(M);
      expect(matmul2(L,L.T)).toBeAllCloseTo(I);
      expect(matmul2(L.T,L)).toBeAllCloseTo(I); }
    { const I = eye(N);
      expect(matmul2(R,R.T)).toBeAllCloseTo(I);
      expect(matmul2(R.T,R)).toBeAllCloseTo(I); }

    const LAR = matmul(L,A,R);
    expect(LAR.dtype).toBe('float64');

    const I =       Int32Array.from({length: 3*M}, (_,i) => i),
          F = new Float64Array( M*(M+2) + 2*M + (M+1)*(M+1) );
    I[M-1] = mid;
    I[mid] = M-1;
    I.subarray(0,M-1).sort((i,j) => A(j,j) - A(i,i));

    const [B_off,V_off] = function(){
      switch(Math.random()*2 | 0) {
        case 0: return [M*(M+2)              , M*(M+2) + 2*M];
        case 1: return [M*(M+2) + (M+1)*(M+1), M*(M+2)      ];
        default: throw new Error('Unexpected fall-through');
      }
    }();
    
    for( let i=0; i < M; i++ )
    { const  j = I[i];
      F[B_off + 2*i  ] = A(  j,j);
      F[B_off + 2*i+1] = A(mid,j);
    }
    F[B_off + 2*M-2] = 0;

    for( let i=0; i <= M; i++ )
    for( let j=0; j <= M; j++ )
      F[V_off + (M+1)*i+j] = R(j,i);

    expect(
      I.slice(0,M).sort()
    ).toEqual(
      Int32Array.from({length: M}, (_,i) => i)
    );

    const U = L.mapElems();
    expect(U.dtype).toBe('float64');

    _svd_dc_neves(N, N, U.data,0, F,B_off,V_off, I);
    Object.freeze(U.data.buffer);
    Object.freeze(F     .buffer);

    const V = new NDArray( Int32Array.of(M+1,M+1), F.slice(V_off,
                                                           V_off + (M+1)*(M+1)) );
    Object.freeze(V.data.buffer);

    // check orthogonality
    { const I = eye(M);
      expect(matmul2(U,U.T)).toBeAllCloseTo(I);
      expect(matmul2(U.T,U)).toBeAllCloseTo(I); }
    { const I = eye(N);
      expect(matmul2(V,V.T)).toBeAllCloseTo(I);
      expect(matmul2(V.T,V)).toBeAllCloseTo(I); }

    // check decomposition
    const S  = tabulate([M,N], (i,j) => F[B_off + 2*i]*(i===j)),
         USV = matmul(U,S,V.T);

    expect(USV).toBeAllCloseTo(LAR);

    const sv = F.subarray(B_off,B_off+2*M).filter((_,i) => i%2 === 0);

    // check order of singular values
    for( let i=1; i < M; i++ )
      expect( sv[i-1] ).not.toBeLessThan( sv[i] );
    expect( sv[M-1] ).not.toBeLessThan(0);
  })


  forEachItemIn(
    function*(){
      for( let run=0; run < 128; run++ )
      for( let N=2; N < 48; N++ )
      {
        const M = N-1;

        const rng = () => (Math.random()*2-1) * (Math.random() < 0.875);

        const A = tabulate([M,N], 'float64', () => 0);
        for( let i=0; i < M; i++ ) {
          A.set( [i,i  ], rng() );
          A.set( [i,i+1], rng() );
        }

        Object.freeze(A.data.buffer);
        yield A;
      }
    }()
  ).it('_svd_dc_bidiag works for random examples', A => {
    const [M,N] = A.shape.slice(-2);

    const U = tabulate([M,M], 'float64', () => 0),
          I = new   Int32Array(M*3),
          F = new Float64Array(M*(M+2) + 2*M + (M+1)*(M+1));

    const [B_off,V_off] = function(){
      switch(Math.random()*2 | 0) {
        case 0: return [M*(M+2)              , M*(M+2) + 2*M];
        case 1: return [M*(M+2) + (M+1)*(M+1), M*(M+2)      ];
        default: throw new Error('Unexpected fall-through');
      }
    }();

    for( let i=0; i < M; i++ ) {
      F[B_off + 2*i  ] = A(i,i);
      F[B_off + 2*i+1] = A(i,i+1);
    }

    _svd_dc_bidiag( N,N, U.data,0, F,B_off,V_off, I);

    const S = tabulate([M,N], (i,j) => i !== j ? 0 : F[B_off + 2*i]),
          V = new NDArray( Int32Array.of(M+1,M+1), F.slice(V_off,
                                                           V_off + (M+1)*(M+1)) );

    // check orthogonality
    { const I = eye(M);
      expect(matmul2(U,U.T)).toBeAllCloseTo(I);
      expect(matmul2(U.T,U)).toBeAllCloseTo(I); }
    { const I = eye(N);
      expect(matmul2(V,V.T)).toBeAllCloseTo(I);
      expect(matmul2(V.T,V)).toBeAllCloseTo(I); }

    const USV = matmul(U,S,V.T);

    expect(USV).toBeAllCloseTo(A);
  })


  const test_accuracy = A =>
  {
    expect(A.ndim).toBe(2);

    const [M,N]= A.shape,   L = Math.min(M,N),
       [U,sv,V]= svd_dc(A), D = diag_mat(sv);

    for( let i=1; i < sv.shape[0]; i++ )
    {
      expect(sv(i-1)).not.toBeLessThan(sv(i))
      expect(sv(i  )).not.toBeLessThan(-0.0)
    }

    const a = matmul(U,D,V);

    const A_TOL = eps(A.dtype) * 4 * Math.max(M,N) * norm(A),
          U_TOL = eps(A.dtype) * 4 * M,
          V_TOL = eps(A.dtype) * 4 * N;

    const I = eye( Math.min(M,N) );
    if( M >= N ) {
      expect( matmul2(U.T,U) ).toBeAllCloseTo(I, {rtol:0, atol:U_TOL})
      expect( matmul2(V.T,V) ).toBeAllCloseTo(I, {rtol:0, atol:V_TOL})
    }
    if( M <= N ) {
      expect( matmul2(U,U.T) ).toBeAllCloseTo(I, {rtol:0, atol:U_TOL})
      expect( matmul2(V,V.T) ).toBeAllCloseTo(I, {rtol:0, atol:V_TOL})
    }
    expect(D).toBeDiagonal()

    const        diff = zip_elems([A,a], (x,y) => x-y);
    expect( norm(diff) ).not.toBeGreaterThan(A_TOL);
  };


  for( const dr of [0,1] )
  for( const dc of [0,1] )
  if( 1 !== dr*dc )
  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from;

      function* sizes() {
        const steps_per_binade = 4;

        for( let run=0; run < 16; run++ )
        {
          for( let N=1; N < 16; N++ )
            yield N;

          for( let run=4*steps_per_binade; run <= 8*steps_per_binade; run++ )
            yield Math.round(2**(run/steps_per_binade))
        }
      }

      for( const L of sizes() )
      {
        const M = L+dr,
              N = L+dc;
        const shape = [M,N],
                  A = tabulate(shape, 'float64', () => Math.random() < 0.1 ? 0 : Math.random()*2-1);

        // CREATE SOME RANK DEFICIENCIES
        if( Math.random() < 0.25 ) {
          for( let i=0; i < M; i++ )
            if( Math.random() < 0.01 ) {
              const l = randInt(0,M),
                scale = Math.random()*4 - 2;
              for( let j=0; j < N; j++ ) A.set( [i,j], scale*A(l,j) );
            }
        }
        Object.freeze(A.data.buffer)
        yield A
      }
    }()
  ).it(
    `svd_dc is accurate for random matrices with aspect ratio [N+${dr},N+${dc}]`,
    test_accuracy
  );


  for( const [rng,suffix] of [
    [() =>                           Math.random()*8 - 4, ''                      ],
    [() => Math.random() < 0.1 ? 0 : Math.random()*8 - 4, ' with occasional zeros']
  ])
    forEachItemIn(
      function*(){
        const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from;

        function* shapes() {
          const steps_per_binade = 4;

          for( let M=8; M > 1; M-- )
          for( let N=8; N > 1; N-- )
            yield [M,N];

          for( let run=1024; run-- > 0; )
          {
            const M = randInt(1,64),
                  N = randInt(1,64);
            yield [M,N];
          }
        }

        for( const [M,N] of shapes() )
        {
          const A = tabulate([M,N], 'float64', rng);
          Object.freeze(A.data.buffer)
          yield A
        }
      }()
    ).it(
      'svd_dc is accurate for random matrices'+suffix,
      test_accuracy
    );


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from;

      function* shapes() {
        const steps_per_binade = 4;

        for( let M=8; M > 1; M-- )
        for( let N=8; N > 1; N-- )
          yield [M,N];

        for( let run=1024; run-- > 0; )
        {
          const M = randInt(1,64),
                N = randInt(1,64);
          yield [M,N];
        }
      }

      for( const [M,N] of shapes() )
      {
        const                         sparseness = Math.random(),
          rng = () => Math.random() < sparseness ? 0 : Math.random()*8 - 4;

        const A = tabulate([M,N], 'float64', rng);
        Object.freeze(A.data.buffer)
        yield A
      }
    }()
  ).it(
    'svd_dc is accurate for random sparse matrices',
    test_accuracy
  );


// FIXME make the following test pass:
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
//    'svd_dc is accurate for random matrices with zero rows',
//    test_accuracy
//  );


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from;

      function* shapes() {
        for( let M=8; M > 1; M-- )
        for( let N=8; N > 1; N-- )
          yield [M,N];

        for( let run=1024; run-- > 0; )
        {
          const M = randInt(1,64),
                N = randInt(1,64);
          yield [M,N];
        }
      }

      for( const [M,N] of shapes() )
        yield _rand_cols0(M,N);
    }()
  ).it(
    'svd_dc is accurate for random matrices with zero columns',
    test_accuracy
  );
})
