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
import {eps} from '../dt'
import {array, NDArray} from '../nd_array'
import {tabulate} from '../tabulate'
import {_rand_rankdef} from '../_test_data_generators'
import {zip_elems} from '../zip_elems'
import math from '../math'

import {diag_mat} from './diag'
import {matmul, matmul2} from './matmul'
import {eye} from './eye'
import {solve} from './solve'
import {norm} from './norm'
import {rand_ortho} from './rand_ortho'

import {svd_rank,
        svd_decomp,
        svd_solve,
        svd_lstsq} from './svd'
import {svd_dc                } from './svd_dc'
import {svd_jac_2sided        } from './svd_jac_2sided'
import {svd_jac_2sided_blocked} from './svd_jac_2sided_blocked'
import {svd_jac_classic       } from './svd_jac_classic'

import {_svd_jac_post} from './_svd_jac_utils'
import {_transpose_inplace} from './transpose_inplace'


const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from


describe('_svd_jac_utils', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  forEachItemIn(
    function*(){
      for( let run=256; run-- > 0; )
      {
        const N = randInt(1,256);
        const U = rand_ortho('float64',N),
              V = rand_ortho('float64',N),
              S = tabulate([N,N], 'float64', (i,j) => (i===j)*(Math.random()*8-4));
        Object.freeze(U.data.buffer);
        Object.freeze(V.data.buffer);
        Object.freeze(S.data.buffer);
        yield [U,S,V];
      }
    }()
  ).it('_svd_jac_post(...) positivizes and sorts singular values.', ([U,S,V]) => {
    const N = U.shape[0];

    const A = matmul(U,S,V);

    const ord =       Int32Array.from({length: N}, (_,i) => i),
           sv = new Float64Array(N);

    // shuffle ord
    for( let i=ord.length; i > 1; )
    {
      const j = randInt(0,i--);
      [ord[i],ord[j]] =
      [ord[j],ord[i]];
    }

    _transpose_inplace(N, U.data,0);
    _svd_jac_post( N, U.data,S.data,V.data,0, sv,0, ord);

    S = diag_mat(sv);

    for( let i=0; i < N; i++ ) expect(sv[i  ]).not.toBeLessThan(0);
    for( let i=1; i < N; i++ ) expect(sv[i-1]).not.toBeLessThan(sv[i]);

    const a = matmul(U,S,V);

    expect(a).toBeAllCloseTo(A);
  });
})


describe('svd', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
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

        const M = randInt(1,24); shapes[0].push(M,M)
        const N = randInt(1,24); shapes[1].push(M,N)

        yield shapes.map(
          s => tabulate(s,'float64', () => Math.random()*2-1)
        )
      }
    }()
  ).it('svd_solve solves random square examples', ([A,y]) => {
    const [U,sv,V] = svd_decomp(A),
                x  = svd_solve(U,sv,V, y)
  
    expect( matmul2(A,x) ).toBeAllCloseTo(y)
  })


  it('svd_lstsq solves int32 examples', () => {
    for( const [dt1,dt2] of [
      [  'int32','float64'],
      ['float64',  'int32'],
      [  'int32',  'int32']
    ])
    {
      const A = array(dt1,
        [[1,1],
         [1,2],
         [1,3],
         [1,4],
         [1,5]]
      )
      const y = array(dt2, [[2,3,4,5,6]]).T
      Object.freeze(A.data.buffer)
      Object.freeze(y.data.buffer)

      const [U,sv,V] = svd_decomp(A),
                  x  = svd_lstsq(U,sv,V, y)

      expect(x).toBeAllCloseTo([[1],
                                [1]])
    }
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

        const M = randInt(1,24); shapes[0].push(M,M)
        const N = randInt(1,24); shapes[1].push(M,N)

        yield shapes.map(
          s => tabulate(s,'float64', () => Math.random()*2-1)
        )
      }
    }()
  ).it('svd_lstsq solves random square examples', ([A,y]) => {
    const [U,sv,V] = svd_decomp(A),
                x  = svd_lstsq(U,sv,V, y)
  
    expect( matmul2(A,x) ).toBeAllCloseTo(y)
  })


  forEachItemIn(
    function*(){
      for( let run=1024; run-- > 0; )
      {
        let ndim = randInt(0,4),
          shapes = [ Array.from({length: ndim}, () => randInt(1,8)) ]
        shapes.splice( randInt(0,2), 0, shapes[0].slice( randInt(0,ndim) ) )

        for( let d=ndim; d > 0; d-- )
        for( let i=randInt(0,2); i-- > 0; ) { // <- 50/50 chance to collapse dimension to test broadcasting
          const    shape = shapes[randInt(0,2)],
               j = shape.length - d
          if(0<=j) shape[j] = 1
        }

        let M = randInt(1,24),
            N = randInt(1,24)
        if( M < N ) {
          const L=M; M=N; N=L
        }

        const J = randInt(1,32)
        shapes[0].push(M,N)
        shapes[1].push(M,J)

        yield shapes.map(
          s => tabulate(s,'float64', () => Math.random()*2-1)
        )
      }
    }()
  ).it('svd_lstsq solves least squares of random over-determined examples', ([A,y]) => {
    const [U,sv,V]= svd_decomp(A),
                x = svd_lstsq(U,sv,V, y),
               Ax = matmul2(A,x)

    // every least square solution satisfies the normal equaltion Aᵀ(Ax - y) = 0
    expect( matmul2(A.T, zip_elems([Ax,y], math.sub) ) ).toBeAllCloseTo(0)
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
  ).it('svd_lstsq computes one solution of random under-determined examples', ([A,y]) => {
    const [U,sv,V]= svd_decomp(A),
                x = svd_lstsq(U,sv,V, y);

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
  ).it('svd_lstsq solves least squares for random rank-deficient examples', ([A,y]) => {
    const [U,sv,V]= svd_decomp(A),
                x = svd_lstsq(U,sv,V, y),
               Ax = matmul2(A,x)
  
    expect( matmul2(A.T, zip_elems([Ax,y], math.sub) ) ).toBeAllCloseTo(0)
  })

  const svd_decomps = {
    svd_decomp,
    svd_dc,
    svd_jac_2sided,
    svd_jac_2sided_blocked,
    svd_jac_classic
  }

  for( const [svd_name,svd_deco] of Object.entries(svd_decomps) )
    it(`${svd_name} decomposes and int32 input correctly`, () => {
      const A = array('int32',
        [[1,1],
         [1,2],
         [1,3],
         [1,4],
         [1,5]]
      )
      Object.freeze(A.data.buffer)

      const [U,sv,V] = svd_deco(A)

      expect(sv.dtype).toMatch(/^float/)
      expect(U.dtype).toBe(sv.dtype)
      expect(V.dtype).toBe(sv.dtype)

      const D = diag_mat(sv)

      for( let i=1; i < A.shape[1]; i++ )
        expect( sv(i-1) ).not.toBeLessThan( sv(i) );

      const a = matmul(U,D,V);

      const U_TOL = eps(A.dtype) * A.shape[0] * 4,
            V_TOL = eps(A.dtype) * A.shape[1] * 4;

      const I = eye(2)
      expect( matmul2(U.T,U) ).toBeAllCloseTo(I, {rtol:0, atol:U_TOL})
      expect( matmul2(V.T,V) ).toBeAllCloseTo(I, {rtol:0, atol:V_TOL})
      expect( matmul2(V,V.T) ).toBeAllCloseTo(I, {rtol:0, atol:V_TOL})

      expect(D).toBeDiagonal()
      expect(a).toBeAllCloseTo(A)
    })

  for( const [svd_name,svd_deco] of Object.entries(svd_decomps) )
    forEachItemIn(
      function*(){
        for( let run=1024; run-- > 0; )
        {
          const ndim = randInt(0,4),
               shape = Array.from({ length: ndim }, () => randInt(1,8) )
          shape.push(
            randInt(1,24),
            randInt(1,24)
          )
          const A = tabulate(shape, 'float64', () => Math.random() < 0.1 ? 0 : Math.random()*2-1),
             [M,N]= shape.slice(-2)

          // CREATE SOME RANK DEFICIENCIES
          if( Math.random() < 0.5 ) {
            const a = A.reshape(-1,M,N);
            for( let k=a.shape[0]; k-- > 0; )
            for( let i=0; i < M; i++ )
              if( Math.random() < 0.1 ) {
                const l = randInt(0,M),
                  scale = Math.random()*4 - 2;
                for( let j=0; j < N; j++ ) a.set( [k,i,j], scale*a(k,l,j) );
              }
          }
          Object.freeze(A.data.buffer)
          yield A
        }
      }()
    ).it(`${svd_name} works on generated examples`, A => {
      const [M,N]= A.shape.slice(-2), L = Math.min(M,N),
         [U,SV,V]= svd_deco(A),       D = diag_mat(SV)

      for( const sv of SV.reshape(-1,L) )
      for( let i=1; i < sv.shape[0]; i++ )
      {
        expect(sv(i-1)).not.toBeLessThan(sv(i))
        expect(sv(i  )).not.toBeLessThan(-0.0)
      }

      const a = matmul(U,D,V);

      const U_TOL = eps(A.dtype) * M*4,
            V_TOL = eps(A.dtype) * N*4;

      if( M >= N ) {
        const I = eye(N)
        expect( matmul2(U.T,U) ).toBeAllCloseTo(I, {rtol:0, atol:U_TOL})
        expect( matmul2(V.T,V) ).toBeAllCloseTo(I, {rtol:0, atol:V_TOL})
        expect( matmul2(V,V.T) ).toBeAllCloseTo(I, {rtol:0, atol:V_TOL})
      }
      else {
        const I = eye(M)
        expect( matmul2(U.T,U) ).toBeAllCloseTo(I, {rtol:0, atol:U_TOL})
        expect( matmul2(U,U.T) ).toBeAllCloseTo(I, {rtol:0, atol:U_TOL})
        expect( matmul2(V,V.T) ).toBeAllCloseTo(I, {rtol:0, atol:V_TOL})
      }
      expect(D).toBeDiagonal()
      expect(a).toBeAllCloseTo(A)
    })


  for( const [svd_name,svd_deco] of Object.entries(svd_decomps) )
    forEachItemIn(
      function*(){
        function* sizes() {
          for( let N=1; N < 16; N++ )
            yield N;

          const steps_per_binade = 4;

          for( let run=4*steps_per_binade; run <= 8*steps_per_binade; run++ )
            yield Math.round(2**(run/steps_per_binade))
        }

        for( const N of sizes() )
        {
          const shape = [N,N],
                    A = tabulate(shape, 'float64', () => Math.random() < 0.1 ? 0 : Math.random()*2-1);

          // CREATE SOME RANK DEFICIENCIES
          if( Math.random() < 0.25 ) {
            for( let i=0; i < N; i++ )
              if( Math.random() < 0.01 ) {
                const l = randInt(0,N),
                  scale = Math.random()*4 - 2;
                for( let j=0; j < N; j++ ) A.set( [i,j], scale*A(l,j) );
              }
          }
          Object.freeze(A.data.buffer)
          yield A
        }
      }()
    ).it(`${svd_name} is accurate for generated square matrices`, A => {
      const [M,N]= A.shape,     L = Math.min(M,N),
         [U,sv,V]= svd_deco(A), D = diag_mat(sv);

      for( let i=1; i < sv.shape[0]; i++ )
      {
        expect(sv(i-1)).not.toBeLessThan(sv(i))
        expect(sv(i  )).not.toBeLessThan(-0.0)
      }

      const a = matmul(U,D,V);

      const A_TOL = eps(A.dtype) * Math.max(M,N)*8 * norm(A),
            U_TOL = eps(A.dtype) * M*4,
            V_TOL = eps(A.dtype) * N*4;

      if( M >= N ) {
        const I = eye(N)
        expect( matmul2(U.T,U) ).toBeAllCloseTo(I, {rtol:0, atol:U_TOL})
        expect( matmul2(V.T,V) ).toBeAllCloseTo(I, {rtol:0, atol:V_TOL})
        expect( matmul2(V,V.T) ).toBeAllCloseTo(I, {rtol:0, atol:V_TOL})
      }
      else {
        const I = eye(M)
        expect( matmul2(U.T,U) ).toBeAllCloseTo(I, {rtol:0, atol:U_TOL})
        expect( matmul2(U,U.T) ).toBeAllCloseTo(I, {rtol:0, atol:U_TOL})
        expect( matmul2(V,V.T) ).toBeAllCloseTo(I, {rtol:0, atol:V_TOL})
      }
      expect(D).toBeDiagonal()
      expect( norm(zip_elems([A,a], (x,y) => x-y)) ).not.toBeGreaterThan(A_TOL);
    })


  for( const [svd_name,svd_deco] of Object.entries(svd_decomps) )
    forEachItemIn(
      function*(){
        for( let run=1024; run-- > 0; )
        {
          const N = randInt(1,24),
             ndim = randInt(0,4),
            shape = Array.from({ length: ndim }, () => randInt(1,8) );
          shape.push(N);

          const S = tabulate(shape, 'float64', () => Math.random() < 0.1 ? 0 : Math.random()*2-1),
                A = diag_mat(S);

          for( let i=S.data.length; i-- > 0; )
            S.data[i] = Math.abs(S.data[i]);

          for( let i=S.data.length; (i -= N) >= 0; )
            S.data.subarray(i,i+N).sort( (x,y) => y-x );

          Object.freeze(S.data.buffer)
          Object.freeze(A.data.buffer)
          yield [S,A];
        }
      }()
    ).it(`${svd_name} works on batches of generated diagonal matrices`, ([S,A]) => {
      const  [N] = A.shape.slice(-1),
        [U,SV,V] = svd_deco(A);

      const I = eye(N),
          tol = svd_name.startsWith('svd_jac') ? {rtol:0, atol:0} : {};

      expect(SV).toBeAllCloseTo(S, tol);
      expect( matmul2(U,U.T) ).toBeAllCloseTo(I, tol);
      expect( matmul2(U.T,U) ).toBeAllCloseTo(I, tol);
      expect( matmul2(V,V.T) ).toBeAllCloseTo(I, tol);
      expect( matmul2(V.T,V) ).toBeAllCloseTo(I, tol);
      expect( matmul(U,diag_mat(SV),V) ).toBeAllCloseTo(A, tol);
    })


  for( const [svd_name,svd_deco] of Object.entries(svd_decomps) )
    forEachItemIn(
      function*(){
        for( let run=1024; run-- > 0; )
        {
          const M = randInt(1,32),
                N = randInt(1,32),
            shape = Array.from({length: randInt(0,5) }, () => randInt(1,8));

          yield _rand_rankdef(...shape, M,N);
        }
      }()
    ).it(`${svd_name} + svd_rank works on random rank-deficient examples`, ([A,R]) => {
      const [U,sv,V] = svd_deco(A),
                  r  = svd_rank(sv);

      expect(r.shape).toEqual(A.shape.slice(0,-2))
      expect(R.shape).toEqual(A.shape.slice(0,-2))
      expect(r.dtype).toBe('int32')
      expect(R.dtype).toBe('int32')
      expect(r).toBeAllCloseTo(R, {rtol:0, atol:0})
    })
})
