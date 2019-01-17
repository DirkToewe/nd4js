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
import math from '../math'
import {eye} from './eye'
import {qr_decomp} from './qr'

import {svd_rank,
        svd_decomp,
        svd_solve,
        svd_lstsq} from './svd'
import {svd_jac_1sided } from './svd_jac_1sided'
import {svd_jac_2sided } from './svd_jac_2sided'
import {svd_jac_classic} from './svd_jac_classic'


describe('svd', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

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


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

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
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

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
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

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
            N = randInt(1,24)
        if( M > N ) {
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
  ).it('svd_lstsq computes one solution of random under-determined examples', ([A,y]) => {
    const [U,sv,V]= svd_decomp(A),
                x = svd_lstsq(U,sv,V, y)
  
    expect( matmul2(A,x) ).toBeAllCloseTo(y)
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
            J = randInt(1,32),
            L = Math.min(M,N)

        const A = tabulate(shapes[0], () => {
          const r  = randInt(0,L),
             [Q,_] = qr_decomp( tabulate([M,N], 'float64', () => Math.random()*2-1) ),
                R  = tabulate([L,N], 'float64', (i,j) => {
                       if(r <= i || j < i) return 0
                       if(j == i) return (0.5 + Math.random()) * (randInt(0,2)*2 - 1)
                       return Math.random()*0.5 - 0.25
                     })
          return matmul2(Q,R)
        })

        shapes[0].push(M,N)
        shapes[1].push(M,J)

        yield [
          new NDArray(
              Int32Array.from(shapes[0]),
            Float64Array.from( function*(){ for( const a of A.data ) yield* a.data }() )
          ),
          tabulate(shapes[1],'float64', () => Math.random()*2-1)
        ]
      }
    }()
  ).it('svd_lstsq solves least squares for random rank-deficient examples', ([A,y]) => {
//    console.log('SHAPES:', [...A.shape], [...y.shape])
    const [U,sv,V]= svd_decomp(A),
                x = svd_lstsq(U,sv,V, y),
               Ax = matmul2(A,x)
  
    expect( matmul2(A.T, zip_elems([Ax,y], math.sub) ) ).toBeAllCloseTo(0)
  })


  const svd_decomps = {
    svd_decomp,
    svd_jac_1sided,
    svd_jac_2sided,
    svd_jac_classic
  }

  for( const [svd_name,svd_deco] of Object.entries(svd_decomps) )
    forEachItemIn(
      function*(){
        const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

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

      if( M >= N ) {
        const I = eye(N)
        expect( matmul2(U.T,U) ).toBeAllCloseTo(I)
        expect( matmul2(V.T,V) ).toBeAllCloseTo(I)
        expect( matmul2(V,V.T) ).toBeAllCloseTo(I)
      }
      else {
        const I = eye(M)
        expect( matmul2(U.T,U) ).toBeAllCloseTo(I)
        expect( matmul2(U,U.T) ).toBeAllCloseTo(I)
        expect( matmul2(V,V.T) ).toBeAllCloseTo(I)
      }
      expect(D).toBeDiagonal()
      expect(a).toBeAllCloseTo(A)
    })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

      for( let run=1024; run-- > 0; )
      {
        let M = randInt(1,32),
            N = randInt(1,32)
        if( N > M ) [M,N] = [N,M]
        const L = Math.min(M,N),
          ndim = randInt(0,7)
        let r_shape = Array.from({ length: ndim-2 }, () => randInt(1,8) ),
            A_shape = r_shape.concat([M,N])
        r_shape = Int32Array.from(r_shape)
        A_shape = Int32Array.from(A_shape)

        const A = tabulate(A_shape, 'float64', () => Math.random() < 0.1 ? 0 : Math.random()*2-1)
        Object.freeze(A.data.buffer)
        const [Q,_]= qr_decomp(A),
               r = tabulate(    r_shape, 'int32', () => Math.random()*L),
               R = tabulate([...r_shape,L,N], 'float64', (...idx) => {
                 const j = idx.pop(),
                       i = idx.pop(), R = r(...idx)
                 if(R <= i || j < i) return 0
                 if(j == i) return (0.5 + Math.random()) * (randInt(0,2)*2 - 1)
                 return Math.random()*0.5 - 0.25
               })
        yield [r,matmul2(Q,R)]
      }
    }()
  ).it('svd_rank works on random examples', ([R,A]) => {
    const [U,sv,V] = svd_decomp(A),
             r     = svd_rank(sv)

    expect(r.shape).toEqual(A.shape.slice(0,-2))
    expect(R.shape).toEqual(A.shape.slice(0,-2))
    expect(r.dtype).toBe('int32')
    expect(R.dtype).toBe('int32')
    expect(r).toBeAllCloseTo(R, {rtol:0, atol:0})
  })
})
