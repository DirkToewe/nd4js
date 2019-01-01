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
import {tabulate} from '../tabulate'
import {zip_elems} from '../zip_elems'
import {matmul2} from './matmul'
import {qr_decomp} from './qr'
import {rrqr_decomp,
        rrqr_decomp_full,
        rrqr_rank,
        rrqr_solve,
        rrqr_lstsq} from './rrqr'
import {eye} from './eye'
import math from '../math'


describe('rrqr', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

      for( let run=1024; run-- > 0; )
      {
        const ndim = randInt(2,5),
             shape = Int32Array.from({ length: ndim }, () => randInt(1,24) )
        const A = tabulate(shape, 'float64', () => Math.random() < 0.1 ? 0 : Math.random()*2-1)
        Object.freeze(A.data.buffer)
        yield A
      }
    }()
  ).it('rrqr_decomp_full works on random examples', A => {
    const [M,N] = A.shape.slice(-2),
        [Q,R,P] = rrqr_decomp_full(A),
             a  = matmul2(Q,R)

    Object.freeze(Q.data.buffer); expect(Q.dtype).toBe('float64')
    Object.freeze(R.data.buffer); expect(R.dtype).toBe('float64')
    Object.freeze(P.data.buffer); expect(P.dtype).toBe(  'int32')
    Object.freeze(a.data.buffer)

    expect( Q.shape ).toEqual( Int32Array.of(...A.shape.slice(0,-2), M, M) )
    expect( R.shape ).toEqual( Int32Array.of(...A.shape.slice(0,-2), M, N) )
    expect( P.shape ).toEqual( Int32Array.of(...A.shape.slice(0,-2),    N) )

    const I = eye(M)
    Object.freeze(I.data.buffer)
    expect( matmul2(Q.T,Q) ).toBeAllCloseTo(I)
    expect( matmul2(Q,Q.T) ).toBeAllCloseTo(I)

    A = tabulate(A.shape, (...idx) => {
      idx[A.ndim-1] = P(...idx.slice(0,-2), idx[A.ndim-1])
      return A(...idx)
    })

    expect(R).toBeUpperTriangular()
    expect(a).toBeAllCloseTo(A)

    for( let off=0; off < R.data.length; off += M*N )
      for( let i=Math.min(M,N); --i > 0; )
      {
        const R_II = math.abs(R.data[off + N* i   + i   ]),
              R_ii = math.abs(R.data[off + N*(i-1)+(i-1)])
        expect(R_ii).not.toBeLessThan(R_II)
      }
  })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

      for( let run=1024; run-- > 0; )
      {
        const ndim = randInt(2,5),
             shape = Int32Array.from({ length: ndim }, () => randInt(1,24) )
        const A = tabulate(shape, 'float64', () => Math.random() < 0.1 ? 0 : Math.random()*2-1)
        Object.freeze(A.data.buffer)
        yield A
      }
    }()
  ).it('rrqr_decomp works on random examples', A => {
    const [M,N] = A.shape.slice(-2),
        [Q,R,P] = rrqr_decomp(A),
             L  = Math.min(M,N),
             a  = matmul2(Q,R)
    Object.freeze(Q.data.buffer); expect(Q.dtype).toBe('float64')
    Object.freeze(R.data.buffer); expect(R.dtype).toBe('float64')
    Object.freeze(P.data.buffer); expect(P.dtype).toBe(  'int32')
    Object.freeze(a.data.buffer)

    expect( Q.shape ).toEqual( Int32Array.of(...A.shape.slice(0,-2), M, L) )
    expect( R.shape ).toEqual( Int32Array.of(...A.shape.slice(0,-2), L, N) )
    expect( P.shape ).toEqual( Int32Array.of(...A.shape.slice(0,-2),    N) )

    A = tabulate(A.shape, (...idx) => {
      idx[A.ndim-1] = P(...idx.slice(0,-2), idx[A.ndim-1])
      return A(...idx)
    })

    expect(R).toBeUpperTriangular()
    expect(a).toBeAllCloseTo(A)

    const I = eye(L)
    Object.freeze(I.data.buffer)
                 expect( matmul2(Q.T,Q) ).toBeAllCloseTo(I)
    if( M <= N ) expect( matmul2(Q,Q.T) ).toBeAllCloseTo(I)

    for( let off=0; off < R.data.length; off += L*N )
      for( let i=L; --i > 0; )
      {
        const R_II = math.abs(R.data[off + N* i   + i   ]),
              R_ii = math.abs(R.data[off + N*(i-1)+(i-1)])
        expect(R_ii).not.toBeLessThan(R_II)
      }
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

        const M = randInt(1,32); shapes[0].push(M,M)
        const N = randInt(1,32); shapes[1].push(M,N)

        yield shapes.map(
          s => tabulate(s,'float64', () => Math.random()*2-1)
        )
      }
    }()
  ).it('rrqr_solve solves random square examples', ([QR,y]) => {
    const [Q,R,P] = rrqr_decomp(QR),
             x    = rrqr_solve(Q,R,P, y)
  
    expect( matmul2(QR,x) ).toBeAllCloseTo(y)
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

        const M = randInt(1,32); shapes[0].push(M,M)
        const N = randInt(1,32); shapes[1].push(M,N)

        yield shapes.map(
          s => tabulate(s,'float64', () => Math.random()*2-1)
        )
      }
    }()
  ).it('rrqr_lstsq solves random square examples', ([QR,y]) => {
    const [Q,R,P] = rrqr_decomp(QR),
             x    = rrqr_lstsq(Q,R,P, y)
  
    expect( matmul2(QR,x) ).toBeAllCloseTo(y)
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
        for( let i=randInt(0,2); i-- > 0; ) { // <- 50/50 chance to collapse dimension to test broadcasting
          const    shape = shapes[randInt(0,2)],
               j = shape.length - d
          if(0<=j) shape[j] = 1
        }

        let M = randInt(1,32),
            N = randInt(1,32)
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
  ).it('rrqr_lstsq+rrqr_decomp solves least squares of random over-determined examples', ([A,y]) => {
    const [Q,R,P]= rrqr_decomp(A),
             x   = rrqr_lstsq(Q,R,P, y),
            Ax   = matmul2(A,x)

    // every least square solution satisfies the normal equaltion Aᵀ(Ax - y) = 0
    expect( matmul2(A.T, zip_elems([Ax,y], math.sub) ) ).toBeAllCloseTo(0)
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
        for( let i=randInt(0,2); i-- > 0; ) { // <- 50/50 chance to collapse dimension to test broadcasting
          const    shape = shapes[randInt(0,2)],
               j = shape.length - d
          if(0<=j) shape[j] = 1
        }

        let M = randInt(1,32),
            N = randInt(1,32)
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
  ).it('rrqr_lstsq+rrqr_decomp_full solves least squares of random over-determined examples', ([A,y]) => {
    const [Q,R,P]= rrqr_decomp_full(A),
             x   = rrqr_lstsq(Q,R,P, y),
            Ax   = matmul2(A,x)

    // every least square solution satisfies the normal equaltion Aᵀ(Ax - y) = 0
    expect( matmul2(A.T, zip_elems([Ax,y], math.sub) ) ).toBeAllCloseTo(0)
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
            N = randInt(1,32)
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
  ).it('rrqr_lstsq computes one solution of random under-determined examples', ([A,y]) => {
    const [Q,R,P]= rrqr_decomp_full(A),
             x   = rrqr_lstsq(Q,R,P, y)
  
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
  ).it('rrqr_lstsq solves least squares for random rank-deficient examples', ([A,y]) => {
//    console.log('SHAPES:', [...A.shape], [...y.shape])
    const [Q,R,P]= rrqr_decomp_full(A),
             x   = rrqr_lstsq(Q,R,P, y),
            Ax   = matmul2(A,x)
  
    expect( matmul2(A.T, zip_elems([Ax,y], math.sub) ) ).toBeAllCloseTo(0)
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
  ).it('rrqr_rank works on random examples', ([R,A]) => {
    const [Q,T,P] = rrqr_decomp(A),
             r    = rrqr_rank(T)

    expect(r.shape).toEqual(A.shape.slice(0,-2))
    expect(R.shape).toEqual(A.shape.slice(0,-2))
    expect(r.dtype).toBe('int32')
    expect(R.dtype).toBe('int32')
    expect(r).toBeAllCloseTo(R, {rtol:0, atol:0})
  })
})
