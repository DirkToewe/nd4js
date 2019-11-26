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
import {array, NDArray} from '../nd_array'
import {tabulate} from '../tabulate'
import {zip_elems} from '../zip_elems'
import {matmul, matmul2} from './matmul'
import math from '../math'
import {qr_decomp} from './qr'
import {lstsq} from './lstsq'


describe('lstsq', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  it('solves int32 examples', () => {
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

      const x = lstsq(A,y)

      expect(x).toBeAllCloseTo([[1],
                                [1]])
    }
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
  ).it('solves random square examples', ([A,y]) => {
    const x = lstsq(A,y)
  
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
  ).it('solves least squares of random over-determined examples', ([A,y]) => {
    const x = lstsq(A,y),
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
  ).it('computes one solution of random under-determined examples', ([A,y]) => {
    const x = lstsq(A,y)
  
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
  ).it('solves least squares for random rank-deficient examples', ([A,y]) => {
//    console.log('SHAPES:', [...A.shape], [...y.shape])
    const x   = lstsq(A,y),
         Ax   = matmul2(A,x)
  
    expect( matmul2(A.T, zip_elems([Ax,y], math.sub) ) ).toBeAllCloseTo(0)
  })
})
