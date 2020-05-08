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

import {lstsq} from './lstsq'
import {matmul,
        matmul2} from './matmul'
import {solve} from './solve'
import {_rand_rankdef} from '../_test_data_generators'


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
  
    expect( matmul2(A,x) ).toBeAllCloseTo(y, {atol: 1e-7})
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
    expect( matmul2(A.T, zip_elems([Ax,y], (x,y) => x-y) ) ).toBeAllCloseTo(0)
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
  ).it('solves least squares for random rank-deficient examples', ([A,y]) => {
    const x   = lstsq(A,y),
         Ax   = matmul2(A,x)
  
    expect( matmul2(A.T, zip_elems([Ax,y], (x,y) => x-y) ) ).toBeAllCloseTo(0)
  })
})
