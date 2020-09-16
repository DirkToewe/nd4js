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

import {matmul, matmul2} from './matmul'
import {cholesky_decomp, cholesky_solve} from './cholesky'


describe('cholesky', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from

      for( let run=2048; run-- > 0; )
      {
        let ndim = randInt(0,4),
          shapes = [ Array.from({length: ndim}, () => randInt(1,4)) ]
        shapes.splice( randInt(0,2), 0, shapes[0].slice( randInt(0,ndim) ) )

        for( let d=ndim; d > 0; d-- )
        for( let i=randInt(0,2); i-- > 0; ) {
          const    shape = shapes[randInt(0,2)],
               j = shape.length - d
          if(0<=j) shape[j] = 1
        }

        const M = randInt(1,32); shapes[0].push(M,M)
        const N = randInt(1,32); shapes[1].push(M,N)

        const y = tabulate(shapes[1],'float64', () => Math.random()*2-1),
              L = tabulate(shapes[0],'float64', (...indices) => {
                const [i,j] = indices.slice(-2);
                return i===j ? Math.random()*1.6 + 0.4
                     : i < j ? 0
                     :         Math.random()*1.2 - 0.6;
              })

        yield [L,y]
      }
    }()
  ).it('cholesky_solve works on random examples', ([L,y]) => {
    const x = cholesky_solve(L,y),
          Y = matmul(L,L.T,x)

    expect(Y).toBeAllCloseTo(y)
  });


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from

      for( let run=2048; run-- > 0; )
      {
        const shape = Int32Array.from({ length: randInt(2,6) }, () => randInt(1,4) );
        shape[shape.length-2] =
        shape[shape.length-1] = randInt(1,32);
        yield tabulate(shape,'float64',(...indices) => {
          const [i,j] = indices.slice(-2);
          if( i < j ) return 0;
          if( i===j ) return Math.random()*1.5 + 0.5;
                      return Math.random()*2.2 - 1.1;
        })
      }
    }()
  ).it('cholesky_decomp works on random examples', L => {
    expect(L).toBeLowerTriangular()

    const LLT = matmul2(L,L.T),
          l = cholesky_decomp(LLT)

    expect(l.shape).toEqual(L.shape)
    expect(l).toBeLowerTriangular()
    expect(l).toBeAllCloseTo(L)  
  });


  forEachItemIn(
    function*(){

      for( let N=1; N < 1024; N = Math.max(N+1,N*1.5) | 0 )
      {
        const A = new NDArray(
            Int32Array.of(N,N),
          Float64Array.from({length: N*N}, () => Math.random()*8 - 4)
        );
        const         LLT = matmul2(A,A.T)
        Object.freeze(LLT.data.buffer)
        Object.freeze(LLT)
        yield         LLT
      }
    }()
  ).it('cholesky_decomp works on random large examples', LLT => {
    const l = cholesky_decomp(LLT),
        llt = matmul2(l,l.T);

    expect(l.shape).toEqual(LLT.shape);
    expect(l).toBeLowerTriangular();
    expect(llt).toBeAllCloseTo(LLT);  
  });
})
