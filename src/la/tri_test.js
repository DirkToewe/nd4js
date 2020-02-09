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
import {matmul2} from './matmul'
import {tril,triu, tril_solve, triu_solve} from './tri'


describe('tri', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  forEachItemIn(
    function*(){
      function* shapes() {
        for( let i=1; i <= 6; i++ )
        for( let j=1; j <= 6; j++ ) { yield [i,j]
        for( let k=1; k <= 6; k++ ) { yield [i,j,k]
        for( let l=1; l <= 6; l++ ) { yield [i,j,k,l] }}}
      }
      for( const shape of shapes() )
      {
        const A = tabulate( shape, 'int32', (...idx) => idx.reduce((flat,s) => 10*flat + s+1, 0) ),
           [M,N]= shape.slice(-2)
        Object.freeze(A.data.buffer)
        for( let k=-M; ++k < N; )
          yield [A,k]
      }
    }()
  ).it('tril works on generated examples', ([A,k]) => {
    const L = tril(A,k)

    for( const [idx,L_idx] of L.elems() )
    {
      const [i,j] = idx.slice(-2)
      expect(L_idx).toBe( i < j-k ? 0 : A(...idx) )
    }
  })

  forEachItemIn(
    function*(){
      function* shapes() {
        for( let i=1; i <= 6; i++ )
        for( let j=1; j <= 6; j++ ) { yield [i,j]
        for( let k=1; k <= 6; k++ ) { yield [i,j,k]
        for( let l=1; l <= 6; l++ ) { yield [i,j,k,l] }}}
      }
      for( const shape of shapes() )
      {
        const A = tabulate( shape, 'int32', (...idx) => idx.reduce((flat,s) => 10*flat + s+1, 0) ),
           [M,N]= shape.slice(-2)
        Object.freeze(A.data.buffer)
        for( let k=-M; ++k < N; )
          yield [A,k]
      }
    }()
  ).it('triu works on generated examples', ([A,k]) => {
    const U = triu(A,k)

    for( const [idx,U_idx] of U.elems() )
    {
      const [i,j] = idx.slice(-2)
      expect(U_idx).toBe( i > j-k ? 0 : A(...idx) )
    }
  })

  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from

      for( let run=1024; run-- > 0; )
      {
        let ndim = randInt(2,6),
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

        const y = tabulate(shapes[1],'float64', () => Math.random()*2-1),
              L = tabulate(shapes[0],'float64', (...indices) => {
                const [i,j] = indices.slice(-2);
                return i===j ? Math.random()*1 + 0.5
                     : i < j ? 0
                     :         Math.random()*2e-1 - 1e-1;
              })

        yield [L,y]
      }
    }()
  ).it('tril_solve works on random examples', ([L,Y]) => {
    expect(L).toBeLowerTriangular()

    const x = tril_solve(L,Y),
          y = matmul2(L,x)

    expect(y).toBeAllCloseTo(Y)
  })

  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from

      for( let run=1024; run-- > 0; )
      {
        let ndim = randInt(2,6),
          shapes = [ Array.from({length: ndim}, () => randInt(1,8)) ]
        shapes.splice( randInt(0,2), 0, shapes[0].slice( randInt(0,ndim) ) )

        for( let d=ndim; d > 0; d-- )
        for( let i=randInt(0,2); i-- > 0; ) {
          const    shape = shapes[randInt(0,2)],
               j = shape.length - d
          if(0<=j) shape[j] = 1
        }

        const M = randInt(1,24); shapes[0].push(M,M);
        const N = randInt(1,24); shapes[1].push(M,N);

        const y = tabulate(shapes[1],'float64', () => Math.random()*2-1),
              U = tabulate(shapes[0],'float64', (...indices) => {
                const [i,j] = indices.slice(-2);
                return i===j ? Math.random()*1 + 0.5
                     : i > j ? 0
                     :         Math.random()*2e-1 - 1e-1;
              })

        yield [U,y]
      }
    }()
  ).it('triu_solve works on random examples', ([U,Y]) => {
    expect(U).toBeUpperTriangular()

    const x = triu_solve(U,Y),
          y = matmul2(U,x)

    expect(y).toBeAllCloseTo(Y)
  })
})
