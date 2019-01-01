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
import {zip_elems} from '../zip_elems'
import {tabulate} from '../tabulate'
import {matmul, matmul2} from './matmul'
import {triu} from './tri'
import {lu_decomp, lu_solve} from './lu'
import {array} from '../nd_array'
import math from '../math'


describe('lu', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  function test_lu_decomp(A)
  {
    const [LU,P]= lu_decomp(A),
            U   = triu(LU),
           L    = LU.mapElems('float64', (LU_ij,...indices) => {
             const [i,j] = indices.slice(-2);
             return i===j ? 1 :
                    i < j ? 0 : LU_ij
           })

     expect(LU.dtype).toBe('float64')
     expect( P.dtype).toBe(  'int32')
  
     const absMax = L.reduceElems((x,y) => Math.max(
       math.abs(x),
       math.abs(y)
     ))
     expect(absMax).toBe(1)
  
     A = tabulate(A.shape, 'float64', (...idx) => {
       idx[idx.length-2] = P(...idx.slice(0,-1) )
       return A(...idx)
     })

     expect(LU.shape).toEqual(A.shape)
     expect( P.shape).toEqual(A.shape.slice(0,-1))
     expect( matmul2(L,U) ).toBeAllCloseTo(A)
  }

  it('lu_decomp works on example of shape [3,3]', () => {
    const A = array(
      [[2, -1,  6],
       [1, -2,  3],
       [0,  2,-10]]
    )
    test_lu_decomp(A)
  })

  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from

      for( let run=1024; run-- > 0; )
      {
        let A_shape = Int32Array.from({ length: randInt(2,5) }, () => randInt(1,24))
        A_shape[A_shape.length-2] = A_shape[A_shape.length-1]
  
        yield tabulate(A_shape, 'float64', () => Math.random()*2 - 1)
      }
    }()
  ).it('lu_decomp works on random examples', test_lu_decomp)

  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from

      for( let run=1024; run-- > 0; )
      {
        const [
          LU_shape,
           P_shape,
           y_shape
        ] = function(){
          const shape = Array.from({ length: randInt(2,6) }, () => randInt(1,8) ),
                shapes = [
                  shape.slice( randInt(0, shape.length) ),
                  shape.slice( randInt(0, shape.length) )
                ]
          shapes.splice( randInt(0,3), 0, shape )
  
          for( let d=shape.length; d > 0; d-- )
          for( let i=randInt(0,3); i-- > 0; )
          {
            const shape = shapes[randInt(0,3)],
                  j = shape.length - d
            if(0<=j)  shape[j] = 1
          }

          return shapes
        }()

        const M = randInt(1,24),
              N = randInt(1,24)
        LU_shape.push(M,M)
         P_shape.push(M)
         y_shape.push(M,N)

        let
          y = tabulate(y_shape, 'float64', () => Math.random()*2-1),
          P = tabulate(P_shape.slice(0,-1), 'int32', () => {
            const idx = Int32Array.from({ length: P_shape[P_shape.length-1] }, (_,i) => i );
            // SHUFFLE
            for( let i=idx.length; --i > 0; ) {
              const j = randInt(0,i+1),
                idx_i = idx[i];
                        idx[i] = idx[j];
                                 idx[j] = idx_i
            }
            return idx;
          }),
          LU = tabulate(LU_shape, 'float64', (...indices) => {
            const [i,j] = indices.slice(-2);
            return i==j ? (Math.random()*1    + 0.5) * (randInt(0,2)*2 - 1)
                        :  Math.random()*2e-1 - 1e-1
          })
        P = tabulate(P_shape, 'int32', (...idx) => P(...idx.slice(0,-1))[idx[idx.length-1]])

        yield [LU,P,y]
      }
    }()
  ).it('lu_solve works on random examples', ([LU,P,Y]) => {
    const L = LU.mapElems('float64', (LU_ij,...indices) => {
            const [i,j] = indices.slice(-2);
            return i===j ? 1 :
                   i < j ? 0 : LU_ij
          }),
          U = triu(LU),
          x = lu_solve(LU,P,Y),
          y = matmul(L,U,x)

    // BROADCAST TO COMMON SHAPE
    Y = zip_elems([Y,P.sliceElems('...','new')], 'float64', y => y)
    P = zip_elems([P,Y.sliceElems('...', 0   )],   'int32', P => P)

    // PERMUTE Y
    Y = tabulate( Y.shape, 'float64', (...idx) => {
      idx[idx.length-2] = P(...idx.slice(0,-1) )
      return Y(...idx)
    })

    expect(y).toBeAllCloseTo(Y)
  })
})
