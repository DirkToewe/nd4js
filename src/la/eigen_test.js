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
import {array} from '../nd_array'
import {tabulate} from '../tabulate'
import {matmul2} from './matmul'
import math from '../math'
import {zip_elems} from '../zip_elems'
import {eigen,
        eigenvals,
        eigen_balance_pre,
        eigen_balance_post} from './eigen'


describe('eigen', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  for( const P of [2,Infinity])
    it(`eigen_balance_pre (p=${isFinite(P) ? P : '∞'}) works on example of shape [2,2]`, () => {
      const A = array([[1,1],
                       [3,2]]);
      for( const P of [Infinity])
      {
        const [_,B]= eigen_balance_pre(A,Infinity),
            absMax = (x,y) => Math.max( Math.abs(x), Math.abs(y) ),
            before = A.reduceElems(absMax),
            after  = B.reduceElems(absMax)
        expect(after).toBeLessThan(before)
      }
    })


  for( const P of [2,Infinity])
    it(`eigen_balance_pre (p=${isFinite(P) ? P : '∞'}) works on example of shape [3,3]`, () => {
      const A = array([[ 1, 1,-3],
                       [ 4, 2, 2],
                       [ 1, 0, 2]]);
      for( const P of [Infinity])
      {
        const [_,B]= eigen_balance_pre(A,Infinity),
            absMax = (x,y) => Math.max( Math.abs(x), Math.abs(y) ),
            before = A.reduceElems(absMax),
            after  = B.reduceElems(absMax)
        expect(after).toBeLessThan(before)
      }
    })


  for( const P of [2,Infinity] )
    forEachItemIn(
      function*(){
        const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from
  
        for( let run=1024; run-- > 0; )
        {
          const N = randInt(1,32),
            shape = Array.from({length: randInt(0,3)}, () => randInt(1,8) )
          shape.push(N,N)
    
          let A = tabulate(shape, 'float64', (...idx) => {
            const [i,j] = idx.slice(-2)
            if( Math.random() < 0.1 )
              return 0
            return Math.random()*2 - 1
          })
          if( Math.random < 0.5 )
            A = A.T
          Object.freeze(A.data.buffer)
          yield A
        }
      }()
    ).it(`eigen_balance_pre (p=${isFinite(P) ? P : '∞'}) works on random examples`, A => {

      if( P !== 2 && isFinite(P) )
        throw new Error('Assertion failed.')
      const norm = P===2
        ? math.hypot
        : (x,y) => math.max( math.abs(x), math.abs(y) )

      const [D,S] = eigen_balance_pre(A,P)
      expect(A.shape).toEqual(S.shape)

      const before= A.reduceElems(norm),
            after = S.reduceElems(norm)
      expect(after).not.toBeGreaterThan(before)
  
      const a = zip_elems(
        [D.reshape(...D.shape,1), S,
         D.reshape(...D.shape.slice(0,-1),1,-1)],
        'float64',
        (D_i,S_ij,D_j) => D_i * S_ij / D_j
      );
  
      expect(a).toBeAllCloseTo(A,{rtol:0, atol:0})
    })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

      for( let run=1024; run-- > 0; )
      {
        const N = randInt(1,16),
          shape = Array.from({length: randInt(0,3)}, () => randInt(1,8) )
        shape.push(N,N)
  
        let A = tabulate(shape, 'float64', (...idx) => {
          const [i,j] = idx.slice(-2)
          if( Math.random() < 0.1 )
            return 0
          return Math.random()*2 - 1
        })
        if( Math.random < 0.5 )
          A = A.T
        Object.freeze(A.data.buffer)
        yield A
      }
    }()
  ).it('eigen_balance_post works on random examples', A => {
    const [D,B]= eigen_balance_pre(A,2),
          [Λ,U]= eigen(B),
             V = eigen_balance_post(D,U)

    expect(V.shape).toEqual( A.shape )
    expect(Λ.shape).toEqual( A.shape.slice(0,-1) )

    // ASSERT THAT THE EIGENVECTORS ARE NORMALIZED
    expect(
      V.   mapElems(            'float64', math.abs  )
       .reduceElems(/*axis=*/-2,'float64', math.hypot)
    ).toBeAllCloseTo(1)

    const AV = matmul2(A,V);

    const λ = zip_elems([AV,V], 'complex128', (x,y) => math.mul(x, math.conj(y)) ).reduceElems(-2, 'complex128', math.add)
    expect(Λ).toBeAllCloseTo(λ)

    const λv = zip_elems([V,Λ.reshape(...Λ.shape.slice(0,-1),1,-1)], 'complex128', math.mul)
    expect(AV).toBeAllCloseTo(λv)
  })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

      for( let run=1024; run-- > 0; )
      {
        const N = randInt(1,16),
          shape = Array.from({length: randInt(0,3)}, () => randInt(1,8) )
        shape.push(N,N)
  
        let A = tabulate(shape, 'float64', (...idx) => {
          const [i,j] = idx.slice(-2)
          if( Math.random() < 0.1 )
            return 0
          return Math.random()*2 - 1
        })
        if( Math.random < 0.5 )
          A = A.T
        Object.freeze(A.data.buffer)
        yield A
      }
    }()
  ).it('eigen works on random examples', A => {
    const [Λ,V] = eigen(A)

    expect(V.shape).toEqual( A.shape )
    expect(Λ.shape).toEqual( A.shape.slice(0,-1) )

    // ASSERT THAT THE EIGENVECTORS ARE NORMALIZED
    expect(
      V.   mapElems(   'float64', math.abs  )
       .reduceElems(-2,'float64', math.hypot)
    ).toBeAllCloseTo(1)

    const AV = matmul2(A,V);

    const λ = zip_elems([AV,V], 'complex128', (x,y) => math.mul(x, math.conj(y)) ).reduceElems(-2, 'complex128', math.add)
    expect(Λ).toBeAllCloseTo(λ)

    const λv = zip_elems([V,Λ.reshape(...Λ.shape.slice(0,-1),1,-1)], 'complex128', math.mul)
    expect(AV).toBeAllCloseTo(λv)
  })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

      for( let run=1024; run-- > 0; )
      {
        const N = randInt(1,16),
          shape = Array.from({length: randInt(0,3)}, () => randInt(1,8) )
        shape.push(N,N)
  
        let A = tabulate(shape, 'float64', (...idx) => {
          const [i,j] = idx.slice(-2)
          if( Math.random() < 0.1 )
            return 0
          return Math.random()*2 - 1
        })
        if( Math.random < 0.5 )
          A = A.T
        Object.freeze(A.data.buffer)
        yield A
      }
    }()
  ).it('eigenvals works on random examples', A => {
    const [Λ] = eigen    (A),
           L  = eigenvals(A)

    expect(L).toBeAllCloseTo(Λ)
  })
})
