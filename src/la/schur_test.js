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
import {schur_decomp,
        schur_eigen,
        schur_eigenvals} from './schur'
import {matmul, matmul2} from './matmul'
import {eye} from './eye'
import math from '../math'
import {zip_elems} from '../zip_elems'
import {rank} from './rank'
import {qr_decomp,
        qr_lstsq} from './qr'
import {NDArray} from '../nd_array'
import {tril_solve} from './tri'


describe('schur', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

      for( let run=1024; run-- > 0; )
      {
        const N = randInt(1,16),
          shape = Array.from({length: randInt(0,3)}, () => randInt(1,8) )
        shape.push(N,N)
  
        const [Q,_]= qr_decomp( tabulate(shape, 'float64', () => Math.random()*2-1) ),
                 R = tabulate(shape, 'float64', (...idx) => {
                   const j=idx.pop(),
                         i=idx.pop()
                   if( i > j ) return 0
                   if( i===j ) return (Math.random()*1 + 0.5) * (Math.random() < 0.5 ? -1 : +1)
                   return Math.random()*0.4 - 0.2
                 })

        const Λ = new NDArray(
          Int32Array.from(shape),
          Float64Array.from(
            function*(){
              for( let k=shape.slice(0,-2).reduce(math.mul,1); k-- > 0; )
              {
                const ev = Float64Array.from(
                  {length: randInt(1,N)},
                  () => (Math.random()*1 + 0.5) * (Math.random() < 0.5 ? -1 : +1)
                )

                for( let i=0; i < N; i++ )
                for( let j=0; j < N; j++ )
                  if( i===j ) yield ev[Math.trunc(Math.random()*ev.length)]
                  else        yield 0
              }
            }()
          )
        )

        // T = R @ Λ @ inv(R)
        const T = tril_solve(R.T, matmul2(R,Λ).T).T

        yield [Q,T]
      }
    }()
  ).it('schur_eigen works on random examples with repeated eigenvalues and full set of eigenvectors', ([Q,T]) => {
    let [N] = Q.shape.slice(-1),
      [Λ,V] = schur_eigen(Q,T)

    expect( V.mapElems('float64',c=>c.im) ).toBeAllCloseTo(0)
    V = V.mapElems('float64',c=>c.re)

    expect( rank(V) ).toBeAllCloseTo(N, {rtol:0,atol:0})
  })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

      for( let run=1024; run-- > 0; )
      {
        const N = randInt(1,32),
           ndim = randInt(2,5),
          shape = Array.from({ length: ndim-2 }, () => randInt(1,8) )
        shape.push(N,N)

        const A = tabulate(shape, 'float64', () => Math.random() < 0.1 ? 0 : Math.random()*2-1)
        Object.freeze(A.data.buffer)
        yield A
      }
    }()
  ).it('schur_decomp works on random examples', A => {
    const [N] = A.shape.slice(-1),       
        [Q,T] = schur_decomp(A),
           a  = matmul(Q, T, Q.T);

    expect(T).toBeUpperHessenberg()

    for( const t of T.reshape(-1,N,N) )
    for( let i=2; i < N; i++ )
      expect( t(i,i-1)*t(i-1,i-2) ).toBe(0)
  
    expect(Q.shape).toEqual(A.shape)
    expect(T.shape).toEqual(A.shape)

    const I = eye(N)
    expect( matmul2(Q,Q.T) ).toBeAllCloseTo(I)
    expect( matmul2(Q.T,Q) ).toBeAllCloseTo(I)
  
    expect(a).toBeAllCloseTo(A)
  })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

      for( let run=1024; run-- > 0; )
      {
        const N = randInt(1,16),
           ndim = randInt(2,5),
          shape = Array.from({ length: ndim-2 }, () => randInt(1,8) )
        shape.push(N,N)

        const A = tabulate(shape, 'float64', () => Math.random() < 0.1 ? 0 : Math.random()*2-1)
        Object.freeze(A.data.buffer)
        yield A
      }
    }()
  ).it('schur_eigen works on random examples', A => {
    const [Q,T] = schur_decomp(A),
          [Λ,V] = schur_eigen(Q,T)

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
           ndim = randInt(2,5),
          shape = Array.from({ length: ndim-2 }, () => randInt(1,8) )
        shape.push(N,N)

        const A = tabulate(shape, 'float64', () => Math.random() < 0.1 ? 0 : Math.random()*2-1)
        Object.freeze(A.data.buffer)
        yield A
      }
    }()
  ).it('schur_eigenvals works on random examples', A => {
    const [Q,T]= schur_decomp(A),
          [Λ,_]= schur_eigen(Q,T),
          EV   = schur_eigenvals(T)

    expect(EV.shape).toEqual(Λ.shape)
    expect(EV).toBeAllCloseTo(Λ)
  })
})
