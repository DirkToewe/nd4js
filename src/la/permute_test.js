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
import {permute_rows,
        permute_cols,
      unpermute_rows,
      unpermute_cols} from './permute'


describe('permute', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  it('permute_rows works on example pair of shapes [3,2] and [3]', () => {
    const A = [[11,12],
               [21,22],
               [31,32]]

    const permutations = [
      [0,1,2],
      [0,2,1],
      [1,0,2],
      [1,2,0],
      [2,0,1],
      [2,1,0]
    ];
    for( const P of permutations )
    {
      const B = permute_rows(A,P)
      expect(B.dtype).toBe('int32')
      expect(B.shape).toEqual( Int32Array.of(3,2) )
      expect(B).toBeAllCloseTo( A.map( (_,i) => A[P[i]] ) )
    }
  })


  it('permute_rows works on example pair of shapes [2,3,2] and [3]', () => {
    const A = [[[111,112],
                [121,122],
                [131,132]],
               [[211,212],
                [221,222],
                [231,232]]]

    const permutations = [
      [0,1,2],
      [0,2,1],
      [1,0,2],
      [1,2,0],
      [2,0,1],
      [2,1,0]
    ]

    for( const P of permutations )
    {
      const B = permute_rows(A,P)
      expect(B.dtype).toBe('int32')
      expect(B.shape).toEqual( Int32Array.of(2,3,2) )
      expect(B).toBeAllCloseTo(
        A.map( A => A.map( (_,i) => A[P[i]] ) )
      )
    }
  })


  it('permute_rows works on example pair of shapes [3,2] and [2,3]', () => {
    const A = [[11,12],
               [21,22],
               [31,32]]

    const permutations = [
      [0,1,2],
      [0,2,1],
      [1,0,2],
      [1,2,0],
      [2,0,1],
      [2,1,0]
    ];
    for( const P of permutations )
    for( const Q of permutations )
    {
      const B = permute_rows(A,[P,Q])
      expect(B.dtype).toBe('int32')
      expect(B.shape).toEqual( Int32Array.of(2,3,2) )
      expect(B).toBeAllCloseTo([
        A.map( (_,i) => A[P[i]] ),
        A.map( (_,i) => A[Q[i]] )
      ])
    }
  })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

      for( let run=1024; run-- > 0; )
      {
        let  ndim = randInt(0,5),
          shape_A = Int32Array.from({length: ndim}, () => randInt(1,8)),
          shape_P = shape_A.slice(),
              cut = randInt(0,ndim+1)

        for( let i=shape_A.length; i-- > cut; )
          switch( randInt(0,4) )
          {
            case 0 : shape_A[i] = 1; break
            case 1 : shape_P[i] = 1; break
            default:
          }

        if( Math.random() < 0.5 ) shape_A = shape_A.slice(cut)
        else                      shape_P = shape_P.slice(cut)

        const M = randInt(1,24),
              N = randInt(1,24)

        let A = tabulate( shape_A, () => tabulate([M,N], 'float64', () => Math.random()*2-1) ),
            P = tabulate( shape_P, () => {
              const P = Int32Array.from({length: M}, (_,i) => i)
              for( let i=M; i-- > 0; )
              {
                const P_i = P[i],    j = randInt(0,i+1)
                            P[i] = P[j]
                                   P[j] = P_i
              }
              return P
            }),
            B = zip_elems([A,P], (A,P) => tabulate([M,N], (i,j) => A(P[i],j)) )

        A = new NDArray(  Int32Array.of(...A.shape,M,N), Float64Array.from(function*(){ for( const a of A.data ) yield* a.data }())  )
        P = new NDArray(  Int32Array.of(...P.shape,M  ),   Int32Array.from(function*(){ for( const p of P.data ) yield* p      }())  )
        B = new NDArray(  Int32Array.of(...B.shape,M,N), Float64Array.from(function*(){ for( const b of B.data ) yield* b.data }())  )

        yield [A,P,B]
      }
    }()
  ).it('permute_rows works on random examples', ([A,P,B]) => {
    const b = permute_rows(A,P)
    expect(b.shape).toEqual(B.shape)
    expect(b).toBeAllCloseTo(B, {rtol:0, atol:0})
  })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

      for( let run=1024; run-- > 0; )
      {
        let  ndim = randInt(0,5),
          shape_A = Int32Array.from({length: ndim}, () => randInt(1,8)),
          shape_P = shape_A.slice(),
              cut = randInt(0,ndim+1)

        for( let i=shape_A.length; i-- > cut; )
          switch( randInt(0,4) )
          {
            case 0 : shape_A[i] = 1; break
            case 1 : shape_P[i] = 1; break
            default:
          }

        if( Math.random() < 0.5 ) shape_A = shape_A.slice(cut)
        else                      shape_P = shape_P.slice(cut)

        const M = randInt(1,24),
              N = randInt(1,24)

        let A = tabulate( shape_A, () => tabulate([M,N], () => Math.random()*2-1) ),
            P = tabulate( shape_P, () => {
              const P = Int32Array.from({length: N}, (_,i) => i)
              for( let i=N; i-- > 0; )
              {
                const P_i = P[i],    j = randInt(0,i+1)
                            P[i] = P[j]
                                   P[j] = P_i
              }
              return P
            }),
            B = zip_elems([A,P], (A,P) => tabulate([M,N], (i,j) => A(i,P[j])) )

        A = new NDArray(  Int32Array.of(...A.shape,M,N), Float64Array.from(function*(){ for( const a of A.data ) yield* a.data }())  )
        P = new NDArray(  Int32Array.of(...P.shape,  N),   Int32Array.from(function*(){ for( const p of P.data ) yield* p      }())  )
        B = new NDArray(  Int32Array.of(...B.shape,M,N), Float64Array.from(function*(){ for( const b of B.data ) yield* b.data }())  )

        yield [A,P,B]
      }
    }()
  ).it('permute_cols works on random examples', ([A,P,B]) => {
    const b = permute_cols(A,P)
    expect(b.shape).toEqual(B.shape)
    expect(b).toBeAllCloseTo(B, {rtol:0, atol:0})
  })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

      for( let run=1024; run-- > 0; )
      {
        let  ndim = randInt(0,5),
          shape_A = Int32Array.from({length: ndim}, () => randInt(1,8)),
          shape_P = shape_A.slice(),
              cut = randInt(0,ndim+1)

        for( let i=shape_A.length; i-- > cut; )
          switch( randInt(0,4) )
          {
            case 0 : shape_A[i] = 1; break
            case 1 : shape_P[i] = 1; break
            default:
          }

        if( Math.random() < 0.5 ) shape_A = shape_A.slice(cut)
        else                      shape_P = shape_P.slice(cut)

        const M = randInt(1,24),
              N = randInt(1,24)

        let A = tabulate( shape_A, () => tabulate([M,N], () => Math.random()*2-1) ),
            P = tabulate( shape_P, () => {
              const P = Int32Array.from({length: M}, (_,i) => i)
              for( let i=M; i-- > 0; )
              {
                const P_i = P[i],    j = randInt(0,i+1)
                            P[i] = P[j]
                                   P[j] = P_i
              }
              return P
            })

        A = new NDArray(  Int32Array.of(...A.shape,M,N), Float64Array.from(function*(){ for( const a of A.data ) yield* a.data }())  )
        P = new NDArray(  Int32Array.of(...P.shape,M  ),   Int32Array.from(function*(){ for( const p of P.data ) yield* p      }())  )

        yield [A,P]
      }
    }()
  ).it('unpermute_rows works on random examples', ([A,P]) => {
    {
      const a = unpermute_rows(permute_rows(A,P),P)
      expect(a).toBeAllCloseTo(A, {rtol:0, atol:0})
    }{
      const a = permute_rows(unpermute_rows(A,P),P)
      expect(a).toBeAllCloseTo(A, {rtol:0, atol:0})
    }
  })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

      for( let run=1024; run-- > 0; )
      {
        let  ndim = randInt(0,5),
          shape_A = Int32Array.from({length: ndim}, () => randInt(1,8)),
          shape_P = shape_A.slice(),
              cut = randInt(0,ndim+1)

        for( let i=shape_A.length; i-- > cut; )
          switch( randInt(0,4) )
          {
            case 0 : shape_A[i] = 1; break
            case 1 : shape_P[i] = 1; break
            default:
          }

        if( Math.random() < 0.5 ) shape_A = shape_A.slice(cut)
        else                      shape_P = shape_P.slice(cut)

        const M = randInt(1,24),
              N = randInt(1,24)

        let A = tabulate( shape_A, () => tabulate([M,N], () => Math.random()*2-1) ),
            P = tabulate( shape_P, () => {
              const P = Int32Array.from({length: N}, (_,i) => i)
              for( let i=N; i-- > 0; )
              {
                const P_i = P[i],    j = randInt(0,i+1)
                            P[i] = P[j]
                                   P[j] = P_i
              }
              return P
            })

        A = new NDArray(  Int32Array.of(...A.shape,M,N), Float64Array.from(function*(){ for( const a of A.data ) yield* a.data }())  )
        P = new NDArray(  Int32Array.of(...P.shape,  N),   Int32Array.from(function*(){ for( const p of P.data ) yield* p      }())  )

        yield [A,P]
      }
    }()
  ).it('unpermute_cols works on random examples', ([A,P]) => {
    {
      const a = unpermute_cols(permute_cols(A,P),P)
      expect(a).toBeAllCloseTo(A, {rtol:0, atol:0})
    }{
      const a = permute_cols(unpermute_cols(A,P),P)
      expect(a).toBeAllCloseTo(A, {rtol:0, atol:0})
    }
  })
})
