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

import {NDArray} from './nd_array'
import {zip_elems} from './zip_elems'
import {forEachItemIn} from './jasmine_utils'


describe('zip_elems', () => {
  it('works on an example pair of shapes [2,2] and [2,2]', () => {
    const
      a = new NDArray(Int32Array.of(2,2), [1,2,3,4]),
      b = new NDArray(Int32Array.of(2,2), [5,6,7,8]),
      c = zip_elems([a,b], (a_ij,b_ij,i,j) => {
        expect(a_ij).toBe( a(i,j) )
        expect(b_ij).toBe( b(i,j) )
        return 10*a_ij + b_ij
      })
    expect(c.shape).toEqual( Int32Array.of(2,2) )
    expect(c.data ).toEqual([
      15, 26,
      37, 48
    ])
  })

  it('works on an example pair of shapes [3,1] and [1,4]', () => {
    const
      a = new NDArray(Int32Array.of(3,1), [1,2,3]),
      b = new NDArray(Int32Array.of(1,4), [4,5,6,7]),
      c = zip_elems([a,b], (a_i0,b_0j,i,j) => {
        expect(a_i0).toBe( a(i,0) )
        expect(b_0j).toBe( b(0,j) )
        return 10*a_i0 + b_0j
      })
    expect(c.shape).toEqual( Int32Array.of(3,4) )
    expect(c.data ).toEqual([
      14, 15, 16, 17,
      24, 25, 26, 27,
      34, 35, 36, 37,
    ])
  })

  forEachItemIn(
    function*(){
      const shapePairs = (...shape) => {
        shape = Int32Array.from(shape)

        function* shapePairs(d, shapeA, shapeB)
        {
          if( shape.length === d )
            yield Object.freeze([
              Int32Array.from(shapeA),
              Int32Array.from(shapeB)
            ])
          else {
            yield* shapePairs(d+1, [...shapeA,shape[d]], [...shapeB,shape[d]])
            yield* shapePairs(d+1, [...shapeA,      1 ], [...shapeB,shape[d]])
            yield* shapePairs(d+1, [...shapeA,shape[d]], [...shapeB,      1 ])
            yield* shapePairs(d+1, [...shapeA,      1 ], [...shapeB,      1 ])
            if( shapeA.length === 0 ) yield* shapePairs(d+1, [], [...shapeB,      1 ])
            if( shapeB.length === 0 ) yield* shapePairs(d+1, [...shapeA,      1 ], [])
          }
        }

        return shapePairs(0,[],[])
      }

      for( let i=2; i <= 4; i++ ) { yield* shapePairs(i)
      for( let j=2; j <= 4; j++ ) { yield* shapePairs(i,j)
      for( let k=2; k <= 4; k++ ) { yield* shapePairs(i,j,k)
      for( let l=2; l <= 4; l++ ) { yield* shapePairs(i,j,k,l) }}}}
    }()
  ).it('works on generated example pairs', ([shapeA,shapeB]) => {
    const shape = Int32Array.from(
      { length: Math.max(shapeA.length,shapeB.length)},
      (_,i) => Math.max(
        shapeA[shapeA.length-1-i] || 1,
        shapeB[shapeB.length-1-i] || 1
      )
    ).reverse()
    Object.freeze(shape.buffer)

    function toA(indices){ return shapeA.map( (shp,i) => indices[i+shape.length-shapeA.length] % shp ) }
    function toB(indices){ return shapeB.map( (shp,i) => indices[i+shape.length-shapeB.length] % shp ) }

    const a = new NDArray(shapeA, Int32Array.from({ length: shapeA.reduce((a,b) => a*b, 1) }, (_,i) => (i+1)*1000 ) ),
          b = new NDArray(shapeB, Int32Array.from({ length: shapeB.reduce((a,b) => a*b, 1) }, (_,i) => (i+1)      ) )
    Object.freeze(a.data.buffer)
    Object.freeze(b.data.buffer)

    function testWith( dtype, mapper )
    {
      const c = zip_elems([a,b], dtype, mapper)

      expect(c.shape).toEqual(shape)
      for( const [indices,cVal] of c.elems() )
        expect(cVal).toBe(
            a(...toA(indices))
          + b(...toB(indices))
        )
    }

    const mapper = (aVal,bVal,...indices) => {
      expect(aVal).toBe( a(...toA(indices)) )
      expect(bVal).toBe( b(...toB(indices)) )
      return aVal + bVal
    }

    testWith(mapper)
    for( const dtype of [undefined,'int32','float32','float64','complex128','object'] )
      testWith
    
    expect( () => zip_elems([a,b], 'not_a_type', mapper) ).toThrow()
  })

  forEachItemIn(
    function*(){
      yield [[],[]]
      for( let i=1; i <= 6; i++ ) {
      for( let j=1; j <= 6; j++ ) { yield [[i],[j]]
      for( let k=1; k <= 6; k++ ) { yield [[i,j],[k]];     yield [[i],[j,k]];
      for( let l=1; l <= 6; l++ ) { yield [[i],[j,k,l]];   yield [[i,j],[k,l]];   yield [[i,j,k],[l]]
      for( let m=1; m <= 6; m++ ) { yield [[i],[j,k,l,m]]; yield [[i,j],[k,l,m]]; yield [[i,j,k],[l,m]]; yield [[i,j,k,l],[m]]}}}}}
    }()
  ).it('fails on generated example pairs if not broadcast-compatible', ([shapeA, shapeB]) => {
    let broadcastCompatible = true
    for( let i=shapeA.length, j=shapeB.length; i-- > 0 && j-- > 0; )
      broadcastCompatible &= shapeA[i] === shapeB[j] || shapeA[i] === 1 || shapeB[j] === 1
    if( ! broadcastCompatible )
    {
      shapeA = Int32Array.from(shapeA);
      shapeB = Int32Array.from(shapeB);
      const
        a = new NDArray(shapeA, Int32Array.from({ length: shapeA.reduce((a,b) => a*b, 1) }, (_,i) => (i+1)*1000 ) ),
        b = new NDArray(shapeB, Int32Array.from({ length: shapeB.reduce((a,b) => a*b, 1) }, (_,i) => (i+1)      ) )
      expect( () => nd.Array.from([a,b]                                           ) ).toThrow()
      expect( () => nd.Array.from([a,b],            (aVal,bVal,...indices) => null) ).toThrow()
      expect( () => nd.Array.from([a,b], 'float32', (aVal,bVal,...indices) => null) ).toThrow()
      expect( () => nd.Array.from([a,b], 'float32'                                ) ).toThrow()
    }
  })
})
