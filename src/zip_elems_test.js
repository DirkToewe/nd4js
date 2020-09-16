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

import {Complex128Array} from './dt';
import {forEachItemIn} from './jasmine_utils'
import {NDArray} from './nd_array'
import {zip_elems} from './zip_elems'

import {is_array} from './arrays/is_array'
import {Comparator as ArrayComparator} from './arrays/comparator';


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


  const ARR_TYPES = Object.freeze([
    [     'int32',      Int32Array],
    [   'float64',    Float64Array],
    ['complex128', Complex128Array],
    [    'object',           Array]
  ]);

  for( const [A_str, A_type] of ARR_TYPES )
  for( const [B_str, B_type] of ARR_TYPES )
  for( const dtype of [
    [],
    [undefined],
    [     'int32'],
    [   'float64'],
    ['complex128'],
    ['object']
  ])
    forEachItemIn(
      function*(){

        function shapePairs(...shape)
        {
          Object.freeze(shape);
          if( ! shape.every(s => s > 0) )
            throw new Error('Assertion failed.');
          const                     len = shape.length,
            shapeA = new Int32Array(len),
            shapeB = new Int32Array(len);

          function* shapePairs(d)
          {
            if( shape.length === d ) {
              for( let m=-1; m++ < len && (m===0 || shapeA[m-1]===1); )
              for( let n=-1; n++ < len && (n===0 || shapeB[n-1]===1); ) yield [shapeA.slice(m), shapeB.slice(n)];
            }
            else {
              shapeA[d] = shapeB[d] = 1; yield* shapePairs(d+1);
              for( let n=1; n++ < shape[d]; ) {
                shapeA[d] = 1; shapeB[d] = n; yield* shapePairs(d+1);
                shapeA[d] = n; shapeB[d] = 1; yield* shapePairs(d+1);
                shapeA[d] = n; shapeB[d] = n; yield* shapePairs(d+1);
              }
            }
          }

          return shapePairs(0);
        }

        yield *shapePairs(4,4,3,2);
      }()
    ).it(`works on generated example pairs of types (${A_str.padStart(10)},${B_str.padStart(10)}) -> ${dtype.map(x => `${x}`).join('').padStart(10)}`, ([shapeA,shapeB]) => {
      // console.log([...shapeA],[...shapeB])
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

      const a = new NDArray(shapeA, A_type.from({ length: shapeA.reduce((a,b) => a*b, 1) }, (_,i) => (i+1)*1000 ) ),
            b = new NDArray(shapeB, B_type.from({ length: shapeB.reduce((a,b) => a*b, 1) }, (_,i) => (i+1)      ) );
      Object.freeze(a.data.buffer);
      Object.freeze(b.data.buffer);

      const mapper = (aVal,bVal,...indices) => {
        expect(aVal).toEqual( a(...toA(indices)) )
        expect(bVal).toEqual( b(...toB(indices)) )
        return aVal + bVal
      }

      const  c = zip_elems([a,b], ...dtype, mapper)
      expect(c.dtype).toBe( dtype[0] || 'object' );
      expect(c.shape).toEqual(shape);

      for( const [indices,cVal] of c.elems() )
        expect(cVal*1).toBe(
            a(...toA(indices))
          + b(...toB(indices))
        );
      
      expect( () => zip_elems([a,b], 'not_a_type', mapper) ).toThrow();
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
