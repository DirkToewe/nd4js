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

import {array, NDArray} from './nd_array'
import {forEachItemIn} from './jasmine_utils'


describe('NDArray.dtype', () => {
  it('returns the correct dtype', () => {
    expect( new NDArray(Int32Array.of(3  ),Float64Array.of(1,2,3)       ).dtype ).toBe('float64')
    expect( new NDArray(Int32Array.of(2,3),Float32Array.of(1,2,3,4,5,6) ).dtype ).toBe('float32')
    expect( new NDArray(Int32Array.of(1,3),  Int32Array.of(1,2,3)       ).dtype ).toBe(  'int32')
    expect( new NDArray(Int32Array.of(3,3),         [1,2,3,4,5,6,7,8,9] ).dtype ).toBe( 'object')
  })
})


describe('nd.array', () => {
  const test = shape => {
    shape = Int32Array.from(shape)
    function arr(d, ...indices)
    {
      if( shape.length === d )
        return indices.reduce( (a,b) => 10*a+b )
      else
        return Array.from({ length: shape[d] }, (_,i) => arr(d+1, ...indices, i) )
    }
    arr = array( arr(0) )
    expect(shape).toEqual(arr.shape)
    function test(d, ...indices)
    {
      if( shape.length === d )
        expect( arr(...indices) ).toBe( indices.reduce((a,b) => 10*a+b) )
      else for( let i=0; i < shape[d]; i++ )
        test(d+1, ...indices, i)
    }
    test(0)
  }

  forEachItemIn(
    function*(){
      for( let i=1; i <= 8; i++ ) { yield [i]
      for( let j=1; j <= 8; j++ ) { yield [i,j]
      for( let k=1; k <= 8; k++ ) { yield [i,j,k]
      for( let l=1; l <= 8; l++ ) { yield [i,j,k,l] }}}}
    }()
  ).it('works on generated examples of nested JS arrays', test)
})


describe('NDArray.toString()', () => {
  it('works on [2,3] example', () => {
    const a = array([
      [1,2,3],
      [4,5,6]
    ])
    expect( a.toString().replace(/\s/g,'') ).toBe('[[1,2,3],[4,5,6]]')
  })
  it('works on [2,2,2] example', () => {
    const a = array([
      [[1,2],
       [3,4]],
      [[5,6],
       [7,8]]
    ])
    expect( a.toString().replace(/\s/g,'') ).toBe('[[[1,2],[3,4]],[[5,6],[7,8]]]')
  })
})


describe('NDArray.toNestedArray', () => {
  forEachItemIn(
    function*(){
      function* shapes(){
        yield Int32Array.of()
        for( let l=1; l < 9; l++ ) { yield Int32Array.of(l)
        for( let m=1; m < 9; m++ ) { yield Int32Array.of(l,m)
        for( let n=1; n < 9; n++ ) { yield Int32Array.of(l,m,n) }}}
      }

      for( const shape of shapes() )
      {
        const length = shape.reduce( (len,d) => len*d, 1 ),
              A = new NDArray( shape, Int32Array.from({length}, (_,i) => i+1) )
        yield Object.freeze(A)
      }
    }()
  ).it('works on random examples.', A => {
    const B = array(A.toNestedArray())

    expect(B.shape).toEqual(A.shape)
    expect(B.dtype).toEqual(A.dtype)
    expect(B.data ).toEqual(A.data )
  })
})


describe('NDArray(...indices)', () => {
  it('works on example of shape [2,3]', () => {
    const arr = new NDArray(Int32Array.of(2,3), [
      1,2,3,
      4,5,6
    ])
    function test()
    {
      for( let [i,j,a_ij] of [
        [0,0, 1], [0,1, 2], [0,2, 3],
        [1,0, 4], [1,1, 5], [1,2, 6]
      ])
      {
        expect( arr(i,  j  ) ).toBe(a_ij)
        expect( arr(i-2,j  ) ).toBe(a_ij)
        expect( arr(i,  j-3) ).toBe(a_ij)
        expect( arr(i-2,j-3) ).toBe(a_ij)
      }
      for( let j=-10; j <= 10; j++ ) {
      for( let i=  2; i <= 10; i++ ) expect( () => arr(i,j) ).toThrow()
      for( let i=-10; i <  -2; i++ ) expect( () => arr(i,j) ).toThrow()
      }
      for( let i=-10; i <= 10; i++ ) {
      for( let j=  3; j <= 10; j++ ) expect( () => arr(i,j) ).toThrow()
      for( let j=-10; j <  -3; j++ ) expect( () => arr(i,j) ).toThrow()
      }
    }
    test(); arr.data =   Int32Array.from(arr.data)
    test(); arr.data = Float64Array.from(arr.data)
    test(); arr.data = Float32Array.from(arr.data)
    test()
  })

  it('works on example of shape [4]', () => {
    const arr = new NDArray(Int32Array.of(4), [1,2,3,4])
    function test()
    {
      for( let [i,a_i] of [ [0,1], [1,2], [2,3], [3,4] ])
      {
        expect( arr(i  ) ).toBe(a_i)
        expect( arr(i-4) ).toBe(a_i)
      }
      for( let i=   4; i <= 100; i++ ) expect( () => arr(i) ).toThrow()
      for( let i=-100; i <   -4; i++ ) expect( () => arr(i) ).toThrow()
    }
    test(); arr.data =   Int32Array.from(arr.data)
    test(); arr.data = Float64Array.from(arr.data)
    test(); arr.data = Float32Array.from(arr.data)
    test()
  })

  forEachItemIn(
    function*(){ 
      for( const T of [Array, Int32Array, Float32Array, Float64Array] )
      {
        yield [T,[]]
        for( let i=1; i <= 4; i++ ) { yield [T,Int32Array.of(i)]
        for( let j=1; j <= 4; j++ ) { yield [T,Int32Array.of(i,j)]
        for( let k=1; k <= 4; k++ ) { yield [T,Int32Array.of(i,j,k)] }}}
      }
    }()
  ).it('works on generated examples', ([ArrayType,shape]) => {
    expect( shape.every( d => 0 < d && d < 10 ) ).toBe(true)

    const len = shape.length,
          arr = new NDArray( Int32Array.from(shape), ArrayType.from(
            { length: shape.reduce( (len,d) => len*d, 1 ) },
            (_,i) => i
          )),
          visited = Uint8Array.from(arr.data, () => 0 )

    function test(d, ...multi_idx)
    {
      if( shape.length == d )
      {
        let idx = 0,
            stride = 1,
            out_of_bounds = false
        for( let i=len; i-- > 0; stride *= shape[i] )
        {
          idx += multi_idx[i] * stride
          if(0 > multi_idx[i])
            idx  +=  shape[i] * stride
          out_of_bounds |= (
               multi_idx[i] < -shape[i] 
            || multi_idx[i] >= shape[i]
          )
        }
        if( out_of_bounds )
          expect( () => arr(...multi_idx) ).toThrow()
        else {
          expect( 4*(1<<len) ).not.toBeLessThan( ++visited[idx] )
          expect( arr(...multi_idx) ).toBe(idx)
        }
      }
      else for( let i=-8; i < +8; i++ )
        test(d+1, ...multi_idx, i)
    }
    test(0)
    expect( visited.every( x => x === 1<<len ) ).toBe(true)
  })
})


describe('NDArray.elems', () => {
  forEachItemIn(
    function*(){
      yield []
      for( let i=1; i <= 6; i++ ) { yield [i]
      for( let j=1; j <= 6; j++ ) { yield [i,j]
      for( let k=1; k <= 6; k++ ) { yield [i,j,k] }}}
    }()
  ).it('works on generated examples', shape => {
    shape = Int32Array.from(shape)
    const len = shape.reduce( (len,d) => len*d, 1 ),
          a = new NDArray(shape, Int32Array.from({length: len}, (_,i) => i)),
          visited = new Int32Array(len)

    for( const [idx,val] of a.elems() )
    {
      const i = a(...idx)
      expect(val).toBe(i)
      expect(++visited[i]).toBe(1)
    }

    expect( visited.every( n => 1===n ) ).toBe(true)
  })
})


describe('NDArray.forElems', () => {
  forEachItemIn(
    function*(){
      yield []
      for( let i=1; i <= 6; i++ ) { yield [i]
      for( let j=1; j <= 6; j++ ) { yield [i,j]
      for( let k=1; k <= 6; k++ ) { yield [i,j,k] }}}
    }()
  ).it('works on generated examples', shape => {
    shape = Int32Array.from(shape)
    const len = shape.reduce( (len,d) => len*d, 1 ),
          a = new NDArray(shape, Int32Array.from({length: len}, (_,i) => i)),
          visited = new Int32Array(len)

    a.forElems( (val,...idx) => {
      const i = a(...idx)
      expect(val).toBe(i)
      expect(++visited[i]).toBe(1)
    })

    expect( visited.every( n => 1===n ) ).toBe(true)
  })  
})


describe('NDArray.sliceElems', () => {
  it('works on example of shape [3,4]', () => {
    const a = new NDArray(Int32Array.of(3,4), [
      11, 12, 13, 14,
      21, 22, 23, 24,
      31, 32, 33, 34
    ])
    function test()
    {
      let b
      
      for( const row of [0,1,2] )
      {
        b = a.sliceElems(row)
        expect(a.dtype).toBe(b.dtype)
        expect(b.shape).toEqual( Int32Array.of(4) )
        expect(b.data ).toEqual(
          a.data.constructor.of( a(row,0), a(row,1), a(row,2), a(row,3) )
        )
      }

      for( const rest of ['...', [], [,,], [,,1], [,4,], [,4,1], [0,,], [0,,1], [0,4,], [0,4,1]])
        for( const row of [0,1,2] )
        {
          b = a.sliceElems(row,rest)
          expect(a.dtype).toBe(b.dtype)
          expect(b.shape).toEqual( Int32Array.of(4) )
          expect(b.data ).toEqual(
            a.data.constructor.of( a(row,0), a(row,1), a(row,2), a(row,3) )
          )
        }

      for( const rest of ['...',[], [,,], [,,1], [,3,], [,3,1], [0,,], [0,,1], [0,3,], [0,3,1]])
        for( const col of [0,1,2,3] )
        {
          b = a.sliceElems(rest,col)
          expect(a.dtype).toBe(b.dtype)
          expect(b.shape).toEqual( Int32Array.of(3) )
          expect(b.data ).toEqual(
            a.data.constructor.of( a(0,col), a(1,col), a(2,col) )
          )
        }

      b = a.sliceElems([1,,],[1,3])
      expect(a.dtype).toBe(b.dtype)
      expect(b.shape).toEqual( Int32Array.of(2,2) )
      expect(b.data ).toEqual(
        a.data.constructor.of(22,23,32,33)
      )

      b = a.sliceElems([,,2],[1,3])
      expect(a.dtype).toBe(b.dtype)
      expect(b.shape).toEqual( Int32Array.of(2,2) )
      expect(b.data ).toEqual(
        a.data.constructor.of(12,13,32,33)
      )

      b = a.sliceElems([,,2],[1,,2])
      expect(a.dtype).toBe(b.dtype)
      expect(b.shape).toEqual( Int32Array.of(2,2) )
      expect(b.data ).toEqual(
        a.data.constructor.of(12,14,32,34)
      )
    }
    test(); a.data =  Int32Array.from(a.data)
    test(); a.data =Float64Array.from(a.data)
    test(); a.data =Float32Array.from(a.data)
    test()
  })
})


describe('NDArray.reduceElems', () => {
  it('works on example of shape [2,3]', () => {
    const a = array([
      [1,2,3],
      [4,5,6]
    ])
    expect( a.reduceElems( (a,b) => a+b ) ).toEqual(21)
    let
    b = a.reduceElems([0,1], (a,b) => a+b)
    expect(b.shape).toEqual( Int32Array.of() )
    expect(b.data ).toEqual([21])

    b = a.reduceElems([], (a,b) => a+b )
    expect(a.shape).toEqual(b.shape)
    expect(b.data ).toEqual([...a.data])

    for( const axis of [0,[0]] )
    {
      b = a.reduceElems(axis, (a,b) => a+b )
      expect(b.shape).toEqual( Int32Array.of(3) )
      expect(b.data ).toEqual([5,7,9])
    }

    for( const axis of [1,[1]] )
    {
      b = a.reduceElems(axis, (a,b) => a+b )
      expect(b.shape).toEqual( Int32Array.of(2) )
      expect(b.data ).toEqual([6,15])
    }
  })

  it('works on example of shape [2,2,2]', () => {
    const a = array([
      [[1,2], [3,4]],
      [[5,6], [7,8]]
    ])
    expect( a.reduceElems( (a,b) => a+b ) ).toEqual(36)
    let
    b = a.reduceElems([0,1,2], (a,b) => a+b)
    expect(b.shape).toEqual( Int32Array.of() )
    expect(b.data ).toEqual([36])


    b = a.reduceElems([], (a,b) => a+b )
    expect(b.shape).toEqual(a.shape)
    expect(b.data ).toEqual([...a.data])


    b = a.reduceElems(0, (a,b) => a+b)
    expect(b.shape).toEqual( Int32Array.of(2,2) )
    expect(b.data ).toEqual([6,8,10,12])

    b = a.reduceElems(1, (a,b) => a+b)
    expect(b.shape).toEqual( Int32Array.of(2,2) )
    expect(b.data ).toEqual([4,6,12,14])

    b = a.reduceElems(2, (a,b) => a+b)
    expect(b.shape).toEqual( Int32Array.of(2,2) )
    expect(b.data ).toEqual([3,7,11,15])


    b = a.reduceElems([0,1], (a,b) => a+b)
    expect(b.shape).toEqual( Int32Array.of(2) )
    expect(b.data ).toEqual([16,20])

    b = a.reduceElems([0,2], (a,b) => a+b)
    expect(b.shape).toEqual( Int32Array.of(2) )
    expect(b.data ).toEqual([14,22])

    b = a.reduceElems([1,2], (a,b) => a+b)
    expect(b.shape).toEqual( Int32Array.of(2) )
    expect(b.data ).toEqual([10,26])
  })
})


describe('NDArray.map', () => {
  forEachItemIn(
    function*(){
      yield []
      for( let i=1; i <= 8; i++ ) { yield Int32Array.of(i)
      for( let j=1; j <= 8; j++ ) { yield Int32Array.of(i,j)
      for( let k=1; k <= 8; k++ ) { yield Int32Array.of(i,j,k) }}}
    }()
  ).it('works on generated examples', shape => {
    const a = new NDArray(
      Int32Array.from(shape),
      Int32Array.from({ length: shape.reduce( (s,t) => s*t, 1 ) }, (_,i) => i+1 )
    )
    const b = a.mapElems( (a_idx,...idx) => {
      expect(a_idx).toBe( a(...idx) )
      return [a_idx,idx]
    })
    expect(b.shape).toEqual(a.shape)
    for( const [a_idx,idx] of b.data )
      expect(a_idx).toBe( a(...idx) )
  })
})


describe('NDArray.transpose', () => {
  it('works on example of shape [3]', () => { 
    const arr = array([1,2,3])
    for( const axes of [ [], [0] ] )
    {
      const arr_T = arr.transpose(...axes)
      expect(arr_T.shape).toEqual(arr.shape)
      expect(arr_T.data ).toEqual(arr.data )
      expect(arr_T.data).not.toBe(arr.data )
    }
  })

  it('works on example of shape [2,3]', () => { 
    const arr = array([
      [1,2,3],
      [4,5,6]
    ])
    for( const axes of [ [], [1], [1,0] ] )
    {
      const arr_T = arr.transpose(...axes)
      expect(arr_T.shape).toEqual( Int32Array.of(3,2) )
      expect(arr_T.data ).toEqual( Int32Array.of(1,4, 2,5, 3,6) )
    }
    for( const axes of [ [0], [0,1] ] )
    {
      const arr_T = arr.transpose(...axes)
      expect(arr_T.shape).toEqual(arr.shape)
      expect(arr_T.data ).toEqual(arr.data )
      expect(arr_T.data).not.toBe(arr.data )
    }
  })

  it('works on example of shape [2,3,4]', () => { 
    const arr = array(
      [[[ 1, 2, 3, 4],
        [ 5, 6, 7, 8],
        [ 9,10,11,12]],
       [[13,14,15,16],
        [17,18,19,20],
        [21,22,23,24]]]
    )
    for( const axes of [ [], [0,2], [0,2,1] ] )
    {
      const arr_T = arr.transpose(...axes)
      expect(arr_T.shape).toEqual( Int32Array.of(2,4,3) )
      expect(arr_T.data ).toEqual( Int32Array.of(
         1, 5, 9,
         2, 6,10,
         3, 7,11,
         4, 8,12,

        13,17,21,
        14,18,22,
        15,19,23,
        16,20,24
      ))
    }
    for( const axes of [ [1], [1,0], [1,0,2] ] )
    {
      const arr_T = arr.transpose(...axes)
      expect(arr_T.shape).toEqual( Int32Array.of(3,2,4) )
      expect(arr_T.data ).toEqual( Int32Array.of(
         1, 2, 3, 4,
        13,14,15,16,

         5, 6, 7, 8,
        17,18,19,20,

         9,10,11,12,
        21,22,23,24
      ))
    }
    for( const axes of [ [2], [2,0], [2,0,1] ] )
    {
      const arr_T = arr.transpose(...axes)
      expect(arr_T.shape).toEqual( Int32Array.of(4,2,3) )
      expect(arr_T.data ).toEqual( Int32Array.of(
         1, 5, 9,
        13,17,21,
        
         2, 6,10,
        14,18,22,

         3, 7,11,
        15,19,23,

         4, 8,12,
        16,20,24
      ))
    }
  })
})
