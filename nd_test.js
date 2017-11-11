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
{
  'use strict'

  const nd = require('./nd.js')

   //
  // SETUP HOMEBREW TEST FRAMEWORK
 //
  function test( name, procedure )
  {
    console.log('Testing '+name+'...')

      const assert = (assertion,err_msg) => {
        if( ! assertion )
          throw new Error( null  ==  err_msg
            ? 'Assertion failed.'
            : 'Assertion failed: ' + err_msg + '.'
          )
      }
      assert.err = (procedure,err_msg) => {
        try {
          procedure()
        }
        catch(err) { return }
        throw new Error( null == err_msg
          ? 'Error expected but none thrown.'
          : 'Error expected but none thrown: ' + err_msg + '.'
        )
      }
       //
      // EQUALITY ASSERTIONS
     //
      assert.eq = (expect, actual, err_msg) => {
        if( expect != actual )
          throw new Error( null  ==  err_msg
            ?             '"'+expect+'" expected but "'+actual+'" encountered.'
            : err_msg + '. "'+expect+'" expected but "'+actual+'" encountered.'
          )
      }
       //
      // STRICT EQUALITY ASSERTIONS
     //
      assert.strict_eq = (expect, actual, err_msg) => {
        if( expect !== actual )
          throw new Error( null  ==  err_msg
            ?             '"'+expect+'" expected but "'+actual+'" encountered.'
            : err_msg + '. "'+expect+'" expected but "'+actual+'" encountered.'
          )
      }
       //
      // ARRAY EQUALITY ASSERTIONS
     //
      assert.array_eq = (expect, actual, err_msg) => {
        if( expect.length !== actual.length || expect.some( (x,i) => x != actual[i] ) )
          throw new Error( null  ==  err_msg
            ?             '['+expect+'] expected but ['+actual+'] encountered.'
            : err_msg + '. ['+expect+'] expected but ['+actual+'] encountered.'
          )
      }
       //
      // ARRAY STRICT EQUALITY ASSERTIONS
     //
      assert.array_strict_eq = (expect, actual, err_msg) => {
        if( expect.length !== actual.length || expect.some( (x,i) => x !== actual[i] ) )
          throw new Error( null  ==  err_msg
            ?             '['+expect+'] expected but ['+actual+'] encountered.'
            : err_msg + '. ['+expect+'] expected but ['+actual+'] encountered.'
          )
      }

    console.time('Time')
    procedure(assert)
    console.timeEnd('Time')

    console.log('Finished testing "'+name+'" successfully.\n')
  }

   //
  // FACTORY METHODS
 //
  test('nd.array', assert => {
    function test(...shape)
    {
      function array(d, ...indices)
      {
        if( shape.length === d )
          return indices.reduce( (a,b) => 10*a+b )
        else
          return Array.from({ length: shape[d] }, (_,i) => array(d+1, ...indices, i) )
      }
      array = nd.array( array(0) )
      assert.array_eq(shape, array.shape)
      function test(d, ...indices)
      {
        if( shape.length === d )
          assert.strict_eq(
            indices.reduce( (a,b) => 10*a+b ),
            array(...indices)
          )
        else for( let i=0; i < shape[d]; i++ )
          test(d+1, ...indices, i)
      }
      test(0)
    }
    for( let i=1; i < 8; i++ ) {test(i)
    for( let j=1; j < 8; j++ ) {test(i,j)
    for( let k=1; k < 8; k++ ) {test(i,j,k)
    for( let l=1; l < 8; l++ ) {test(i,j,k,l)}}}}
  })

  test('nd.Array.from#1', assert => {
    const
      a = new nd.Array([2,2], [1,2,3,4]),
      b = new nd.Array([2,2], [5,6,7,8])
      c = nd.Array.from([a,b], (a_ij,b_ij,i,j) => {
        assert.strict_eq(a(i,j), a_ij)
        assert.strict_eq(b(i,j), b_ij)
        return 10*a_ij + b_ij
      })
    assert.array_strict_eq([2,2], c.shape)
    assert.array_strict_eq([
      15, 26,
      37, 48
    ], c.data)
  })

  test('nd.Array.from#2', assert => {
    const
      a = new nd.Array([3,1], [1,2,3]),
      b = new nd.Array([1,4], [4,5,6,7])
      c = nd.Array.from([a,b], (a_i0,b_j,i,j) => {
        assert.strict_eq(a(i,0), a_i0)
        assert.strict_eq(b(0,j),   b_j)
        return 10*a_i0 + b_j
      })
    assert.array_strict_eq([3,4], c.shape)
    assert.array_strict_eq([
      14, 15, 16, 17,
      24, 25, 26, 27,
      34, 35, 36, 37,
    ], c.data)
  })

  test('nd.Array.from#3', assert => {
    function test(...shape)
    {
      function test(d, shapeA, shapeB)
      {
        if( shape.length === d )
        {
          function toA(indices){ return shapeA.map( (shp,i) => indices[i+shape.length-shapeA.length] % shp ) }
          function toB(indices){ return shapeB.map( (shp,i) => indices[i+shape.length-shapeB.length] % shp ) }
          const
            a = new nd.Array(shapeA, Int32Array.from({ length: shapeA.reduce((a,b) => a*b, 1) }, (_,i) => (i+1)*1000 ) ),
            b = new nd.Array(shapeB, Int32Array.from({ length: shapeB.reduce((a,b) => a*b, 1) }, (_,i) => (i+1)      ) )
          const
            c =     nd.Array.from([a,b], (aVal,bVal,...indices) => {
              assert.strict_eq(a(...toA(indices)), aVal)
              assert.strict_eq(b(...toB(indices)), bVal)
              return a+b
            })
          assert.err( () => nd.Array.from([a,b]) )
          assert.err( () => nd.Array.from([a,b], 'int32') )
          assert.err( () => nd.Array.from([a,b], 'not_a_type', (aVal,bVal,...indices) => {}) )
        }
        else {
          test(d+1, [...shapeA,shape[d]], [...shapeB,shape[d]])
          test(d+1, [...shapeA,      1 ], [...shapeB,shape[d]])
          test(d+1, [...shapeA,shape[d]], [...shapeB,      1 ])
          test(d+1, [...shapeA,      1 ], [...shapeB,      1 ])
          if( shapeA.length === 0 ) test(d+1, [], [...shapeB,      1 ])
          if( shapeB.length === 0 ) test(d+1, [...shapeA,      1 ], [])
        }
      }
      test(0,[],[])
    }
    for( let i=2; i <= 4; i++ ) {test(i)
    for( let j=2; j <= 4; j++ ) {test(i,j)
    for( let k=2; k <= 4; k++ ) {test(i,j,k)
    for( let l=2; l <= 4; l++ ) {test(i,j,k,l)}}}}
  })

  test('nd.Array.from#4', assert => {
    function test(shapeA, shapeB)
    {
      let broadcastCompatible = true
      for( let i=shapeA.length, j=shapeB.length; i-- > 0 && j-- > 0; )
        broadcastCompatible &= shapeA[i] === shapeB[j] || shapeA[i] === 1 || shapeB[j] === 1
      if( ! broadcastCompatible ) {
        const
          a = new nd.Array(shapeA, Int32Array.from({ length: shapeA.reduce((a,b) => a*b, 1) }, (_,i) => (i+1)*1000 ) ),
          b = new nd.Array(shapeB, Int32Array.from({ length: shapeB.reduce((a,b) => a*b, 1) }, (_,i) => (i+1)      ) )
        assert.err( () => nd.Array.from([a,b]                                           ) )
        assert.err( () => nd.Array.from([a,b],            (aVal,bVal,...indices) => null) )
        assert.err( () => nd.Array.from([a,b], 'float32', (aVal,bVal,...indices) => null) )
        assert.err( () => nd.Array.from([a,b], 'float32'                                ) )
      }
    }
    for( let i=1; i <= 6; i++ ) {
    for( let j=1; j <= 6; j++ ) {test([i],[j])
    for( let k=1; k <= 6; k++ ) {test([i,j],[k]);     test([i],[j,k]);
    for( let l=1; l <= 6; l++ ) {test([i],[j,k,l]);   test([i,j],[k,l]);   test([i,j,k],[l])
    for( let m=1; m <= 6; m++ ) {test([i],[j,k,l,m]); test([i,j],[k,l,m]); test([i,j,k],[l,m]); test([i,j,k,l],[m])}}}}}
  })

  test('nd.tabulate', assert => {
    function test(...shape)
    {
      arr = nd.tabulate(shape, (...indices) => indices.reduce((a,b) => 10*a+b) )
      assert.array_eq(shape, arr.shape)
      function test(d,...indices)
      {
        if( d == shape.length )
          assert.strict_eq( indices.reduce((a,b) => 10*a+b), arr(...indices) )
        else for( let i=0; i < shape[d]; i++ )
          test(d+1, ...indices, i)
      }
      test(0)
    }
    for( let i=1; i < 10; i++ ) {test(i)
    for( let j=1; j < 10; j++ ) {test(i,j)
    for( let k=1; k < 10; k++ ) {test(i,j,k)
    for( let l=1; l < 10; l++ ) {test(i,j,k,l)}}}}
  })

   //
  // GENERAL
 //
  test('nd.Array.toString#1', assert => {
    const a = nd.array([
      [1,2,3],
      [4,5,6]
    ])
    assert.strict_eq('[[1,2,3],[4,5,6]]', a.toString().replace(/\s/g,'') )
  })

  test('nd.Array.toString#2', assert => {
    const a = nd.array([
      [[1,2],[3,4]],
      [[5,6],[7,8]]
    ])
    assert.strict_eq('[[[1,2],[3,4]],[[5,6],[7,8]]]', a.toString().replace(/\s/g,'') )
  })

  test('nd.Array(...indices)#1', assert => {
    const arr = new nd.Array([2,3], [
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
        assert.strict_eq( a_ij, arr(i,  j  ) )
        assert.strict_eq( a_ij, arr(i-2,j  ) )
        assert.strict_eq( a_ij, arr(i,  j-3) )
        assert.strict_eq( a_ij, arr(i-2,j-3) )
      }
      for( let j=-100; j <= 100; j++ ) {
      for( let i=   2; i <= 100; i++ ) assert.err( () => arr(i,j) )
      for( let i=-100; i <   -2; i++ ) assert.err( () => arr(i,j) )
      }
      for( let i=-100; i <= 100; i++ ) {
      for( let j=   3; j <= 100; j++ ) assert.err( () => arr(i,j) )
      for( let j=-100; j <   -3; j++ ) assert.err( () => arr(i,j) )
      }
    }
    test(); arr.data =   Int32Array.from(arr.data)
    test(); arr.data = Float64Array.from(arr.data)
    test(); arr.data = Float32Array.from(arr.data)
    test()
  })

  test('nd.Array(...indices)#2', assert => {
    const arr = new nd.Array([4], [1,2,3,4])
    function test()
    {
      for( let [i,a_i] of [ [0,1], [1,2], [2,3], [3,4] ])
      {
        assert.strict_eq( a_i, arr(i  ) )
        assert.strict_eq( a_i, arr(i-4) )
      }
      for( let i=   4; i <= 100; i++ ) assert.err( () => arr(i) )
      for( let i=-100; i <   -4; i++ ) assert.err( () => arr(i) )
    }
    test(); arr.data =   Int32Array.from(arr.data)
    test(); arr.data = Float64Array.from(arr.data)
    test(); arr.data = Float32Array.from(arr.data)
    test()
  })

  test('nd.Array(...indices)#3', assert => {
    function test(...shape) {
      assert( shape.every( d => 0 < d && d < 10 ) )
      const
        len = shape.length,
        arr = new nd.Array( shape, Array.from(
          { length: shape.reduce( (len,d) => len*d, 1 ) },
          (_,i) => i
        )),
        visited = Uint8Array.from(arr.data, () => 0 )
      function test(d, ...multi_idx)
      {
        if( shape.length == d )
        {
          let
            idx = 0,
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
            assert.err( () => arr(...multi_idx) )
          else {
            assert( 4*(1<<len) >= ++visited[idx] )
            assert.strict_eq( idx, arr(...multi_idx) )
          }
        }
        else for( let i=-8; i < +8; i++ )
          test(d+1, ...multi_idx, i)
      }
      test(0); arr.data =   Int32Array.from(arr.data)
      test(0); arr.data = Float64Array.from(arr.data)
      test(0); arr.data = Float32Array.from(arr.data)
      test(0)
      assert( visited.every( x => x === 4*(1<<len) ) )
    }
    test()
    for( let i=1; i <= 4; i++ ) {test(i)
    for( let j=1; j <= 4; j++ ) {test(i,j)
    for( let k=1; k <= 4; k++ ) {test(i,j,k)}}}
  })

   //
  // TRANSFORMATIONS
 //
  test('nd.Array.slice#1', assert => {
    const a = new nd.Array([3,4],[
      11, 12, 13, 14,
      21, 22, 23, 24,
      31, 32, 33, 34
    ])
    function test()
    {
      let b
      
      for( row of [0,1,2] )
      {
        b = a.slice(row)
        assert.strict_eq(a.dtype, b.dtype)
        assert.array_eq([4], b.shape)
        assert.array_eq([a(row,0),a(row,1),a(row,2),a(row,3)], b.data)
      }
  
      for( rest of ['...', [], [,,], [,,1], [,4,], [,4,1], [0,,], [0,,1], [0,4,], [0,4,1]])
        for( row of [0,1,2] )
        {
          b = a.slice(row,rest)
          assert.strict_eq(a.dtype, b.dtype)
          assert.array_eq([4], b.shape)
          assert.array_eq([a(row,0),a(row,1),a(row,2),a(row,3)], b.data)
        }
  
      for( rest of ['...',[], [,,], [,,1], [,3,], [,3,1], [0,,], [0,,1], [0,3,], [0,3,1]])
        for( col of [0,1,2,3] )
        {
          b = a.slice(rest,col)
          assert.strict_eq(a.dtype, b.dtype)
          assert.array_eq([3], b.shape)
          assert.array_eq([a(0,col),a(1,col),a(2,col)], b.data)
        }
  
      b = a.slice([1,,],[1,3])
      assert.strict_eq(a.dtype, b.dtype)
      assert.array_eq([2,2], b.shape)
      assert.array_eq([22,23,32,33], b.data)

      b = a.slice([,,2],[1,3])
      assert.strict_eq(a.dtype, b.dtype)
      assert.array_eq([2,2], b.shape)
      assert.array_eq([12,13,32,33], b.data)

      b = a.slice([,,2],[1,,2])
      assert.strict_eq(a.dtype, b.dtype)
      assert.array_eq([2,2], b.shape)
      assert.array_eq([12,14,32,34], b.data)
    }
    test(); a.data =  Int32Array.from(a.data)
    test(); a.data =Float64Array.from(a.data)
    test(); a.data =Float32Array.from(a.data)
    test()
  })

  test('nd.Array.reduce#1', assert => {
    const a = nd.array([
      [[1,2], [3,4]],
      [[5,6], [7,8]]
    ])
    assert.strict_eq(36, a.reduce( (a,b) => a+b ) )
    let
    b = a.reduce([0,1,2], (a,b) => a+b)
    assert.array_strict_eq([], b.shape)
    assert.array_strict_eq([36], b.data)

    b = a.reduce([], (a,b) => a+b )
    assert.array_strict_eq(a.shape, b.shape)
    assert.array_strict_eq(a.data,  b.data )

    b = a.reduce(0, (a,b) => a+b)
    assert.array_strict_eq([2,2], b.shape)
    assert.array_strict_eq([6,8,10,12], b.data)
    b = a.reduce(1, (a,b) => a+b)
    assert.array_strict_eq([2,2], b.shape)
    assert.array_strict_eq([4,6,12,14], b.data)
    b = a.reduce(2, (a,b) => a+b)
    assert.array_strict_eq([2,2], b.shape)
    assert.array_strict_eq([3,7,11,15], b.data)

    b = a.reduce([0,1], (a,b) => a+b)
    assert.array_strict_eq([2], b.shape)
    assert.array_strict_eq([16,20], b.data)
    b = a.reduce([0,2], (a,b) => a+b)
    assert.array_strict_eq([2], b.shape)
    assert.array_strict_eq([14,22], b.data)
    b = a.reduce([1,2], (a,b) => a+b)
    assert.array_strict_eq([2], b.shape)
    assert.array_strict_eq([10,26], b.data)
  })

  test('nd.Array.reduce#2', assert => {
    const a = nd.array([
      [1,2,3],
      [4,5,6]
    ])
    assert.strict_eq(21, a.reduce( (a,b) => a+b ) )
    let
    b = a.reduce([0,1], (a,b) => a+b)
    assert.array_strict_eq([], b.shape)
    assert.array_strict_eq([21], b.data)

    b = a.reduce([], (a,b) => a+b )
    assert.array_strict_eq(a.shape, b.shape)
    assert.array_strict_eq(a.data,  b.data )

    b = a.reduce(0, (a,b) => a+b)
    assert.array_strict_eq([3], b.shape)
    assert.array_strict_eq([5,7,9], b.data)
    b = a.reduce(1, (a,b) => a+b)
    assert.array_strict_eq([2], b.shape)
    assert.array_strict_eq([6,15], b.data)
  })

  test('nd.concat#1', assert => {
    const
      a = nd.array([
        [11,12,13],
        [21,22,23]
      ]),
      b = nd.array([
        [31,32,33]
      ]),
      c = nd.array([
        [41,42,43],
        [51,52,53],
        [61,62,63],
      ])
    let
    d = nd.concat(0,[a,b,c])
    assert.array_strict_eq([6,3], d.shape)
    assert.array_strict_eq([11,12,13, 21,22,23, 31,32,33, 41,42,43, 51,52,53, 61,62,63], d.data)
  })

  test('nd.stack#1', assert => {
    const
      a = nd.array([
        [11,12,13],
        [21,22,23]
      ]),
      b = nd.array([
        [31,32,33],
        [41,42,43]
      ]),
      c = nd.array([
        [51,52,53],
        [61,62,63]
      ])
    let
    d = nd.stack(0,[a,b,c])
    assert.array_strict_eq([3,2,3], d.shape)
    assert.array_strict_eq([11,12,13, 21,22,23, 31,32,33, 41,42,43, 51,52,53, 61,62,63], d.data)

    d = nd.stack(1,[a,b,c])
    assert.array_strict_eq([2,3,3], d.shape)
    assert.array_strict_eq([11,12,13, 31,32,33, 51,52,53, 21,22,23, 41,42,43, 61,62,63], d.data)

    d = nd.stack(2,[a,b,c])
    assert.array_strict_eq([2,3,3], d.shape)
    assert.array_strict_eq([11,31,51, 12,32,52, 13,33,53, 21,41,61, 22,42,62, 23,43,63], d.data)
  })

   //
  // DATA TYPE OPERATIONS
 //
  test('nd.Array.dtype', assert => {
    assert.strict_eq('float64', new nd.Array([3],  Float64Array.of(1,2,3)       ).dtype )
    assert.strict_eq('float32', new nd.Array([2,3],Float32Array.of(1,2,3,4,5,6) ).dtype )
    assert.strict_eq(  'int32', new nd.Array([1,3],  Int32Array.of(1,2,3,4)     ).dtype )
    assert.strict_eq( 'object', new nd.Array([3,3],         [1,2,3,4,5,6,7,8,9] ).dtype )
  })

  test('nd._check_dtype', assert => {
    nd._check_dtype(  'int32')
    nd._check_dtype('float32')
    nd._check_dtype('float64')
    nd._check_dtype( 'object')
    assert.err( () => nd._check_dtype(12) )
    assert.err( () => nd._check_dtype('not_a_type') )
  })

  test('nd.dtypes', assert => {
    assert( 'object' in nd.dtypes )
    assert('float64' in nd.dtypes )
    assert('float32' in nd.dtypes )
    assert(  'int32' in nd.dtypes )
    assert( Object.keys(nd.dtypes).length === 4 )
    assert.strict_eq(       Array, nd.dtypes.object )
    assert.strict_eq(Float64Array, nd.dtypes.float64)
    assert.strict_eq(Float32Array, nd.dtypes.float32)
    assert.strict_eq(  Int32Array, nd.dtypes.  int32)
  })

  test('nd.dtypeof', assert => {
    assert.strict_eq(  'int32', nd.dtypeof( 1337) )
    assert.strict_eq('float64', nd.dtypeof(1.337) )
    assert.strict_eq( 'object', nd.dtypeof('str') )
    assert.strict_eq( 'object', nd.dtypeof({x:2}) )
    assert.strict_eq( 'object', nd.dtypeof([1,2]) )
  })

  test('nd.super_dtype', assert => {
    assert.strict_eq(  'int32', nd.super_dtype('int32', 'int32') )

    assert.strict_eq('float32', nd.super_dtype('int32',  'float32') )
    assert.strict_eq('float32', nd.super_dtype('float32',  'int32') )
    assert.strict_eq('float32', nd.super_dtype('float32','float32') )

    assert.strict_eq('float64', nd.super_dtype('int32',  'float64') )
    assert.strict_eq('float64', nd.super_dtype('float32','float64') )
    assert.strict_eq('float64', nd.super_dtype('float64',  'int32') )
    assert.strict_eq('float64', nd.super_dtype('float64','float32') )
    assert.strict_eq('float64', nd.super_dtype('float64','float64') )

    assert.strict_eq('object', nd.super_dtype(  'int32','object') )
    assert.strict_eq('object', nd.super_dtype('float32','object') )
    assert.strict_eq('object', nd.super_dtype('float64','object') )
    assert.strict_eq('object', nd.super_dtype('object',  'int32') )
    assert.strict_eq('object', nd.super_dtype('object','float32') )
    assert.strict_eq('object', nd.super_dtype('object','float64') )
    assert.strict_eq('object', nd.super_dtype('object', 'object') )
  })

  test('nd.is_subdtype', assert => {
    assert( nd.is_subdtype(  'int32','object') )
    assert( nd.is_subdtype('float32','object') )
    assert( nd.is_subdtype('float64','object') )
    assert( nd.is_subdtype( 'object','object') )

    assert( nd.is_subdtype(  'int32','float64') )
    assert( nd.is_subdtype('float32','float64') )
    assert( nd.is_subdtype('float64','float64') )
    assert(!nd.is_subdtype( 'object','float64') )

    assert( nd.is_subdtype(  'int32','float32') )
    assert( nd.is_subdtype('float32','float32') )
    assert(!nd.is_subdtype('float64','float32') )
    assert(!nd.is_subdtype( 'object','float32') )

    assert( nd.is_subdtype(  'int32','int32') )
    assert(!nd.is_subdtype('float32','int32') )
    assert(!nd.is_subdtype('float64','int32') )
    assert(!nd.is_subdtype( 'object','int32') )
  })
}