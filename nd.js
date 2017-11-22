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
let nd = {}
try {
  module.exports = nd
}
catch(err){
  console.log(err)
}
{
  'use strict'

/******************
 * IMPLEMENTATION *
 ******************/
  nd.ellipsis= '...'
  nd.newaxis = 'new'

  nd.Array = class NDArray extends Function
  {
     //
    // FACTORY METHODS
   //
    constructor(shape, data)
    {
      if( shape.some(s => s < 1) )
        throw new Error('Every entry of shape must be >= 1.')
      function self(...indices) { return self.data[self._flat_idx(indices)] }
      Object.setPrototypeOf(self, NDArray.prototype)
      self.shape = shape
      self.data  = data
      return self
    }

    static from( ndarray, dtype, mapper )
    {
      // PREPROCESS ARGUMENTS
      if( null == mapper )
        if( dtype instanceof Function ) { mapper = dtype; dtype = undefined }

      if( ! (ndarray instanceof Array) )
        ndarray = [ndarray]
      else if( null == mapper && ndarray.length > 1 )
        throw new Error('A mapping function is required, when more than 1 ndarray is provided.')

      if( null == dtype ) dtype = null == mapper ? ndarray[0].dtype : 'object'
      if( null == mapper) mapper = x => x

      nd._check_dtype(dtype)

      const
        ndim = ndarray.reduce( (ndim,arr) => Math.max(ndim,arr.ndim), 0),
        shape = Int32Array.from({ length: ndim }, () => 1 )

      // FIND COMMON (BROADCASTED) SHAPE
      for( let arr of ndarray )
        for( let i=ndim, j=arr.ndim; i-- > 0 && j-- > 0; )
          if( 1 === shape[i] )
            shape[i] = arr.shape[j]
          else if( shape[i] != arr.shape[j] && arr.shape[j] != 1 )
            throw new Error('Shapes are not broadcast-compatible.')

      // GENERATE RESULT DATA
      const
        multi_idx = new Int32Array(ndim), // <- index in result
        data = new nd.dtypes[dtype]( shape.reduce((a,b) => a*b, 1) ),

        values  = new      Array(ndarray.length), // <- cache of ndarray[indices]
        indices = new Int32Array(ndarray.length), // <- indices in ndarray(s)
        strides = new Int32Array(ndarray.length)

      let flat_idx = 0

      function write(d) {
        if( d === ndim ) {
          strides.fill(1)
          ndarray.forEach( (arr,i) => values[i] = arr.data[indices[i]++] )
          data[flat_idx++] = mapper(...values, ...multi_idx)
        }
        else for( multi_idx[d] = 0;; )
        {
          write(d+1)
          if( ++multi_idx[d] >= shape[d] )
            return ndarray.forEach( (arr,i) => strides[i] *= (arr.shape[d-ndim+arr.ndim] || 1) )
          ndarray.forEach( (arr,i) => {
            if( ! (arr.shape[d-ndim+arr.ndim] > 1) ) // <- handles undefined (index out of bounds)
              indices[i] -= strides[i]
          })
        }
      }
      write(0)

      return new NDArray(shape,data)
    }

     //
    // GENERAL
   //
    get ndim() { return this.shape.length }

    get dtype() {
      if( this.data instanceof Float64Array ) return 'float64'
      if( this.data instanceof Float32Array ) return 'float32'
      if( this.data instanceof   Int32Array ) return   'int32'
      if( this.data instanceof        Array ) return  'object'
      throw new Error('Data type could not determined.')
    }

    set(indices, value) {
      this.data[this._flat_idx(indices)] = value
    }

    _flat_idx(indices)
    {
      const shape = this.shape

      if( 'number' === typeof indices )
        indices = [indices]

      if(indices.length != shape.length)
        throw new Error(   shape.length  + ' indices expected, but ' + indices.length + ' given.')
      
      let
        flat_idx = 0,
        stride = 1
      for( let i=shape.length; i-- > 0; stride *= shape[i] )
      {
        let idx = indices[i]
        if( idx < 0 )  idx += shape[i]
        if( idx < 0 || idx >= shape[i] ) throw new Error('Multi-index out of bounds.')
        flat_idx  +=   idx * stride
      }
      return flat_idx;
    }

    toString( max_len )
    {
      if( null==max_len ) max_len = 10
      const
        last   = max_len >>> 1,
        first  = max_len - last,
        shape  = this.shape,
        data   = this.data,
        strides= new Int32Array(shape.length)

      strides[strides.length-1] = 1
      for( let i=strides.length; --i > 0; )
        strides[i-1] = strides[i] * shape[i]

      function str( indent, d, idx )
      {
        if( shape.length === d )
        {
          let result = ''+data[idx]
          if( shape.length > 1 )
            while( result.length < 10 )
              result = ' '+result
          return result
        }
        else
        {
          indent += ' '
          let
            prefix = '[',
            infix  = ', ',
            suffix = ']'
          if( d < shape.length-1 )
            infix = ',\n'+indent
          function* elems() {
            if( shape[d] > max_len )
            {
              for( let j=0;             j < first;    j++ ) yield str(indent, d+1, idx+j*strides[d])
              yield ' ...' + (shape[d]-max_len) + ' more...'
              for( let j=shape[d]-last; j < shape[d]; j++ ) yield str(indent, d+1, idx+j*strides[d])
            }
            else for( let j=0;          j < shape[d]; j++ ) yield str(indent, d+1, idx+j*strides[d])
          }
          elems = [...elems()]
          return prefix + elems.join(infix) + suffix
        }
      }
      return str('',0,0)
    }

     //
    // ITERATION
   //
    *[Symbol.iterator]() {
      const
        shape = this.shape.slice(1),
        stride = shape.reduce((a,b) => a*b, 1)
      for( let i=0; i < this.data.length; i += stride )
        yield new nd.Array( shape, this.data.slice(i,i+stride) )
    }

    forEach( consumer ) {
      const
        len = this.shape[0],
        shape = this.shape.slice(1),
        stride = shape.reduce((a,b) => a*b, 1)
      for( let i=0, idx=0; i < len; i++, idx += stride )
        consumer( new nd.Array( shape, this.data.slice(idx,idx+stride) ), i )
    }
    
    *entries() {
      const
        shape= this.shape,
        data = this.data,
        multi_idx = new Int32Array(shape.length) // <- index in result
      let flat_idx = 0

      function* entries(d) {
        if( d === shape.length )
          yield [ multi_idx.slice(), data[flat_idx++] ]
        else for(
          multi_idx[d] = 0;
          multi_idx[d] < shape[d];
          multi_idx[d]++
        )
          yield* entries(d+1)
      }
      yield* entries(0)
    }

    forEntries( consumer ) {
      const
        shape= this.shape,
        data = this.data,
        multi_idx = new Int32Array(shape.length) // <- index in result
      let flat_idx = 0
  
      function forEntries(d) {
        if( d === shape.length )
          consumer( data[flat_idx++], ...multi_idx )
        else for(
          multi_idx[d] = 0;
          multi_idx[d] < shape[d];
          multi_idx[d]++
        )
          forEntries(d+1)
      }
      forEntries(0)
    }

     //
    // TRANSFORMATION
   //
    valueOf() {
      if( this.shape.length === 0 )
        return this.data[0]
      return this
    }

    map( dtype, mapper ) { return NDArray.from(this, dtype, mapper) }

    reshape( ...shape )
    {
      const len = this.data.length;
      let
        rest = 1,
        infer = -1;

      for( let i=shape.length; --i >= 0; )
        if( shape[i] < 0 ) {
          if( shape[i] !== -1 ) throw new Error('Invalid dimension: ' + shape[i])
          if( infer    !== -1 ) throw new Error('At most on dimension may be -1.')
          infer = i;
        }
        else
          rest *= shape[i];

      if( infer !== -1 )
      {
        if( 0 !== len % rest ) throw new Error('Shape cannot be inferred.')
        shape[infer] = len / rest
      }

      return new NDArray(shape, this.data)
    }

    reduce( axes, dtype, reducer )
    {
      if( null == reducer )
      {
        if( null    ==  dtype      ) return this.data.reduce(axes)
        if('string' === typeof axes) return this.data.reduce(dtype)
        reducer = dtype; dtype = undefined
      }
      if( null == dtype ) dtype = 'object'

      if( ! nd.is_subdtype(this.dtype, dtype) )
        throw new Error('New dtype must be a super-dtype.')

      const oldNDim = this.ndim

      if( 'number' === typeof axes ) axes = [axes]
      axes = new Set( function*() {
        for( let ax of axes ) {
          if( 0 > ax )  ax += oldNDim
          if( 0 > ax || ax >= oldNDim ) throw new Error('Reduction axis '+ax+' out of bounds.')
          yield ax
        }
      }() )

      const
        oldShape= this.shape,
        newShape= oldShape.filter( (size,i) => ! axes.has(i) ),
        newData = new nd.dtypes[dtype]( newShape.reduce((a,b) => a*b, 1) ),
        oldData = this.data
      let
        newIdx = 0,
        oldIdx = 0

      function fill(d, reduce)
      {
        if( oldShape.length === d )
        {
          if(reduce) newData[newIdx]=reducer(oldData[oldIdx++], newData[newIdx])
          else       newData[newIdx]=        oldData[oldIdx++]
          ++newIdx
        }
        else if( axes.has(d) ) {
          const idx = newIdx
          for( let j=oldShape[d]; j-- > 0; reduce=true )
          {
            newIdx = idx
            fill(d+1, reduce)
          }
        }
        else for( let j=oldShape[d]; j-- > 0; )
          fill(d+1, reduce)
      }
      fill(0,false)

      if( newIdx !=   newData.length ) throw new Error([newIdx,   newData.length])
      if( oldIdx != this.data.length ) throw new Error([oldIdx, this.data.length])

      return new NDArray(newShape,newData)
    }

    slice(...slices)
    {
      const
        oldNdim = this.ndim,
        oldShape= this.shape
      let
        nEllipses = 0,
        nNew = 0,
        nReduced = 0, // <- dimensions of this matrix that to not appear in the sliced one
        nSliced = 0
      // some basic acCOUNTing
      for( const slc of slices )
      {
             if( slc === nd.ellipsis) ++nEllipses;
        else if( slc === nd.newaxis ) ++nNew;
        else if( (typeof slc) === 'number' )
          if( ! (slc % 1 == 0) )
            throw new Error('Only integral indices allowed.')
          else
            ++nReduced;
        else {
          const len = slc.length;
          if( (typeof len) !== 'number' || ! (len % 1 == 0) || len > 3 || len==1 )
            throw new Error('Illegal argument(s).')
          ++nSliced;
        }
      }
      if( nEllipses >  1 ) throw new Error('Highlander! There can be only one ... (ellipsis).')
      if( nEllipses == 0 )
        slices.push(nd.ellipsis)
      if( nReduced + nSliced > oldNdim )
        throw new Error('Cannot sliced more dimensions than available.')

      // fill in omitted slice arguments and replace ellipsis
      let off = 0;
      const
        strides = [],
        newShape = Int32Array.from(
          function* () {
            let
              d = oldNdim,
              stride = 1;
            for( let i=slices.length; --i >= 0; )
            {
              let slc = slices[i]
              switch(slc)
              {
                case nd.newaxis:
                  strides.push(0)
                  yield 1
                  break
                case nd.ellipsis:
                  for( let i = nReduced + nSliced; i < oldNdim; i++ )
                  {
                    strides.push(stride)
                    stride *= oldShape[--d]
                    yield oldShape[d]
                  }
                  break
                default:
                  --d;
                  if( (typeof slc) === 'number' )
                  {
                    if( 0 > slc )  slc += oldShape[d]
                    if( 0 > slc || slc >= oldShape[d] ) throw new Error('Index out of bounds.')
                    off += stride * slc
                  }
                  else
                  {
                    let [start=undefined,end=undefined,step=1] = slc;

                    if( end   < 0 )  end   += oldShape[d]
                    if( start < 0 )  start += oldShape[d]
                    if( 0 > start || start >= oldShape[d] ) throw new Error('Slice start out of bounds.')

                    if( step == 0 ) throw new Error('Stride/step cannot be zero.')
                    if( step <  0 ) {
                      if(start== null) start= oldShape[d]-1
                      if(end  == null) end  = -1
                      if(  -1 > end || end >= oldShape[d] ) throw new Error('Slice end out of bounds.')
                    } else {
                      if(start== null) start= 0
                      if(end  == null) end  = oldShape[d]
                      if(   0 > end || end >  oldShape[d] ) throw new Error('Slice end out of bounds.')
                    }
                    strides.push(stride * step)
                    off += stride * start
                    yield 1 + Math.trunc( (-Math.sign(step) + end - start) / step )
                  }
                  stride *= oldShape[d]
              }
            }
          }()
        );

      const
        oldData = this.data,
        newData = new nd.dtypes[this.dtype]( newShape.reduce((a,b) => a*b, 1) )
      let newIdx = 0
      function fill(d, oldIdx)
      {
        if( 0 === d )
          newData[newIdx++] = oldData[oldIdx]
        else for( let i=newShape[--d]; i-- > 0; oldIdx += strides[d] )
          fill(d, oldIdx)
      }
      fill(newShape.length, off)

      newShape.reverse()

      return new NDArray(newShape,newData)
    }

//    convolve(options, conv_op) // <- IDEA: img.convolve({shape=[4,4,3], paddings=['same','same','none']}, (values, ...indices) => values.sum() )
  }

  nd.array = (dtype, content) =>
  {
    if( null == content ){ content = dtype; dtype = undefined }

    if( ('object' == typeof content) && 'shape' in content && 'data' in content )
    {
      let data = content.data
      if( null != dtype && ! (data instanceof nd.dtypes[dtype]) )
        data = nd.dtypes[dtype].from(data)
      return new nd.Array(content.shape, data)
    }

    function* shape(content)
    {
      if( 'object' === typeof content && 'length' in content )
      {
        const len = content.length;
        if( (typeof len) !== 'number' || ! (len % 1 === 0) || len === 0 )
          throw 'Illegal argument(s).';
        yield len
        yield* shape(content[0])
      }
    }
    shape = Int32Array.from( shape(content) )

    if( null == dtype )
      dtype = function dtype(d,content)
      {
        if( d === shape.length )
          return nd.dtypeof(content)
        else {
          let
            j = shape[d]-1,
            dt = dtype(d+1,content[j])
          for( ; j >= 0; j-- )
            dt = nd.super_dtype(dt, dtype(d+1,content[j]) )
          return dt;
        }
      }(0,content)

    const data = new nd.dtypes[dtype]( shape.reduce((a,b) => a*b, 1) )
    let idx = 0

    function fill(d,content)
    {
      if( d === shape.length )
        data[idx++] = content
      else {
        if( content.length !== shape[d] )
          throw new Error('Nested array not axis-aligned.')
        for( let j=0; j < shape[d]; j++ )
          fill(d+1, content[j])
      }
    }
    fill(0,content)

    return new nd.Array(shape,data)
  }

  nd.tabulate = ( shape, dtype, idx2val ) =>
  {
    if( null == idx2val){ idx2val = dtype; dtype = undefined }
    if( null == idx2val) throw new Error('Missing argument: idx2val: Function.')
    if( null == dtype ) dtype = 'object'

    const
      multi_idx = new Int32Array(shape.length), // <- index in result
      data = new nd.dtypes[dtype]( shape.reduce((a,b) => a*b, 1) )

    let flat_idx = 0

    function write(d) {
      if( d === shape.length )
        data[flat_idx++] = idx2val(...multi_idx)
      else for(
        multi_idx[d] = 0;
        multi_idx[d] < shape[d];
        multi_idx[d]++
      )
        write(d+1)
    }
    write(0)

    return new nd.Array(shape,data)
  }

  nd.concat = (axis, dtype, ndarrays) =>
  {
    if( null == ndarrays )
    {
      if( null == dtype ){ ndarrays = axis; axis = undefined }
      else {
        ndarrays = dtype; dtype = undefined
        if( 'string' === typeof axis ) { dtype = axis; axis = undefined }
      }
    }
    if( ! ('length' in ndarrays) ) ndarrays = [...ndarrays]
    if( null == axis ) axis = 0
    if( null == dtype) dtype = ndarrays.map( a => a.dtype ).reduce(nd.super_dtype)

    if( 0 > axis )  axis += ndarrays[0].shape.length
    if( 0 > axis || axis >= ndarrays[0].shape.length ) throw new Error('Axis out of bounds.')

    const newShape = Int32Array.from(ndarrays[0].shape)
    for( let i=ndarrays.length; --i > 0; )
    {
      const shape = ndarrays[i].shape

      if(   newShape.length != shape.length )
        throw new Error('All shapes must have the same length.')
      if( ! newShape.every( (len,d) => newShape[d] === shape[d] || axis === d ) )
        throw new Error('Shape along all axes but the concatentation axis must match.')

      newShape[axis] += shape[axis]
    }

    const
      rest = newShape.slice(axis+1).reduce((a,b) => a*b, 1),
      indices = new Int32Array(ndarrays.length),
      newData = new nd.dtypes[dtype]( newShape.reduce((a,b) => a*b, 1) )
    let newIdx = 0

    function fill(d)
    {
      if( d === axis )
        ndarrays.forEach( (arr,i) => {
          for( let j=arr.shape[d]*rest; j-- > 0; )
            newData[newIdx++] = arr.data[indices[i]++]
        })
      else for( let i=newShape[d]; i-- > 0; )
        fill(d+1)
    }
    fill(0)

    return new nd.Array(newShape, newData)
  }

  nd.stack = (axis, dtype, ndarrays) =>
  {
    if( null == ndarrays )
    {
      if( null == dtype ){ ndarrays = axis; axis = undefined }
      else {
        ndarrays = dtype; dtype = undefined
        if( 'string' === typeof axis ) { dtype = axis; axis = undefined }
      }
    }
    if( ! ('length' in ndarrays) ) ndarrays = [...ndarrays]
    if( null == axis ) axis = 0

    if( 0 > axis )  axis += ndarrays[0].shape.length+1
    if( 0 > axis || axis >  ndarrays[0].shape.length ) throw new Error('Axis out of bounds.')

    ndarrays = ndarrays.map( arr => arr.reshape(
      ...arr.shape.slice(0,axis),
      1,
      ...arr.shape.slice(axis)
    ))
    return nd.concat(axis,dtype,ndarrays)
  }

/*
//  nd.str = ()

//  nd.indices = ()

//  nd.filter = ()

   //
  // MODIFICATION
 //
  nd.splice = ()

   //
  // FACTORY METHODS
 //
  nd.array = ()

  nd.scalar = ()

  nd.full = ()

  nd.tabulate = ()

   //
  // TRANSFORMATION
 //
  nd.map = ()

  nd.stack = ()

  nd.reshape = ()
  
//  nd.transpose = ()

  nd.conv = ()

  nd.einsum = ()
*/
   //
  // DATA TYPE OPERATIONS
 //
  nd.dtypes = {
      'int32':   Int32Array,
    'float32': Float32Array,
    'float64': Float64Array,
     'object':        Array
  }

  nd.dtypeof = value => {
    if( 'number' === typeof value )
    {
      if( value % 1 == 0 )
      {
        if(    value <= ~(1 << 31)
            && value >=  (1 << 31) )
          return 'int32'
        return 'object'
      }
      return 'float64'
    }
    return 'object'
  }

  nd.super_dtype = (dtype1,dtype2) => {
    nd._check_dtype(dtype1)
    nd._check_dtype(dtype2)
    if( dtype1 === 'object'  || dtype2 === 'object'  ) return 'object'
    if( dtype1 === 'float64' || dtype2 === 'float64' ) return 'float64'
    if( dtype1 === 'float32' || dtype2 === 'float32' ) return 'float32'
    return 'int32'
  }

  nd.is_subdtype = (sub_dtype, sup_dtype) => {
    nd._check_dtype(sub_dtype)
    nd._check_dtype(sup_dtype)
    const rank = {
        'int32': 0,
      'float32': 1,
      'float64': 2,
       'object': 3
    }
    return rank[sub_dtype] <= rank[sup_dtype]
  }

  nd._check_dtype = dtype => {
    if( ! nd.dtypes.hasOwnProperty(dtype) )
      throw "Invalid dtype '" + dtype + "'. Must be one of {'" + Object.getOwnPropertyNames(nd.dtypes).join("', '") + "'}."
  }

   //
  // LINEAR ALGEBRA
 //
/*
  const la = {}
  nd.la = la
  {
//    la.cholesky_decomp()

//    la.cholesky_solve()

//    la.lu_decomp()

//    la.lu_solve()

//    la.solve()
  }
*/

/*****************
 * DOCUMENTATION *
 *****************/
  nd.createDocHTML = () => {
    const
      jsdom= require('jsdom'),
      dom  = new jsdom.JSDOM('<!DOCTYPE html>')
      doc  = dom.window.document,
      meta = doc.createElement('meta'),
      title= doc.createElement('title'),
      ul   = doc.createElement('ul'),
      h1   = doc.createElement('h1')
    h1   .innerHTML = 'nd.js Documentation'
    title.innerHTML = 'nd.js Documentation'
    meta.setAttribute('charset', 'utf-8')
    doc.lang = 'en'
    doc.head.appendChild(meta)
    doc.head.appendChild(title)
    doc.body.style = 'font-family:monospace'
    doc.body.appendChild(h1)
    doc.body.appendChild(ul)

    const ND = {}
    Object.assign(ND,nd)
    delete ND.Array

    function addSection( key, val )
    {
      const
        h2 = doc.createElement('h2'),
        pre= doc.createElement('pre'),
        li = doc.createElement('li'),
        a  = doc.createElement('a')

      h2.id = key.replace(' ','_').replace('[','').replace(']','')
      a.innerHTML = key
      a.setAttribute('href', '#'+h2.id)
      ul.appendChild(li)
      li.appendChild(a)

      h2 .innerHTML ='nd.'+key
      pre.innerHTML = nd.help_string(val)
      doc.body.appendChild(h2)
      doc.body.appendChild(pre)
    }

    function addClassDoc( name, clazz )
    {
      addSection( name, clazz )

      for( const key of Object.getOwnPropertyNames(clazz) )
        if( ! key.startsWith('__') )
          try {
            addSection(name+'.'+key+' [STATIC]', clazz[key] )
          }
          catch(err) {}

      for( const key of Object.getOwnPropertyNames(clazz.prototype) )
        if( ! key.startsWith('__') )
          try {
            addSection(name+'.'+key, clazz.prototype[key] )
          }
          catch(err) {}
    }

    addClassDoc('Array', nd.Array)

    for( const [key,val] of Object.entries(ND) )
      addSection(key,val)

    return dom
  }

  nd.help = obj => console.log( nd.help_string(obj) )

  nd.help_string = obj => {
    switch(obj)
    {
      case nd.dtypes  : return `\
A mapping of all supported data types and their respective
(primitive) array types.
`
      case nd.ellipsis: return `\
The "symbol" used in nd.Array.slice() to indicate that the
missing slice indices are to be filled in in its stead.
`
      case nd.newaxis : return `\
The "symbol" used in nd.Array.slice() to create an new axis
of size 1.
`
      case null:
      case undefined:
        return '<Doc N/A>'
      default:
        if( obj.hasOwnProperty('__doc__') )
          return obj.__doc__
        else
          return '<Doc N/A>'
    }
  }



  //
 // FACTORY METHODS
//
nd.Array.constructor.__doc__ = `\
Creates a new nd.Array with given shape and data array. This
constructor may not create any protection copies of the given
data but use it directly instead. This method is intended
for internal use mostly. More user-friendly an convenient
factory methods are available like nd.array, nd.tabulate or
nd.Array.from.

Parameters
----------
shape: int[]
  The shape of new nd.Array. May be used directly (no protection
  copy may be perfomed).
data: dtype[]
  The flat data array to be used in this nd.Array. May be used
  directly (no protection copy may be performed).

Returns
-------
ndarray: nd.Array
  A new nd.Array.

Examples
--------
>>> a = new nd.Array([2,3,2], [111,112,121,122,131,132,211,212,221,222,231,232])
>>> console.log( a.toString() )
  [[[111,112],
    [121,122],
    [131,132]],
   [[211,212],
    [221,222],
    [231,232]]]
`



nd.Array.from.__doc__ = `\
Creates a new nd.Array from one or more nd.Arrays, using a function to
map the values of the nd.Array(s) to the values of the new nd.Array.

This can be used to perform unary (sin, cos, exp, ...) and binary
(+, -, *, ...) operations.

This method supports NumPy-style broadcasting of the nd.Array arguments.

Parameters
----------
ndarray: nd.Array or nd.Array[]
  The nd.Array or list of nd.Arrays from which the new nd.Array is built.
dtype: String
  [OPTIONAL] The data type of the new nd.Array.
mapper: (...values, ...indices) => dtype
  A function used to determine the entry values of the new nd.Array from
  the entry values of the old nd.Array(s). May be omitted if only a single
  nd.Array is given by ndarray, in which case the entry values are copied.

Returns
-------
ndarray: nd.Array
  The newly created nd.Array, where:
  ndarray(i0,i1,...,i[n]) = mapper(ndarray[0](i0,i1,...,i[n]), ndarray[1](i0,i1,...,i[n]), ..., ndarray[m](i0,i1,...,i[n]), i0,i1,...,i[n]).

Examples
--------
>>> let a = nd.array([1,2,3,4])
>>> let b = nd.Array.from(a, x => x*x)
>>> console.log( b.toString() )
  [1, 4, 9, 16]


>>> let c = nd.array([[1],[2]])
>>> let d = nd.array([1,2,3,4])
>>> let e = nd.Array.from([c,d], (x,y) => 10*x+y )
>>> console.log( e.toString() )
  [[11,12,13,14],
   [21,22,23,24]]
`



nd.array.__doc__ = `\
Tries to heuristically create an nd.Array from the input. This
method is usually used to create nd.Arrays nested arrays, which
allows for a creation of nd.Array which is well readable by
humans.

  * If the input is an object containing a 'shape' and 'data'
    property, a new nd.Array is created using said shape and data
    directly (without protection copy), unless a conversion is
    necessary due to the dtype.
  * If the input is a (possibly nested) JavaScript array-like, it
    is converted to an nd.Array, copying the array data in the
    process. The nesting defines the shape of the resulting
    nd.Array. 
  * Otherwise the data is interpreted as a scalar value and
    is put into an nd.Array of shape [].

Parameters
----------
dtype: String
  [Optional] The dtype of the returned nd.Array.
content: { length: int[], data: dtype[] } or dtype or dtype[] or dtype[][] or ...
  The content of the returned nd.Array.

Returns
-------
ndarray: nd.Array
  An nd.Array representation of the content.

Examples
--------
>>> let content = { shape: [2,3], data: [1,2,3,4,5,6] }
>>> let a = nd.array(content)
>>> console.log( a.toString() )
  [[1,2,3],
   [4,5,6]]
>>> content.data[3] = 7
>>> console.log( a.toString() )
  [[1,2,3],
   [7,5,6]]


>>> let content = { shape: [2,3], data: [1,2,3,4,5,6] }
>>> let a = nd.array('int32',content)
>>> console.log( a.toString() )
  [[1,2,3],
   [4,5,6]]
>>> content.data[3] = 7
>>> console.log( a.toString() )
  [[1,2,3],
   [4,5,6]]


>>> let a = nd.array([
...   [1,2,3,4],
...   [5,6,7,8]
... ])
>>> console.log( a.toString() )
  [[1,2,3,4],
   [5,6,7,8]]


>>> console.log( nd.array(12) )
  { [NDArray: self] shape: [], data: Int32Array[12] }
`



nd.tabulate.__doc__ = `\
Creates a new nd.Array of the given shape, calling a function for
each entry index and using the returned value as entry value of
the new nd.Array.

Parameters
----------
shape: int[]
  The shape of the new nd.Array.
dtype: String
  [OPTIONAL] The data type of the newly created nd.Array.
idx2val: (...int) => dtype
  The function returning the entry value for a given entry index as input.

Returns
-------
ndarray: nd.Array
  The newly tabulated nd.Array, where:
  ndarray(i0,i1,...,i[n]) = idx2val(i0,i1,...i[n])

Examples
--------
>>> let a = nd.tabulate([3,2], (i,j) => 10*(i+1) + (j+1) )
>>> console.log( a.toString() )
  [[11,12],
   [21,22],
   [31,32]]
`



  //
 // GENERAL
//
nd.Array.__doc__ = `\
An n-dimensional, axis-aligned grid of entries, very similar
to NumPy ndarrays in Python.

Attributes
----------
shape: int[]
  The shape of the nd.Array.
data : dtype[]
  The flat data array backing this nd.Array. May be either
  a standard JavaScript Array, or one of JavaScripts primitive
  array types. See nd.dtypes for a list of supported data types.
ndim: int
  The number of dimensions in this nd.Array, i.e. the number of
  indices used to address the entries. Equivalent to shape.length.
  
dtype: String
  The data type of this nd.Array. Specifies which type of values
  may be stored in the nd.Array. See nd.dtypes for the available
  data types and their respective (primitive) array types.
`


nd.Array.prototype.call.__doc__ = `\
Returns the value of the entry specified by the given index.

Parameters
----------
this: nd.Array
  The array whose entry is being read.
indices: ...int
  The multi-index of the entry whose value is to be returned.

Returns
-------
value: dtype

Examples
--------
>>> let a = nd.array([[1,2,3],[4,5,6]])
>>> console.log( a(0,0), a(0,1), a(0,2), a(1,0), a(1,1), a(1,2) )
  1 2 3 4 5 6
`
nd.Array.prototype.apply = nd.Array.prototype.call



nd.Array.prototype.set.__doc__ = `\
Sets the value of the specified entry to the given value.

Parameters
----------
indices: int[]
  The multi index of the entry that is to be changes.
value: dtype
  The new value for the entry at indices.

Examples
--------
>>> let a = nd.array([ [0,0,0], [0,0,0], [0,0,0] )
>>> a.set([0,0], 1)
>>> a.set([0,1], 2)
>>> a.set([1,1], 3)
>>> a.set([2,2], 4)
>>> console.log( a.toString() )
  [[1,2,0],
   [0,3,0],
   [0,0,4]]
`



nd.Array.prototype.toString.__doc__ = `\
Returns a readable string representation of this nd.Array.

Parameters
----------
max_len: int
  [OPTIONAL] The maximum number of elements along each to be represented
  in the String. If there is more elements along an axis than that, an
  ellipsis (...) is inserted instead of the central elements.

Returns
-------
repr: String
  A readable string representation of this nd.Arrray.

Examples
--------
>>> let a = nd.tabulate([6,6], (i,j) => 10*(i+1) + (j+1) )
>>> console.log( a.toString() )
  [[11,12,13,14,15,16],
   [21,22,23,24,25,26],
   [31,32,33,34,35,36],
   [41,42,43,44,45,46],
   [51,52,53,54,55,56],
   [61,62,63,64,65,66]]

>>> console.log( a.toString(4) )
  [[11, 12, ...2 more..., 15, 16],
   [21, 22, ...2 more..., 25, 26],
    ...2 more...,
   [51, 52, ...2 more..., 55, 56],
   [61, 62, ...2 more..., 65, 66]]
`



  //
 // ITERATIONS
//
nd.Array.prototype[Symbol.iterator].__doc__ =`\
Returns an iterator that iterates through the slices of this nd.Array along
the first axis. Equivalent to:

function*() {
  for( let i=0; i < this.shape[0]; i++ )
    yield this.slice(i)
}

Returns
-------
iter: Iterator
  An iterator that yields the slices of this nd.Array along its first axis.

Examples
--------
>>> a = nd.array([[1,2],[3,4],[5,6]])
>>> for( let row of a )
...   console.log( row.toString() )
  [1,2]
  [3,4]
  [5,6]
`



nd.Array.prototype.forEach.__doc__ = `\
Calls the given callback for each slice of this nd.Array along the
first axis. Equivalent to:

for( let i=0; i < this.shape[0]; i++ )
  consumer( this.slice(i), i )

Parameters
----------
consumer: (slice: nd.Array, index: int) => ()
  The callback to be called for each slice of this nd.Array along the first axis.

Examples
--------
>>> let a = nd.array([[1,2,3],[4,5,6]])
>>> a.forEach( (row,i) => console.log(\`\${i} -> \${row}\`) )
  0 -> [1,2,3]
  1 -> [4,5,6]
`



nd.Array.prototype.entries.__doc__ = `\
Returns an iterator of all multiindex-value pairs of entries in
this nd.Array.

Returns
-------
iter: *[[...int], dtype]
  An iterator of multiindex-value pairs, one for each entry in this
  nd.Array.

Examples
--------
>>> a = nd.array([[1,2,3],[4,5,6]])
>>> for( let [[i,j],a_ij] of a.entries() )
...   console.log(\`a[\${i},\${j}] = \${val}\`)
  a[0,0] = 1
  a[0,1] = 2
  a[0,2] = 3
  a[1,0] = 4
  a[1,1] = 5
  a[1,2] = 6
`



nd.Array.prototype.forEntries.__doc__ = `\
Calls the given callback for each entry in this nd.Array.

Parameters
----------
consumer: (val: dtype, indices: ...int) => ()

Examples
--------
>>> let a = nd.array([[1,2,3],[4,5,6]])
>>> a.forEntries( (a_ij,i,j) => console.log(\`a[\${i},\${j}] = \${a_ij}\`) )
  a[0,0] = 1
  a[0,1] = 2
  a[0,2] = 3
  a[1,0] = 4
  a[1,1] = 5
  a[1,2] = 6
`



  //
 // TRANSFORMATION
//
nd.Array.prototype.valueOf.__doc__ = `
If this nd.Array is scalar (shape == []), the only entry's value is returned,
otherwise this nd.Array itself is returned.

This method rarely has to be called explicitly but is used by JavaScript
whenever it deems it appropriate.

Returns
-------
value: dtype or nd.Array
  Returns the only entry if this nd.Array is scalar. Returns this otherwise.

Example
-------
>>> a = nd.array(12)
>>> console.log( a+1 )
  13
`



nd.Array.prototype.map.__doc__ = `\
Creates a new nd.Array by applying the given mapping function to each entry
of this array and writing the results back into a new array.

Equivalent to nd.Array.from(this, dtype, mapper).

Parameters
----------
dtype: String
  [OPTIONAL] The data type of the newly created nd.Array.
mapper: (value, ...indices) => dtype
  The function used to map the old nd.Array's entries to the values of the
  new nd.Array.

Returns
-------
mapped: nd.Array
  The newly created nd.Array, where:
  mapped(i0,i1,...,i[n]) = mapper(this(i0,i1,...,i[n]), i0,i1,...,i[n]).

Examples
--------
>>> let a = nd.array([[1,2],[3,4]])
>>> let b = a.map( x => x*x )
>>> console.log( b.toString() )
  [[1, 4],
   [9,16]]
`



nd.Array.prototype.reshape.__doc__ = `\
Returns view of this nd.Array with a different shape. This is similar
to NumPy's reshape with a 'C'-order.

Parameters
----------
shape: ...int
  The shape of the view. May contain a single -1 entry, in which case
  the respective axis' size is inferred.

Returns
-------
reshaped: nd.Array
  A reshaped view of this nd.Array.

Examples
--------
>>> let a = nd.array([1,2,3,4,5,6,7,8,9])
>>> let b = a.reshape(3,-1)
>>> console.log( b.toString() )
  [[1,2,3],
   [4,5,6],
   [7,8,9]]
`



nd.Array.prototype.reduce.__doc__ = `\
Uses the given binary operator to reduce the entries of of this nd.Array
along the specified axes. If no axes are specified all entries are reduce
to a single value and said value is returned instead of an nd.Array.

Parameters
----------
axes: int or int[]
  [OPTIONAL] The axes along which the array is to be reduced. If not defined
  the nd.Array is reduced along all axes and a scalar value is returned instead
  of an nd.Array.
dtype: String
  The data type of the reduced nd.Array. Has to be a super-dtype of the original
  nd.Array.
reducer: (dtype,dtype) => dtype
  The function used to reduce the entries along the axes.

Returns
-------
reduced: nd.Array or dtype
  The reduced nd.Array if axes were specified or the reduced value id no
  axes were specified.

Examples
--------
>>> let a = nd.array([
...   [1,2,3],
...   [4,5,6]
... ])
>>> console.log( a.reduce( (x,y) => x+y ) )
  21


>>> console.log( a.reduce( [0,1], (x,y) => x+y ) )
  { [NDArray: self] shape: [], data: [21] }

>>> console.log( a.reduce( 0, (x,y) => x+y ).toString() )
  [5,7,9]

>>> console.log( a.reduce( [1], (x,y) => x+y ).toString() )
  [6,15]
`



nd.Array.prototype.slice.__doc__ = `\
Extracts a sub-region specified by combination of indices, ranges, newaxis symbols.

Parameters
----------
slices: int or nd.newaxis or nd.ellipsis or [start,stop,step]
  The slices to be taken along each axis, e.g 3 would only take the fourth entries
  along an axis, [,,] would take all elements along and axis, [1,,] would take all
  but the first entries along an axis, [,,3] would take ever third element along
  an axis, [,-1,] would take all but the last entries along an axis, nd.newaxis
  would insert a new axis of size 1.

  nd.ellipsis can be used to fill up with [,,] for the remaining axes.

Returns
-------
sliced: nd.Array
  A sliced subregion of this nd.Array.

Examples
--------
>>> let a = nd.array([
...   [11,12,13,14],
...   [21,22,23,24],
...   [31,32,33,34]
... ])
>>> console.log( a.slice('...', -2).toString() )
  [13, 23, 33]

>>> console.log( a.slice([,,2]).toString() )
  [[11,12,13,14],
   [31,32,33,34]]

>>> console.log( a.slice(2,[1,3,]).toString() )
  [32,33]
`



nd.stack.__doc__ = `\
Arranges a list of nd.Arrays into a new nd.Array. The nd.Arrays are stacked
along a newly axis, inserted at the specified index. All stacked nd.Arrays
must have the same shape.

Parameters
----------
axis: int
  [OPTIONAL] The index at which the new index is to be inserted. Default value: 0.
dtype: String
  [OPTIONAL] The type of the stacked nd.Array. Should be a super-dtype of all nd.Arrays'
  dtypes.
ndarrays: nd.Array[]
  The nd.Arrays that are to be stacked. All nd.Arrays must have the same shape.

Returns
-------
stacked: nd.Array

Examples
--------
>>> let a = nd.array([1,2,3])
>>> let b = nd.array([4,5,6])
>>> console.log( nd.stack([a,b]).toString )
  [[1,2,3],
   [4,5,6]]

>>> console.log( nd.stack(1,'float64',[a,b]).toString )
  [[1.0, 4.0],
   [2.0, 5.0],
   [3.0, 6.0]]
`



nd.concat.__doc__ = `\
Arranges a list of nd.Arrays into a new nd.Array. The nd.Arrays are concatenated
along an existing specified axis. Aside from said axis, the shape of all nd.Arrays
must be the same.

Parameters
----------
axis: int
  [OPTIONAL] The axis along which the nd.Arrays are concatenated. Default value: 0.
dtype: String
  [OPTIONAL] The type of the concatenated nd.Array. Should be a super-dtype of all nd.Arrays'
  dtypes.
ndarrays: nd.Array[]
  The nd.Arrays that are to be concatenated. All nd.Arrays must have the same shape.

Returns
-------
concatenated: nd.Array


Examples
--------
>>> let a = nd.array([[1,2],[3,4]])
>>> let b = nd.array([[5,6],[7,8]])
>>> console.log( nd.concat([a,b]).toString() )
  [[1,2],
   [3,4],
   [5,6],
   [7,8]]
   
>>> console.log( nd.concat(1,[a,b]).toString() )
  [[1,2,5,6],
   [3,4,7,8]]
`
}


