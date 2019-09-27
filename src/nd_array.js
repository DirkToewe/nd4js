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

import {zip_elems} from './zip_elems'
import {ARRAY_TYPES, dtypeof, is_subdtype, super_dtype, _check_dtype} from './dt'


export function array(dtype, content)
{
  if( null == content ){ content = dtype; dtype = undefined }

  if( ('object' == typeof content) && 'shape' in content && 'data' in content )
  {
    let data = content.data
    if( null != dtype && ! (data instanceof ARRAY_TYPES[dtype]) )
      data = ARRAY_TYPES[dtype].from(data)
    return new NDArray( content.shape, data.slice() )
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
        return dtypeof(content)
      else {
        let
          j = shape[d]-1,
          dt = dtype(d+1,content[j])
        for( ; j >= 0; j-- )
          dt = super_dtype(dt, dtype(d+1,content[j]) )
        return dt;
      }
    }(0,content)

  _check_dtype(dtype)
  const data = new ARRAY_TYPES[dtype]( shape.reduce((a,b) => a*b, 1) )
  let idx = 0;

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
  fill(0,content);

  return new NDArray(shape,data)
}

export function asarray( arrayLike ) // <- TODO: add dtype
{
  if( arrayLike instanceof NDArray )
    return arrayLike
  return array(arrayLike)
}

export class NDArray extends Function
{
  static get name() { return 'nd.Array'; }

   //
  // FACTORY METHODS
 //
  constructor(shape, data)
  {
    if(             ! (shape instanceof Int32Array)    ) throw new Error('Shape must be Int32Array.');
    if(                shape.some( s => s < 1 )        ) throw new Error(`Invalid shape: ${shape}.`)
    if( data.length != shape.reduce( (x,y) => x*y, 1 ) ) throw new Error(`Shape [${shape}] does not match array length of ${data.length}.`)
    const self = (...indices) => self.data[self._flat_idx(indices)];
    Object.setPrototypeOf(self, NDArray.prototype)
    Object.freeze(shape.buffer);
    self.shape = shape;
    self.data  = data;
    return self;
  }

   //
  // GENERAL
 //
  get ndim() { return this.shape.length }

  get dtype() {
    for( const [dtype,ArrayType] of Object.entries(ARRAY_TYPES) )
      if( this.data instanceof ArrayType ) return dtype
    throw new Error('Data type could not determined.')
  }

//    TODO implements this in combination with a PROXY
//    get length() {
//      return this.shape[0];
//    }

  set(indices, value) {
    this.data[this._flat_idx(indices)] = value
  }

  modify( indices, modifier ) {
    const i = this._flat_idx(indices);
    this.data[i] = modifier(this.data[i],...indices);
  }

  _flat_idx(indices)
  {
    const shape = this.shape

    if( indices.length != shape.length ) throw new Error(`Multi-index [${indices}] does not have expected length of ${shape.length}.`);
    
    let
      flat_idx = 0,
      stride = 1;
    for( let i=shape.length; i-- > 0; stride *= shape[i] )
    {
      let idx = indices[i];
      if( idx % 1 != 0 ) throw new Error(`Multi-index [${indices}] contains non-integer entries.`);
      if( idx < 0 )  idx += shape[i]
      if( idx < 0 || idx >= shape[i] ) throw new Error(`Multi-index [${indices}] out of bounds [${shape}].`);
      flat_idx  +=   idx * stride;
    }
    return flat_idx;
  }

  [Symbol.for('nodejs.util.inspect.custom')]( depth, options ) {
    return this.toString();
  }

  toString( max_len_fn = (d,shape) => d >= shape.length-2 ? 6 : 6 )
  {
    if( ! (max_len_fn instanceof Function) ) {
      const max_len = max_len_fn*1
      max_len_fn = () => max_len
    }
    const
      shape  = this.shape,
      data   = this.data,
      strides= new Int32Array(shape.length)

    strides[strides.length-1] = 1
    for( let i=strides.length; --i > 0; )
      strides[i-1] = strides[i] * shape[i]

    /** Collects the String representations of all (displayed) entries.
     */
    function* entries( d, idx )
    {
      const max_len = max_len_fn(d,shape),
            last   = max_len >>> 1,
            first  = max_len - last
      if( shape.length === d ) {
        const data_idx = data[idx]
        switch(data_idx) {
          case null     : yield 'null'     ; break
          case undefined: yield 'undefined'; break
          default       : yield data_idx.toString()
        }
      }
      else if( shape[d] > max_len+1 ) {
        for( let j=0;             j < first;    j++ ) yield* entries(d+1, idx+j*strides[d])
        for( let j=shape[d]-last; j < shape[d]; j++ ) yield* entries(d+1, idx+j*strides[d])
      }
      else for( let j=0;          j < shape[d]; j++ ) yield* entries(d+1, idx+j*strides[d])
    }
    entries = [...entries(0,0) ];

    // pad all entries to the same string length
    if( this.ndim > 1 ) {
      const padLen = entries.reduce( (a,b) => Math.max(a,b.length), 0 );
      entries = entries.map( e => e.padStart(padLen) );
    }

    let iEntries = 0;
    function* str( indent, d, idx )
    {
      if( shape.length === d ) { yield entries[iEntries++]; return; }
      indent += ' '
      let
        prefix = '[ ',
        infix  = ', ',
        suffix = ' ]'
      if( d < shape.length-1 ) {
        infix = ',\n'+indent;
        prefix = '[';
        suffix = ']';
      }
      yield prefix;

      const max_len = max_len_fn(d,shape),
            last   = max_len >>> 1,
            first  = max_len - last

      if( shape[d] > max_len+1 ) {
        for( let j=0;             j < first;    j++ ) {                      yield* str(indent, d+1, idx+j*strides[d]); yield infix; }
                                                                             yield `...${shape[d]-max_len} more...`
        for( let j=shape[d]-last; j < shape[d]; j++ ) {         yield infix; yield* str(indent, d+1, idx+j*strides[d]) }
      }
      else for( let j=0;          j < shape[d]; j++ ) { if(j>0) yield infix; yield* str(indent, d+1, idx+j*strides[d]) }
      yield suffix;
    }
    return [...str('',0,0) ].join('');
  }

  toNestedArray() {
    const  data = this.data,
          shape = this.shape;
    let flat_idx = 0
  
    const toNested = d => d === shape.length
      ? data[flat_idx++]
      : Array.from({length: shape[d]}, () => toNested(d+1));

    return toNested(0)
  }

   //
  // ITERATION
 //
  *[Symbol.iterator]() {
    const shape = this.shape.slice(1),
         stride = shape.reduce((a,b) => a*b, 1)
    for( let i=0; i < this.data.length; )
      yield new NDArray( shape, this.data.slice(i,i+=stride) )
  }

  forEach( consumer ) {
    const
      len = this.shape[0],
      shape = this.shape.slice(1),
      stride = shape.reduce((a,b) => a*b, 1)
    for( let i=0, idx=0; i < len; i++, idx += stride )
      consumer( new NDArray( shape, this.data.slice(idx,idx+stride) ), i )
  }
  
  *elems() {
    const
      shape= this.shape,
      data = this.data,
      multi_idx = new Int32Array(shape.length) // <- index in result
    let flat_idx = 0

    function* elems(d) {
      if( d === shape.length )
        yield [ multi_idx.slice(), data[flat_idx++] ]
      else for(
        multi_idx[d] = 0;
        multi_idx[d] < shape[d];
        multi_idx[d]++
      )
        yield* elems(d+1)
    }
    yield* elems(0)
  }

  forElems( consumer ) {
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
    return this;
  }

  mapElems( dtype, mapper ) {
    if( null == mapper ) {
      if( dtype == null ) return new NDArray(this.shape, this.data.slice());
      if( dtype instanceof Function ) { mapper = dtype; dtype = undefined }
    }
    if( null == mapper ) mapper = x => x;
    return zip_elems([this], dtype, mapper)
  }

  get T() {
    if( this.ndim == 0 ) return this.sliceElems() // <- TODO maybe [[]]
    if( this.ndim == 1 ) return this.sliceElems([],'new')
    return this.transpose()
  }

  get H() {
    const result = this.T;
    if( this.dtype != 'object' &&
        this.dtype != 'complex128' ) return result;
    return result.mapElems(math.conj);
  }

  transpose( ...axes ) // <- TODO allow for ellipsis '...' as input
  {
    axes = [...axes]
    const
      newShape = Int32Array.from(this.shape),
      ndim = newShape.length

    let strides = new Int32Array(ndim)
    if( ndim > 0 )
      strides[ndim-1] = 1
    for( let i=ndim; --i > 0; )
      strides[i-1] = newShape[i]*strides[i]

    // BY DEFAULT THE LAST 2 AXES ARE SWAPPED
    if( axes.length == 0 ) {
      if( ndim <= 1 ) return this.sliceElems()
      newShape[ndim-2] = this.shape[ndim-1]
      newShape[ndim-1] = this.shape[ndim-2]
      let tmp = strides[ndim-2]; strides[ndim-2] = strides[ndim-1]; strides[ndim-1] = tmp
    }
    // 
    else {
      const
        set = new Set(axes),
        _strides = strides; strides = new Int32Array(ndim)

      if( set.size != axes.length ) throw new Error('Duplicate axes are not allowed.')

      for( let i=axes.length; i-- > 0; )
      {
        const j = axes[i];
        if( 0 > j || j >= ndim ) throw new Error('Axis out of bounds.')
        newShape[i] = this.shape[j]
        strides[i] = _strides[j]
      }

      // COMPLETE WITH REMAINING/MISSING/IMPLIED INDICES
      for( let i=0, j=set.size; i < ndim; i++ )
        if( ! set.has(i) )
        {
          newShape[j] = this.shape[i]
          strides[j++] = _strides[i]
        }
    }

    const
      oldData =   this.data,
      newData = new oldData.__proto__.constructor(oldData.length)

    let newI = oldData.length

    function copy( d, oldI )
    {
      if( ndim > d )
        for( let i = newShape[d]; i-- > 0; oldI -= strides[d] )
          copy(d+1,oldI)
      else newData[--newI] = oldData[oldI]
    }
    copy(0,newI-1)

    return new NDArray(newShape,newData)
  }

  reshape( ...shape )
  {
    shape = Int32Array.from(shape);
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

  reduceElems( axes, dtype, reducer )
  {
    if( null == reducer )
    {
      if( null    ==  dtype      ) return this.data.reduce(axes)
      if('string' === typeof axes) return this.data.reduce(dtype)
      reducer = dtype; dtype = undefined
    }
    if( null == dtype ) dtype = 'object'

    if( ! is_subdtype(this.dtype, dtype) )
      throw new Error('New dtype must be a super-dtype.')

    const oldNDim = this.ndim

    if( 'number' === typeof axes ) axes = [axes]
    if( axes instanceof NDArray ) {
      if( ! is_subdtype(axes.dtype,'int32') ) throw new Error(`Invalid dtype ${axes.dtype} for axes.`) 
      if( axes.ndim === 1 )
        axes = axes.data
      else
        throw new Error('Only 1D nd.Array allowed for axes.')
    }
    axes = new Set( function*() {
      for( let ax of axes ) {
        if( 0 > ax )  ax += oldNDim
        if( 0 > ax || ax >= oldNDim ) throw new Error('Reduction axis '+ax+' out of bounds.')
        yield +ax
      }
    }() )

    const
      oldShape= this.shape,
      newShape= oldShape.filter( (size,i) => ! axes.has(i) ),
      newData = new ARRAY_TYPES[dtype]( newShape.reduce((a,b) => a*b, 1) ),
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
        for( let j=oldShape[d]; j-- > 0; reduce=true ) // <- copy over very first value and only then start reduction
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

  sliceElems(...slices)
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
           if( slc === '...') ++nEllipses;
      else if( slc === 'new') ++nNew;
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
      slices.push('...')
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
              case 'new':
                strides.push(0)
                yield 1
                break
              case '...':
                for( let i = nReduced + nSliced; i < oldNdim; i++ )
                {
                  strides.push(stride)
                  stride *= oldShape[--d]
                  yield oldShape[d]
                }
                break
              default:
                --d;
                if( 'number' === typeof slc )
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
      newData = new ARRAY_TYPES[this.dtype]( newShape.reduce((a,b) => a*b, 1) )
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
