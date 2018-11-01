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
let nd = {}
try {
  module.exports = nd;
}
catch(err){
  console.log(err);
}

nd.version = '0.1.0';

{
  nd.Complex = class Complex
  {
    constructor( re, im ) {
      if( null == re ) re = 0.0;
      if( null == im ) im = 0.0;
      if( re.re != null ||
          re.im != null ) {
        im = re.im || 0;
        re = re.re || 0;
      }
      this.re = re * 1;
      this.im = im * 1;
      // TODO: remove these assertions
      if( isNaN(re) ) throw new Error(     'Real value is NaN.');
      if( isNaN(im) ) throw new Error('Imaginary value is NaN.');
      Object.seal(this);
    }

    add( re, im ) {
      if( null == re ) re = 0.0;
      if( null == im ) im = 0.0;
      if( re.re != null ||
          re.im != null ) {
        im = re.im || 0;
        re = re.re || 0;
      }
      return new Complex(
        this.re + re,
        this.im + im
      );
    }

    sub( re, im ) {
      if( null == re ) re = 0.0;
      if( null == im ) im = 0.0;
      if( re.re != null ||
          re.im != null ) {
        im = re.im || 0;
        re = re.re || 0;
      }
      return new Complex(
        this.re - re,
        this.im - im
      );
    }

    mul( re, im ) {
      if( null == re ) re = 0.0;
      if( null == im ) im = 0.0;
      if( re.re != null ||
          re.im != null ) {
        im = re.im || 0;
        re = re.re || 0;
      }
      return new Complex(
        this.re*re - this.im*im,
        this.re*im + this.im*re
      );
    }

    div( re, im ) {
      if( null == re ) re = 0.0;
      if( null == im ) im = 0.0;
      if( re.re != null ||
          re.im != null ) {
        im = re.im || 0;
        re = re.re || 0;
      }
      if( im == 0 ) return new Complex(this.re/re, this.im/re);

//      return new Complex(
//        (this.re*re + this.im*im) / (re*re + im*im),
//        (this.im*re - this.re*im) / (re*re + im*im)
//      );

      if( Math.abs(re) >= Math.abs(im) )
      {
        const R = im / re;
        return new Complex(
          (this.re + this.im*R) / (re + im*R),
          (this.im - this.re*R) / (re + im*R)
        );
      } else {
        const R = re / im;
        return new Complex(
          (this.re*R + this.im) / (re*R + im),
          (this.im*R - this.re) / (re*R + im)
        );        
      }
    }

    abs() {
      return Math.hypot(this.re, this.im);
    }

    arg() {
      return Math.atan2(this.im, this.re);
    }

    conj() {
      return new Complex(this.re, -this.im );
    }

    toFixed(digits) {
      return new Complex(
        this.re.toFixed(digits),
        this.im.toFixed(digits)
      );
    }

    sqrt() {
      // https://en.wikipedia.org/wiki/Square_root#Algebraic_formula
      const abs = this.abs();
      return new Complex(
        math.sqrt( math.mul( math.add(abs,this.re), 0.5 ) ),
        math.sqrt( math.mul( math.sub(abs,this.re), 0.5 ) )
      );
    }

    // https://nodejs.org/api/util.html#util_util_inspect_custom
    [Symbol.for('nodejs.util.inspect.custom')]() { return this.toString(); }

    toString() {
      if( this.im == 0 ) return this.re.toString();
      if( this.re == 0 ) return `${this.im}j`;
      if( this.im <  0 ) return `${this.re} - ${-this.im}j`;
      return `${this.re} + ${this.im}j`;
    }

    valueOf() {
      if( this.im == 0 )
        return this.re;
      return this;
    }
  }

/******************
 * IMPLEMENTATION *
 ******************/
  nd.ellipsis= '...'
  nd.newaxis = 'new'
//  nd.all     = '*'

  nd.Array = class NDArray extends Function
  {
//    static get name() { return 'nd.Array'; }

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

    static from( ndarray, dtype, mapper )
    {
      // PREPROCESS ARGUMENTS
      if( null == mapper )
        if( dtype instanceof Function ) { mapper = dtype; dtype = undefined }

      if( ! (ndarray instanceof Array) )
        ndarray = [ndarray]
      else if( null == mapper && ndarray.length > 1 )
        throw new Error('A mapping function is required, when more than 1 ndarray is provided.')

      ndarray = ndarray.map( m => nd.asarray(m) )

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
        strides = new Int32Array(ndarray.length);

      let flat_idx = 0

      function write(d) {
        if( d === ndim ) {
          strides.fill(1)
          for( let i=ndarray.length; i-- > 0; )
            values[i] = ndarray[i].data[indices[i]++]
          data[flat_idx++] = mapper(...values, ...multi_idx)
          return
        }
        for( multi_idx[d] = 0;; )
        {
          write(d+1)
          if( ++multi_idx[d] >= shape[d] )
            break
          for( let i=ndarray.length; i-- > 0; )
            if( ! (ndarray[i].shape[ d - ndim + ndarray[i].ndim ] > 1) ) // <- handles undefined (index out of bounds)
              indices[i] -= strides[i]
        }
        for( let i=ndarray.length; i-- > 0; )
          strides[i] *= ( ndarray[i].shape[ d - ndim + ndarray[i].ndim ] || 1 )
      }
      write(0)

      return new NDArray(shape,data)
    }

     //
    // GENERAL
   //
    get ndim() { return this.shape.length }

    get dtype() {
      for( const [dtype,ArrayType] of Object.entries(nd.dtypes) )
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

      /** Collects the String representations of all (displayed) entries.
       */
      function* entries( d, idx )
      {
        if( shape.length === d ) yield data[idx].toString();
        else if( shape[d] > max_len ) {
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
        if( shape[d] > max_len ) {
          for( let j=0;             j < first;    j++ ) {                      yield* str(indent, d+1, idx+j*strides[d]); yield infix; }
                                                                               yield `...${shape[d]-max_len} more...`
          for( let j=shape[d]-last; j < shape[d]; j++ ) {         yield infix; yield* str(indent, d+1, idx+j*strides[d]) }
        }
        else for( let j=0;          j < shape[d]; j++ ) { if(j>0) yield infix; yield* str(indent, d+1, idx+j*strides[d]) }
        yield suffix;
      }
      return [...str('',0,0) ].join('');
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
    
    *elems() {
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

    mapElems( dtype, mapper ) { return NDArray.from(this, dtype, mapper) }

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

      if( ! nd.is_subdtype(this.dtype, dtype) )
        throw new Error('New dtype must be a super-dtype.')

      const oldNDim = this.ndim

      if( 'number' === typeof axes ) axes = [axes]
      if( axes instanceof nd.Array ) {
        if( ! nd.is_subdtype(axes.dtype,'int32') ) throw new Error(`Invalid dtype ${axes.dtype} for axes.`) 
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

  {
    const
      ch_b64= 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/',
      ch_ws = '\f\n\r\t\v ',
      LUT = new Int32Array(256)
    LUT.fill(-2)

    for( let i=ch_b64.length; i-- > 0; ) LUT[ ch_b64.charCodeAt(i) ] =  i
    for( let i=ch_ws .length; i-- > 0; ) LUT[ ch_ws .charCodeAt(i) ] = -1

    LUT[ '='.charCodeAt(0) ] = 0

    nd.arrayFromB64 = (dtype, shape, content) =>
    {
      shape = Int32Array.from(shape)
  
      const
        result = new nd.dtypes[dtype]( shape.reduce( (m,n) => m*n, 1 ) ),
        buf = new Uint8Array(result.buffer)

      buf[0] = 1
      const bigEndian = ( new Int16Array(result.buffer) )[0] !== 1

      let i=0, j=0
  
      function next6() {
        let result = LUT[content.charCodeAt(j++)]
        // skip whitespaces
        while( result === -1 )
          result = LUT[content.charCodeAt(j++)]

        if( ! (result >= 0) ) // <- handles undefined
          throw new Error(`Illegal Base64 character: '${content.charCodeAt(j-1)}'`)

        return result
      }
  
      while( i < buf.length ) {
        const
          bit0to5  = next6(),
          bit6to11 = next6(),
          bit12to17= next6(),
          bit18to23= next6()
        buf[i++] = ( bit0to5             << 2) | (bit6to11  >> 4)
        buf[i++] = ((bit6to11  & 0b1111) << 4) | (bit12to17 >> 2)
        buf[i++] = ((bit12to17 & 0b0011) << 6) | (bit18to23 & 0b111111)
      }

      if( bigEndian )
        throw new Error('Not yet implemented for BigEndian machines.') // <- FIXME

      return new nd.Array(shape,result)
    }
  }

  nd.asarray = (arrayLike) => // <- TODO: add dtype
  {
    if( arrayLike instanceof nd.Array )
      return arrayLike
    return nd.array(arrayLike)
  }

  nd.array = (dtype, content) =>
  {
    if( null == content ){ content = dtype; dtype = undefined }

    if( ('object' == typeof content) && 'shape' in content && 'data' in content )
    {
      let data = content.data
      if( null != dtype && ! (data instanceof nd.dtypes[dtype]) )
        data = nd.dtypes[dtype].from(data)
      return new nd.Array( content.shape, data.slice() )
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

    return new nd.Array(shape,data)
  }

  nd.tabulate = ( shape, dtype, idx2val ) =>
  {
    shape = Int32Array.from(shape);
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
    if( null == dtype) dtype = nd.super_dtype( ...ndarrays.map( a => a.dtype ) )

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
      indices = Int32Array.from(ndarrays, ndarr => ndarr.data.length),
      newData = new nd.dtypes[dtype]( newShape.reduce((a,b) => a*b, 1) )
    let newIdx = newData.length

    function fill(d)
    {
      if( d === axis )
        for( let i=ndarrays.length; i-- > 0; )
          for( let j=ndarrays[i].shape[d]*rest; j-- > 0; )
            newData[--newIdx] = ndarrays[i].data[--indices[i]]
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
  // UTILITIES
 //
  const util = {}
  nd.util = util
  {
    const ComplexArrayType = FloatArray => {
  
      const HANDLER = {
        get( complexArray, property ) {
          if( typeof(property) !== 'symbol' && property % 1 === 0 )
            return complexArray.get(property);
          return complexArray[property];
        },
        set( complexArray, property, value ) {
          if( typeof(property) !== 'symbol' && property % 1 === 0 )
          {
            complexArray.set(property,value);
            return true; // <- FIXME bounds check?
          }
          complexArray[property] = value;
          return true;
        }
      }
  
      class ComplexArray
      {
        static get BYTES_PER_ELEMENT() {
          return FloatArray.BYTES_PER_ELEMENT*2;
        }
  
        static get name() {
          return 'Complex'+ComplexArray.BYTES_PER_ELEMENT*8+'Array';
        }
  
        constructor( buffer, byteOffset, length )
        {
          if( buffer % 1 === 0 )
          {
            length = buffer;
            this._array = new FloatArray(length*2)
          }
          else {
            if( null != length ) length *= 2;
            this._array = new FloatArray(buffer,byteOffset,length);
          }
          return new Proxy(this,HANDLER)
        }
  
        static from( source, mapFn, thisArg ) {
          if( ! source.hasOwnProperty('length') )
            source = [...source ];
          const result = new ComplexArray(source.length);
          mapFn = mapFn || (x => x);
          for( let i=0; i < source.length; i++ )
            result[i] = mapFn.call(thisArg, source[i], i);
          return result;
        }
  
        static of(...elements )
        {
          const result = new ComplexArray(elements.length)
          for( let i=0; i < elements.length; i++ )
            result[i] = elements[i];
          return result;
        }
  
        get buffer    () { return this._array.buffer;       }
        get byteOffset() { return this._array.byteOffset;   }
        get byteLength() { return this._array.byteLength;   }
        get length    () { return this._array.length >>> 1; }
  
        *[Symbol.iterator]() {
          for( let i=0; i < this.length; i++ )
            yield this[i];
        }
  
        *keys() {
          for( let i=0; i < this.length; i++ )
            yield i;
        }
  
        *values() {
          yield *this[Symbol.iterator]();
        }
  
        *entries() {
          for( let i=0; i < this.length; i++ )
            yield [i, this[i]];
        }
  
        forEach( callback, thisArg ) {
          for( let i=0; i < source.length; i++ )
            callback.call(thisArg, source[i], i);
        }
  
        map( mapFn, thisArg ) {
          return ComplexArray.from(this,mapFn,thisArg);
        }
  
        reduce( reduceFn, initialValue ) {
          let i=0;
          if( null == initialValue ) {
            if( this.length == 0 )
              throw new TypeError('TypeError: Reduce of empty array with no initial value.');
            initialValue = this[i++];
          }
          for( ; i < this.length; i++ )
            initialValue = reduceFn(initialValue,this[i],i,this);
          return initialValue;
        }
  
        slice(begin,end) {
          if( null != begin ) begin *= 2;
          if( null != end   ) end   *= 2;
          return new ComplexArray(this._array.slice(begin,end).buffer);
        }
  
        subarray(begin,end) {
          if( null != begin ) begin *= 2;
          if( null != end   ) end   *= 2;
          const sub = this._array.subarray(begin,end);
          return new ComplexArray(sub.buffer,sub.byteOffset,sub.length/2)
        }
  
        join( separator ) {
          if( null == separator ) separator = ', ';
          let result = '';
          for( let i=0; i < this.length; i++ ) {
            if( i > 0 ) result += separator;
            result += this[i].toString();
          }
          return result;
        }
  
        fill( value, start, end )
        {
          if( ! (value instanceof nd.Complex) ) value = new nd.Complex(value);
  
          if( null == start ) start = 0;
          if( null == end   ) end   = this.length;
          if( 0 > end   ) end   += this.length;
          if( 0 > start ) start += this.length;
  
          for( let i=start; i < end; i++ )
            this[i] = value;
  
          return this;
        }
  
        get( index ) {
          return Object.freeze(
            new nd.Complex(
              this._array[2*index+0],
              this._array[2*index+1]
            )
          )
        }
  
        set( index, value ) {
          if( ! (value instanceof nd.Complex) )
          {
            if( value % 1 === 0 )
            {
              this._array[2*index+0] = value;
              this._array[2*index+1] = 0;
              return;
            }
            else value = new nd.Complex(value);
          }
          this._array[2*index+0] = value.re;
          this._array[2*index+1] = value.im;
        }
  
        [Symbol.for('nodejs.util.inspect.custom')]( depth, options ) {
          return this.toString(options.maxArrayLength);
        }
  
        toString( max_len ) {
          if( null == max_len ) max_len = 16;
          let prefix = ComplexArray.name + ' ['
          if( this.length > max_len )
          {
            const last   = max_len >>> 1,
                  first  = max_len - last - 1;
            for( let i=0; i < first; i++ ) {
              prefix += this[i].toString() + ', '
            }
            prefix += `...${this.length - last - first}more...`
            for( let i=this.length-last; i < this.length; i++ ) {
              prefix += ', ' + this[i].toString();
            }
            return prefix + ']'
          }
          return prefix + this.join(', ') + ']'
        }
      }
  
      return ComplexArray;
    }
    util.Complex128Array = ComplexArrayType(Float64Array);
//    util.KahanSum = class KahanSum
//    {
//      constructor( initVal ) {
//        this[':='](initVal);
//      }
//      [':=']( value ) {
//        this.sum = initVal;
//        this.rem = 0.0;
//      }
//      ['+=']( summand ) {
//        throw new Error('Not yet implemented.');
//      }
//      ['-=']( subtrahend ) {
//        throw new Error('Not yet implemented.');
//      }
//      valueOf() {
//        return this.sum;
//      }
//    }
  }

   //
  // DATA TYPE OPERATIONS
 //
  nd.dtypes = {
         'int32':              Int32Array,
       'float32':            Float32Array,
       'float64':            Float64Array,
    'complex128': nd.util.Complex128Array,
       'object' :                   Array
  }

  nd.eps = {
    'float32': 1.1920928955078125e-7,
    'float64': 2.2204460492503130808472633361816E-16
  }

  nd.dtypeof = value => {
    if( value % 1 === 0 )
    {
      if(    value <= ~(1 << 31)
          && value >=  (1 << 31) )
        return 'int32'
      return 'object'
    }
    if( value * 1 == value ) return 'float64'
    if( value instanceof nd.Complex ) return 'complex128';
    return 'object'
  }

  nd.super_dtype = (...dtypes) => dtypes.reduce( (dtype1,dtype2) => {
    nd._check_dtype(dtype1)
    nd._check_dtype(dtype2)
    if( dtype1 ===     'object' || dtype2 ===     'object' ) return 'object'
    if( dtype1 === 'complex128' || dtype2 === 'complex128' ) return 'complex128'
    if( dtype1 ===    'float64' || dtype2 ===    'float64' ) return 'float64'
    if( dtype1 ===    'float32' || dtype2 ===    'float32' ) return 'float32'
    return 'int32'
  });

  nd.is_subdtype = (sub_dtype, sup_dtype) => {
    nd._check_dtype(sub_dtype)
    nd._check_dtype(sup_dtype)
    const rank = {
           'int32': 0,
         'float32': 1,
         'float64': 2,
      'complex128': 3,
          'object': 4
    }
    return rank[sub_dtype] <= rank[sup_dtype]
  }

  nd._check_dtype = dtype => {
    if( ! nd.dtypes.hasOwnProperty(dtype) )
      throw new Error("Invalid dtype '" + dtype + "'. Must be one of {'" + Object.getOwnPropertyNames(nd.dtypes).join("', '") + "'}.")
  }

   //
  // SCALAR MATH
 //
  const math = {}
  nd.math = math
  {
    math.nextUp = x => {
      // FIXME implement this completely
      // https://gist.github.com/Yaffle/4654250
      const
        EPSILON   = Number.EPSILON,
        MIN_VALUE = Number.MIN_VALUE;
      let y = x * (x < 0 ? 1 - EPSILON/2 : 1 + EPSILON);
      if( ! isFinite(x) ) throw new Error('Assertion failed!');
      if( ! isFinite(y) ) throw new Error('Assertion failed!');
      if (y === x) y += MIN_VALUE;
      let b = x + (y - x) / 2; if (x < b && b < y) { y = b; }
      let c =     (y + x) / 2; if (x < c && c < y) { y = c; }
      return y === 0 ? -0 : y;
    };

    math.add = (x,y) => {
      if( x instanceof nd.Complex ) return x.add(y);
      if( y instanceof nd.Complex ) return y.add(x);
      return x + y;
    }

    math.sub = (x,y) => {
      if( x instanceof nd.Complex ) return                x .sub(y);
      if( y instanceof nd.Complex ) return new nd.Complex(x).sub(y);
      return x - y;
    }

    math.mul = (x,y) => {
      if( x instanceof nd.Complex ) return x.mul(y);
      if( y instanceof nd.Complex ) return y.mul(x);
      return x * y;
    }

    math.cast = (x,dtype) => {
      if( dtype ===      'int32' ) return x & 0xFFFFFFFF;
      if( dtype ===    'float32' ) return Math.fround(x);
      if( dtype === 'complex128' ) return x instanceof nd.Complex ? x : new nd.Complex(x);
      return x;
    }

    math.zero = dtype => math.cast(0,dtype);
    math.one  = dtype => math.cast(1,dtype);

    math.div = (x,y) => {
      if( x instanceof nd.Complex ) return                x .div(y);
      if( y instanceof nd.Complex ) return new nd.Complex(x).div(y);
      return x / y;
    }

    math.neg = x => x instanceof nd.Complex ? x.neg() : -x;

    math.abs = x => x instanceof nd.Complex ? x.abs() : Math.abs(x);

    math.sqrt = x => x instanceof nd.Complex
      ?   x.sqrt()
      : ( x >= 0
        ?              Math.sqrt(x)
        : new nd.Complex(x).sqrt()
      );

    math.exp = x => x instanceof nd.Complex ? x.exp() : Math.exp(x);

    math.min = (x,y) => Math.min(x,y);
    math.max = (x,y) => Math.max(x,y);

    math.hypot = (x,y) => Math.hypot(x,y);

    math.atan2 = (x,y) => Math.atan2(x,y);

    math.conj = x => x instanceof nd.Complex ? x.conj() : x;

    math.is_equal = (x,y) => {
      if( x instanceof nd.Complex ) return x.equals(y);
      if( y instanceof nd.Complex ) return y.equals(x);
      return x == y;
    }

    math.is_close = (x,y) => {
      const atol = 1e-8,
            rtol = 1e-5,
             tol = atol + rtol * math.max(
              math.abs(x),
              math.abs(y)
            );
      return math.abs( math.sub(x,y) ) <= tol;
    }

    math.not_close = (x,y) => ! math.is_close(x,y);
  }

   //
  // OPTIMIZATION
 //
  const opt = {}
  nd.opt = opt
  {
    opt.root1d = (F,x_min,x_max) => {
      let f_min = F(x_min),
          f_max = F(x_max);
      if( f_min > f_max ) [f_min,f_max,x_min,x_max] = [f_max,f_min,x_max,x_min];
      if( f_min > 0 ) { if( f_min <= +Number.EPSILON ) return x_min; else throw new Error('Assertion failed.') }
      if( f_max < 0 ) { if( f_max >= -Number.EPSILON ) return x_max; else throw new Error('Assertion failed.') }
      if( f_min > 0 ) { return x_min; }
      if( f_max < 0 ) { return x_max; }

      for(;;) {
        const x = (x_min + x_max) / 2;
        if( x_min == x ||
            x_max == x ) break;
        const f = F(x);
        if( isNaN(f) ) throw new Error('NaN encountered.');
        if( f <= 0 ) x_min = x;
        if( f >= 0 ) x_max = x;
      }

      f_min = F(x_min);
      f_max = F(x_max);
      return Math.abs(f_min) <= Math.abs(f_max) ? x_min : x_max;
    }
  }


   //
  // LINEAR ALGEBRA
 //
  const la = {}
  nd.la = la
  {
    la.eye = (...shape) => { // TODO: Add dtype
      const dtype = shape[0] in nd.dtypes ? shape.shift() : 'float64';

      if( shape.length <  1 ) throw new Error('Size parameter missing.');
      if( shape.length == 1 ) shape.push( shape[shape.length-1] );

      return nd.tabulate(shape, dtype, (...idx) => {
        const [i,j] = idx.slice(-2);
        return i==j ? 1 : 0;
      });
    };

    la.tril = (m,k) => {
      if( undefined == k ) k = 0;
      m = nd.asarray(m);
      if( m.ndim < 2 ) throw new Error('Input must be at least 2D.');
      return m.mapElems( m.dtype, (x,...indices) => {
        const [i,j] = indices.slice(-2);
        return i < j-k ? 0 : x
      });
    };


    la.triu = (m,k) => {
      if( undefined == k ) k = 0;
      m = nd.asarray(m);
      if( m.ndim < 2 ) throw new Error('Input must be at least 2D.');
      return m.mapElems( m.dtype, (x,...indices) => {
        const [i,j] = indices.slice(-2);
        return i > j-k ? 0 : x
      });
    };


    la.diag_mat = (diag) => {
      diag = nd.asarray(diag);
      const
        shape = diag.shape,
        N = shape[shape.length-1],
        d = diag.data,
        D = new nd.dtypes[diag.dtype]( shape.reduce( (a,b) => a*b, 1 )*N );
      diag = undefined;

      if( N <= 0 ) throw new Error('Assertion Failed!');

      for( let d_off=0,
               D_off=0;; )
      {
        const d_end = d_off+N;
        while(true) {
          D[D_off] = d[d_off];
          if( ++d_off == d_end ) break;
          D_off += N+1;
        }
        if( ++D_off == D.length ) break;
      }
      return new nd.Array(Int32Array.of(...shape,N), D)
    };


    la.diag = (A, offset) => {
      A = nd.asarray(A);
      if( null == offset ) offset = 0;
      const
        [N,M] = A.shape.slice(-2),
        shape = A.shape.slice(0,-1),
        DTypeArray = nd.dtypes[A.dtype];
      A = A.data;

      if( offset < -N ) throw new Error(`Offset ${offset} out of shape [${A.shape}].`);
      if( offset > +M ) throw new Error(`Offset ${offset} out of shape [${A.shape}].`);

      const L = Math.min(
        N, N+offset,
        M, M-offset
      );
      shape[shape.length-1] = L;

      const
        d = new DTypeArray( shape.reduce( (a,b) => a*b, 1 ) ),
        a_off = Math.max( 0, offset, -offset*M );

      for( let d_off=0,
               A_off=0; A_off < A.length;
               d_off += L,
               A_off += N*M
      )
      {
        for( let i=0; i < L; i++ )
          d[d_off+i] = A[A_off + a_off + i*(M+1)]
      }
      return new nd.Array(shape, d)
    };


    la._transpose_inplace = A => {
      const [N,M] = A.shape.slice(-2);
      if( N != M ) throw new Error('In-place transposition is only supported for square matrices.');
      A = A.data;
      for( let off=0; off < A.length; off += N*N )
      for( let i=0  ; i < N-1; i++ )
      for( let j=1+i; j < N  ; j++ ) {
        const
          ij = off + N*i+j,
          ji = off + N*j+i,
          A_ij = A[ij]; A[ij] = A[ji]; A[ji] = A_ij;
      }
    }


    la.matmul = (...matrices) => {
      matrices = matrices.map(nd.asarray)
      if( matrices.length == 1 ) return matrices[0];
      if( matrices.length == 2 ) return la.matmul2(...matrices);

      /** Returns the number of floating point operations necessary to
       *  matrix multiply two arrays of the given shapes.
       */
      function nOps( shapeA, shapeB ) {
        const [I,K] = shapeA.slice(-2),
               J    = shapeB[shapeB.length-1];
        if( shapeB[shapeB.length-2] != K )
          throw new Error('Shape mismatch.');

        const
          ndim = Math.max(shapeA.length, shapeB.length),
          shape = Int32Array.from({ length: ndim }, () => 1 );
        shape[ndim-2] = I;
        shape[ndim-1] = J;

        // FIND COMMON (BROADCASTED) SHAPE
        for( let shp of [shapeA,shapeB] )
          for( let i=ndim-2, j=shp.length-2; i-- > 0 && j-- > 0; )
            if( 1 === shape[i] )
              shape[i] = shp[j];
            else if( shape[i] != shp[j] && shp[j] != 1 )
              throw new Error('Shapes are not broadcast-compatible.');

        return [ shape.reduce( (a,b) => a*b, 1 )*K, shape ]
      };
      // https://en.wikipedia.org/wiki/Matrix_chain_multiplication
      // https://www.geeksforgeeks.org/matrix-chain-multiplication-dp-8/
      
      // opMap[i][j] (j >= i) caches the optimal number of FLOPs required to multiply matrices[i] up to (including) matrices[j]
      const opMap = Array.from({ length: matrices.length }, () => [] );

      // initialized opMap
      for( let i=0; i < matrices.length; i++ ) opMap[i][i] = [ 0, matrices[i].shape ]
      // compute remaining opMap
      for( let len=2; len <= matrices.length;   len++ )
      for( let  i =0;  i  <= matrices.length-len; i++ )
      {
        let
          minFlops = Infinity,
          minShape;
        for( let j=1; j < len; j++ )
        {
          const [lFlops,lShape] = opMap[i  ][i+ j -1]
          const [rFlops,rShape] = opMap[i+j][i+len-1]
          let [flops,shape] = nOps(lShape,rShape)
          flops += lFlops + rFlops;
          if( flops < minFlops ) {
            minFlops = flops;
            minShape = shape; // <- the shape should always be the same so this is not strictly necessary
          }
        }
        if( minShape === undefined ) throw new Error('Integer overflow (too many FLOPs).');
        opMap[i][i+len-1] = [ minFlops, minShape ]
      }

      // compute the result using the minimal number of FLOPs using opMap
      function product(from, to)
      {
        if( from == to ) return matrices[from];
        let
          minFlops = Infinity,
          minIdx;
        for( let i=from; i < to; i++ )
        {
          const [lFlops,lShape] = opMap[from][i]
          const [rFlops,rShape] = opMap[i+1][to]
          let [flops,_] = nOps(lShape,rShape)
          flops += lFlops + rFlops;
          if( flops < minFlops ) {
            minFlops = flops;
            minIdx   = i;
          }
        }
        return la.matmul2(
          product(from, minIdx),
          product(minIdx+1, to)
        );
      }

      return product(0, matrices.length-1)
    }


    la.matmul2 = (a,b) => {
      a = nd.asarray(a)
      b = nd.asarray(b)
      if( a.ndim < 2 ) throw new Error('A must be at least 2D.');
      if( b.ndim < 2 ) throw new Error('B must be at least 2D.');

      const
        [I,K] = a.shape.slice(-2),
         J    = b.shape[b.ndim-1];
      if( b.shape[b.ndim-2] != K )
        throw new Error('The last dimension of A and the 2nd to last dimension of B do not match.');

      const
        ndim = Math.max(a.ndim, b.ndim),
        shape = Int32Array.from({ length: ndim }, () => 1 );
      shape[ndim-2] = I;
      shape[ndim-1] = J;

      // FIND COMMON (BROADCASTED) SHAPE
      for( let arr of [a,b] )
        for( let i=ndim-2, j=arr.ndim-2; i-- > 0 && j-- > 0; )
          if( 1 === shape[i] )
            shape[i] = arr.shape[j];
          else if( shape[i] != arr.shape[j] && arr.shape[j] != 1 )
            throw new Error('Shapes are not broadcast-compatible.');

      // GENERATE RESULT DATA
      const DTypeArray = nd.dtypes[nd.super_dtype(a.dtype,b.dtype)];
      const C = new DTypeArray( shape.reduce((a,b) => a*b, 1) );
      C.fill(0.0);

      const
        A = a.data,
        B = b.data;
      let
        aOff = 0, aStride = 1,
        bOff = 0, bStride = 1,
        cOff = 0;

      const loops = new Int32Array(ndim-2);
      for( let d=0; d >= 0; )
        if( d === ndim-2 ) {
          aStride = I*K;
          bStride = K*J;
          for( const cEnd = cOff + I*J; cOff < cEnd; cOff += J, bOff -= bStride ) 
          for( const aEnd = aOff + K  ; aOff < aEnd; cOff -= J, aOff++ ) 
          for( const bEnd = bOff +   J; bOff < bEnd; cOff += 1, bOff++ )
            C[cOff] = math.add( C[cOff], math.mul(A[aOff],B[bOff]) );
          bOff += bStride;
          d -= 1;
        }
        else
        {
          if( loops[d]++ > 0 ) {
            if( loops[d] > shape[d] ) {
              aStride *= a.shape[ d - ndim + a.ndim ] || 1;
              bStride *= b.shape[ d - ndim + b.ndim ] || 1;
              loops[d--] = 0;
              continue;
            }
            if( ! (a.shape[ d - ndim + a.ndim ] > 1) ) aOff -= aStride;
            if( ! (b.shape[ d - ndim + b.ndim ] > 1) ) bOff -= bStride;
          }
          ++d;
        }

      return new nd.Array(shape,C);
    }

    la.cholesky_decomp = S =>
    {
      S = nd.asarray(S);
      const
        shape = S.shape,
        [N,M] =   shape.slice(-2),
        L = Float64Array.from(S.data);
      S = undefined;

      if( N != M )
        throw new Error('Last two dimensions must be quadratic.')

      for( let off=0; off < L.length; off += N*N )
        // https://de.wikipedia.org/wiki/Cholesky-Zerlegung
        for( let i=0; i<N; i++ )
        for( let j=0; j<N; j++ )
          if( i < j ) {
            L[off+N*i+j] = 0;
          } else {
            let
              sum = L[off+N*i+j],
              cor = 0.0;
            for( let k=0; k<j; k++ ) {
              // https://en.wikipedia.org/wiki/Kahan_summation_algorithm
              const
                r = cor + L[off+N*i+k] * L[off+N*j+k],
                s = sum - r;
              cor = (s - sum) + r;
              sum =  s;
            }
            if( i > j ) L[off+N*i+j] = sum / L[off+N*j+j];
            else {      L[off+N*i+i] = Math.sqrt(sum);
              if( isNaN(L[off+N*i+i]) )
                throw new Error('Matrix contains NaNs or is (near) singular.');
            }
          }

      return new nd.Array(shape,L);
    }


    la.cholesky_solve = (L,y) => {
      L = nd.asarray(L);
      y = nd.asarray(y);
      if( L.ndim < 2 ) throw new Error('L must be at least 2D.');
      if( y.ndim < 2 ) throw new Error('y must be at least 2D.');

      const
        [N,M] = L.shape.slice(-2),
        [I,J] = y.shape.slice(-2);
      if( N != M ) throw new Error('Last two dimensions of L must be quadratic.')
      if( I != M ) throw new Error("L and y don't match.");

      const
        ndim = Math.max(L.ndim, y.ndim),
        shape = Int32Array.from({ length: ndim }, () => 1 );
      shape[ndim-2] = I;
      shape[ndim-1] = J;

      // FIND COMMON (BROADCASTED) SHAPE
      for( let arr of [L,y] )
        for( let i=ndim-2, j=arr.ndim-2; i-- > 0 && j-- > 0; )
          if( 1 === shape[i] )
            shape[i] = arr.shape[j];
          else if( shape[i] != arr.shape[j] && arr.shape[j] != 1 )
            throw new Error('Shapes are not broadcast-compatible.');

      // GENERATE RESULT DATA
      const
        dtype = nd.super_dtype(L.dtype,y.dtype),
        x_dat = new nd.dtypes[dtype]( shape.reduce((a,b) => a*b, 1) ),
        L_dat = L.data,
        y_dat = y.data;
      let
        L_off = 0, L_stride = 1,
        y_off = 0, y_stride = 1,
        x_off = 0;

      function solv(d) {
        if( d === ndim-2 ) {
          L_stride = N*N;
          y_stride = N*J;

          // COPYING y
          for( let i=0; i < y_stride; i++ )
            x_dat[x_off+i] = y_dat[y_off+i]

          // FORWARD SUBSTITUTION
          for( let i=0; i < I; i++ )
          for( let j=0; j < J; j++ )
          {
            for( let k=0; k < i; k++ )
              x_dat[x_off+i*J+j] -= L_dat[L_off+N*i+k] * x_dat[x_off+k*J+j]
            x_dat[x_off+i*J+j] /= L_dat[L_off+N*i+i]
          }

          // BACKWARD SUBSTITUTION
          for( let i=I; i-- > 0; )
          for( let j=J; j-- > 0; )
          {
            x_dat[x_off+i*J+j] /= L_dat[L_off+N*i+i]
            for( let k=i; k-- > 0; )
              x_dat[x_off+k*J+j] -= L_dat[L_off+N*i+k] * x_dat[x_off+i*J+j]
          }

          L_off += L_stride;
          y_off += y_stride;
          x_off += y_stride;

          return;
        }
        for( let l=shape[d]; ; l-- ) {
          solv(d+1);
          if( l == 1 ) break;
          if( ! (L.shape[ d - ndim + L.ndim ] > 1) ) L_off -= L_stride;
          if( ! (y.shape[ d - ndim + y.ndim ] > 1) ) y_off -= y_stride;
        }
        L_stride *= L.shape[ d - ndim + L.ndim ] || 1;
        y_stride *= y.shape[ d - ndim + y.ndim ] || 1;
      }
      solv(0);

      return new nd.Array(shape,x_dat);
    }


    la.tril_solve = (L,y) => {
      if( L.ndim < 2 ) throw new Error('L must be at least 2D.');
      if( y.ndim < 2 ) throw new Error('y must be at least 2D.');

      const
        [N,M] = L.shape.slice(-2),
        [I,J] = y.shape.slice(-2);
      if( N != M ) throw new Error('Last two dimensions of L must be quadratic.')
      if( I != M ) throw new Error("L and y don't match.");

      const
        ndim = Math.max(L.ndim, y.ndim),
        shape = Int32Array.from({ length: ndim }, () => 1 );
      shape[ndim-2] = I;
      shape[ndim-1] = J;

      // FIND COMMON (BROADCASTED) SHAPE
      for( let arr of [L,y] )
        for( let i=ndim-2, j=arr.ndim-2; i-- > 0 && j-- > 0; )
          if( 1 === shape[i] )
            shape[i] = arr.shape[j];
          else if( shape[i] != arr.shape[j] && arr.shape[j] != 1 )
            throw new Error('Shapes are not broadcast-compatible.');

      // GENERATE RESULT DATA
      const
        dtype = nd.super_dtype(L.dtype,y.dtype),
        x_dat = new nd.dtypes[dtype]( shape.reduce((a,b) => a*b, 1) ),
        L_dat = L.data,
        y_dat = y.data;
      let
        L_off = 0, L_stride = 1,
        y_off = 0, y_stride = 1,
        x_off = 0;

      function solv(d) {
        if( d === ndim-2 ) {
          L_stride = N*N;
          y_stride = N*J;

          // COPYING y
          for( let i=0; i < y_stride; i++ )
            x_dat[x_off+i] = y_dat[y_off+i]

          // FORWARD SUBSTITUTION
          for( let i=0; i < I; i++ )
          for( let j=0; j < J; j++ )
          {
            for( let k=0; k < i; k++ )
              x_dat[x_off+i*J+j] -= L_dat[L_off+N*i+k] * x_dat[x_off+k*J+j]
            x_dat[x_off+i*J+j] /= L_dat[L_off+N*i+i]
          }

          L_off += L_stride;
          y_off += y_stride;
          x_off += y_stride;

          return;
        }
        for( let l=shape[d]; ; l-- ) {
          solv(d+1);
          if( l == 1 ) break;
          if( ! (L.shape[ d - ndim + L.ndim ] > 1) ) L_off -= L_stride;
          if( ! (y.shape[ d - ndim + y.ndim ] > 1) ) y_off -= y_stride;
        }
        L_stride *= L.shape[ d - ndim + L.ndim ] || 1;
        y_stride *= y.shape[ d - ndim + y.ndim ] || 1;
      }
      solv(0);

      return new nd.Array(shape,x_dat);
    }


    la.triu_solve = (U,y) => {
      if( U.ndim < 2 ) throw new Error('U must be at least 2D.');
      if( y.ndim < 2 ) throw new Error('y must be at least 2D.');

      const
        [K,N] = U.shape.slice(-2),
        [I,J] = y.shape.slice(-2);
      if( K != N ) throw new Error('Last two dimensions of U must be quadratic.')
      if( I != N ) throw new Error("U and y don't match.");

      const
        ndim = Math.max(U.ndim, y.ndim),
        shape = Int32Array.from({ length: ndim }, () => 1 );
      shape[ndim-2] = I;
      shape[ndim-1] = J;

      // FIND COMMON (BROADCASTED) SHAPE
      for( let arr of [U,y] )
        for( let i=ndim-2, j=arr.ndim-2; i-- > 0 && j-- > 0; )
          if( 1 === shape[i] )
            shape[i] = arr.shape[j];
          else if( shape[i] != arr.shape[j] && arr.shape[j] != 1 )
            throw new Error('Shapes are not broadcast-compatible.');

      // GENERATE RESULT DATA
      const x_dat = new nd.dtypes[nd.super_dtype(U.dtype,y.dtype)]( shape.reduce((a,b) => a*b, 1) );

      const
        U_dat = U.data,
        y_dat = y.data;
      let
        U_off = 0, U_stride = 1,
        y_off = 0, y_stride = 1,
        x_off = 0;

      function solv(d) {
        if( d === ndim-2 ) {
          U_stride = N*N;
          y_stride = N*J;

          // COPYING y
          for( let i=0; i < y_stride; i++ )
            x_dat[x_off+i] = y_dat[y_off+i]

          // BACKWARD SUBSTITUTION
          for( let i=I; i-- > 0; )
          for( let j=J; j-- > 0; )
          {
            for( let k=K; --k > i; )
              x_dat[x_off+i*J+j] -= U_dat[U_off+N*i+k] * x_dat[x_off+k*J+j]
            x_dat[x_off+i*J+j] /= U_dat[U_off+N*i+i]
          }

          U_off += U_stride;
          y_off += y_stride;
          x_off += y_stride;

          return;
        }
        for( let l=shape[d]; ; l-- ) {
          solv(d+1);
          if( l == 1 ) break;
          if( ! (U.shape[ d - ndim + U.ndim ] > 1) ) U_off -= U_stride;
          if( ! (y.shape[ d - ndim + y.ndim ] > 1) ) y_off -= y_stride;
        }
        U_stride *= U.shape[ d - ndim + U.ndim ] || 1;
        y_stride *= y.shape[ d - ndim + y.ndim ] || 1;
      }
      solv(0);

      return new nd.Array(shape,x_dat);
    }


    la.lu_decomp = A =>
    {
      const
        [N,M] = A.shape.slice(-2),
         LU_dat=    Float64Array.from(A.data),
          P_dat = new Int32Array(LU_dat.length/M);

      if( N != M )
        throw new Error('Last two dimensions must be quadratic.')

      for(
        let LU_off=0,
             P_off=0;
        LU_off < LU_dat.length;
        LU_off += N*N,
         P_off += N
      )
      {
        for( let i=0; i < N; i++ ) P_dat[P_off+i] = i;
        for( let i=0; i < N; i++ )
        {
          const row_i = LU_off + i*N;
          // ROW PIVOTING
          {
            let p=i;
            for( let j=i+1; j < N; j++ )
              if( Math.abs( LU_dat[LU_off + N*j+i] )
                > Math.abs( LU_dat[LU_off + N*p+i] ) )
                p=j;

            if( i != p )
            {
              const   P_p = P_dat[P_off+i]; P_dat[P_off+i] = P_dat[P_off+p]; P_dat[P_off+p] = P_p; // KEEP TRACK OF ROW SWAPS
              const row_p = LU_off + p*N;
              // SWAP ROWS
              for( let j=0; j < N; j++ ) {
                const tmp = LU_dat[row_i+j]; LU_dat[row_i+j] = LU_dat[row_p+j]; LU_dat[row_p+j] = tmp;
              }
            }
          }
          // ELIMINATE ELEMENTS BELOW PIVOT
          for( let j=i+1; j < N; j++ )
          {
            const
              row_j = LU_off + j*N,
              scale = LU_dat[row_j+i] / LU_dat[row_i+i];
            LU_dat[row_j+i] = scale;
            for( let k=i+1; k < N; k++ )
              LU_dat[row_j+k] -= scale * LU_dat[row_i+k];
          }
        }
      }

      return [
        new nd.Array(A.shape,           LU_dat),
        new nd.Array(A.shape.slice(0,-1),P_dat)
      ];
    }


    la.lu_solve = (LU,P,y) => {
      if( undefined == y ) { y=P; [LU,P] = LU; }
      if( LU.ndim < 2 ) throw new Error('LU must be at least 2D.');
      if( P .ndim < 1 ) throw new Error( 'P must be at least 1D.');
      if( y .ndim < 2 ) throw new Error( 'y must be at least 2D.');

      const
        [N,M] = LU.shape.slice(-2),
        [I,J] =  y.shape.slice(-2);
      if( M != N ) throw new Error('Last two dimensions of LU must be quadratic.')
      if( M != I ) throw new Error("LU and y don't match.");
      if( M != P.shape.slice(-1) )
        throw new Error("LU and P don't match.")

      const
        ndim = Math.max(LU.ndim, P.ndim+1, y.ndim),
        shape = Int32Array.from({ length: ndim }, () => 1 );
      shape[ndim-2] = I;
      shape[ndim-1] = J;

      // FIND COMMON (BROADCASTED) SHAPE
      for( let arr of [LU,y] )
        for( let i=ndim-2, j=arr.ndim-2; i-- > 0 && j-- > 0; )
          if( 1 === shape[i] )
            shape[i] = arr.shape[j];
          else if( shape[i] != arr.shape[j] && arr.shape[j] != 1 )
            throw new Error('LU and y are not broadcast-compatible.');

      for( let i=ndim-2, j=P.ndim-1; i-- > 0 && j-- > 0; )
        if( 1 === shape[i] )
          shape[i] = P.shape[j]
        else if( shape[i] != P.shape[j] && P.shape[j] != 1 )
          throw new Error('P is not broadcast-compatible.');

      // GENERATE RESULT DATA
      const x_dat = new nd.dtypes[nd.super_dtype(LU.dtype,P.dtype,y.dtype)]( shape.reduce((a,b) => a*b, 1) );

      const
        LU_dat = LU.data,
         P_dat =  P.data,
         y_dat =  y.data;
      let
        LU_off = 0, LU_stride = 1,
         P_off = 0,  P_stride = 1,
         y_off = 0,  y_stride = 1,
         x_off = 0;

      function solv(d) {
        if( d === ndim-2 ) {
          LU_stride = N*N;
           P_stride = N;
           y_stride = N*J;

          // COPYING PERMUTED y
          for( let i=0; i < I; i++ )
          {
            const
              row_x = x_off + J * i,
              row_y = y_off + J * P_dat[P_off+i];
            for( let j=0; j < J; j++ )
              x_dat[row_x+j] = y_dat[row_y+j];
          }

          // FORWARD SUBSTITUTION
          for( let i=0; i < I; i++ )
          for( let j=0; j < J; j++ )
          for( let k=0; k < i; k++ )
            x_dat[x_off+i*J+j] -= LU_dat[LU_off+N*i+k] * x_dat[x_off+k*J+j]

          // BACKWARD SUBSTITUTION
          for( let i=I; i-- > 0; )
          for( let j=J; j-- > 0; )
          {
            for( let k=N; --k > i; )
              x_dat[x_off+i*J+j] -= LU_dat[LU_off+N*i+k] * x_dat[x_off+k*J+j]
            x_dat[x_off+i*J+j] /= LU_dat[LU_off+N*i+i]
          }

          LU_off += LU_stride;
           P_off +=  P_stride;
           y_off +=  y_stride;
           x_off +=  y_stride;

          return;
        }
        for( let l=shape[d]; ; l-- ) {
          solv(d+1);
          if( l == 1 ) break;
          if( ! (LU.shape[ d - ndim + LU.ndim   ] > 1) ) LU_off -= LU_stride;
          if( ! ( P.shape[ d - ndim +  P.ndim+1 ] > 1) )  P_off -=  P_stride;
          if( ! ( y.shape[ d - ndim +  y.ndim   ] > 1) )  y_off -=  y_stride;
        }
        LU_stride *= LU.shape[ d - ndim + LU.ndim   ] || 1;
         P_stride *=  P.shape[ d - ndim +  P.ndim+1 ] || 1;
         y_stride *=  y.shape[ d - ndim +  y.ndim   ] || 1;
      }
      solv(0);

      return new nd.Array(shape,x_dat);
    }


    la.qr_decomp_full = A => {
      A = nd.asarray(A);
      if( A.ndim < 2 ) throw new Error('A must be at least 2D.');
      const
        DTypeArray = nd.dtypes[ nd.super_dtype(A.dtype,'float64') ], // <- ensure at least double precision
        R_shape = A.shape,
        Q_shape = Int32Array.from(R_shape),
        [N,M]   = R_shape.slice(-2),
         R = DTypeArray.from(A.data);
      A = undefined
      Q_shape[Q_shape.length-1] = N;
      const Q = new DTypeArray(R.length/M*N);
      Q.fill(0); // <- in case of an object array

      for(
        let Q_off=0,
            R_off=0;
        Q_off < Q.length;
        Q_off += N*N,
        R_off += N*M
      )
      {
        // INIT Q TO IDENTITY MATRIX
        for( let i=0; i < N; i++ ) Q[Q_off+N*i+i] = 1.0;

        for( let i=1; i < N; i++ ) { const I = Math.min(i,M);
        for( let j=0; j < I; j++ )
        {
          // USE GIVENS ROTATION TO ELIMINATE ELEMENT R_ji
          const R_ij = R[R_off+M*i+j]; if( R_ij == 0.0 ) continue;
          const R_jj = R[R_off+M*j+j],
                       norm = Math.hypot(R_jj,R_ij),
            c = R_jj / norm,
            s = R_ij / norm;
          R[R_off+M*i+j] = 0;
          R[R_off+M*j+j] = norm;
          // ROTATE ROW i AND j IN R
          for( let k=j; ++k < M; ) {
            const ik = R_off+M*i+k, R_ik = R[ik],
                  jk = R_off+M*j+k, R_jk = R[jk];
            R[ik] = c*R_ik - s*R_jk;
            R[jk] = s*R_ik + c*R_jk;
          }
          // ROTATE COL i AND j IN Q (Q TRANSPOSED FOR CACHE LOCALITY REASONS) 
          for( let k=0; k <= i; k++ ) {
            const ik = Q_off+N*i+k, Q_ik = Q[ik],
                  jk = Q_off+N*j+k, Q_jk = Q[jk];
            Q[ik] = c*Q_ik - s*Q_jk;
            Q[jk] = s*Q_ik + c*Q_jk;
          }
        }}
        // TRANSPOSE Q (Q TRANSPOSED FOR CACHE LOCALITY REASONS)
        for( let i=0;   i < N; i++ )
        for( let j=i+1; j < N; j++ ) {
          const
            ij = Q_off+N*i+j,
            ji = Q_off+N*j+i,
            Q_ij = Q[ij]; Q[ij] = Q[ji]; Q[ji] = Q_ij;
        }
      }

      return [
        new nd.Array(Q_shape, Q),
        new nd.Array(R_shape, R)
      ];
    }

    la.solve = (A,y) => la.qr_solve(la.qr_decomp(A),y) // <- FIXME: does not work for rank-deficient matrices. Maybe split la.solve into la.solve and la.lstsq would be safer
    la.qr_decomp = A => {
      A = nd.asarray(A);
      if( A.ndim < 2 ) throw new Error('A must be at least 2D.');
      const
        DTypeArray = nd.dtypes[ nd.super_dtype(A.dtype,'float64') ], // <- ensure at least double precision
        Q_shape = A.shape,
        R_shape = Int32Array.from(Q_shape),
        [N,M] = Q_shape.slice(-2);
      R_shape[R_shape.length-2] = M;

      if( N <= M ) return la.qr_decomp_full(A);

      const Q = DTypeArray.from(A.data); // <- we could encourage GC by setting `A = undefined` after this line
      A = undefined
      const R = new DTypeArray(Q.length/N*M),
           cs = new DTypeArray( N*M - (M*(M+1) >>> 1) ),
            r = function(){
      	      try      { return    cs.subarray(M); }
      	      catch(e) { return new DTypeArray(M); }
            }();  // <- additional space to temp. store rows of R not contained in the result

      for(
        let R_off=0,
            Q_off=0; Q_off < Q.length; Q_off += N*M,
                                       R_off += M*M
      )
      {
        let csi=0;

        // COMPUTE R (inside of Q)
        for( let i=1; i < N; i++ ) { const I = Math.min(i,M);
        for( let j=0; j < I; j++ )
        { // USE GIVENS ROTATION TO ELIMINATE ELEMENT R_ji
          const R_ij = Q[Q_off + M*i+j]; if( R_ij === 0.0 ) { cs[csi++] = 0.0; continue; }
          const R_jj = Q[Q_off + M*j+j];
          let          norm = Math.hypot(R_jj,R_ij),
            c = R_jj / norm,
            s = R_ij / norm;
          if( c < 0 ) {
               c *= -1;
               s *= -1;
            norm *= -1;
          }
          cs[csi++] = s;
          Q[Q_off + M*i+j] = 0;
          Q[Q_off + M*j+j] = norm;
          // ROTATE ROW i AND j IN R
          for( let k=j; ++k < M; ) {
            const ik = Q_off + M*i+k, R_ik = Q[ik],
                  jk = Q_off + M*j+k, R_jk = Q[jk];
            Q[ik] = c*R_ik - s*R_jk;
            Q[jk] = s*R_ik + c*R_jk;
          }
        }}

        if( csi != cs.length ) throw new Error('Assertion failed!');

        // MOVE R FROM Q -> R AND INIT Q TO I
        for( let i=0; i < M; i++ )
        for( let j=i; j < M; j++ ) {
          R[R_off + M*i+j] = Q[Q_off + M*i+j];
                             Q[Q_off + M*i+j] = i !== j ? 0.0 : 1.0;
        }

        // COMPUTE Q
        for( let i=N; --i > 0; ) { const I = Math.min(i,M);
        for( let j=I; j-- > 0; )
        {
          const s = cs[--csi]; if( 0.0 === s ) continue;
          const c = Math.sqrt( (1-s)*(1+s) );
          // ROTATE ROW i AND j IN Q
          for( let k=j; k < M; k++ ) {
            const ik = Q_off + M*i+k, R_ik = Q[ik],
                  jk = Q_off + M*j+k, R_jk = Q[jk];
            Q[ik] = s*R_jk + c*R_ik;
            Q[jk] = c*R_jk - s*R_ik;
          }
        }}

        if( csi != 0 ) throw new Error('Assertion failed!');
      }

      return [
        new nd.Array(Q_shape, Q),
        new nd.Array(R_shape, R)
      ];
    }


    la.qr_solve = (Q,R,y) =>
    {
      if( undefined == y ) { y=R; [Q,R] = Q; }
      if( Q.ndim < 2 ) throw new Error('Q must be at least 2D.');
      if( R.ndim < 2 ) throw new Error('R must be at least 2D.');
      if( y.ndim < 2 ) throw new Error('y must be at least 2D.');

      const
        [N,M] = Q.shape.slice(-2),
        [I]   = R.shape.slice(-1),
        [J]   = y.shape.slice(-1);
      if( N != y.shape[y.ndim-2] ) throw new Error("Q and y don't match.")
      if( M != R.shape[R.ndim-2] ) throw new Error("Q and R don't match.")

      const
        ndim = Math.max(Q.ndim, R.ndim, y.ndim),
        shape = Int32Array.from({ length: ndim }, () => 1 );
      shape[ndim-2] = M;
      shape[ndim-1] = J;

      // FIND COMMON (BROADCASTED) SHAPE
      for( let arr of [Q,R,y] )
        for( let i=ndim-2, j=arr.ndim-2; i-- > 0 && j-- > 0; )
          if( 1 === shape[i] )
            shape[i] = arr.shape[j];
          else if( shape[i] != arr.shape[j] && arr.shape[j] != 1 )
            throw new Error('Q, R, y are not broadcast-compatible.');

      // GENERATE RESULT DATA
      const
        DTypeArray = nd.dtypes[ nd.super_dtype(Q.dtype,R.dtype,y.dtype) ],
        x_dat = new DTypeArray( shape.reduce((a,b) => a*b, 1) ),
        Q_dat = Q.data,
        R_dat = R.data,
        y_dat = y.data;
      let
        Q_off = 0, Q_stride = 1,
        R_off = 0, R_stride = 1,
        y_off = 0, y_stride = 1,
        x_off = 0;

      function solv(d) {
        if( d === ndim-2 ) {
          Q_stride = N*M;
          R_stride = M*I;
          y_stride = N*J;

          // Q.T @ y
          for( let i=0; i < M; i++ )
          for( let j=0; j < J; j++ )
          for( let k=0; k < N; k++ )
            x_dat[x_off+i*J+j] += Q_dat[Q_off+k*M+i] * y_dat[y_off+k*J+j]

          // BACKWARD SUBSTITUTION
          for( let i=M; i-- > 0; )
          for( let j=J; j-- > 0; ) {
            for( let k=I; --k > i; )
              x_dat[x_off+i*J+j] -= R_dat[R_off+I*i+k] * x_dat[x_off+k*J+j]
            x_dat[x_off+i*J+j] /= R_dat[R_off+I*i+i]
          }

          Q_off += Q_stride;
          R_off += R_stride;
          y_off += y_stride;
          x_off += M*J;

          return;
        }
        for( let l=shape[d]; ; l-- ) {
          solv(d+1);
          if( l == 1 ) break;
          if( ! (Q.shape[ d - ndim + Q.ndim ] > 1) ) Q_off -= Q_stride;
          if( ! (R.shape[ d - ndim + R.ndim ] > 1) ) R_off -= R_stride;
          if( ! (y.shape[ d - ndim + y.ndim ] > 1) ) y_off -= y_stride;
        }
        Q_stride *= Q.shape[ d - ndim + Q.ndim ] || 1;
        R_stride *= R.shape[ d - ndim + R.ndim ] || 1;
        y_stride *= y.shape[ d - ndim + y.ndim ] || 1;
      }
      solv(0);

      return new nd.Array(shape,x_dat);
    }


    la.rrqr_decomp = A => {
      // TODO: implement Strong RRQR as well, e.g. as la.rrqr_decomp_strong 
      // SEE: Ming Gu, Stanley C. Eisenstat,
      //     "EFFICIENT ALGORITHMS FOR COMPUTING A STRONG RANK-REVEALING QR FACTORIZATION"
      //      https://math.berkeley.edu/~mgu/MA273/Strong_RRQR.pdf

      A = nd.asarray(A);
      if( A.ndim < 2 ) throw new Error('A must be at least 2D.');
      const
        DTypeArray = nd.dtypes[ nd.super_dtype(A.dtype,'float64') ], // <- ensure at least double precision
        Q_shape =                 A.shape,
        R_shape = Int32Array.from(Q_shape),
        P_shape =                 Q_shape.slice(0,-1),
        [N,M]   =                 Q_shape.slice(  -2);
      R_shape[R_shape.length-2] = M;
      P_shape[P_shape.length-1] = M;

      if( N < M ) throw new Error('Not yet implemented for A.shape[-2] < A.shape[-1].');

      const Q = DTypeArray.from(A.data); // <- we could encourage GC by setting `A = undefined` after this line
      A = undefined
      const
        norm_sum = new DTypeArray(M), // <─┬─ underflow-safe representation of the column norm
        norm_max = new DTypeArray(M), // <─╯  (https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Math/hypot#Polyfill)
        P = new Int32Array(Q.length/N), // <- tracks column permutations
        R = new DTypeArray(Q.length/N*M),
        cs= new DTypeArray( 2*N*M - M*(M+1) ),// <- cache cos() and sin() values to apply M column rotations to Q at once
        q = new DTypeArray(N-M); // <- caches the part of Q that is not contained in the result
      const norm = i => {
        let max = norm_max[i];
        return isFinite(max) ? Math.sqrt(norm_sum[i])*max : max;
      }

      for(
        let Q_off=0,
            R_off=0,
            P_off=0; Q_off < Q.length; Q_off += N*M,
                                       R_off += M*M,
                                       P_off +=   M
      )
      {
        // INIT P
        for( let i=0; i < M; i++ ) P[P_off + i] = i;

        // COMPUTE COLUMN NORM
        // (https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Math/hypot#Polyfill)
        for( let j=0; j < N; j++ )
        for( let k=0; k < M; k++ ) {
          const A_jk = Math.abs(Q[Q_off + M*j+k]);
          if(   A_jk > 0 ) {
            if( A_jk > norm_max[k] ) {
              const scale = norm_max[k] / A_jk; norm_max[k] = A_jk;
              norm_sum[k] *= scale*scale;
            }
            const ratio = A_jk / norm_max[k];
            norm_sum[k] += ratio*ratio;
          }
        }

        let csi = 0;
        // ELIMINATE COLUMN BY COLUMN OF R (WHICH IS CURRENTLY STORED IN Q)
        for( let i=0; i < M; i++ )
        { // FIND PIVOT COLUMN THAT (HOPEFULLY) ENSURES RANK REVEAL (RRQR is inherently not guaranteed to do that)
          let
            p = i,
            norm_p = norm(p);
          for( let j=i; ++j < M; )
            if( norm(j) > norm(p) ) { p = j; norm_p = norm(p); }
          // swap pivot to column i
          if( p != i ) {
            for( let j=0; j < N; j++ ) {
              const ji = Q_off + M*j+i,
                    jp = Q_off + M*j+p, A_ji = Q[ji]; Q[ji] = Q[jp]; Q[jp] = A_ji;
            }
            const P_i = P[P_off+i]; P[P_off+i] = P[P_off+p]; P[P_off+p] = P_i;
          }

          // RESET COLUMN NORM (INDEX i IS SET TO ZERO FOR THE NEXT RRQR)
          norm_sum.fill(0.0, i);
          norm_max.fill(0.0, i);

          // REMOVE ENTRIES BELOW DIAGONAL
          for( let j=i+1; j < N; j++ )
          { // compute Given's rotation 
            const A_ji = Q[Q_off + M*j+i]; if( A_ji == 0.0 ) { cs[csi++]=1.0;
                                                               cs[csi++]=0.0; continue };
            const A_ii = Q[Q_off + M*i+i],
                         norm = Math.hypot(A_ii,A_ji),
              c = A_ii / norm,
              s = A_ji / norm;
            Q[Q_off + M*i+i] = norm;
            Q[Q_off + M*j+i] = 0.0;
            // rotate i and j
            for( let k=i; ++k < M; ) {
              // Given's rotation
              const jk = Q_off + M*j+k; {
              const ik = Q_off + M*i+k, 
                  A_ik = Q[ik],
                  A_jk = Q[jk];
                Q[jk] = c*A_jk - s*A_ik;
                Q[ik] = s*A_jk + c*A_ik;
              }
              // re-compute column norm
              const A_jk = Math.abs(Q[jk]);
              if(   A_jk > 0 ) {
                if( A_jk > norm_max[k]  ) {
                  const scale = norm_max[k] / A_jk; norm_max[k] = A_jk;
                  norm_sum[k] *= scale*scale;
                }
                const ratio = A_jk / norm_max[k];
                norm_sum[k] += ratio*ratio;
              }
            }
            cs[csi++] = c;
            cs[csi++] = s;
          }
        }

        if( csi != cs.length ) throw new Error('Assertion failed!')

        // MOVE R FROM Q -> R
        for( let i=0; i < M; i++ )
        for( let j=i; j < M; j++ ) {
          R[R_off + M*i+j] = Q[Q_off + M*i+j];
                             Q[Q_off + M*i+j] = i !== j ? 0.0 : 1.0;
        }

        // COMPUTE Q
        for( let i=M; i-- > 0; )
        for( let j=N; --j > i; )
        {
          const s = cs[--csi],
                c = cs[--csi];
          for( let k=M; k-- > i;) {
            // Given's rotation
            const jk = Q_off + M*j+k; {
            const ik = Q_off + M*i+k, 
                A_ik = Q[ik],
                A_jk = Q[jk];
              Q[jk] = s*A_ik + c*A_jk;
              Q[ik] = c*A_ik - s*A_jk;
            }
          }
        }

        if( csi != 0 ) throw new Error('Assertion failed!')
      }

      return [
        new nd.Array(Q_shape, Q),
        new nd.Array(R_shape, R),
        new nd.Array(P_shape, P)
      ];
    }


    la.hessenberg_decomp = A => {
      A = nd.asarray(A);
      if( A.ndim < 2 ) throw new Error('A must be at least 2D.');
      const
        DTypeArray = nd.dtypes[ nd.super_dtype(A.dtype,'float64') ], // <- ensure at least double precision
        shape = A.shape,
        N = shape[shape.length-1],
        w   = new DTypeArray(N-1),
        tmp = new DTypeArray(N);
      if( N != shape[shape.length-2] ) throw new Error('A is not square.');

      const H = DTypeArray.from(A.data); A = undefined;
      const U = new DTypeArray(H.length);

      for( let off=0; off < H.length; off += N*N )
      { // INIT P TO IDENTITY
        for( let i=0; i < N; i++ ) U[off + N*i+i] = 1.0;

        // ELIMINATION BEGIN
        for( let i=N; --i > 1; )
        { 
          let row_norm; { // BEGIN CREATE HOUSEHOLDER VECTOR
            let
              norm_sum = 0.0,
              norm_max = 0.0;
            for( let j=i; j-- > 0; ) {
              let H_ij = H[off + N*i+j];
              w[j] = H_ij;
              H_ij = Math.abs(H_ij);
              if(   H_ij > 0 ) {
                if( H_ij > norm_max ) {
                  const scale = norm_max / H_ij; norm_max = H_ij;
                  norm_sum *= scale*scale;
                }
                const ratio = H_ij / norm_max;
                norm_sum += ratio*ratio;
              }
            }
            row_norm = isFinite(norm_max) ? Math.sqrt(norm_sum)*norm_max : norm_max;
            w[i-1] -= row_norm;
            norm_sum = 0.0;
            norm_max = 0.0;
            for( let j=i; j-- > 0; ) {
              const w_j = Math.abs(w[j]);
              if(   w_j > 0 ) {
                if( w_j > norm_max ) {
                  const scale = norm_max / w_j; norm_max = w_j;
                  norm_sum *= scale*scale;
                }
                const ratio = w_j / norm_max;
                norm_sum += ratio*ratio;
              }
            }
            norm_sum = isFinite(norm_max) ? Math.sqrt(norm_sum)*norm_max : norm_max;
            if( norm_sum == 0 ) continue;
            for( let j=i; j-- > 0; ) w[j] /= norm_sum;
          } // END CREATE HOUSEHOLDER VECTOR

          // APPLY w TO RIGHT OF H
          for( let j=N; j-- > 0; )
            if( j != i )
            { let sum=0.0
              for( let k=i; k-- > 0; ) sum += w[k] * H[off + N*j+k];
              sum *= 2;
              for( let k=i; k-- > 0; ) H[off + N*j+k] -= w[k] * sum;
            }
          H.fill( 0.0, off + N*i+0, off + N*i+i-1 )
          H[off + N*i+i-1] = row_norm;

          // APPLY w TO LEFT OF H
          tmp.fill(0.0);
          for( let j=i; j-- > 0; )
          for( let k=N; k-- > 0; ) tmp[k] += w[j] * H[off + N*j+k];
          for( let j=i; j-- > 0; )
          for( let k=N; k-- > 0; ) H[off + N*j+k] -= 2 * w[j] * tmp[k];

          // APPLY w TO RIGHT OF U
          for( let j=N; j-- > 0; )
          { let sum=0.0
            for( let k=i; k-- > 0; ) sum += w[k] * U[off + N*j+k];
            sum *= 2;
            for( let k=i; k-- > 0; ) U[off + N*j+k] -= w[k] * sum;
          }
        }
        // ELIMINATION END
      }

      return [
        new nd.Array(shape,U),
        new nd.Array(shape,H)
      ];
    };


//    la.hessenberg_decomp = A => {
//      A = nd.asarray(A);
//      if( A.ndim < 2 ) throw new Error('A must be at least 2D.');
//      const
//        DTypeArray = nd.dtypes[ nd.super_dtype(A.dtype,'float64') ], // <- ensure at least double precision
//        shape = A.shape,
//        N = shape[shape.length-1];
//      if( N != shape[shape.length-2] ) throw new Error('A is not square.');
//
//      const H = DTypeArray.from(A.data); A = undefined;
//      const U = new DTypeArray(H.length);
//
//      for( let off=0; off < H.length; off += N*N )
//      { // INIT P TO IDENTITY
//        for( let i=0; i < N; i++ ) U[off + N*i+i] = 1.0;
//
//        // { ELIMINATION BEGIN
//        for( let i=N-1; --i > 0; )
//        for( let j=i  ; j-- > 0; )
//        { // ROTATE COLUMS/ROWS j and i
//          // determine angle
//          const H_Ij = H[off + N*(i+1)+j]; if( H_Ij == 0 ) continue;
//          const H_Ii = H[off + N*(i+1)+i],
//                norm = Math.hypot(H_Ii,H_Ij),
//                c = H_Ii / norm,
//                s = H_Ij / norm;
//
//          // rotate columns in H
//          for( let k=i+2; k-- > 0; )
//          {
//            const kj = off + N*k+j, H_kj = H[kj],
//                  ki = off + N*k+i, H_ki = H[ki];
//            H[kj] = c*H_kj - s*H_ki;
//            H[ki] = s*H_kj + c*H_ki;
//          }
//
//          // rotate rows in H
//          for( let k=N; k-- > 0; )
//          {
//            const jk = off + N*j+k, H_jk = H[jk],
//                  ik = off + N*i+k, H_ik = H[ik];
//            H[jk] = c*H_jk - s*H_ik;
//            H[ik] = s*H_jk + c*H_ik;
//          }
//
//          // rotate rows in U (transposed for cache locality reasons)
//          for( let k=N; k-- > 0; ) // <- FIXME are there operations to be saved here?
//          {
//            const jk = off + N*j+k, U_jk = U[jk],
//                  ik = off + N*i+k, U_ik = U[ik];
//            U[jk] = c*U_jk - s*U_ik;
//            U[ik] = s*U_jk + c*U_ik;
//          }
//        }
//        // } ELIMINATION END
//
//        // TRANSPOSE U BACK (transposed for cache locality reasons)
//        for( let i=0; i < N; i++ )
//        for( let j=0; j < i; j++ ) {
//          const ij = off + N*i+j,
//                ji = off + N*j+i, U_ij = U[ij]; U[ij] = U[ji]; U[ji] = U_ij;
//        }
//      }
//
//      return [
//        new nd.Array(shape,U),
//        new nd.Array(shape,H)
//      ]
//    };


    la.bidiag_decomp = A => {
      A = nd.asarray(A);
      if( A.ndim < 2 ) throw new Error('A must be at least 2D.');
      const
        DTypeArray = nd.dtypes[ nd.super_dtype(A.dtype,'float64') ], // <- ensure at least double precision
        U_shape = Int32Array.from(A.shape),
        B_shape = Int32Array.from(U_shape),
        V_shape = Int32Array.from(U_shape),
        [N,M] = A.shape.slice(  -2),
        len   = A.shape.slice(0,-2).reduce( (a,b) => a*b, 1 ),
        I     = Math.min(N,M),   // #rows    of bidiag. matrix B
        J     = N >= M ? I : I+1;// #columns of bidiag. matrix B
      U_shape[U_shape.length-1] = I;
      B_shape[B_shape.length-2] = I;
      B_shape[B_shape.length-1] = J;
      V_shape[V_shape.length-2] = J;

      let U,V;
      // copy A to U or V, whichever is large enough
      if( N >= M ) {
        U  =     DTypeArray.from(A.data); A = U;
        V  = new DTypeArray( len*M*M );
      } else {
        V  = new DTypeArray(len*(N+1)*M);
        A = A.data;
        for( let i=A.length,
                 j=V.length; i-- > 0; ) V[--j] = A[i];
        A = V;
        U  = new DTypeArray(len*N*N);
      }

      const
        B  = new DTypeArray( len*I*J ),                                //<- resulting bidiagonal matrix
        csU= new DTypeArray(         I * (N-1 + N-I  ) ),              //<- caches ALL row    rotations later used to compute U
        csV= new DTypeArray( N < M ? I * (M-2 + M-I-1) : (I-1)*(M-2) ),//<- caches ALL column rotations later used to compute V
        uv = new DTypeArray( Math.max(N-I,M) );                        //<- cached remainder of a row of U or a column of V to compute U/V

      for(
        let A_off= N >= M ? 0 : len*M,
            U_off=0,
            B_off=0,
            V_off=0; B_off < B.length;
            A_off += N*M,
            U_off += N*I,
            B_off += I*J,
            V_off += J*M
      )
      {
        let
          csiU = 0,
          csiV = 0;
        for( let i=0; i < I; i++ )
        {
          // ELIMINATE THE ELEMENTS IN (i+1)-th COLUMN BELOW DIAG (IN B)
          for( let j=i+1; j < N; j++ )
          {
            const     A_ji = A[A_off + M*j+i]; if( 0.0 == A_ji ) { csU[csiU++]=1; csU[csiU++]=0; continue; }
            const     A_ii = A[A_off + M*i+i],
                             norm = Math.hypot(A_ii,A_ji),
                  c = A_ii / norm,
                  s = A_ji / norm;
            A[A_off + M*i+i] = norm;
            A[A_off + M*j+i] = 0;
            // ROTATE ROWS IN A
            for( let k=i; ++k < M; ) {
              const
                ik = A_off + M*i+k, A_ik = A[ik],
                jk = A_off + M*j+k, A_jk = A[jk];
              A[ik] = A_jk*s + A_ik*c;
              A[jk] = A_jk*c - A_ik*s;
            }
            csU[csiU++] = c;
            csU[csiU++] = s;
          }
          if( i+3 > M ) continue;
          // ELIMINATE (i+1)-th ROW RIGHT OF BIDIAG
          for( let j=i+2; j < M; j++ ) {
            const ii = A_off + M*i+(i+1), A_ii = A[ii],
                  ij = A_off + M*i+ j   , A_ij = A[ij], norm = Math.hypot(A_ii,A_ij);
            A[ii] = norm;
            A[ij] = 0;
            csV[csiV++] = norm == 0 ? 1 : A_ii / norm;
            csV[csiV++] = norm == 0 ? 0 : A_ij / norm;
          }
          // APPLY ABOVE ELIMINATION TO REMAINING ROWS RIGHT OF BIDIAG (ROW-BY-ROW IMPROVES CACHE LOCALITY)
          for( let k=i+1; k < N; k++ ) { csiV -= 2*(M-i-2);
          for( let j=i+2; j < M; j++ ) {
              const
                c = csV[csiV++],
                s = csV[csiV++],
                ki = A_off + M*k+(i+1), A_ki = A[ki],
                kj = A_off + M*k+ j   , A_kj = A[kj];
              A[ki] = A_kj*s + A_ki*c;
              A[kj] = A_kj*c - A_ki*s;
          }}
        }

        if( csiU != csU.length ) throw new Error('Assertion failed: (csiU='+csiU+") != (csU.length="+csU.length+")" );
        if( csiV != csV.length ) throw new Error('Assertion failed: (csiV='+csiV+") != (csV.length="+csV.length+")" );

        // MOVE A -> B
        for( let i=0; i < I; i++ )
        for( let j=0; j < J; j++ )
          B[B_off + J*i+j] = A[A_off + M*i+j];

         //
        // COMPUTE U
       //
        // ROTATE U
        for( let k = 0; k < N; k++ )
        { csiU = 0
          for( let i=0; i < N-I; i++ )          uv[i] = k != (i+I) ? 0.0 : 1.0;
          for( let i=0; i < I  ; i++ ) U[U_off+I*k+i] = k !=  i    ? 0.0 : 1.0;
          for( let i=0; i < I  ; i++ ) {
            const skip = Math.max(1,k-i);
            csiU += 2*(skip-1)
            for( let j = i+skip; j < I; j++ ) {
              const
                c = csU[csiU++],
                s = csU[csiU++],
                ki = U_off + I*k+i, U_ki = U[ki],
                kj = U_off + I*k+j, U_kj = U[kj];
              U[ki] = U_kj*s + U_ki*c;
              U[kj] = U_kj*c - U_ki*s;              
            }
            for( let j = Math.max(0,i+skip-I); j < (N-I); j++ ) {
              const
                c = csU[csiU++],
                s = csU[csiU++], ki = U_off + I *k+i, U_ki = U[ki];
              U[ki] = uv[j]*s + U_ki*c;
              uv[j] = uv[j]*c - U_ki*s;
            }
          }
        }
         //
        // COMPUTE V
       //
        // ROTATE V
        for( let k = 0; k < M; k++ )
        { csiV = 0
          // uv CACHES (k+1)-th COLUMN OF V
          uv.fill(0.0, 0,M)
          uv[k] = 1
          for( let i=0; i < I; i++ )
          {
            const skip = Math.max(2,k-i);
            csiV += 2*(skip-2);
            for( let j=i+skip; j < M; j++ ) {
              const
                c = csV[csiV++],
                s = csV[csiV++],
                v_i = uv[i+1],
                v_j = uv[j  ];
              uv[i+1] = v_j*s + v_i*c;
              uv[j  ] = v_j*c - v_i*s;
            }
          }
          // WRITE (k+1)-th COLUMN
          for( let i=0; i < J; i++ ) V[V_off + M*i+k] = uv[i];
        }
      }

      return [
        new nd.Array(U_shape, U),
        new nd.Array(B_shape, B),
        new nd.Array(V_shape, V)
      ];
    }


    la._svd_decomp_jac_1sided = A => {
      // SEE:
      //   http://www.netlib.org/lapack/lawnspdf/lawn169.pdf
      //   http://www.netlib.org/lapack/lawnspdf/lawn170.pdf
      A = nd.asarray(A);
      const
         U_shape =                 A.shape,
         V_shape = Int32Array.from(U_shape),
         sv_shape=                 U_shape.slice(0,-1),
        [N,M]    =                 U_shape.slice(  -2);

      if( N < M ) {
        const [U,sv,V] = la._svd_decomp_jac_1sided(A.T);
        // TODO: transpose V in-place
        return [V.T,sv,U.T];
      }

      const
        DTypeArray = nd.dtypes[nd.super_dtype(A.dtype,'float64')],
        TOL = Number.EPSILON * N,// <- FIXME what about float32 ?s
        sqrt2 = Math.sqrt(2);
      sv_shape[sv_shape.length-1] = M;
       V_shape[ V_shape.length-2] = M;

      // TODO: check if QR preconditioning is sensible for N==M
      let [Q,R,P] = la.rrqr_decomp( N >= M ? A : A.T );
      A = undefined;

      Q = Q.data;
      R = R.data;
      P = P.data;

      const
        U  = new DTypeArray(M*M), // <- temporarily stores U
        sv = new DTypeArray(Q.length/N),
        tmp= new DTypeArray(M); // <- temporary storage for matrix multiplication

      function* indices_241(M) {
        // SEE: W. F. Mascarenhas "On the Convergence of the Jacobi Method for Arbitrary Orderings."
        // ┏    ╷    ╷         ┓
        // ┃ B1 ┊ B3 ┊         ┃
        // ┃┄┄┄┄┼┄┄┄┄┤   B7    ┃
        // ┃    ┊ B2 ┊         ┃
        // ┃    └┄┄┄┄┼┄┄┄┄┬┄┄┄┄┃ = RRᵀ
        // ┃         ┊ B4 ┊ B6 ┃
        // ┃         └┄┄┄┄┼┄┄┄┄┃
        // ┃              ┊ B5 ┃
        // ┗              ╵    ┛
        // ORDER: {1,2,3,4,5,6,1,2,4,5,7}
        const
          m = M >>> 2,
          IJ = new Int32Array(2);
        for( let k=0; k < 2; k++ )
        {
          // B1
          for( IJ[0]=0      ; IJ[0] <   m; IJ[0]++ )
          for( IJ[1]=1+IJ[0]; IJ[1] <   m; IJ[1]++ ) yield IJ;
          // B2
          for( IJ[0]=  m    ; IJ[0] < 2*m; IJ[0]++ )
          for( IJ[1]=1+IJ[0]; IJ[1] < 2*m; IJ[1]++ ) yield IJ;
          // B3
          if( 0==k )
          for( IJ[0]=0      ; IJ[0] <   m; IJ[0]++ )
          for( IJ[1]=  m    ; IJ[1] < 2*m; IJ[1]++ ) yield IJ;
          // B4
          for( IJ[0]=2*m    ; IJ[0] < 3*m; IJ[0]++ )
          for( IJ[1]=1+IJ[0]; IJ[1] < 3*m; IJ[1]++ ) yield IJ;
          // B5
          for( IJ[0]=3*m    ; IJ[0] <   M; IJ[0]++ )
          for( IJ[1]=1+IJ[0]; IJ[1] <   M; IJ[1]++ ) yield IJ;
          // B6
          if( 0==k )
          for( IJ[0]=2*m    ; IJ[0] < 3*m; IJ[0]++ )
          for( IJ[1]=3*m    ; IJ[1] <   M; IJ[1]++ ) yield IJ;
        }
        // B7
        for( IJ[0]=  0; IJ[0] < 2*m; IJ[0]++ )
        for( IJ[1]=2*m; IJ[1] <   M; IJ[1]++ ) yield IJ;
      }

      for( let Q_off    = 0,
               R_off    = 0,
               P_sv_off = 0; R_off < R.length;
               Q_off    += N*M,
               R_off    += M*M,
               P_sv_off += M
      )
      {
        // [STEP 1] determine numerical rank of R
        const RANK = function(){
          let rankTol = Number.EPSILON * function(){
            let
              sum = 0.0,
              max = 0.0;
            for( let i=M*M; i-- > 0; )
            {
              const elem = Math.abs( R[R_off + i] );
              if(   elem != 0 ) { // <- handles NaN (by making the result NaN)
                if( elem > max ) {
                  const scale = max / elem; max = elem;
                  sum *= scale*scale;
                }
                const ratio = elem / max;
                sum += ratio*ratio;
              }
            }
            return isFinite(max) ? Math.sqrt(sum)*max : max;
          }();
          let RANK = M;
          while( Math.abs( R[R_off + (M+1)*(--RANK)] ) <= rankTol );
          return RANK+1;
        }();

        // [STEP 2] eliminate the right-most, rank-deficient columns of R, then only perform the Jacobi SVD on the upper left square part
        for( let i=M   ; i-- > RANK; ) { // <- columns that are to be eliminated
        for( let j=RANK; j-- > 0   ; ) {
        for( let k=RANK; --k > j   ; ) {
            const
              s    = U[        M*i+k],
              c    = R[R_off + M*i+k],
              R_ji = R[R_off + M*j+i],
              R_jk = R[R_off + M*j+k];
            R[R_off + M*j+i] = c*R_ji - s*R_jk;
            R[R_off + M*j+k] = s*R_ji + c*R_jk;
        }
          const
            R_ji = R[R_off + M*j+i],
            R_jj = R[R_off + M*j+j],
            norm = Math.hypot(R_ji,R_jj),
            c = R_jj / norm,
            s = R_ji / norm;
          R[R_off + M*j+i] = c*R_ji - s*R_jj;
          R[R_off + M*j+j] = s*R_ji + c*R_jj;
          R[R_off + M*i+j] = c;
          U[        M*i+j] = s;
          if( ! (Math.abs(R[R_off + M*j+i]) < 1e-8) ) throw new Error('Assertion failed.');
          R[R_off + M*j+i] = 0;
        }}

        // [STEP 3]COPY R -> U
        for( let i=0   ; i < RANK; i++ )
        for( let j=0   ; j < RANK; j++ ) U[M*i+j] = R[R_off + M*i+j];

        // [STEP 4] PERFORM JACOBI SWEEPS
        let sweeps = 0;
        for( let maxErr=Infinity; maxErr > TOL*TOL; )
        {
          if( ++sweeps > 128 ) throw new Error('Too many iterations.')
          maxErr=0.0;
          for( const [i,j] of indices_241(RANK) ) // "2.41"-order with locally super-quadratic convergence (in theory)
//          for( let i=0  ; i < RANK; i++ ) // ◀─┬─ alternative, simpler (De Rijk) order
//          for( let j=1+i; j < RANK; j++ ) // ◀─╯
          {
            if( j <= i ) throw new Error('Assertion failed');
            let R_ii = 0.0, // ◀─╮
                R_ij = 0.0, // ◀─┼─ TODO: make this underflow-saf(ish)er
                R_jj = 0.0; // ◀─╯ 
            { // compute dot products
              const    IK = R_off + M*i + RANK;
              for( let ik = R_off + M*i,
                       jk = R_off + M*j; ik < IK; )
              {
                const R_ik = R[ik++]; R_ii += R_ik*R_ik;
                const R_jk = R[jk++]; R_ij += R_ik*R_jk;
                                      R_jj += R_jk*R_jk;
              }
            }
            if( R_ii <= 0 ) throw new Error('Unexpected underflow.');
            if( R_jj <= 0 ) throw new Error('Unexpected underflow.');
            maxErr = Math.max( maxErr, (R_ij / R_ii) * (R_ij / R_jj) );// <- FIXME: check whether underflow makes a problem
            if( maxErr/TOL <= TOL ) continue;

            // determine rotation angle using binary search. TODO: use sth. faster than binary search
            let
              s = nd.opt.root1d(
                s => {
                  const cs = Math.sqrt( (1+s)*(1-s) ) * s;
                  return cs*R_ii + R_ij*(1 + sqrt2*s)*(1 - sqrt2*s) - cs*R_jj; // <- theoretically better approximation (underflow-safe?)
                }, /*s_min=*/-sqrt2/2, /*s_max=*/+sqrt2/2
              ),
              c = Math.sqrt( (1+s)*(1-s) );
            

            // CHOOSE (ANGLE) SOLUTION THAT KEEPS MAIN DIAGONAL IN ORDER, i.e.:
            // c*(c*R_ii - s*R_ij) - s*(c*R_ij - s*R_jj) = r_ii < r_jj = s*(s*R_ii + c*R_ij) + c*(s*R_ij + c*R_jj)
            if( R_ii*(1 + sqrt2*s)*(1 - sqrt2*s) < R_jj + 2*c*s*R_ij ) { const C = c; c = -s; s = C; }

            { // perform rotation
              const    ik0= R_off + M*i;
              for( let ik = R_off + M*i+RANK,
                       jk = R_off + M*j+RANK; ik > ik0; )
              { const R_ik = R[--ik],
                      R_jk = R[--jk];
                R[ik] = c*R_ik - s*R_jk;
                R[jk] = c*R_jk + s*R_ik;
              }
            }
          }
        }

        // [STEP 5] COMPUTE SINGULAR VALUES AND NORMALIZE ROWS OF V
        // TODO: V might also have to be reorthogonalized
        for( let i=0; i < RANK; i++ )
        {
          let sum = 0.0,
              max = 0.0;
          for( let j=0; j < RANK; j++ )
          {
            const R_ij = Math.abs(R[R_off + M*i+j]);
            if(   R_ij != 0 ) { // <- handles NaN (by making the result NaN)
              if( R_ij > max ) {
                const scale = max / R_ij; max = R_ij;
                sum *= scale*scale;
              }
              const ratio = R_ij / max;
              sum += ratio*ratio;
            }
          }
          const norm = isFinite(max) ? Math.sqrt(sum)*max : max;
          if( norm == 0 ) throw new Error('Unhandled underflow. Ask the developer nicely and bribe him with a couple of beers, then he may address this issue.');
          for( let j=RANK; j-- > 0; ) R[R_off + M*i+j] /= norm;
          sv[P_sv_off + i] = norm;
        }

        // [STEP 6] COMPUTE U ( Using U = R @ inv(V) @ inv(SV) = R @ V.T @ diag(1/sv) )
        //   ( should be numerically stable since V is orthogonal, i.e. well-conditioned ... right? )
        for( let i=0; i < RANK; i++ ) {
          for( let j=0; j < RANK; j++ ) { tmp[j] = U[M*i+j]; U[M*i+j]=0; }
  
          for( let j=0; j < RANK; j++ ) {
            for( let k=0; k < RANK; k++ )
              U[M*i+j] += R[R_off + M*j+k]*tmp[k];
            U[M*i+j] /= sv[P_sv_off + j];
          }
        }

        // [STEP 7] UNDO [STEP 2] 
        for( let i=RANK; i < M   ; i++ ) {
        for( let j=0   ; j < i   ; j++ ) {
        for( let k=0   ; k < RANK; k++ ) {
            const
              s    = U[        M*i+k],
              c    = R[R_off + M*i+k],
              R_ji = R[R_off + M*j+i],
              R_jk = R[R_off + M*j+k];
            R[R_off + M*j+i] =  c*R_ji + s*R_jk;
            R[R_off + M*j+k] = -s*R_ji + c*R_jk;
        }}
          // LAST REMAINING ROW
          let R_ii = 1.0;
          for( let k=0; k < RANK; k++ ) {
            const
              s = U[        M*i+k],
              c = R[R_off + M*i+k];
                  R[R_off + M*i+k] = -s*R_ii; R_ii *= c;
          }
          R[R_off + M*i+i] = R_ii;
        }

        // [STEP 8] MULTIPLY Q = Q @ U
        for( let i=0; i < N; i++ )
        {
          for( let j=0; j < RANK; j++ ) { tmp[j] = Q[Q_off + M*i+j]; Q[Q_off + M*i+j]=0; }

          for( let k=0; k < RANK; k++ )
          for( let j=0; j < RANK; j++ ) Q[Q_off + M*i+j] += tmp[k] * U[M*k+j];
        }

        // [STEP 9] APPLY P TO R (COLUMN PERMUTATIONS)
        for( let i=0; i < M; i++ )
        {
          for( let j=0; j   < M; j++ ) tmp[P[P_sv_off + j]] = R[R_off + M*i+j];
          for( let j=M; j-- > 0;     )                        R[R_off + M*i+j] = tmp[j];
        }

        // TODO SORT SV???
      }

      return [
        new nd.Array( U_shape, Q ),
        new nd.Array(sv_shape, sv),
        new nd.Array( V_shape, R )
      ];
    };


//    la._svd_decomp_jac_1sided = A => {
//      // SEE:
//      //   http://www.netlib.org/lapack/lawnspdf/lawn169.pdf
//      //   http://www.netlib.org/lapack/lawnspdf/lawn170.pdf
//      A = nd.asarray(A);
//      const
//         U_shape =                 A.shape,
//         V_shape = Int32Array.from(U_shape),
//         sv_shape=                 U_shape.slice(0,-1),
//        [N,M]    =                 U_shape.slice(  -2);
//
//      if( N < M ) {
//        const [U,sv,V] = la._svd_decomp_jac_1sided(A.T);
//        // TODO: transpose V in-place
//        return [V.T,sv,U.T];
//      }
//
//      const
//        DTypeArray = nd.dtypes[nd.super_dtype(A.dtype,'float64')],
//        TOL = Number.EPSILON * N; // <- FIXME what about float32 ?
//      sv_shape[sv_shape.length-1] = M;
//       V_shape[ V_shape.length-2] = M;
//
//      // TODO: check if QR preconditioning is sensible for N==M
//      let [Q,R,P] = la.rrqr_decomp( N >= M ? A : A.T ); // <- TODO: determine rank and perform second QR decomp.
//      A = undefined;
//
//      Q = Q.data;
//      R = R.data;
//      P = P.data;
//
//      const
//        U  = new DTypeArray(M*M), // <- temporarily stores U
//        sv = new DTypeArray(Q.length/N),
//        tmp= new DTypeArray(M); // <- temporary storage for matrix multiplication
//
//      function* indices_241() {
//        // SEE: W. F. Mascarenhas "On the Convergence of the Jacobi Method for Arbitrary Orderings."
//        // ┏    ╷    ╷         ┓
//        // ┃ B1 ┊ B3 ┊         ┃
//        // ┃┄┄┄┄┼┄┄┄┄┤   B7    ┃
//        // ┃    ┊ B2 ┊         ┃
//        // ┃    └┄┄┄┄┼┄┄┄┄┬┄┄┄┄┃ = RRᵀ
//        // ┃         ┊ B4 ┊ B6 ┃
//        // ┃         └┄┄┄┄┼┄┄┄┄┃
//        // ┃              ┊ B5 ┃
//        // ┗              ╵    ┛
//        // ORDER: {1,2,3,4,5,6,1,2,4,5,7}
//        const
//          m = M >>> 2,
//          IJ = new Int32Array(2);
//        for( let k=0; k < 2; k++ )
//        {
//          // B1
//          for( IJ[0]=0      ; IJ[0] <   m; IJ[0]++ )
//          for( IJ[1]=1+IJ[0]; IJ[1] <   m; IJ[1]++ ) yield IJ;
//          // B2
//          for( IJ[0]=  m    ; IJ[0] < 2*m; IJ[0]++ )
//          for( IJ[1]=1+IJ[0]; IJ[1] < 2*m; IJ[1]++ ) yield IJ;
//          // B3
//          if( 0==k )
//          for( IJ[0]=0      ; IJ[0] <   m; IJ[0]++ )
//          for( IJ[1]=  m    ; IJ[1] < 2*m; IJ[1]++ ) yield IJ;
//          // B5
//          for( IJ[0]=3*m    ; IJ[0] <   M; IJ[0]++ )
//          for( IJ[1]=1+IJ[0]; IJ[1] <   M; IJ[1]++ ) yield IJ;
//          // B4
//          for( IJ[0]=2*m    ; IJ[0] < 3*m; IJ[0]++ )
//          for( IJ[1]=1+IJ[0]; IJ[1] < 3*m; IJ[1]++ ) yield IJ;
//          // B6
//          if( 0==k )
//          for( IJ[0]=2*m    ; IJ[0] < 3*m; IJ[0]++ )
//          for( IJ[1]=3*m    ; IJ[1] <   M; IJ[1]++ ) yield IJ;
//        }
//        // B7
//        for( IJ[0]=  0; IJ[0] < 2*m; IJ[0]++ )
//        for( IJ[1]=2*m; IJ[1] <   M; IJ[1]++ ) yield IJ;
//      }
//
//      for( let Q_off    = 0,
//               R_off    = 0,
//               P_sv_off = 0; R_off < R.length;
//               Q_off    += N*M,
//               R_off    += M*M,
//               P_sv_off += M
//      )
//      {
//        // determine the rank
//        let RANK = M;
//        while( R[R_off + (M+1)*(--RANK)] == 0 );
//        ++RANK;
//        console.log(`RANK: ${RANK}/${M}`);
//
//        // TODO eliminate the right-most, rank-deficient columns of R, then only perform the Jacobi SVD on the upper left square part
//
//        // COPY R -> U
//        for( let i=M*M; i-- > 0; ) U[i] = R[R_off + i];
//
//        // PERFORM JACOBI SWEEPS
//        let sweeps = 0;
//        for( let maxErr=Infinity; maxErr/TOL > TOL; )
//        {
//          if( ++sweeps > 128 ) throw new Error('Too many iterations.')
//          maxErr=0.0;
//          for( const [i,j] of indices_241() ) // "2.41"-order with locally super-quadratic convergence (in theory)
////          for( let i=0  ; i < M; i++ ) // ◀─┬─ alternative, simpler (De Rijk) order
////          for( let j=1+i; j < M; j++ ) // ◀─╯
//          {
//            if( j <= i ) throw new Error('Assertion failed');
//            let R_ii = 0.0,
//                R_ij = 0.0,
//                R_jj = 0.0;
//            { // compute dot products
//              const    IK = R_off + M*(i+1);
//              for( let ik = R_off + M* i,
//                       jk = R_off + M* j; ik < IK; )
//              {
//                const R_ik = R[ik++]; R_ii += R_ik*R_ik;
//                const R_jk = R[jk++]; R_ij += R_ik*R_jk;
//                                      R_jj += R_jk*R_jk;
//              }
//            }
//            if( R_ii <= 0 ) throw new Error('Assertion failed!');
//            if( R_jj <= 0 ) throw new Error('Assertion failed!');
//            maxErr = Math.max( maxErr, (R_ij / R_ii) * (R_ij / R_jj) );// <- FIXME: check whether underflow makes a problem
//
//            // determine rotation angle using binary search. TODO: use faster approximation than binary search
//            let c,s; {
//              let s_min = -Math.sqrt(2)/2,
//                  s_max = +Math.sqrt(2)/2;
//              const
//                sign = R_ii > R_jj ? +1.0 : -1.0, // <- re-orient the function to be montonously increasing
//                F = s => {
//                  const c = Math.sqrt( (1+s)*(1-s) );
////                  return sign * ( s*c*(R_ii - R_jj) + 2*R_ij*(0.5 + s)*(0.5 - s) ); // <- theoretically better approximation (underflow-safe)
//                  return sign * ( // <- approximates actually performed numerical operations
//                      s*(c*R_ii - s*R_ij)
//                    + c*(c*R_ij - s*R_jj)
//                  )
//                };
//              // begin binary search
//              for(;;) {
//                const s = (s_min + s_max)/2;
//                if( s_min == s ||
//                    s_max == s ) break;
//                const f = F(s);
//                if( isNaN(f) ) throw new Error('NaN encountered.');
//                if( f <= 0 ) s_min = s;
//                if( f >= 0 ) s_max = s;
//              }
//              // choose between α_min and α_max
//              s = Math.abs( F(s_min) ) < Math.abs( F(s_max) ) ? s_min : s_max,
//              c = Math.sqrt( (1+s)*(1-s) );
//            }
//
//            { // always choose the rotation that keeps the main diagonal in order
//              const r_ii = c*(c*R_ii - s*R_ij) - s*(c*R_ij - s*R_jj),
//                    r_jj = s*(s*R_ii + c*R_ij) + c*(s*R_ij + c*R_jj);
//              if( r_ii < r_jj ) { const C = c; c = -s; s = C; }
//            }
//
//            { // perform rotation
//              const    ik0= R_off + M* i;
//              for( let ik = R_off + M*(i+1),
//                       jk = R_off + M*(j+1); ik > ik0; )
//              {
//                const R_ik = R[--ik],
//                      R_jk = R[--jk];
//                R[ik] = c*R_ik - s*R_jk;
//                R[jk] = c*R_jk + s*R_ik;
//              }
//            }
//          }
//        }
////        console.log(`${sweeps} sweeps.`)
//
//        // COMPUTE SINGULAR VALUES AND NORMALIZE ROWS OF V
//        for( let i=0; i < M; i++ )
//        {
//          let sum = 0.0,
//              max = 0.0;
//          for( let j=0; j < M; j++ )
//          {
//            const R_ij = Math.abs(R[R_off + M*i+j]);
//            if(   R_ij != 0 ) { // <- handles NaN (by making the result NaN)
//              if( R_ij > max ) {
//                const scale = max / R_ij; max = R_ij;
//                sum *= scale*scale;
//              }
//              const ratio = R_ij / max;
//              sum += ratio*ratio;
//            }
//          }
//          const norm = isFinite(max) ? Math.sqrt(sum)*max : max;
//          for( let j=M; j-- > 0; ) R[R_off + M*i+j] /= norm;
//          sv[P_sv_off + i] = norm;
//        }
//
//        // COMPUTE U ( Using U = R @ inv(V) @ inv(SV) = R @ V.T @ diag(1/sv) )
//        //   ( should be numerically stable since V is orthogonal, i.e. well-conditioned ... right? )
//        for( let i=0; i < M; i++ ) {
//          for( let j=0; j < M; j++ ) { tmp[j] = U[M*i+j]; U[M*i+j]=0; }
//  
//          for( let j=0; j < M; j++ ) {
//            for( let k=0; k < M; k++ )
//              U[M*i+j] += R[R_off + M*j+k]*tmp[k];
//            U[M*i+j] /= sv[P_sv_off + j];
//          }
//        }
//
//        // MULTIPLY Q = Q @ U
//        for( let i=0; i < N; i++ )
//        {
//          for( let j=0; j < M; j++ ) { tmp[j] = Q[Q_off + M*i+j]; Q[Q_off + M*i+j]=0; }
//
//          for( let k=0; k < M; k++ )
//          for( let j=0; j < M; j++ ) Q[Q_off + M*i+j] += tmp[k] * U[M*k+j];
//        }
//
//        // APPLY P TO R (COLUMN PERMUTATIONS)
//        for( let i=0; i < M; i++ )
//        {
//          for( let j=0; j   < M; j++ ) tmp[P[P_sv_off + j]] = R[R_off + M*i+j];
//          for( let j=M; j-- > 0;     )                        R[R_off + M*i+j] = tmp[j];
//        }
//      }
//
//      return [
//        new nd.Array( U_shape, Q ),
//        new nd.Array(sv_shape, sv),
//        new nd.Array( V_shape, R )
//      ];
//    };


    la._svd_decomp_jac_classic = A => {
      A = nd.asarray(A);
      const
        shape = A.shape,
        N = shape[shape.length-2],
        TOL = Number.EPSILON * N;

      // QR DECOMPOSE RECTANGULAR MATRICES AND WORK ON (QUADRATIC) R
      {
        const M = shape[shape.length-1];
        // if A is not square use QR Decomposition
        if( N > M ) {
          const [Q,R]    = la.  qr_decomp(A);
          const [U,sv,V] = la._svd_decomp_jac_classic(R);
          return [la.matmul2(Q,U), sv, V];
        }
        if( N < M ) {
          const [Q,R]    = la.  qr_decomp(A.T);
          const [U,sv,V] = la._svd_decomp_jac_classic(R);
          la._transpose_inplace(V)
          return [V, sv, la.matmul2(Q,U).T];
        }
      }
      // ALLOCATE RESULT DATA
      const DTypeArray = nd.dtypes[A.dtype],
            U = DTypeArray.from(A.data); A = undefined; // <- potentially allow GC
      const D = new DTypeArray(N*N), // <- tempory storage for decomposition
            V = new DTypeArray(U.length),
            sv= new DTypeArray(U.length/N);

      if( 1 >  N ) throw new Error('Assertion failed.');
      if( 1 == N ) {
        for( let i=U.length; i-- > 0; )
          if( U[i] < +0.0 ) {
              U[i] *= -1.0;
             sv[i]  = -1.0;
          }
          else sv[i] = +1.0;
        return [
          new nd.Array(shape,            sv ),
          new nd.Array(shape.slice(0,-1), U ),
          new nd.Array(shape,      V.fill(1))
        ];
      }

       //
      // BUILD TRIANGLE TREE
     //
      // the size of each level of the tree
      const treeSizes = Int32Array.from( function*(){
        for( let n=N; n > 1; ) {
          n = n+1 >>> 1;
          yield n;
        }
      }() );
      const treeData = new DTypeArray( treeSizes.reduce( (len,N) => len + (N*N+N >>> 1), 0 ) );

      /** updates the triangle tree to reflect a changed row in D.
       */
      const update_row = row =>
      {
        row >>>= 1;
        let col=-1;
        // build bottom tree level
        {
          const hyp = (i,j) => {
            i += 2*row;
            j += 2*col; const D_ij = D[N*i+j],
                              D_ji = D[N*j+i]; return D_ij*D_ij + D_ji*D_ji;
          };
          const rowOff = row*row+row >>> 1;
          for( col=0; col < row; col++ )
                 treeData[rowOff + col] = Math.max( hyp(0,0), hyp(0,1) );
          if( 2*row < N-1 ) {
              for( col=0; col < row; col++ )
                 treeData[rowOff + col] = Math.max( hyp(1,0), hyp(1,1), treeData[(row*row+row >>> 1) + col] );
                 treeData[rowOff + col] = hyp(1,0);
          } else treeData[rowOff + col] = 0;
        }
        // build remaining tree levels
        let OFF = -1;
        const val = (i,j) => {
          i += 2*row;
          j += 2*col; return treeData[OFF + (i*i+i >>> 1) + j];
        }
        for( let off=0, h=1; h < treeSizes.length; h++ )
        {
          const N = treeSizes[h-1];
          OFF = off; off += N*N+N >>> 1;
          row >>>= 1;
          for( col=row; col >= 0; col-- )
          {
            let max = val(0,0);
            if( row*2 < N-1 ) max = Math.max( max, val(1,0), val(1,1) );
            if( row != col  ) max = Math.max( max, val(0,1) );
            treeData[off + (row*row+row >>> 1) + col] = max;
          }
        }
      }

      /** updates the triangle tree to reflect a changed row in D.
       */
      const update_col = col =>
      {
        col >>>= 1;
        let row=-1;
        // build bottom tree level
        {
          const hyp = (i,j) => {
            i += 2*row;
            j += 2*col; const D_ij = D[N*i+j],
                              D_ji = D[N*j+i]; return D_ij*D_ij + D_ji*D_ji;
          };
          for( row=treeSizes[0]; row-- > col; ) {
            let max = 0;
            if( 2*row < N-1 ) { const h10 = hyp(1,0); if( h10 > max ) max = h10; }
            if( 2*col < N-1 ) { const h01 = hyp(0,1); if( h01 > max ) max = h01; }
            if( row != col  ){{ const h00 = hyp(0,0); if( h00 > max ) max = h00; }
            if( 2*row < N-1 ) { const h11 = hyp(1,1); if( h11 > max ) max = h11; }}
            treeData[(row*row+row >>> 1) + col] = max;
          }
        }
        // build remaining tree levels
        let OFF = -1;
        const val = (i,j) => {
          i += 2*row;
          j += 2*col; return treeData[OFF + (i*i+i >>> 1) + j];
        }
        for( let off=0, h=1; h < treeSizes.length; h++ )
        {
          const
            N = treeSizes[h-1],
            n = treeSizes[h];
          OFF = off; off += N*N+N >>> 1;
          col >>>= 1;
          for( row=col; row < n; row++ )
          {
            let max = val(0,0);
            if( row*2 < N-1 ) max = Math.max( max, val(1,0), val(1,1) );
            if( row != col  ) max = Math.max( max, val(0,1) );
            treeData[off + (row*row+row >>> 1) + col] = max;
          }
        }
      }

      let UV_off=0;
      for( let sv_off=0; sv_off < sv.length; UV_off += N*N, sv_off += N )
      {
        // MOVE FROM U TO D
        for( let i=0; i < N; i++ )
        for( let j=0; j < N; j++ ) { D[N*i+j] = U[UV_off + N*i+j]; U[UV_off + N*i+j] = i != j ? 0 : 1 };
        // INIT V TO IDENTITY
        for( let i=0; i < N; i++ )
        for( let j=0; j < N; j++ ) V[UV_off + N*i+j] = i != j ? 0 : 1;
        // INIT TRIANGLE TREE
        for( let i=0; i < N; i += 2 ) update_row(i);

         //
        // (CLASSICAL) JACOBI SVD ITERATIONS
       //
        let nIter = 0;
        for( let s=-1,
                 t=-1;; nIter++ )
        {
          // FIND THE OFF DIAGONAL PAIR WITH THE LARGEST HYPOTHENUSE
          let k,l; {
            let i=0,
                j=0;
            for( let off=treeData .length,
                     h  =treeSizes.length; h-- > 0; )
            {
              const val = (k,l) => {
                k += i;
                l += j; return treeData[off + (k*k+k >>> 1) + l];
              };
              const n = treeSizes[h];
              off -= n*n+n >>> 1;       
              let max = val(0,0); k=i;
                                  l=j;
              if( i+1 < n ) { const v10 = val(1,0); if( v10 > max ) { max = v10; k=i+1; l=j  ; }
                              const v11 = val(1,1); if( v11 > max ) { max = v11; k=i+1; l=j+1; } }
              if( i != j )  { const v01 = val(0,1); if( v01 > max ) { max = v01; k=i;   l=j+1; } }
              i = 2*k;
              j = 2*l;
            }
            let max = -Infinity;
            const hyp = (s,t) => {
              s += i;
              t += j; const D_st = D[N*s+t],
                            D_ts = D[N*t+s]; return D_st*D_st + D_ts*D_ts;
            };
            if( i+1 < N ) { const h10 = hyp(1,0); if( h10 > max ) { max=h10; k=i+1; l=j  ; } }
            if( j+1 < N ) { const h01 = hyp(0,1); if( h01 > max ) { max=h01; k=i  ; l=j+1; } }
            if( i != j  ) { const h00 = hyp(0,0); if( h00 > max ) { max=h00; k=i  ; l=j  ; }
            if( i+1 < N ) { const h11 = hyp(1,1); if( h11 > max ) {          k=i+1; l=j+1; } } }
          };
          if( l >  k ) throw new Error('Assertion failed.')
          if( s == k &&
              t == l ) break; // <- DRY
          s=k;
          t=l;
//          {
//            // CHECK THAT THIS IS TRULY THE MAXIMUM
//            const hyp = (i,j) => {
//              D_ij = D[N*i+j],
//              D_ji = D[N*j+i]; return D_ij*D_ij + D_ji*D_ji;
//            };
//            for( let i=1; i < N; i++ )
//            for( let j=0; j < i; j++ )
//              if( ! (hyp(i,j) <= hyp(k,l)) ) throw new Error('Assertion failed');
//          }
          const
            D_kk = D[N*k+k], D_kl = D[N*k+l],
            D_lk = D[N*l+k], D_ll = D[N*l+l];
          if( ! ( Math.hypot(D_kl,D_lk) / Math.sqrt(Math.abs(D_kk)) / Math.sqrt(Math.abs(D_ll)) > TOL ) )
            break; // <- TODO check if really a good stopping criterion (there may be smaller off-diag. entries larger relative to their respective diag. entries)

          // determine rotation angles
          // ┌                ┐ ┌            ┐ ┌                ┐   ┌        ┐
          // │ cos(β) -sin(β) │ │ D_kk  D_kl │ │ cos(α) -sin(α) │ ! │ d1   0 │
          // │                │ │            │ │                │ = │        │
          // │ sin(β)  cos(β) │ │ D_lk  D_ll │ │ sin(α)  cos(α) │   │  0  d2 │
          // └                ┘ └            ┘ └                ┘   └        ┘
          //
          // such that: |d1| >= |d2|
          //
          // => 0 = cos(α)*{D_kl*cos(β) - D_ll*sin(β)} + sin(α)*{D_kk*cos(β) - D_lk*sin(β)} = ½⋅{(D_ll - D_kk)⋅sin(α-β) - (D_kk + D_ll)⋅sin(α+β) + (D_kl + D_lk)⋅cos(α-β) + (D_kl - D_lk)⋅cos(α+β)}
          //    0 = cos(α)*{D_kk*sin(β) + D_lk*cos(β)} - sin(α)*{D_kl*sin(β) + D_ll*cos(β)} = ½⋅{(D_ll - D_kk)⋅sin(α-β) + (D_kk + D_ll)⋅sin(α+β) + (D_kl + D_lk)⋅cos(α-β) - (D_kl - D_lk)⋅cos(α+β)}
          //   d1 = cos(α)*{D_kk*cos(β) - D_lk*sin(β)} - sin(α)*{D_kl*cos(β) - D_ll*sin(β)}
          //   d2 = cos(α)*{D_kl*sin(β) + D_ll*cos(β)} + sin(α)*{D_kk*sin(β) + D_lk*cos(β)}
          //
          // => 0 = {D_kl⋅cos(α) - D_kk⋅sin(α)}⋅cos(β) + {D_lk⋅sin(α) - D_ll⋅cos(α)}⋅sin(β)
          //    0 = {D_kl⋅sin(α) + D_kk⋅cos(α)}⋅sin(β) + {D_lk⋅cos(α) + D_ll⋅sin(α)}⋅cos(β)
          //   d1 = {D_kl⋅sin(α) + D_kk⋅cos(α)}⋅cos(β) + {D_lk⋅cos(α) - D_ll⋅sin(α)}⋅sin(β)
          //   d2 = {D_kl⋅cos(α) - D_kk⋅sin(α)}⋅sin(β) + {D_lk⋅sin(α) + D_ll⋅cos(α)}⋅cos(β)
          //
          // => 0 = (D_kk+D_ll)⋅sin(α-β) + (D_kl-D_lk)⋅cos(α-β)  
          //    0 = (D_kk-D_ll)⋅sin(α+β) + (D_kl+D_lk)⋅cos(α+β)  
          const [cα,sα,cβ,sβ] = function() {
            let Cα,Sα,Cβ,Sβ; {
              const m = Math.atan2(D_lk - D_kl, D_ll + D_kk),// = α - β
                    p = Math.atan2(D_lk + D_kl, D_ll - D_kk),// = α + β
                    α = (p+m)/2,
                    β = (p-m)/2;
              Cα = Math.cos(α); Sα = Math.sin(α);
              Cβ = Math.cos(β); Sβ = Math.sin(β);
            }

            // tan is 180°-periodical so lets try all possible solutions
            for( const [cα,sα,cβ,sβ] of [
              [ Cα, Sα,   Cβ, Sβ], // [p,m]+[  0°,  0°]
              [-Sα, Cα,   Sβ,-Cβ], // [p,m]+[  0°,180°]
              [-Sα, Cα,  -Sβ, Cβ], // [p,m]+[180°,  0°]
              [-Cα,-Sα,   Cβ, Sβ], // [p,m]+[180°,180°]
            ]) {
              const d00 = (D_kl*sα + D_kk*cα)*cβ + (D_lk*cα - D_ll*sα)*sβ,
                    d11 = (D_kl*cα - D_kk*sα)*sβ + (D_lk*sα + D_ll*cα)*cβ;
              // ROTATE IN A WAY THAT ENSURES DESCENDING ORDER
              if( d11 >= Math.abs(d00) ) {
                if( d00 >= -0.0 )
                  return [cα,sα,cβ,sβ];
                Cα = cα; Sα = sα;
                Cβ = cβ; Sβ = sβ;
              }
            }
            return [Cα,Sα,Cβ,Sβ];
          }();

          // ROTATE COLUMNS IN D
          for( let i=N; i-- > 0; ) {
            const D_ik = D[N*i+k],
                  D_il = D[N*i+l];
            D[N*i+k] = D_ik*cα - D_il*sα;
            D[N*i+l] = D_il*cα + D_ik*sα;
          }

          // ROTATE ROWS IN D
          for( let i=N; i-- > 0; ) {
            const D_ki = D[N*k+i],
                  D_li = D[N*l+i];
            D[N*k+i] = D_ki*cβ - D_li*sβ;
            D[N*l+i] = D_li*cβ + D_ki*sβ;
          }

          // ENTRIES (k,l) AND (l,k) ARE REMAINDERS (CANCELLATION ERROR) FROM ELIMINATION => SHOULD BE SAFELY ZEROABLE
          D[N*k+l] = 0.0;
          D[N*l+k] = 0.0;

          // UPDATE TRIANGLE TREE ROWS AND COLUMNS
          update_row(k); update_row(l);
          update_col(k); update_col(l);

          // ROTATE ROWS IN U (TRANSPOSED FOR CACHE LOCALITY REASONS)
          for( let i=0; i < N; i++ ) {
            const ki = UV_off + N*k+i, U_ki = U[ki],
                  li = UV_off + N*l+i, U_li = U[li];
            U[ki] = U_ki*cβ - U_li*sβ;
            U[li] = U_li*cβ + U_ki*sβ;
          }

          // ROTATE ROWS IN V
          for( let i=0; i < N; i++ ) {
            const ki = UV_off + N*k+i, V_ki = V[ki],
                  li = UV_off + N*l+i, V_li = V[li];
            V[ki] = V_ki*cα - V_li*sα;
            V[li] = V_li*cα + V_ki*sα;
          }
        }
//        console.log('Sweep ratio:', nIter / (N*N-N) )

        // MOVE D TO SV
        for( let i=0; i < N; i++ ) sv[sv_off + i] = D[N*i+i];

        // USE INSERTION SORT SV (should take O(N) because it is heavily pre-sorted)
        for( let i=0; i < N; i++ )
        {
          let sv_i = sv[sv_off + i];
          // flip sign
          if( sv_i < 0.0 ) {
            sv_i *= -1;
            sv[sv_off + i] *= -1;
            for( let k=0; k < N; k++ ) U[UV_off + N*i+k] *= -1;
          }
          // merge sort
          for( let j=i; j-- > 0; ) {
            const sv_j = sv[sv_off + j];
            if( sv_j >= sv_i ) break;
            // swap sv
            { const sv_j = sv[sv_off + j]; sv[sv_off + j] = sv[sv_off + j+1]; sv[sv_off + j+1] = sv_j; }
            // swap U and V rows
            for( let k=0; k < N; k++ ) { const U_jk = U[UV_off + N*j+k]; U[UV_off + N*j+k] = U[UV_off + N*(j+1)+k]; U[UV_off + N*(j+1)+k] = U_jk; }
            for( let k=0; k < N; k++ ) { const V_jk = V[UV_off + N*j+k]; V[UV_off + N*j+k] = V[UV_off + N*(j+1)+k]; V[UV_off + N*(j+1)+k] = V_jk; }
          }
        }

        // TRANSPOSE U IN-PLACE (TRANSPOSED FOR CACHE LOCALITY REASONS)
        for( let i=0  ; i < N-1; i++ )
        for( let j=1+i; j < N  ; j++ ) {
          const
            ij = UV_off + N*i+j,
            ji = UV_off + N*j+i, U_ij = U[ij]; U[ij] = U[ji]; U[ji] = U_ij;
        }
      }

      return [
        new nd.Array(shape,             U),
        new nd.Array(shape.slice(0,-1),sv),
        new nd.Array(shape,             V)
      ]
    }


    la._svd_decomp_jac_round_robin = A => {
      A = nd.asarray(A);
      const
        shape = A.shape,
        N = shape[shape.length-2],
        TOL = Number.EPSILON * N * 2;
      // QR DECOMPOSE RECTANGULAR MATRICES AND WORK ON (QUADRATIC) R
      {
        const M = shape[shape.length-1];
        // if A is not square use QR Decomposition
        if( N > M ) {
          const [Q,R]    = la.  qr_decomp(A);
          const [U,sv,V] = la._svd_decomp_jac_round_robin(R);
          return [la.matmul2(Q,U), sv, V];
        }
        if( N < M ) {
          const [Q,R]    = la.  qr_decomp(A.T);
          const [U,sv,V] = la._svd_decomp_jac_round_robin(R);
          la._transpose_inplace(V)
          return [V, sv, la.matmul2(Q,U).T];
        }
      }
      // ALLOCATE RESULT DATA
      const
        DTypeArray = nd.dtypes[ nd.super_dtype(A.dtype,'float64') ],
        U = DTypeArray.from(A.data);
      A = undefined; // <- potentially allow GC

      // the rotations are performed in a Round-Robin/Tournament style fashion
      // think of it as every index playing in tournament exactly once against
      // every other index
      // SEE: https://en.wikipedia.org/wiki/Round-robin_tournament#Scheduling_algorithm
      const
        N_ROUNDS= N - (N+1)%2, 
        N_GAMES = N+1 >>> 1; // <- number of games per round (including 1 Bye if necessary)

      const
        sv = new DTypeArray(U.length/N),
        V  = new DTypeArray(U.length),
        D  = new DTypeArray(N*N),
        csU= new DTypeArray(2*N_GAMES),
        csV= new DTypeArray(2*N_GAMES),
        playerOrder = Int32Array.from({ length: 2*N_GAMES }, (_,i) => i ),
        gameOrder   = Int32Array.from({ length:   N_GAMES }, (_,i) => i );

      if( 1 == N )
        return [
          new nd.Array(shape,           sv.fill(1)),
          new nd.Array(shape.slice(0,-1),U        ),
          new nd.Array(shape,            V.fill(1))
        ];

      // Why does this implementations use the two sided jacobi method?
      // Well there's multiple reasons.
      // The one-sided Jacobi method requires explicit handling of rank
      // deficiencies. This is usually done via RRQR Decomposition. Aside
      // from the fact that tedious tweaking of numeric tolerance is
      // necessary to correctly estimate the rank. RRQR might even fail,
      // especially if not implemented correctly (http://www.netlib.org/lapack/lawnspdf/lawn176.pdf).
      // The Strong RRQR would solve the problem (https://math.berkeley.edu/~mgu/MA273/Strong_RRQR.pdf),
      // but is even harder to implement both perfomantly and correctly.
      //
      // In the later stages of one-sided Jacobi, the dot product that
      // is required to compute an entry of the implied symmetric matrix
      // is one of the most expensive step since most of the entries are
      // already close enough to zero to skip during a sweep. To alleviate
      // this problem some additional housekeeping is required that makes
      // things more complicated. For two-sided jacobi, it's very cheap
      // to skip close-to-zero entries.

      for( let sv_off=0,
               UV_off=0; sv_off < sv.length;
               sv_off += N,
               UV_off += N*N
      )
      {
        // MOVE U -> D
        for( let i=0; i < N; i++ )
        for( let j=0; j < N; j++ ) {
          D[N*i+j] = U[UV_off + N*i+j];
                     U[UV_off + N*i+j] = i != j ? 0.0 : 1.0;
        }
        // INIT V
        for( let i=0; i < N; i++ ) V[UV_off + N*i+i] = 1;

        // PERFORM JACOBI SWEEPS
        let nSweeps = 1;
        for( let DONE=false; ! DONE; nSweeps++ )
        {
          if( nSweeps > 1024*1024*N ) throw new Error(`Iteration limit exceeded (N=${N}).`);
          DONE = true; // <- FIXME this should be checked for every element, right?

          // SHUFFLE PLAYER ORDER
          for( let i=playerOrder.length; i > 0; ) {
            const
              j = Math.floor( Math.random() * i-- ),
              tmp = playerOrder[i]; playerOrder[i] = playerOrder[j]; playerOrder[j] = tmp;
          }

           //
          // FIND THE OPTIMAL ROUND TO PLAY
         //
          let round=-1; {

            let err = -Infinity,
               norm = -Infinity;

            // [ALTERNATIVE 1] FIND THE ROUND WITH THE LARGEST OFF DIAGONAL ENTRY
            for( let rnd =N_ROUNDS; rnd -- > 0; )
            for( let game=N_GAMES ; game-- > 0; ) {
              const i = playerOrder[ game==0 ? 0 : 1 +           (game+rnd-1) % (2*N_GAMES-1) ],
                    j = playerOrder[               1 + (2*N_GAMES-game+rnd-2) % (2*N_GAMES-1) ],
                 D_ii = Math.abs(D[N*i+i]), D_ij = Math.abs(D[N*i+j]),
                 D_ji = Math.abs(D[N*j+i]), D_jj = Math.abs(D[N*j+j]);

              if( ! ( Math.hypot(D_ij,D_ji) / Math.sqrt(D_ii) / Math.sqrt(D_jj) > TOL ) )
                continue;
              DONE = false;

              const e = Math.hypot(D_ij,D_ji);
              if( ! (e <= err) ) {
                err=e;
                round=rnd;
              }
            }
//            console.log('ERR:', err)

//            // [ALTERNATIVE 2] FIND THE ROUND WITH THE LARGEST OFF DIAGONAL NORM
//            for( let rnd =N_ROUNDS; rnd -- > 0; )
//            {
//              let sum = 0.0, max = 0.0, add = elem => {
//                if( elem > 0 ) {
//                  if( elem > max ) {
//                    const scale = max / elem; max = elem; sum *= scale*scale;
//                  } const ratio = elem / max; sum += ratio*ratio;
//                }
//              };
//              for( let game=N_GAMES ; game-- > 0; )
//              {
//                const
//                  i = playerOrder[ game==0 ? 0 : 1 +           (game+rnd-1) % (2*N_GAMES-1) ],
//                  j = playerOrder[               1 + (2*N_GAMES-game+rnd-2) % (2*N_GAMES-1) ],
//                  D_ii = Math.abs(D[N*i+i]), D_ij = Math.abs(D[N*i+j]),
//                  D_ji = Math.abs(D[N*j+i]), D_jj = Math.abs(D[N*j+j]);
//
//                if( ! ( Math.hypot(D_ij,D_ji) / Math.sqrt(D_ii) / Math.sqrt(D_jj) > TOL ) )
//                  continue;
//
//                DONE = false;
//                add(D_ij)
//                add(D_ji)
//              }
//
//              const nrm = isFinite(max) ? Math.sqrt(sum)*max : max;
//              if( nrm > norm ) {
//                 norm = nrm;
//                round = rnd;
//              }
//            }
//            console.log('ERR:', norm)

            if(DONE) break;
          }


          if( round < 0 ) throw new Error('Assertion failed!')
//          for( let round=N_ROUNDS; round-- > 0; ) // <- TOURNAMENT ROUNDS
          {
            let nPlayed = 0;
            // PERFORM ROW ROTATIONS
            for( let game=N_GAMES; game-- > 0; ) { // <- TOURNAMENT GAMES
              const
                i = playerOrder[ game==0 ? 0 : 1 +           (game+round-1) % (2*N_GAMES-1) ],
                j = playerOrder[               1 + (2*N_GAMES-game+round-2) % (2*N_GAMES-1) ];
              if( i==N || j==N ) continue;

              const
                D_ii = D[N*i+i], D_ij = D[N*i+j],
                D_ji = D[N*j+i], D_jj = D[N*j+j];

              if( ! (Math.hypot(D_ij,D_ji) / Math.sqrt(Math.abs(D_ii)) / Math.sqrt(Math.abs(D_jj)) > TOL) ) {
                csU[2*game+0] = NaN; csV[2*game+0] = NaN;
                csU[2*game+1] = NaN; csV[2*game+1] =   0;
                continue;
              }
              gameOrder[nPlayed++] = game;
//              DONE=false;

              // determine rotation angles
              // ┌                ┐ ┌            ┐ ┌                ┐   ┌        ┐
              // │ cos(β) -sin(β) │ │ D_ii  D_ij │ │ cos(α) -sin(α) │ ! │ d1   0 │
              // │                │ │            │ │                │ = │        │
              // │ sin(β)  cos(β) │ │ D_ji  D_JJ │ │ sin(α)  cos(α) │   │  0  d2 │
              // └                ┘ └            ┘ └                ┘   └        ┘
              //
              // such that: |d2| >= |d1| <- FIXME should be other way around
              //
              // => 0 = cos(α)*(D_ij*cos(β) - D_jj*sin(β)) + sin(α)*(D_ii*cos(β) - D_ji*sin(β))
              //    0 = cos(α)*(D_ii*sin(β) + D_ji*cos(β)) - sin(α)*(D_ij*sin(β) + D_jj*cos(β))
              //   d1 = cos(α)*(D_ii*cos(β) - D_ji*sin(β)) - sin(α)*(D_ij*cos(β) - D_jj*sin(β))
              //   d2 = cos(α)*(D_ij*sin(β) + D_jj*cos(β)) + sin(α)*(D_ii*sin(β) + D_ji*cos(β))
              //
              // => 0 = (D_ii+D_jj)⋅sin(α-β) + (D_ij-D_ji)⋅cos(α-β)  
              //    0 = (D_ii-D_jj)⋅sin(α+β) + (D_ij+D_ji)⋅cos(α+β)  
              const [cα,sα,cβ,sβ] = function() {
                let Cα,Sα,Cβ,Sβ; {
                  const m = Math.atan2(D_ji - D_ij, D_jj + D_ii),// = α - β
                        p = Math.atan2(D_ji + D_ij, D_jj - D_ii);// = α + β
                        α = (p+m)/2,
                        β = (p-m)/2;

//                  const F = (α,β) => {
//                    const cα = Math.cos(α), sα = Math.sin(α),
//                          cβ = Math.cos(β), sβ = Math.sin(β);
//                    return Math.max(
//                      Math.abs( cα*(D_ij*cβ - D_jj*sβ) + sα*(D_ii*cβ - D_ji*sβ) ),
//                      Math.abs( cα*(D_ii*sβ + D_ji*cβ) - sα*(D_ij*sβ + D_jj*cβ) )
//                    );
//                  };
//                  // REFINE THE RESULT
//                  let nIter = 0
//                  for( const α_range = new DTypeArray(3),
//                             β_range = new DTypeArray(3);; ) {
//                    nIter++;
//                    if( nIter % 1e6 == 0 ) console.log(`(${α},${β} -> ${F(α,β)}) (${[D_ii,D_ij,D_ji,D_jj]})`);
//
//                    α_range[0] = α; α_range[1] = -nd.math.nextUp(-α); α_range[2] = +nd.math.nextUp(+α);
//                    β_range[0] = β; β_range[1] = -nd.math.nextUp(-β); β_range[2] = +nd.math.nextUp(+β);
//    
//                    let min = +Infinity;
//    
//                    for( const x of α_range )
//                    for( const y of β_range )
//                    {
//                      const F_xy = F(x,y)
//                      if( ! (F_xy >= min) ) {
//                        α = x;
//                        β = y; min = F_xy;
//                      }
//                    }
//                    if( α == α_range[0] &&
//                        β == β_range[0] ) break;
//                  }
//                  console.log(`OPT. DONE IN ${nIter} ITERATIONS! α,β=${α},${β} -> ${F(α,β)}`)

                  Cα = Math.cos(α), Sα = Math.sin(α);
                  Cβ = Math.cos(β), Sβ = Math.sin(β);
                }

                // tan is 180°-periodical so lets try all possible solutions
                for( const [cα,sα,cβ,sβ] of [
                  [ Cα, Sα,   Cβ, Sβ], // [p,m]+[  0°,  0°]
                  [-Sα, Cα,   Sβ,-Cβ], // [p,m]+[  0°,180°]
                  [-Sα, Cα,  -Sβ, Cβ], // [p,m]+[180°,  0°]
                  [-Cα,-Sα,   Cβ, Sβ], // [p,m]+[180°,180°]
                ]) {
                  let
                    d00 = cα*(D_ii*cβ - D_ji*sβ) - sα*(D_ij*cβ - D_jj*sβ),
                    d11 = cα*(D_ij*sβ + D_jj*cβ) + sα*(D_ii*sβ + D_ji*cβ);;
                  if( i > j ) { const tmp=d00; d00=d11; d11=tmp; }
                  // ROTATE IN A WAY THAT ENSURES DESCENDING ORDER
                  if( d00 >= Math.abs(d11) ) {
                    if( d11 >= 0 )
                      return [cα,sα,cβ,sβ];
                    Cα = cα; Sα = sα;
                    Cβ = cβ; Sβ = sβ;
                  }
                }
                return [Cα,Sα,Cβ,Sβ];
              }();

              csU[2*game+0] = cβ; csU[2*game+1] = sβ;
              csV[2*game+0] = cα; csV[2*game+1] = sα;

              // rotate row i and j
              for( let k=0; k < N; k++ ) {
                const
                  D_ik = D[N*i+k],
                  D_jk = D[N*j+k];
                D[N*i+k] = D_ik*cβ - D_jk*sβ;
                D[N*j+k] = D_jk*cβ + D_ik*sβ;
              }
            }

            // PERFORM COLUMN ROTATIONS
            for( let k=0; k < N; k++ )
            for( let iGame=nPlayed; iGame-- > 0; ) { // <- TOURNAMENT GAMES
              const game = gameOrder[iGame],
                       i = playerOrder[ game==0 ? 0 : 1 +           (game+round-1) % (2*N_GAMES-1) ],
                       j = playerOrder[               1 + (2*N_GAMES-game+round-2) % (2*N_GAMES-1) ];
              if( i==N || j==N ) continue;
              const cα = csV[2*game+0],
                    sα = csV[2*game+1],
                  D_ki = D[N*k+i],
                  D_kj = D[N*k+j];
              D[N*k+i] = D_ki*cα - D_kj*sα;
              D[N*k+j] = D_kj*cα + D_ki*sα;
//              if( k==i && Math.abs(D[N*i+j]) > 1e-8 ) throw new Error(`Assertion failed! ${D[N*i+j]}`);
//              if( k==j && Math.abs(D[N*j+i]) > 1e-8 ) throw new Error(`Assertion failed! ${D[N*j+i]}`);
              if( k==i ) D[N*i+j] = 0.0;
              if( k==j ) D[N*j+i] = 0.0;
            }

            // ROTATE U
            for( let k=0; k < N; k++ )
            for( let iGame=nPlayed; iGame-- > 0; ) { // <- TOURNAMENT GAMES
              const game = gameOrder[iGame],
                       i = playerOrder[ game==0 ? 0 : (game+round-1) % (2*N_GAMES-1)  +  1 ],
                       j = playerOrder[     (2*N_GAMES-game+round-2) % (2*N_GAMES-1)  +  1 ];
              if( i==N || j==N ) continue;
              const cβ = csU[2*game+0],
                    sβ = csU[2*game+1],
                  U_ki = U[UV_off + N*k+i],
                  U_kj = U[UV_off + N*k+j];
              U[UV_off + N*k+i] = U_ki*cβ - U_kj*sβ;
              U[UV_off + N*k+j] = U_kj*cβ + U_ki*sβ;
            }
            // ROTATE V
            for( let iGame=nPlayed; iGame-- > 0; ) { // <- TOURNAMENT GAMES
              const game = gameOrder[iGame],
                       i = playerOrder[ game==0 ? 0 : (game+round-1) % (2*N_GAMES-1)  +  1 ],
                       j = playerOrder[     (2*N_GAMES-game+round-2) % (2*N_GAMES-1)  +  1 ];
              if( i==N || j==N ) continue;
              const cα = csV[2*game+0],
                    sα = csV[2*game+1];
              if( sα == 0 ) continue;
              
              // rotate row i and j
              for( let k=0; k < N; k++ ) {
                const V_ik = V[UV_off + N*i+k],
                      V_jk = V[UV_off + N*j+k];
                V[UV_off + N*i+k] = V_ik*cα - V_jk*sα;
                V[UV_off + N*j+k] = V_jk*cα + V_ik*sα;
              }
            }
          }
        }

        // MOVE D -> sv
        for( let i=0; i < N; i++ ) sv[sv_off+i] = D[N*i+i];

        // USE INSERTION SORT SV (should take O(N) because it is heavily pre-sorted)
        for( let i=0; i < N; i++ )
        {
          let sv_i = sv[sv_off + i];
          if( sv_i < 0 ) {
            sv[sv_off + i] *= -1;
            sv_i           *= -1;
            for( k=0; k < N; k++ ) U[UV_off + N*i+k] *= -1;
          }
          for( let j=i; j-- > 0; ) {
            const sv_j = sv[sv_off + j];
            if( sv_j >= sv_i ) break;
            // swap sv
            { const sv_j = sv[sv_off + j]; sv[sv_off + j] = sv[sv_off + j+1]; sv[sv_off + j+1] = sv_j; }
            // swap U and V rows
            for( let k=0; k < N; k++ ) { const U_jk = U[UV_off + N*j+k]; U[UV_off + N*j+k] = U[UV_off + N*(j+1)+k]; U[UV_off + N*(j+1)+k] = U_jk; }
            for( let k=0; k < N; k++ ) { const V_jk = V[UV_off + N*j+k]; V[UV_off + N*j+k] = V[UV_off + N*(j+1)+k]; V[UV_off + N*(j+1)+k] = V_jk; }
          }
        }
      }

      return [
        new nd.Array(shape,             U),
        new nd.Array(shape.slice(0,-1),sv),
        new nd.Array(shape,             V)
      ]
    }


//    la.svd_decomp = A => la.svd_decomp_jac_1sided(A)
    la.svd_decomp = A => la.svd_decomp_jac_classic  (A)
//    la.svd_decomp = A => la.svd_decomp_jac_round_robin(A)


    la.schur_eigenvals = T =>
    {
      const [N] = T.shape.slice(-1),
             FloatArray   = nd.dtypes[ nd.super_dtype(   'float64',T.dtype) ],
             ComplexArray = nd.dtypes[ nd.super_dtype('complex128',T.dtype) ],
             EV_shape = T.shape.slice(0,-1);
      let    EV = new FloatArray(T.data.length/N);
      if( T.shape[T.ndim-2] != N )
        throw new Error('T is not square.');
      T = T.data;

      for( let T_off=0,
              EV_off=0; T_off < T.length; T_off += N*N,
                                         EV_off += N )
      {
        const t = (i,j) => T[T_off + N*i+j]; 

        for( let i=0; i < N-1; )
        {
          const j = i+1;
          if( t(j,i) == 0 ) {
            EV[EV_off + i] = T[T_off + N*i+i];
            EV[EV_off + j] = T[T_off + N*j+j]; ++i;
          } else {
            // FIXME make this work for complex T
            const tr = t(i,i) + t(j,j),
                 det = t(i,i) * t(j,j)  -  t(i,j) * t(j,i);
            let ev1,ev2;
            const sqrt = nd.math.sqrt(tr*tr - 4*det);
            if( tr*tr >= 4*det )
            {
              let sign = tr < 0 ? -1.0 : +1.0;
              ev1  =       0.5 * (tr + sign*sqrt);
              ev2  = det * 2.0 / (tr + sign*sqrt); // <- Citardauq Formula
            } else {
              if( ! (EV instanceof ComplexArray) ) EV = ComplexArray.from(EV);
              ev1 = math.mul( math.add(tr,sqrt), 0.5 );
              ev2 = math.mul( math.sub(tr,sqrt), 0.5 );
            }
            EV[EV_off + i++] = ev1;
            EV[EV_off + i++] = ev2;
          }
        }
      }

      return new nd.Array(EV_shape,EV);
    }


    la.schur_eigen = (Q,T) => {
      if( Q.ndim != T.ndim ) throw new Error('Q.ndim != T.ndim.');
      for( let i=T.ndim; i-- > 0; )
        if( Q.shape[i] != T.shape[i] ) throw new Error('Q.shape != T.shape.')

      const [N] = T.shape.slice(-1),
        ComplexArray = nd.dtypes[ nd.super_dtype('complex128',T.dtype) ],
        V_shape = T.shape,
        Λ_shape = V_shape.slice(0,-1);

      // T is quasi-triangular. For simplification let's assume it's upper triangular:
      // ┌                        ┐
      // │ λ₁ ...                 │
      // │                        │
      // │ 0   λ₂ ...             │
      // │       .                │
      // │ 0   0   .              │
      // │           .            │
      // │ .       .   λₖ ...     │
      // │ .         .   .        │
      // │ .           .   .      │
      // │                   .    │
      // │ 0    ...        0   λₙ │
      // └                        ┘
      //
      // Let's say we want to find an eigenvector for λₖ, we can just solve:
      // ┌                                   ┐ ┌     ┐
      // │ λ₁-λₖ ...                         │ │ x₁  │
      // │                                   │ │     │
      // │ 0   λ₂-λₖ ...                     │ │ x₂  │
      // │          .                        │ │ .   │
      // │ 0     0    .                      │ │ .   │
      // │              .                    │ │ .   │ ! ⇀
      // │ .         .  λₖ₋₁-λₖ ...          │ │ xₖ₋₁│ = 0
      // │                                   │ │     │
      // │ .             .    λₖ-λₖ ...      │ │ 1   │
      // │                        .          │ │     │
      // │ .                 .      .        │ │ 0   │
      // │                            .      │ │ :   │
      // │ 0     .   .   .       0    λₙ-λₖ  │ │ 0   │
      // └                                   ┘ └     ┘
      //
      // Since T is quasi-triangular this is solvable via a modified Backward Substition.
      // For multiple (equivalent) eigenvalues, there is not guaranteed to be more than
      // one linearily independent eigenvector. If that is the case, we at some point
      // arrive at an row in the backward substitution where:
      //
      // (λₛ-λₖ)xₛ = 0·xₛ = 0 - T[s,s+1]·xₛ₊₁ ... - T[s,k-1]*xₖ₋₁ - T[s,k] = -xₛ ≠ 0.
      //
      // We can resolve this by setting:
      //
      // xₜ := 0; for t > s
      // xₛ := 1
      //
      // This will of course than mean that for s and k, we have the same eigenvector.
      //

      if( T.shape[T.ndim-2] != N ) throw new Error('Q is not square.');
      
      const V  =     ComplexArray.from(T.data); T = undefined;
      const Λ  = new ComplexArray(V.length/N),
            // temporary vectors for the eigenvalue computation
            v1 = new ComplexArray(N),
            v2 = new ComplexArray(N);

      const TOL = Number.EPSILON*N;

      /** Computes indices j < J of an eigenvector v using backward substition (see amazing UTF-8 art above).
       */
      function computeVec( λ, v, J )
      {
        const K = Math.min(N,J+2);

        for( let j=J; j-- > 0; )
        {
          for( let k=K; --k > j; )
            v[j] = math.sub(v[j], math.mul(v[k], V[V_off + N*j+k]) );
          if( 0==j || t(j,j-1) == 0 )
          {  //
            // 1x1 BLOCK
           //
            const V_jj = math.sub( V[V_off + N*j+j], λ );
            if( V_jj == 0 ) {
              if( math.is_close(v[j],0) ) { // <- FIXME find a better estimation of near-zeroness (e.g. the norm of the summands?)
                // v is already a valid eigenvalue, let's return it
                v[j] = 0;
                return;
              }
              // v is invalid, let's reset
              v[j] = 1.0;
              v.fill(0.0, j+1,K);
            }
            else v[j] = math.div(v[j],V_jj);
          }
          else
          {  //
            // 2x2 BLOCK
           //
            const i = j-1;
            for( let k=K; --k > j; )
              v[i] = math.sub(v[i], math.mul(v[k], V[V_off + N*i+k]) );
            const T_ii = math.sub( t(i,i), λ ), T_ij = t(i,j),
                  T_jj = math.sub( t(j,j), λ ), T_ji = t(j,i),
                  det = math.sub(
                    math.mul(T_ii,T_jj),
                    math.mul(T_ij,T_ji)
                  );
            if( det == 0 ) throw new Error('Assertion failed.');
            const v_j = math.div( math.sub( math.mul(T_ii,v[j]), math.mul(T_ji,v[i]) ), det );
            const v_i = math.div( math.sub( math.mul(T_jj,v[i]), math.mul(T_ij,v[j]) ), det ); 
            v[i  ] = v_i;
            v[j--] = v_j;
          }
        }
      }

      let V_off=0;
      const t = (i,j) => V[V_off + N*i+j]; 
      for( let Λ_off=0; V_off < V.length; V_off += N*N,
                                          Λ_off += N )
      {

        // COMPUTE EIGENVECTORS (right -> left)
        for( let j=N-1; j >= 0; j-- )
        {
          const i = j-1;
          if( 0==j || t(j,i) == 0 ) {
             //
            // 1x1 BLOCK
           //
            // the eigenvalue is the diagonal value
            const λ = V[V_off + N*j+j];
            Λ[Λ_off + j] = λ;
            // since 0*1 is zero, the eigenequation should be solvable for vec1[j] = 1
            // (unless there is a duplicate eigenvalue with linarily non-independent eigenvectors, but that will be resolved by computeVec)
            v1.fill(0.0, 0,j+2);
            v1[j] = 1;
            computeVec(λ,v1,j);
            // write the solution in the the (j+1)-th column
            for( let k=Math.min(N,j+2); k-- > 0; )
              V[V_off + N*k+j] = v1[k];
          } else {
             //
            // 2x2 BLOCK
           //
            // STEP1: compute eigenpairs of the 2x2 matrix
            const T_ii = t(i,i), T_ij = t(i,j),
                  T_ji = t(j,i), T_jj = t(j,j),
                  tr = math.add( T_ii, T_jj ),
                  det = math.sub(
                    math.mul(T_ii,T_jj),
                    math.mul(T_ij,T_ji)
                  );
            let λ1,λ2;
            const sqrt = math.sqrt( math.sub(
              math.mul(tr,tr),
              math.mul(det,4)
            ));
            if( tr*tr >= 4*det )
            { throw new Error('The given Schur decomposition must not contain 2x2 diagonal blocks with real eigenvalues.');
              let sign = tr < 0 ? -1.0 : +1.0;
              λ1  =       0.5 * (tr + sign*sqrt);
              λ2  = det * 2.0 / (tr + sign*sqrt); // <- Citardauq Formula
            } else {
              λ1 = math.mul( math.add(tr,sqrt), 0.5 );
              λ2 = math.mul( math.sub(tr,sqrt), 0.5 );
            }
            // TODO: the whole following section should be feasible with only a single temporary vector instead of two (vec1,vec2)
            v1.fill(0.0, 0,j+1);
            v2.fill(0.0, 0,j+1);
            // http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/
            if( math.abs(T_ij) >= math.abs(T_ji) )
            {        v1[i] = T_ij; v1[j] = math.sub( λ1, T_ii );
                     v2[i] = T_ij; v2[j] = math.sub( λ2, T_ii );
            } else { v1[j] = T_ji; v1[i] = math.sub( λ1, T_jj );
                     v2[j] = T_ji; v2[i] = math.sub( λ2, T_jj );
            }
            computeVec(λ1,v1,i);
            computeVec(λ2,v2,i);
            for( let k=j+1; k-- > 0; ) {
              V[V_off + N*k+i] = v1[k];
              V[V_off + N*k+j] = v2[k];
            }
            Λ[Λ_off + i  ] = λ1;
            Λ[Λ_off + j--] = λ2;
          }
        }

        // COMPUTE COLUMN NORMS
        const norm_sum = v1,
              norm_max = v2;
        norm_sum.fill(0.0);
        norm_max.fill(0.0);
        for( let i=0; i < N; i++ )
        for( let j=0; j < N; j++ ) {
          const V_ij = math.abs(V[V_off + N*i+j]);
          if(   V_ij > 0 ) {
            if( V_ij > norm_max[j] ) {
              const scale = norm_max[j] / V_ij; norm_max[j] = V_ij;
              norm_sum[j] *= scale*scale;
            }
            const ratio = V_ij / norm_max[j];
            norm_sum[j] += ratio*ratio;
          }
        }
        const norm = norm_sum;
        for( let i=0; i < N; i++ ) {
          const max = norm_max[i];
          norm[i] = isFinite(max) ? Math.sqrt(norm_sum[i])*max : max;
        }

        // NORMALIZE COLUMNS
        for( let i=0; i < N; i++ )
        for( let j=0; j < N; j++ )
          V[V_off + N*i+j] = math.div(V[V_off + N*i+j], norm[j]);
      }

      return [
                          new nd.Array(Λ_shape, Λ),
        nd.la.matmul2( Q, new nd.Array(V_shape, V) ),
      ];
    }


    la.eigen = A => {
      const [Q,T] = la.schur_decomp(A);
      return la.schur_eigen(Q,T);
    }


    la.schur_decomp = A => {
      // HESSENBERG DECOMPOSITION
      const
        N = A.shape[A.ndim-1],
        [Q,H] = la.hessenberg_decomp(A); A = undefined;
      // FRANCIS QR
      la._schur_decomp_qrfrancis_INPLACE(Q,H);
      return [Q,H];
    }


    la._schur_preprocess_balance = (A,p) =>
    {
      // SEE: RODNEY JAMES, JULIEN LANGOU, AND BRADLEY R. LOWERY,
      //      "ON MATRIX BALANCING AND EIGENVECTOR COMPUTATION"
      if( p == null ) p = 2;
      if( p > Number.MAX_VALUE ) return la._schur_preprocess_balance_inf(A);

      const N = A.shape[A.ndim-1],
            DTypeArray = nd.dtypes[A.dtype],
            A_shape = A.shape,
            D_shape = A_shape.slice(0,-1);
      if( ! (p >= 1) ) throw new Error(`Invalid norm p=${p};`);
      if( A.shape[A.ndim-2] != N ) throw new Error('A is not square');
      A = DTypeArray.from(A.data); // <- protection copy

      const D = new DTypeArray(A.length/N),
          TOL = 0.95 ** (1/p);
      D.fill(1.0);

      for( let A_off=A.length,
               D_off=D.length; A_off > 0; ) { A_off -= N*N;
                                              D_off -= N;
        for( let i,   done=false; !done; )
        for(     i=0, done=true; i < N; i++ )
        { let c=0.0, c_max=0.0,
              r=0.0, r_max=0.0;
          // COMPUTE ROW AND COLUMN NORM
          for( let j=0; j < N; j++ )
//            if(true)
            if( i != j )
            {{const A_ij = Math.abs(A[A_off + N*i+j]);
              if(   A_ij > 0 ) {
                if( A_ij > r_max ) {
                  const scale = r_max / A_ij ; r_max = A_ij; r *= scale**p;
                } const ratio =         A_ij / r_max;        r += ratio**p;
              }
            }{const A_ji = Math.abs(A[A_off + N*j+i]);
              if(   A_ji > 0 ) {
                if( A_ji > c_max ) {
                  const scale = c_max / A_ji ; c_max = A_ji; c *= scale**p;
                } const ratio =         A_ji / c_max;        c += ratio**p;
              }
            }}
          r = ! isFinite(r) ? r : r**(1/p) * r_max;
          c = ! isFinite(c) ? c : c**(1/p) * c_max;

          if( r*c == 0.0 ) continue;
          if( ! isFinite(r*c) ) throw new Error('NaN encountered.');
          let old_norm = ( c >= r
            ? ( 1 + (r/c)**p )**(1/p) * c
            : ( 1 + (c/r)**p )**(1/p) * r
          );

          // FIND A POWER OF TWO THAT MORE OR LESS BALANCES r AND c
          let scale = 1.0
          while( r >= c*2 ) { c*=2; r/=2; scale*=2; }
          while( c >= r*2 ) { c/=2; r*=2; scale/=2; }

          let new_norm = ( c >= r
            ? ( 1 + (r/c)**p )**(1/p) * c
            : ( 1 + (c/r)**p )**(1/p) * r
          );

          if( new_norm >= TOL*old_norm ) continue

          done = false;
          D[D_off + i] *= scale;
          for( let j=0; j < N; j++ ) {
            A[A_off + N*i+j] /= scale;
            A[A_off + N*j+i] *= scale;
          }
        }
      }

      return [
        new nd.Array(D_shape, D),
        new nd.Array(A_shape, A)
      ];
    }


    la._schur_preprocess_balance_inf = A =>
    {
      // SEE: RODNEY JAMES, JULIEN LANGOU, AND BRADLEY R. LOWERY,
      //      "ON MATRIX BALANCING AND EIGENVECTOR COMPUTATION"
      const N = A.shape[A.ndim-1],
            DTypeArray = nd.dtypes[A.dtype],
            A_shape = A.shape,
            D_shape = A_shape.slice(0,-1);
      if( A.shape[A.ndim-2] != N ) throw new Error('A is not square');
      A = DTypeArray.from(A.data);

      const D = new DTypeArray(A.length/N);
      D.fill(1.0);

      for( let A_off=A.length,
               D_off=D.length; A_off > 0; ) { A_off -= N*N;
                                              D_off -= N;
        for( let i,   done=false; !done; )
        for(     i=0, done=true; i < N; i++ )
        { let c=0.0,
              r=0.0;
          // COMPUTE ROW AND COLUMN NORM
          for( let j=0; j < N; j++ )
            if( i != j ) {
              const A_ij = Math.abs(A[A_off + N*i+j]); r = Math.max(r,A_ij);
              const A_ji = Math.abs(A[A_off + N*j+i]); c = Math.max(c,A_ji);
            }

          let old_norm = Math.max(c,r);
          if( old_norm == 0.0 ) continue;
          if( ! isFinite(old_norm) ) throw new Error('NaN encountered.');

          // FIND A POWER OF TWO THAT MORE OR LESS BALANCES r AND c
          let scale = 1.0
          while( r >= c*2 ) { c*=2; r/=2; scale*=2; }
          while( c >= r*2 ) { c/=2; r*=2; scale/=2; }

          let new_norm = Math.max(c,r);
          if( new_norm >= old_norm ) continue;

          done = false;
          D[D_off + i] *= scale;
          for( let j=0; j < N; j++ ) {
            A[A_off + N*i+j] /= scale;
            A[A_off + N*j+i] *= scale;
          }
        }
      }

      return [
        new nd.Array(D_shape,D),
        new nd.Array(A_shape,A)
      ];
    }


    la._schur_decomp_unbalance_INPLACE = (P,D,A) => {
      throw new Error('Not yet implemented.')
    }

    /** Takes a Hessenberg Decomposition as input and performs a real Schur Decomposition IN-PLACE.
     *  Does not perform any scaling.
     *  Does not check Hessenberg property.
     */
    la._schur_decomp_qrfrancis_INPLACE = (Q,H) =>
    {
      const N = Q.shape[Q.ndim-1],
            DTypeArray = nd.dtypes[Q.dtype],
            tmp = new DTypeArray(N);
      if( Q.shape[Q.ndim-2] != N ) throw new Error('Q is not square.');
      if( Q.ndim != H.ndim ) throw new Error('Q.ndim != H.ndim.')
      for( let i=Q.ndim; i-- > 0; )
        if( Q.shape[i] != H.shape[i] ) throw new Error('Q.shape != H.shape.')
      Q = Q.data;
      H = H.data;

      const TOL = Number.EPSILON;

      /** This function is recursively called to perform the Francis QR Algorithm, which is
       *  an implicit double-shift version of the QR Algorithm. This function only works on
       *  a subregion of H and an independent Q on each recursive call. That nested Q is
       *  applied to the remaining H and the outer Q once the nested schur decomposition
       *  is finished.
       *
       *  SEE: Prof. Dr. Peter Arbenz
       *       252-0504-00 G
       *       Numerical Methods for Solving Large Scale Eigenvalue Problems
       *       (Spring semester 2018)
       *       Chapter 3: The QR Algorithm
       *       http://people.inf.ethz.ch/arbenz/ewp/Lnotes/2010/chapter3.pdf 
       */
      const francis_qr = ( Q,Q_off,Q_stride, H_off ) =>
      {
        const h = (i,j) => H[H_off + N*i+j],
              is_zero = i => {
                if( Math.abs(h(i,i-1)) > TOL * ( Math.abs(h(i-1,i-1)) + Math.abs(h(i,i)) ) ) // <- goes to else on NaN
                  return false;
                else {
                  H[H_off + N*i+i-1] *= 0.0; // <- Handles NaN. If a value is that small, its digits are likely nonsense (due to cancellation error) so let's set it to zero.
                  return true;
                }
              },
              /** Applies a two-sided given rotation.
               */
              giv = (i,j,c,s) => {
                if( j <= i ) throw new Error('Assertion failed.')
                // ROTATE ROWS IN H
                for( let k=Math.max(0,i-1); k < Q_stride; k++ )
                { const H_i = H[H_off + N*i+k],
                        H_j = H[H_off + N*j+k];
                  H[H_off + N*i+k] = s*H_j + c*H_i;
                  H[H_off + N*j+k] = c*H_j - s*H_i;
                }
                // ROTATE COLUMNS IN H
                for( let k=Math.min(Q_stride,j+2); k-- > 0; )
                { const H_i = H[H_off + N*k+i],
                        H_j = H[H_off + N*k+j];
                  H[H_off + N*k+i] = s*H_j + c*H_i;
                  H[H_off + N*k+j] = c*H_j - s*H_i;
                }
                // ROTATE ROWS IN Q
                for( let k=Q_stride; k-- > 0; )
                { const Q_i = Q[Q_off + Q_stride*i+k],
                        Q_j = Q[Q_off + Q_stride*j+k];
                  Q[Q_off + Q_stride*i+k] = s*Q_j + c*Q_i;
                  Q[Q_off + Q_stride*j+k] = c*Q_j - s*Q_i;
                }
              },
              /** Recursively performs the schur-decomposition of the sub-region [s,e).
               * 
               *  During iteration, the givens rotations are only applied to current
               *  subregion and are accumulated in a temporary matrix. Only after the
               *  Francis QR algorithm is done, the transformations are applied to 
               */
              recurse = (s,e) => {
                if( e-s > Q_stride>>>1 ) throw new Error('Assertion failed.'); // <- assert that memory is bounded by O(n)
                const n = e-s; if( n < 3 ) throw new Error('Assertion failed.');
                const q = new DTypeArray(n*n);
                // INIT q TO IDENTITY
                for( let i=0; i < n; i++ ) q[n*i+i] = 1.0;
    
                // RUN FRANCIS QR ON SUB-PROBLEM
                francis_qr(q,0,n, H_off + N*s+s );
    
                // TRANSPOSE q
                for( let i=0; i < n; i++ )
                for( let j=0; j < i; j++ ) { const q_ij = q[n*i+j]; q[n*i+j] = q[n*j+i]; q[n*j+i] = q_ij; }
    
                // APPLY q TO LEFT OF Q (Q' = q.T @ Q)
                for( let i=0; i < Q_stride; i++ ) // <- each column in Q
                { tmp.fill(0.0, 0,n);
                  for( let j=0; j < n; j++ )
                  for( let k=0; k < n; k++ ) tmp[k] += q[n*j+k] * Q[Q_off + Q_stride*(s+j)+i];
    
                  for( let j=0; j < n; j++ ) Q[Q_off + Q_stride*(s+j)+i] = tmp[j];
                }
                // APPLY q TO LEFT OF H (H' = q.T @ H)
                for( let i=e; i < Q_stride; i++ ) // <- each column in H
                  { tmp.fill(0.0, 0,n);
                    for( let j=0; j < n; j++ )
                    for( let k=0; k < n; k++ ) tmp[k] += q[n*j+k] * H[H_off + N*(s+j)+i];
    
                    for( let j=0; j < n; j++ ) H[H_off + N*(s+j)+i] = tmp[j];
                  }
                // APPLY q TO RIGHT H (H" = H' @ q)
                for( let i=0; i < s; i++ ) // <- each row in H
                  { tmp.fill(0.0, 0,n);
                    for( let j=0; j < n; j++ )
                    for( let k=0; k < n; k++ ) tmp[k] += q[n*j+k] * H[H_off + N*i+(s+j)];
    
                    for( let j=0; j < n; j++ ) H[H_off + N*i+(s+j)] = tmp[j];
                  }
              };
        let stuck_o_meter = 0,
            start = 0,
            end   = Q_stride;

        while(true)
        { // DETECT ZEROS ON THE SUB-DIAGONAL AND SHRINK WORK SIZE ACCORDINGLY
          for( let done=false; !done; )
          {
            if( end-start < 3 ) return;
            if( end-start < Q_stride>>>4 ) return recurse(start,end); // <- ZOOM IN IF SUBPROBLEM IS SMALL ENOUGH

            done = false;
                 if( is_zero(start+1) ) start+=1; else if( is_zero(end-1) ) end-=1;
            else if( is_zero(start+2) ) start+=2; else if( is_zero(end-2) ) end-=2;
            else done = true;

            const mid = start+end >>> 1;
            for( let i=start+2; done && ++i < end-2; )
              if( is_zero(i) )
              { done = false;
                // RUN NESTED FRANCIS QR ON THE SMALLER OF THE TWO 
                if( i > mid ) { recurse(i,  end);   end=i; }
                else          { recurse(start,i); start=i; }
              }

            if( !done ) stuck_o_meter = 0;
          }
          stuck_o_meter += 1;

          // DETERMINE (DOUBLE) SHIFT FROM LOWER RIGHT 2x2 BLOCK
          let i = end-2,
              j = end-1,
              tr = h(i,i) + h(j,j),
              det= h(i,i) * h(j,j)  -  h(i,j) * h(j,i);

          // FOR REAL EIGENVALUES LETS USE THE ONE THAT'S CLOSER TO THE LOWER RIGHT ENTRY
          if( tr*tr > 4*det )
          { let sign = tr >= 0 ? +1.0 : -1.0,
                ev1  =       0.5 * (tr + sign*Math.sqrt(tr*tr - 4*det)),
                ev2  = det * 2.0 / (tr + sign*Math.sqrt(tr*tr - 4*det)); // <- Citardauq Formula
            // use the eigenvalue closer to A[j,j]
            // SEE: Bindel, Fall 2016, Matrix Computations (CS 6210)
            //      Notes for 2016-10-24
            //      https://www.cs.cornell.edu/~bindel/class/cs6210-f16/lec/2016-10-24.pdf
            if( Math.abs(h(j,j) - ev1) > Math.abs(h(j,j) - ev2) )
              ev1 = ev2;
            tr = ev1*2;
            det= ev1*ev1;
          }

          // IF WE'RE STUCK, LET'S WIGGLE LIKE A FISH ON LAND (... well except maybe that fella: https://www.youtube.com/watch?v=fJLCSsnhLFc)
          // SEE: NUMERICAL RECIPES Webnote No. 16, Rev. 1
          //      Description of the QR Algorithm for Hessenberg Matrices
          //      http://numerical.recipes/webnotes/nr3web16.pdf
          if( stuck_o_meter % 16 == 0 ) {
            if( stuck_o_meter > 1024*1024 ) throw new Error('Too many iterations for a single eigenvalue.');
            tr  = Math.abs(h(j,i)) + Math.abs(h(i,end-3))
            det = tr*tr
            tr *= Math.random()*0.5 + 1.25
          }

          // FIRST COLUMN OF DOUBLE SHIFTED MATRIX (H-sI)*(H-conj(s)I) = H² - 2*Re(s)*H + |s|²I
              i = start+0;
              j = start+1;
          let k = start+2,
              a1 = h(i,i)* h(i,i) + h(i,j)*h(j,i) - tr*h(i,i) + det,
              a2 = h(j,i)*(h(i,i) + h(j,j)        - tr),
              a3 = h(j,i)* h(k,j);

          // APPLY "DOUBLE SHIFTED" GIVENS ROTATIONS
          for( let row=2; row-- > 0; j=k, a2=a3 )
            if( a2 != 0 )
            { const      norm = Math.hypot(a1,a2),
                c = a1 / norm,
                s = a2 / norm;
              giv(i,j,c,s); // if( Math.abs(c*a2 - s*a1) > 1e-8 ) throw new Error('Assertion failed.')
              a1 = c*a1 + s*a2
            }

          // REINSTATE HESSENBERG PROPERTY
          for( let col=start; col < end-2; col++ )
          { i = col+1
            const J = Math.min(end,col+4);
            for( j=col+2; j < J; j++ )
            { const H_i = h(i,col),
                    H_j = h(j,col);
              if( H_j == 0 ) continue;
              const norm = Math.hypot(H_i,H_j),
                    c = H_i / norm,
                    s = H_j / norm;
              giv(i,j,c,s); // if( Math.abs( h(j,col) ) > 1e-8 ) throw new Error('Assertion failed!')
              H[H_off + N*j+col] *= 0.0; // <- handles NaN
            }
          }
        }
      }

      // BEGIN SCHUR DECOMPOSING MATRICES
      for( let off=0; off < Q.length; off += N*N )
      { // TRANSPOSE Q
        for( let i=0; i < N; i++ )
        for( let j=0; j < i; j++ ) { const Q_ij = Q[off + N*i+j]; Q[off + N*i+j]=Q[off + N*j+i]; Q[off + N*j+i] = Q_ij; }

        // RUN FRANCIS QR ALGORITHM
        francis_qr(Q,off,N, off);

        // BEGIN RESOLVE REAL-VALUED 2x2 BLOCKS
        for( let j=1; j < N; j++ )
        { const i = j-1;
          if( H[off + N*j+i] != 0 )
          { // The goal is to find a givens rotation that Schur-decomposes a real-eigenvalue 2x2 matrix.
            // ┌                ┐ ┌            ┐ ┌                ┐   ┌      ┐
            // │ cos(α) -sin(α) │ │ H_ii  H_ij │ │ cos(α) -sin(α) │ ! │ λ₁ p │
            // │                │ │            │ │                │ = │      │
            // │ sin(α)  cos(α) │ │ H_ji  H_jj │ │ sin(α)  cos(α) │   │ 0  λ₂│ => 0 == (H_ji⋅cos(α) - H_ii⋅sin(α))⋅cos(α) + (H_jj⋅cos(α) - H_ij⋅sin(α))⋅sin(α)
            // └                ┘ └            ┘ └                ┘   └      ┘ => 0 == (H_jj-H_ii)⋅sin(2⋅α) + (H_ij+H_ji)⋅cos(2⋅α) + H_ji-H_ij =: f(α)
            //
            // A simple and very numerically accurate solution would be Binary Search. In order to do that, we have to bracket a solution. So let's determine the extrema of f(α).
            // 
            // f'(α_max) =!= 0 = 2*(H_jj-H_ii)⋅cos(2⋅α_max) - 2*(H_ij+H_ji)⋅sin(2⋅α_max) => α_max = atan2( H_jj-H_ii, H_ij+H_ji ) / 2 + n*π/2
            const H_ii = H[off+N*i+i], H_ij = H[off+N*i+j],
                  H_ji = H[off+N*j+i], H_jj = H[off+N*j+j], tr = H_ii+H_jj, det= H_ii*H_jj - H_ij*H_ji;

            if( tr*tr < 4*det ) continue;

            const α_min =         Math.atan2(H_jj-H_ii, H_ij+H_ji) / 2,
                  α_max = α_min + Math.PI/2 * (α_min <= 0 ? +1 : -1),
                  α = nd.opt.root1d(
                    α => ( ( H_ji*Math.cos(α) - H_ii*Math.sin(α) ) * Math.cos(α)
                         + ( H_jj*Math.cos(α) - H_ij*Math.sin(α) ) * Math.sin(α) ),
                    α_min, α_max
                  ),
                  c = Math.cos(α),
                  s = Math.sin(α);

            for( let k=i; k < N; k++ ) // ROTATE ROWS IN H
            { const H_i = H[off + N*i+k],
                    H_j = H[off + N*j+k];
              H[off + N*i+k] = s*H_j + c*H_i;
              H[off + N*j+k] = c*H_j - s*H_i;
            }
            for( let k=j+1; k-- > 0; ) // ROTATE COLUMNS IN H
            { const H_i = H[off + N*k+i],
                    H_j = H[off + N*k+j];
              H[off + N*k+i] = s*H_j + c*H_i;
              H[off + N*k+j] = c*H_j - s*H_i;
            }
            for( let k=N; k-- > 0; ) // ROTATE ROWS IN Q
            { const Q_i = Q[off + N*i+k],
                    Q_j = Q[off + N*j+k];
              Q[off + N*i+k] = s*Q_j + c*Q_i;
              Q[off + N*j+k] = c*Q_j - s*Q_i;
            }
            H[off + N*j+i] *= 0.0; // <- Handles NaN. If a value is that small, its digits are likely nonsense (due to cancellation error) so let's set it to zero.
          }
        }
        // END RESOLVE REAL-VALUED 2x2 BLOCKS

        // TRANSPOSE Q BACK
        for( let i=0; i < N; i++ )
        for( let j=0; j < i; j++ ) { const Q_ij = Q[off + N*i+j]; Q[off + N*i+j]=Q[off + N*j+i]; Q[off + N*j+i] = Q_ij; }
      }
      // END SCHUR DECOMPOSING MATRICES
    }// END la._schur_decomp_qrfrancis_INPLACE
  }

/*****************
 * DOCUMENTATION *
 *****************/
  nd.createDocHTML = () => {
    const
      jsdom= require('jsdom'),
      dom  = new jsdom.JSDOM('<!DOCTYPE html>'),
      doc  = dom.window.document,
      meta = doc.createElement('meta'),
      title= doc.createElement('title'),
      ul   = doc.createElement('ul'),
      h1   = doc.createElement('h1');
    h1   .innerHTML = 'nd.js Documentation'
    title.innerHTML = 'nd.js Documentation'
    meta.setAttribute('charset', 'utf-8')
    doc.lang = 'en'
    doc.head.appendChild(meta)
    doc.head.appendChild(title)
    doc.body.style = 'font-family:monospace'
    doc.body.appendChild(h1)
    doc.body.appendChild(ul)

    const ND = {...nd};
    delete ND.Array

    function addSection( key, val )
    {
      const
        h2 = doc.createElement('h2'),
        pre= doc.createElement('pre'),
        li = doc.createElement('li'),
        a  = doc.createElement('a');

      h2.id = key.replace(' ','_').replace('[','').replace(']','')
      a.innerHTML = key
      a.setAttribute('href', '#'+h2.id)
      ul.appendChild(li)
      li.appendChild(a)

      h2 .innerHTML = key
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

    addClassDoc('nd.Array', nd.Array)

    for( const [key,val] of Object.entries(ND) )
      addSection('nd.'+key,val)

    for( const [key,val] of Object.entries(ND.la) )
      if( key != '__doc__' )
        addSection('nd.la.'+key,val)

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
        return obj.__doc__ || '<Doc N/A>';
    }
  }



nd.Complex.__doc__ = `\
A very rudimentary implementation of the Complex number type,
used mainly for the calculation of eigenvalues. Other than that,
a complex dtype is not yet supported.
`

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

WARNING: This method freezes the buffer that's underlying the shape it's given.

Parameters
----------
shape: Int32Array
  The shape of new nd.Array. May be used directly (no protection
  copy may be perfomed). The buffer behind this Int32Array will be frozen.
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
    directly (creating a protection copy of data).
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


nd.asarray.__doc__ = `\
Similar to \`nd.array(content)\` except that input is not copied if
it already is an \`nd.Array\`.

Parameters
----------
content: { length: int[], data: dtype[] } or dtype or dtype[] or dtype[][] or ...
  The content of the returned nd.Array.

Returns
-------
ndarray: nd.Array
  \`ndarray = content instanceof nd.Array ? content : nd.array(content);\`
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



nd.arrayFromB64.__doc__ = `\
Creates an array reading Base64 encoded data. This currently only works on
big endian machines.

Parameters
----------
dtype: String
  The data type of this nd.Array. Specifies which type of values is read.
  See nd.dtypes for the available data types and their respective (primitive)
  array types.
shape: int[]
  The shape of the array. Must match the content length.
content: String
  The flattened nd.Array data encoded in Base64.

Returns
-------
A: nd.Array
  The array initialized with the given Base64 data as content.
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
  The multi index of the entry that is to be changed.
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



nd.Array.prototype.modify.__doc__ = `\
Modifies an array entriy by applying the given functions to it.
This method call is equivalent to:
\`ndarray.set( indices, modifier(ndarray(...indices), ...indices) )\`

Parameters
----------
indices: int[]
  The multi index of the entry that is to modified.
modifier: (dtype, ...int) => dtype
  The function that is applied to the nd.Array entry.
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



nd.Array.prototype.elems.__doc__ = `\
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



nd.Array.prototype.forElems.__doc__ = `\
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



nd.Array.prototype.mapElems.__doc__ = `\
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
>>> let b = a.mapElems( x => x*x )
>>> console.log( b.toString() )
  [[1, 4],
   [9,16]]
`



nd.Array.prototype.transpose.__doc__ = `\
Reorders the axes of this nd.Array. By default the last two axes are
swapped.

Parameters
----------
axes: ...int
  The order in which the axes of the this nd.Array appear in the
  transposed array. If indices are missing from \`axes\`, they are
  appended to \`axes\` in order.

Returns
-------
transposed: nd.Array
  A transposed copy A of this nd.Array, where:
  \`A[i[0], i[1], ...]  =  this[i[axes[0]], i[axes[1], ...]\`

Examples
--------
>>> let a = nd.array([[1,2,3]])
>>> console.log( a.transpose() )
  [[1],
   [2],
   [3]]

>>> let a = nd.array([ [[1,2,3], [4,5,6]], [[7,8,9], [10,11,12]] ])
>>> console.log( a.transpose(2).toString() )
  [[[ 1,  4],
    [ 7, 10]],
   [[ 2,  5],
    [ 8, 11]],
   [[ 3,  6],
    [ 9, 12]]]
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



nd.Array.prototype.reduceElems.__doc__ = `\
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



nd.Array.prototype.sliceElems.__doc__ = `\
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
>>> console.log( a.sliceElems('...', -2).toString() )
  [13, 23, 33]

>>> console.log( a.sliceElems([,,2]).toString() )
  [[11,12,13,14],
   [31,32,33,34]]

>>> console.log( a.sliceElems(2,[1,3,]).toString() )
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



nd.help.__doc__ = `\
Outputs a documentation string for the given method ND.JS method to the console. This method
is intended to be used in interactive mode to explore and learn the ND.JS API.

Parameters
----------
fun: Function
  The function or package ND.JS API method for which the documentation is to be retrieved.

Examples
--------
>>> nd.help(nd.la)
  The Linear Algebra subpackage of ND.JS. Contains methods for the handling and decomposition
  of matrices as well as solvers for linear equations and linear least square systems.

  As convention the last two axes of an \`nd.Array\` are considered as the matrix dimensions.
  Higher dimensional arrays are considered "arrays of matrices" and the broadcasting rules
  apply for the leading dimenions, i.e. for all dimensions of the shape except the last two.
`



nd.help_string.__doc__ = `\
Returns a documentation string for the given method ND.JS method to the console. This method
is used by \`nd.createDocHTML()\` to create the documentation and by \`nd.help()\` to
print help to the console.

Parameters
----------
fun: Function
  The function or package ND.JS API method for which the documentation is to be retrieved.

Returns
-------
help: String
  A documentation String for the given method or package.

Examples
--------
>>> console.log( nd.help(nd.la) )
  "The Linear Algebra subpackage of ND.JS. Contains methods for the handling and decomposition
   of matrices as well as solvers for linear equations and linear least square systems.

   As convention the last two axes of an \`nd.Array\` are considered as the matrix dimensions.
   Higher dimensional arrays are considered "arrays of matrices" and the broadcasting rules
   apply for the leading dimenions, i.e. for all dimensions of the shape except the last two."
`



nd.createDocHTML.__doc__ = `\
Created a JSDom representation of the documentation HTML file.

Returns
-------
doc: JSDOM
  The HTML document containing the HTML documentation file.
`



  //
 // LINEAR ALGEBRA
//
nd.la.__doc__ = `\
The Linear Algebra subpackage of ND.JS. Contains methods for the handling and decomposition
of matrices as well as solvers for linear equations and linear least square systems.

As convention the last two axes of an \`nd.Array\` are considered as the matrix dimensions.
Higher dimensional arrays are considered "arrays of matrices" and the broadcasting rules
apply for the leading dimenions, i.e. for all dimensions of the shape except the last two.
`



nd.la.eye.__doc__ = `\
Returns an nd.Array of identity matrices, where the last two axes are the matrix dimensions.

Parameters
----------
shape: ...int
  The shape of the resulting array. If only one value N is given, an N*N identity
  matrix is returned.

Returns
-------
I: nd.Array
  The array of identity matrices, where the last two axes are the matrix dimensions.

Examples
--------
>>> const I = nd.la.eye(2);
>>> console.log( I.toString() );
  [[1,0],
   [0,1]]

>>> const J = nd.la.eye(2,3,4);
>>> console.log( J.toString() );
  [[[ 1, 0, 0, 0],
    [ 0, 1, 0, 0],
    [ 0, 0, 1, 0]],
   [[ 1, 0, 0, 0],
    [ 0, 1, 0, 0],
    [ 0, 0, 1, 0]]]
`



nd.la.tril.__doc__ = `\
Returns a copy of an array with all entries outside of the lower triangular part set to 0.
The last two axes are considered the matrix dimensions.

Parameters
----------
A: nd.Array[...,N,M]
  The array of matrices whose lower triangular parts are to be returned.
offset: int
  [Optional] The offset from the main diagonal of entries included in the  result.
  If \`offset = +1\` the entries above/right of the main diagonal are not set to zero.
  If \`offset = -2\` the entries on the main diagonal and one below/left of if are set to zero.
  The default value is 0.

Returns
-------
L: nd.Array[...,N,M]
  The lower triangular part of A, where \`L(i,j) == i >= j-offset ? A(i,j) : 0\`.

Examples
--------
>>> const A = nd.tabulate([2,3,3], (k,i,j) => 100*(k+1) + 10*(i+1) + (j+i) );
>>> console.log( A.toString() );
  [[[ 110, 111, 112 ],
    [ 121, 122, 123 ],
    [ 132, 133, 134 ]],
   [[ 210, 211, 212 ],
    [ 221, 222, 223 ],
    [ 232, 233, 234 ]]]

>>> const L1 = nd.la.tril(A);
>>> console.log( L1.toString() );
  [[[ 110,   0,   0 ],
    [ 121, 122,   0 ],
    [ 132, 133, 134 ]],
   [[ 210,   0,   0 ],
    [ 221, 222,   0 ],
    [ 232, 233, 234 ]]]
>>> const L2 = nd.la.tril(A,-1);
>>> console.log( L2.toString() );
  [[[   0,   0,  0 ],
    [ 121,   0,  0 ],
    [ 132, 133,  0 ]],
   [[   0,   0,  0 ],
    [ 221,   0,  0 ],
    [ 232, 233,  0 ]]]
`



nd.la.triu.__doc__ = `\
Returns a copy of an array with all entries outside of the upper triangular part set to 0.
The last two axes are considered the matrix dimensions.

Parameters
----------
A: nd.Array[...,N,M]
  The array of matrices whose upper triangular parts are to be returned.
offset: int
  [Optional] The offset from the main diagonal of entries included in the result.
  If \`offset = +1\` the main diagonal is set to zero.
  If \`offset = -2\` the entries on the main diagonal and one below/left of if are not set to zero.
  The default value is 0.

Returns
-------
U: nd.Array[...,N,M]
  The upper triangular part of A, where \`L(i,j) == i <= j-offset ? A(i,j) : 0\`.

Examples
--------
>>> const A = nd.tabulate([2,3,3], (k,i,j) => 100*(k+1) + 10*(i+1) + (j+i) );
>>> console.log( A.toString() );
  [[[ 110, 111, 112 ],
    [ 121, 122, 123 ],
    [ 132, 133, 134 ]],
   [[ 210, 211, 212 ],
    [ 221, 222, 223 ],
    [ 232, 233, 234 ]]]

>>> const L1 = nd.la.triu(A);
>>> console.log( L1.toString() );
[[[ 110, 111, 112 ],
  [   0, 122, 123 ],
  [   0,   0, 134 ]],
 [[ 210, 211, 212 ],
  [   0, 222, 223 ],
  [   0,   0, 234 ]]]
>>> const L2 = nd.la.triu(A,+1);
>>> console.log( L2.toString() );
[[[ 0, 111, 112 ],
  [ 0,   0, 123 ],
  [ 0,   0,   0 ]],
 [[ 0, 211, 212 ],
  [ 0,   0, 223 ],
  [ 0,   0,   0 ]]]
`



nd.la.diag_mat.__doc__ = `\
Creates a (main) diagonal matrix with the given diagonal values.

Parameters
----------
diag: nd.Array[...,N]
  The diagonal values.

Returns
-------
D: nd.Array[...,N,N]
  Such that \`D(i,j) == i==j ? diag(i) : 0\`.

Examples
--------
>>> const D = nd.la.diag_mat([[1,2,3],[4,5,6]]);
>>> console.log( D.toString() );
  [[[ 1, 0, 0 ],
    [ 0, 2, 0 ],
    [ 0, 0, 3 ]],
   [[ 4, 0, 0 ],
    [ 0, 5, 0 ],
    [ 0, 0, 6 ]]]
`



nd.la.matmul.__doc__ = `\
Computes the matrix product of a series of matrices. The order
of multiplications is optimized to minimize the number of floating
point operations required. The last two axes are considered to be
the matrix dimensions. For the leading axes, broadcasting rules apply.
Along with \`nd.la.matmul2\`, this is one of the only two methods
that currently allows complex inputs.

Parameters
----------
matrices: ...nd.Array
  The matrices that are to be multiplied.

Returns
-------
product: nd.Array
  Where \`product = matrices[0] @ matrices[1] @ ... @ matrices[n-1]\`.

Examples
--------
>>> const v = nd.tabulate([1e6,1], () => 1);
>>> console.log( v.toString() );
  [[ 1 ],
   [ 1 ],
   [ 1 ],
   [ 1 ],
   [ 1 ],
   ...999990 more...,
   [ 1 ],
   [ 1 ],
   [ 1 ],
   [ 1 ],
   [ 1 ]]
>>> const u = nd.la.matmul(v,v.T,v,v.T,v);
>>> console.log( u.toString() );
  [[1000000000000],
   [1000000000000],
   [1000000000000],
   [1000000000000],
   [1000000000000],
   ...999990 more...,
   [1000000000000],
   [1000000000000],
   [1000000000000],
   [1000000000000],
   [1000000000000]]
`


nd.la.matmul2.__doc__ = `\
Computes the matrix product of exactly two matrices. The last two axes
are considered to be the matrix dimensions. For the leading axes,
broadcasting rules apply. Along with \`nd.la.matmul\`, this is one
of the only two methods that currently allows complex inputs.

Parameters
----------
a: nd.Array[...,N,K]
  The left factor matrices.
b: nd.Array[...,K,M]
  The right factor matrices.

Returns
-------
product: nd.Array[...,N,M]
  Where \`product = a @ b\`.

Examples
--------
>>> const A = [
      [2, 0],
      [0,-1]
    ];
>>> const x = nd.array([
      [[1,2]],
      [[3,4]],
      [[5,6]]
    ]).T;
>>> console.log( x.toString() );
  [ [[1],
     [2]],

    [[3],
     [4]],

    [[5],
     [6]] ]
>>> const y = nd.la.matmul2(A,x);
>>> console.log( y.toString() );
  [ [[ 2 ],
     [-2 ]],
  
    [[ 6 ],
     [-4 ]],
  
    [[10 ],
     [-6 ]] ]
`



nd.la.cholesky_decomp.__doc__ = `\
Computes the Cholesky Decomposition of an \`nd.Array\` of symmetric, positive
definite (square) matrices. The implementation assume the matrices to be symmetric
and only looks at the lower triangular values. The last two axes of the \`nd.Array\`
are considered to be the matrix dimensions.

Parameters
----------
S: nd.Array[...,N,N]
  An \`nd.Array\` of symmetric, positive definite (square) matrices. In other words
  S is symmetric and has only positive eigenvalues. 

Returns
-------
L: nd.Array[...,N,N]
  A lower triangular matrix, such that \`S = L @ L.T\`.

Examples
--------
>>> const S = [
      [ 25, -50],
      [-50, 101]
    ];
>>> const L = nd.la.cholesky_decomp(S);
>>> console.log( nd.la.matmul2(L,L.T).toString() );
  [[ 25, -50],
   [-50, 101]]
`



nd.la.cholesky_solve.__doc__ = `\
Given a cholesky decomposition and the right hand side of a linear equation system,
this method computes the result of said system.

Parameters
----------
L: nd.Array[...,N,N]
  An array of lower triangular (square) matrices. L is assumed to be lower triangular
  without looking at the upper triangular values.
y: nd.Array[...,N,M]
  The right hand side of the linear equations system.

Returns
-------
x: nd.Array[...,N,M]
  Such that \`L @ L.T @ x == y\`.

Examples
--------
>>> const S = [
      [ 25, -50],
      [-50, 101]
    ];
>>> const y = [
      [[1],
       [2]],
      [[3],
       [4]]
    ];
>>> const L = nd.la.cholesky_decomp(S);
>>> const x = nd.la.cholesky_solve(L,y);
>>> console.log( nd.la.matmul2(S,x) );
 [[[1.0],
   [2.0]],
  [[3.0],
   [4.0]]]
`



nd.la.tril_solve.__doc__ = `\
Given the lower triangular (square) matrix and the right hand side of a linear
equation system, this method computes the result.

Parameters
----------
L: nd.Array[...,N,N]
  The lower triangular matrix of the linear equations system. L is assumed to
  be lower triangular without ever looking at the values in the upper triangular
  region.
y: nd.Array[...,N,M]
  The right hand side matrix of the linear equation system.

Returns
-------
x: nd.Array[...,N,M]
  Such that \`L @ x = y\`.
`



nd.la.triu_solve.__doc__ = `\
Given the upper triangular (square) matrix and the right hand side of a linear
equation system, this method computes the result.

Parameters
----------
U: nd.Array[...,N,N]
  The upper triangular matrix of the linear equations system. U is assumed to
  be upper triangular without ever looking at the values in the lower triangular
  region.
y: nd.Array[...,N,M]
  The right hand side matrix of the linear equation system.

Returns
-------
x: nd.Array[...,N,M]
  Such that \`U @ x = y\`.
`



nd.la.lu_decomp.__doc__ = `\
Given an \`nd.Array\` of square matrices, this method Computes the
LU(P) decomposition with Column Pivotization.

Parameters
----------
A: nd.Array[...,N,N]

Returns
-------
LU: nd.Array[...,N,N]
  A matrix containing both the lower and upper triangular part of the
  LU decomposition. The diagonal of ther lower triangular matrix L
  is not contained in LU. It only containes ones and is implied.
P : nd.Array[...,N]
  The order in which the row indices of A appear in the LU decomposition, i.e.:
  \`A(P[i],j) == (L @ U)[i,j]\`.
`



nd.la.lu_solve.__doc__ = `\
Given the LU(P) decomposition and the right hand side of a
Linear Equations System, this method computes the result of
said system.

Parameters
----------
LU: nd.Array[...,N,N]
P : nd.Array[...,N]
y : nd.Array[...,N,M]

Returns
-------
x : nd.Array[...,N,M]
  The solution of the Linear Equation System, such that:
  \`y[P[i],:] == (L @ U @ x)[i,:] \`
`



nd.la.qr_decomp_full.__doc__ = `\
Computes the (full) QR Decomposition of a matrix. The QR
Decomposition can be used to solve both Linear Equations
and Linear Least Square problems with high numeric accuracy.
Under normal circumstances, the incomplete QR Decomposition
(\`nd.la.qr_decomp\`) is to be preferred over this method
as it may be significantly more memory efficient.

Parameters
----------
A: nd.Array[...,N,M]
  The matrix for which the QR Decomposition is computed.

Returns
-------
Q: nd.Array[...,N,N]
  An orthogonal square matrix, i.e. \`Q @ Q.T == Q.T @ Q == nd.la.eye(N)\`.
R: nd.Array[...,N,M]
  An upper triangular matrix, such that: \`A = Q @ R\`
`



nd.la.hessenberg_decomp.__doc__ = `\
Computes the (Upper) Hessenberg Decomposition of a matrix. It
is worth mentioning that the Hessenberg Decomposition of a
symmetric matrix is tridiagonal.

Parameters
----------
A: nd.Array[...,N,N]
  The matrix for which the Hessemberg Decomposition is computed.

Returns
-------
U: nd.Array[...,N,N]
  An orthogonal matrix, i.e. \`U @ U.T == U.T @ U == nd.la.eye(N)\`.
H: nd.Array[...,N,N]
  An upper Hessemberg matrix, i.e. \`nd.la.tril(H,-2) == 0\`, such that \`U @ H @ U.T == A\`.
`


nd.la.eigen.__doc__ = `\
Returns the eigenpairs of a real non-symmetric square matrix.

Parameters
----------
A: nd.Array[...,N,N]
  The non-symmetric real squre matrix for which the eigenvalues are computed.

Returns
-------
Λ: nd.Array[...,N]
  A matrix containing all eigenvalues.
V: nd.Array[...,N,N]
  A matrix containing the eigenvectors corresponding to Λ as columns, i.e.:
  \`Λ[i]*V[j,i] == (A @ V)[j,i]\`.
  The columns are normalized using the 2-norm. 
`



nd.la.schur_decomp.__doc__ = `\
Computes the (real) Schur Decomposition of a matrix. The
Schur Decomposition is used to compute the Eigenvalues and
Eigenvectors.

Parameters
----------
A: nd.Array[...,N,N]

Returns
-------
Q: nd.Array[...,N,N]
  An orthogonal square matrix.
T: nd.Array[...,N,N]
  An upper quasi-triangular matrix, such that \`Q @ T @ Q.T == A\`.
  The diagonal of T consists of 2x2 blocks that represent complex
  eigenpairs and 1x1 matrix that represent a real eigenvalue.
`



nd.la.qr_decomp.__doc__ = `\
Computes the (economic) QR Decomposition of a matrix. The
QR Decomposition can be used to solve both Linear Equations
and Linear Least Square problems with high numeric accuracy.
Under normal circumstances, this method is to
be preferred over the full QR Decomposition (\`nd.la.qr_decomp_full\`)
as it may be significantly more memory efficient.

Parameters
----------
A: nd.Array[...,N,M]
  The matrix for which the QR Decomposition is computed.

Returns
-------
Q: nd.Array[...,N,min(N,M)]
  An orthogonal rectangular matrix, i.e. \`Q.T @ Q == nd.la.eye(min(N,M))\`.
R: nd.Array[...,min(N,M),M]
  An upper triangular matrix, such that: \`A = Q @ R\`
`



nd.la.rrqr_decomp.__doc__ = `\
Computes the economic Rank-Revealing QR Decomposition of a matrix.

Parameters
----------
A: nd.Array[...,N,M]
  The matrix for which the QR Decomposition is computed.

Returns
-------
Q: nd.Array[...,N,min(N,M)]
  An orthogonal rectangular matrix, i.e. \`Q.T @ Q == nd.la.eye(min(N,M))\`.
R: nd.Array[...,min(N,M),M]
  An upper triangular matrix, where \`R[i,i] >= R[j,j]\` if and only if \`i <= j\`.
P: nd.Array[...,M]
  The permuted column indices, such that:
  \`(Q @ R)[:,j] == A[:,P[j]]\`
`



nd.la.qr_solve.__doc__ = `\
Given the QR Decompostion and the right hand side of a Linear
Equation System, this method solves it. If the Linear Equation
System is overdetermined, the Linear Least Square solution is
computed.

Parameters
----------
Q: nd.Array[...,N,K]
  The orthognal rectangular matrix of the QR Decomposition.
R: nd.Array[...,K,M]
  The upper triangular matrix of the QR Decomposition.
y: nd.Array[...,N,L]
  The right hand side of the Linear Equation System.

Returns
-------
x: nd.Array[...,N,L]
  Such that \`‖(Q @ R @ x) - y‖₂\` is minimal.
`



nd.la.bidiag_decomp.__doc__ = `\
Computes the Upper Bidiagonal Decomposition of a matrix. Bidiagonal
Decomposition is a common preprocessing step for the Singular Value
Decomposition.

Parameters
----------
A: nd.Array[...,N,M]
  The matrix for which the Upper Bidiagonal Decomposition is
  computed.

Returns
-------
U: nd.Array[..., N, min(N,M)]
  An orthognal rectangular matrix.
B: nd.Array[..., min(N,M), N >= M ? M : N+1 ]
  An upper bidiagonal rectangular matrix.
V: nd.Array[..., N >= M ? M : N+1, M]
  An orthogonal rectangular matrix. Such that \`A = U @ B @ V @\`
`



nd.la.solve.__doc__ = `\
Solves a Linear Equations System (LES) or a Linear Least Squares Problem (LSTSQ).

Parameters
----------
A: nd.Array[...,N,K]
  The matrix of the LES/LSTSQ.
y: nd.Array[...,K,M]
  The right hand side of the LES/LSTSQ.

Returns
-------
x: nd.Array[...,N,M]
  The solution of the LES/LSTSQ such that \`‖(A @ x) - y‖₂\` is minimal.
`



nd.la.svd_decomp.__doc__ = `\
Computes the Singular Value Decomposition (SVD) of a matrix.

Parameters
----------
A: nd.Array[...,N,M]
  The matrix for which the SVD is computed.

Returns
-------
U : nd.Array[...,N,min(N,M)]
  An orthogonal, rectangular matrix.
sv: nd.Array[..., min(N,M) ]
  The singular values of A sorted in descending order.
V : nd.Array[...,min(N,M),M]
  An orthogonal, rectangular matrix, such that \`A == U @ diag(sv) @ V\`.
`
}
