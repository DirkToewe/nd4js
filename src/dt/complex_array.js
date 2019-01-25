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

import {Complex} from './complex';
import 'util';

function createComplexArrayType( FloatArray )
{
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
  };

  class ComplexArray
  {
    static get BYTES_PER_ELEMENT() {
      return FloatArray.BYTES_PER_ELEMENT*2;
    }
  
    static get name() {
      return `Complex${ComplexArray.BYTES_PER_ELEMENT*8}Array`;
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

      // fast-copy complex arrays
      if( null == mapFn )
      {
        if( source instanceof ComplexArray )
          return new ComplexArray( FloatArray.from(source._array).buffer, 0, source.length );
        if( source instanceof   Int32Array ||
            source instanceof Float32Array ||
            source instanceof Float64Array )
        {
          const array = new FloatArray(source.length*2);
          for( let i=source.length; i-- > 0; )
            array[i<<1] = source[i];
          return new ComplexArray(array.buffer, 0, source.length);
        }
      }

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
      if( ! (value instanceof Complex) ) value = new Complex(value);
  
      if( null == start ) start = 0;
      if( null == end   ) end   = this.length;
      if( 0 > end   ) end   += this.length;
      if( 0 > start ) start += this.length;

      for( let i=this._array.length; (i-=2) >= 0; ) {
        this._array[i+1] = value.im;
        this._array[i+0] = value.re;
      }
  
      return this;
    }
  
    get( index ) {
      return new Complex(
        this._array[2*index+0],
        this._array[2*index+1]
      )
    }
  
    set( index, value ) {
      index *= 2;

      if( typeof(value) === 'number' || ! isNaN(value) )
      {
        this._array[index+0] = value;
        this._array[index+1] = isNaN(value) ? NaN : 0;
        return;
      }

      if( ! ('re' in value) ||
          ! ('im' in value) )
        value = new Complex(value);

      this._array[index+0] = value.re;
      this._array[index+1] = value.im;
    }

    get [Symbol.toStringTag]() {
      throw new Error('not yet implemented.');
    }
  
    [Symbol.for('nodejs.util.inspect.custom')]( options ) {
      return 'Complex128Array ' + util.inspect( Array.from(this), options );
    }
  
    toString( max_len ) {
      return this.join(',');
    }
  };

  return ComplexArray;
}

export const Complex128Array = createComplexArrayType(Float64Array)
