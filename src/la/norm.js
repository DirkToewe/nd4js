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

export class FrobeniusNorm
{
  constructor() {
    this.sum = 0.0;
    this.max = 0.0;
    Object.seal(this);
  }

  reset() {
    this.sum = this.max = 0;
  }

  include( x ) {
    x = Math.abs(x);
    if( x !== 0 ) {
      if(         this.max < x ) {
        const s = this.max / x; this.sum *= s*s;
                  this.max = x;
      }
      x /= this.max;
      this.sum += x*x;
    }
  }

  resultIncl( x ) {
    x = Math.abs(x);
    let {sum,max} = this;
    if( x !== 0 ) {
      if(         max < x ) {
        const s = max / x; sum *= s*s;
                  max = x;
      }      x /= max;
      sum += x*x;
    }
    return isFinite(max) ? Math.sqrt(sum)*max : max;
  }

  get result() {
    return isFinite(this.max) ? Math.sqrt(this.sum)*this.max : this.max;
  }
}

export function norm( A, ord='fro', axis=undefined )
{
  let norm;
  switch(ord) {
    case 'fro': norm = new FrobeniusNorm(); break;
    default   : throw new Error(`norm(A,ord,axis): Unsupported ord: ${ord}.`);
  }
  if( null == axis ) {
    for( const x of A.data )
      norm.include(x);
    return norm.result;
  }
  throw new Error('norm(A,ord,axis): axis argument not yet supported.');
}

//    if( null == reducer )
//    {
//      if( null    ==  dtype      ) return this.data.reduce(axes)
//      if('string' === typeof axes) return this.data.reduce(dtype)
//      reducer = dtype; dtype = undefined
//    }
//    if( null == dtype ) dtype = 'object'
//
//    if( ! is_subdtype(this.dtype, dtype) )
//      throw new Error('New dtype must be a super-dtype.')
//
//    const oldNDim = this.ndim
//
//    if( 'number' === typeof axes ) axes = [axes]
//    if( axes instanceof NDArray ) {
//      if( ! is_subdtype(axes.dtype,'int32') ) throw new Error(`Invalid dtype ${axes.dtype} for axes.`) 
//      if( axes.ndim === 1 )
//        axes = axes.data
//      else
//        throw new Error('Only 1D nd.Array allowed for axes.')
//    }
//    axes = new Set( function*() {
//      for( let ax of axes ) {
//        if( 0 > ax )  ax += oldNDim
//        if( 0 > ax || ax >= oldNDim ) throw new Error('Reduction axis '+ax+' out of bounds.')
//        yield +ax
//      }
//    }() )
//
//    const
//      oldShape= this.shape,
//      newShape= oldShape.filter( (size,i) => ! axes.has(i) ),
//      newData = new ARRAY_TYPES[dtype]( newShape.reduce((a,b) => a*b, 1) ),
//      oldData = this.data
//    let
//      newIdx = 0,
//      oldIdx = 0
//
//    function fill(d, reduce)
//    {
//      if( oldShape.length === d )
//      {
//        if(reduce) newData[newIdx]=reducer(oldData[oldIdx++], newData[newIdx])
//        else       newData[newIdx]=        oldData[oldIdx++]
//        ++newIdx
//      }
//      else if( axes.has(d) ) {
//        const idx = newIdx
//        for( let j=oldShape[d]; j-- > 0; reduce=true ) // <- copy over very first value and only then start reduction
//        {
//          newIdx = idx
//          fill(d+1, reduce)
//        }
//      }
//      else for( let j=oldShape[d]; j-- > 0; )
//        fill(d+1, reduce)
//    }
//    fill(0,false)
//
//    if( newIdx !=   newData.length ) throw new Error([newIdx,   newData.length])
//    if( oldIdx != this.data.length ) throw new Error([oldIdx, this.data.length])
//
//    return new NDArray(newShape,newData)
