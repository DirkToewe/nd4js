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

import {NDArray, asarray} from './nd_array'
import math from './math'


function unravel( flat, shape ) {
  const result = Int32Array.from(shape)
  for( let i=result.length; i-- > 0; ) {
    result[i] =       flat % shape[i]
    flat = Math.trunc(flat / shape[i])
  }
  return result
}

const toBeBand = (label, lower, upper) => (util, customEq) => ({
  compare(act, {tol=0}={})
  {
    const LEN = act.data.length,
        [M,N] = act.shape.slice(-2)

    for( let flat=0; flat < LEN; flat++ )
    {
      const i = Math.trunc(flat / N) % M,
            j =            flat % N
  
      if( (i-j > lower  ||
           j-i > upper) && ! (math.abs(act.data[flat]) <= tol) )
        return {
          pass: false,
          message: `Expected\n${act}\nto be ${label}, but actual(${unravel(flat,act.shape)}) = ${act.data[flat]} ≉ 0.`
        }
    }
  
    return { pass: true }
  }
})

function toBeAll(act,exp, predicate, str)
{
  act = asarray(act)
  exp = asarray(exp)
  if( Object.is(act,exp) )
    console.warn('Actual and expected are identical.')

  const ndim = Math.max(
    act.ndim,
    exp.ndim
  )

  const shape = Array.from({length: ndim}, () => 1)
  for( let a = act.ndim,
           e = exp.ndim,
           i =     ndim; i-- > 0; )
  {
    const as = act.shape[--a] || 1,
          es = exp.shape[--e] || 1
         if( as === 1 ) shape[i] = es
    else if( es !== 1 && es !== as )
      return {
        pass: false,
        message: `Expected shape [${act.shape}] to be broadcast-compatible to [${exp.shape}].`
      }
    else shape[i] = as
  }

  let flatAct=0, strideAct,
      flatExp=0, strideExp;

  function visit(d)
  {
    if( ndim === d ) {
      if( ! predicate( act.data[flatAct], exp.data[flatExp] ) )
        throw `A:\n${act}\n` +
              `B:\n${exp}\n` +
              `Expected A to ${str} B, but:\n` +
              `A(${unravel(flatAct,act.shape)}) = ${act.data[flatAct]}\n` +
              `B(${unravel(flatExp,exp.shape)}) = ${exp.data[flatExp]}`;
      flatAct += strideAct = 1
      flatExp += strideExp = 1
      return
    }
    for( let i=0;; )
    {
      visit(d+1)
      if( ++i === shape[d] ) break
      if( ! (act.shape[d - ndim + act.ndim] > 1) ) flatAct -= strideAct
      if( ! (exp.shape[d - ndim + exp.ndim] > 1) ) flatExp -= strideExp
    }
    strideAct *= (act.shape[d - ndim + act.ndim] || 1)
    strideExp *= (exp.shape[d - ndim + exp.ndim] || 1)
  }

  try {
    visit(0)
    return { pass: true }
  }
  catch(message) {
    return { pass: false, message }
  }
}

export const CUSTOM_MATCHERS = {
  toBeAllCloseTo: (util, customEq) => ({
    compare(act, exp, { rtol=1e-5, atol=1e-8 } = {})
    {
      const is_close = (x,y) => {
        const tol = atol + rtol * math.max(
          math.abs(x),
          math.abs(y)
        );
        return math.abs( math.sub(x,y) ) <= tol
      }

      return toBeAll(act,exp, is_close, 'be all close to');
    }
  }),

  toBeAllLessOrClose: (util, customEq) => ({
    compare(act, exp, { rtol=1e-5, atol=1e-8 } = {})
    {
      const clothless = (x,y) => {
        const tol = atol + rtol * math.max(
          math.abs(x),
          math.abs(y)
        );
        return x - y <= tol
      }

      return toBeAll(act,exp, clothless, 'be all less than or close to');
    }
  }),

  toBeCloseTo: (util, customEq) => ({
    compare(act, exp, { rtol=1e-5, atol=1e-8 } = {})
    {
      const tol = atol + rtol * math.max(
        math.abs(act),
        math.abs(exp)
      )
      const result = { pass: math.abs( math.sub(act,exp) ) <= tol }
      if( ! result.pass )
        result.message = `Expected actual=${act} to be close to expected=${exp}.`
      return result
    }
  }),

  toBeDiagonal       : toBeBand(        'diagonal', 0,       0),
  toBeLowerBidiagonal: toBeBand('lower bidiagonal',+1,       0),
  toBeUpperBidiagonal: toBeBand('upper bidiagonal', 0,      +1),
  toBeTridiagonal    : toBeBand(     'tridiagonal',+1,      +1),
  toBeLowerTriangular: toBeBand('lower triangular', Infinity,0),
  toBeUpperTriangular: toBeBand('upper triangular', 0,Infinity),
  toBeLowerHessenberg: toBeBand('lower hessenberg',Infinity,+1),
  toBeUpperHessenberg: toBeBand('upper hessenberg',+1,Infinity)
}

function toStr(item)
{
  if( ARRAY_TYPES.some( Type => item instanceof Type ) )
  {
    let lf = '',
        infix = ', ',
        result
    if( item.length > 9 ) {
      const head = Array.from(item.slice(0,4),toStr),
            tail = Array.from(item.slice( -4),toStr)
      
      if( [...head, ...tail].some(s => /\n/.test(s)) ) {
        infix = '\n,\n'
        lf = '\n'
      }
      result = `${head.join(infix)},${lf || ' '}... ${item.length-8} more...${lf || ' '}, ${tail.join(infix)}`
    }
    else {
      result = Array.from(item,toStr)
      if( result.some(s => /\n/.test(s)) ) {
        infix = '\n,\n'
        lf = '\n'
      }
      result = result.join(infix)
    }
    return item instanceof Array
      ?                            `[${lf}${result}${lf}]`
      : `${item.constructor.name}.of(${lf}${result}${lf})`
  }
  return 'string' === typeof item
    ? `"${item}"`
    :  `${item}`
}

const ARRAY_TYPES = [Array, Int8Array,  Uint8Array, Uint8ClampedArray,   
              Float32Array, Int16Array, Uint16Array, 
              Float64Array, Int32Array, Uint32Array]

export const forEachItemIn = items => ({
  it: (description, specFn, timeout) => {
    if( 'string' !== typeof description )
      throw new Error('description must be string.')

    let item, i;

    const spec = it(description, () => new Promise( (resolve,reject) => {
      const iter = items[Symbol.iterator]()
      i = -1

      const intervalID = setInterval( () => {
        try {
          const t0 = performance.now()

          for(;;)
          {
            const { value, done } = iter.next()
            ++i

            if(done) {
              if( 0 === i )
                throw new Error('forEachItemIn: items may not be empty.')
//              console.log(`${i} tests ↴`)
              clearInterval(intervalID)
              return resolve()
            }

            specFn(item=value)

            if( 250 < performance.now() - t0 )
              return
          }
        }
        catch(err) {
          clearInterval(intervalID)
          return reject(err)
        }
      }, 0)
    }), timeout)


    if( spec !== undefined )
      spec.addExpectationResult = (passed, data, isError) => {
        if( passed )
          return
        try {
          if( !passed )
          {
            function msg( msg )
            {
              let str = toStr(item)
              if( /\n/.test(str) )
                str = '\n' + str
              return `For item[${i}] = ${str}:\n${msg}`
            }

            if( 'message' in data ) data      .message = msg(data      .message)
            else                    data.error.message = msg(data.error.message)
          }
        }
        catch(err) {
          console.error(err)
        }
        finally {
          Object.getPrototypeOf(spec).addExpectationResult.call(spec, passed, data, isError)
        }
      };


    return spec;
  }
});
