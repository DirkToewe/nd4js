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

import {IS_LITTLE_ENDIAN} from '../io'


const val =   Float64Array.of(NaN),
     bits = new Int32Array(val.buffer);


/*DEBUG*/ if( bits.length !== 2 ) throw new Error('Assertion failed.');


export function nextUp(x)
{
  x *= 1;

  if( 0 === x )
    return +Number.MIN_VALUE;
  if( ! (x < +Infinity) ) // <- handles NaN and Infinity
    return x;

  val[0] = x;

  const i = 1-IS_LITTLE_ENDIAN,
        j = 1*IS_LITTLE_ENDIAN;

  if( x > 0 )
  {
    bits[j] += -1 === bits[i];
    bits[i] +=  1;
  }
  else {
    bits[j] -= 0 === bits[i];
    bits[i] -= 1;
  }

  return val[0];
};


export function nextDown(x)
{
  x *= 1;

  if( 0 === x )
    return -Number.MIN_VALUE;
  if( ! (x > -Infinity) ) // <- handles NaN and Infinity
    return x;

  val[0] = x;

  const i = 1-IS_LITTLE_ENDIAN,
        j = 1*IS_LITTLE_ENDIAN;

  if( x > 0 )
  {
    bits[j] -= 0 === bits[i];
    bits[i] -= 1;
  }
  else {
    bits[j] += -1 === bits[i];
    bits[i] +=  1;
  }

  return val[0];
}


export function midl( x, y )
{
  x *= 1;
  y *= 1;
  if( isNaN(x) ) throw new Error('mid(x,y): x must be number.');
  if( isNaN(y) ) throw new Error('mid(x,y): x must be number.');

  if( Math.sign(x)*y < 0 ) // <- check for opposite signs
  {
    const mid = (x+y) / 2; // <- avoids underflow if e.g x=+MAX_VALUE, y=-MAX_VALUE

    if( ! (mid < Math.max(x,y)) ) throw new Error('Assertion failed.');
    if( ! (mid > Math.min(x,y)) ) throw new Error('Assertion failed.');

    return mid;
  }
  else
  {
    const [a,b] = Math.abs(x) <= Math.abs(y) ? [x,y] : [y,x];

    // Returns the mid point between two floats (x,y) or x if there is no mid point
    const mid = a + (b-a)*0.5; // <- avoids underflow if e.g x=y=MAX_VALUE

    if( ! (mid <= Math.max(x,y)) ) throw new Error('Assertion failed.');
    if( ! (mid >= Math.min(x,y)) ) throw new Error('Assertion failed.');

    if( mid === y ) return x;

    return mid;
  }
}
