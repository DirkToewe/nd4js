'use strict';

/* This file is part of ND4JS.
 *
 * ND4JS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ND4JS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ND4JS. If not, see <http://www.gnu.org/licenses/>.
 */

import {is_array} from "./is_array";


export class IndexOutOfBoundsError extends Error
{}

const ARRAY_BOUNDS_CHECKER = {
  get( array, key, proxy ) {
    if( typeof(key) !== 'symbol' && key%1 === 0 ) {
      key |= 0;
      if( !(0 <= key && key < array.length) ) throw new IndexOutOfBoundsError();
    }
    let val = array[key];
    if( val instanceof Function )
      return function( ...args ) {
        return val.apply(this===proxy ? array : this, args);
      };
    return val;
  },
  set( array, key, val ) {
    if( typeof(key) !== 'symbol' && key%1 === 0 ) {
      key |= 0;
      if( !(0 <= key && key < array.length) ) throw new IndexOutOfBoundsError();
    }
    array[key] = val;
    return true;
  }
};

export function checked_array(arr)
{
  if( ! is_array(arr) )
    throw new Error('Assertion failed.');
  return new Proxy(arr,ARRAY_BOUNDS_CHECKER);
}
