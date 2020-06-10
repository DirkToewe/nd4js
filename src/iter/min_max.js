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


export function argmin( seq, compare_fn = (x,y) => (x>y) - (x<y) )
{
  if( !(compare_fn instanceof Function) ) throw new Error('argmin(seq,compare_fn): compare_fn must be a function.');

  if( Symbol.iterator in seq )
    seq = seq[Symbol.iterator]();
  else
    console.warn('argmin(seq,compare_fn): seq[Symbol.iterator] does not exist.');

  let val = function(){
    const first = seq.next();
    if( first.done )
      throw new Error('argmin(seq,compare_fn): seq must be an iterable containing at least one item.');
    return first.value;
  }();

  let idx = 0;

  for( let i=1;; i++ ) {
    const next = seq.next();
    if(   next.done ) return idx;
    if( compare_fn(val,next.value) > 0 ) {
      idx = i;     val=next.value;
    }
  }
}


export function argmax( seq, compare_fn = (x,y) => (x>y) - (x<y) )
{
  if( !(compare_fn instanceof Function) ) throw new Error('argmax(seq,compare_fn): compare_fn must be a function.');

  return argmin( seq, (x,y) => -compare_fn(x,y) );
}


export function min( seq, compare_fn = (x,y) => (x>y) - (x<y) )
{
  if( !(compare_fn instanceof Function) ) throw new Error('argmin(seq,compare_fn): compare_fn must be a function.');

  if( Symbol.iterator in seq )
    seq = seq[Symbol.iterator]();
  else
    console.warn('argmin(seq,compare_fn): seq[Symbol.iterator] does not exist.');

  let val = function(){
    const first = seq.next();
    if( first.done )
      throw new Error('argmin(seq,compare_fn): seq must be an iterable containing at least one item.');
    return first.value;
  }();

  for(;;) {
    const next = seq.next();
    if(   next.done ) return val;
    if( compare_fn(val, next.value) > 0 )
                   val= next.value;
  }
}


export function max( seq, compare_fn = (x,y) => (x>y) - (x<y) )
{
  if( !(compare_fn instanceof Function) ) throw new Error('argmax(seq,compare_fn): compare_fn must be a function.');

  return min( seq, (x,y) => -compare_fn(x,y) );
}
