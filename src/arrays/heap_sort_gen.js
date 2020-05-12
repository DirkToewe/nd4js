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


export function* heap_sort_gen( items, compare_fn = (x,y) => (x>y)-(x<y) ) // <- TODO: move to src/_array_utils.js
{
  if(  items.length%1 !== 0 ) throw new Error('Assertion failed.');
  if(!(items.length   >=  0)) throw new Error('Assertion failed.');
  if(!(compare_fn instanceof Function)) throw new Error('Assertion failed.');

  let i = 0;
  const len = items.length;

  const less = (i,j) => compare_fn(items[i], items[j]) < 0,
        swap = (i,j) => {
          const x = items[i];
                    items[i] = items[j];
                               items[j] = x;
        };

  const parent = child  => len-1 - (len-2-child  >> 1),
        child  = parent => len-2 - (len-1-parent << 1);

  const siftDown = p => {
    let c = child(p);
    if( c > i ) {
      if( ! less(c,c-1) ) --c;
      if(   less(c,p  ) ) { swap(c,p); siftDown(c); }
    }
    else if( c===i && less(c,p) ) swap(c,p);
  };

  // HEAPIFY
  for( let j = parent(0); j < len; j++ )
    siftDown(j);

  // EXTRACT MINIMA FROM HEAP
  while( i < len ) {
    swap(i,  len-1); yield items[i++];
    siftDown(len-1); // <- reinstate heap property
  }
}