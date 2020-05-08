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


const ARITY = 16;


const child = parent => ARITY*parent + 1;


const parent = child => (child-1) / ARITY | 0;


export class NAryHeap
{
  constructor()
  {
    this._heap = [];
  }

  get size() {
    return this._heap.length;
  }

  add( item )
  {
    if( ! ('key' in item) )
      throw new Error('Assertion failed.');

    const heap = this._heap;
    heap.push(item);

    // SIFT UP
    let i = heap.length-1;
    for(;;)
    {
      const p = parent(i);
      if( i===0 || item.key >= heap[p].key )
        break;
      heap[i] = heap[p]; i=p;
    }
    heap[i] = item;
  }

  get min() {
    return this._heap[0];
  }

  popMin()
  {
    const heap = this._heap,
        result = heap[0],
        filler = heap.pop(); // <- fill the root with the rightmost entry (then we sift it down again)

    // SIFT DOWN
    if( 0 < heap.length )
    {
      for( let i=0;; )
      {
        const p = i;
        heap[i] = filler;
        let c = child(i);

        // find largest child that is larger than the filler
        const      C = Math.min(c+ARITY, heap.length);
        for( ; c < C; c++ )
          if( heap[c].key < heap[i].key )
            i = c;

        if( i === p )
          break;

        heap[p] = heap[i]; // <- move smallest value to root/parent
      }
    }

    return result;
  }
}
