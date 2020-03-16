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

export function _min1d_interp_quad( x1, x2, f1, f2, g1 )
{
  x1*=1; if( ! isFinite(x1) ) throw new Error('Assertion failed.');
  x2*=1; if( ! isFinite(x2) ) throw new Error('Assertion failed.');
  f1*=1; if( ! isFinite(f1) ) throw new Error('Assertion failed.');
  f2*=1; if( ! isFinite(f2) ) throw new Error('Assertion failed.');
  g1*=1; if( ! isFinite(g1) ) throw new Error('Assertion failed.');

  if( x1===x2 ) throw new Error('Assertion failed.');

  // polynomial: p(x) = a + b*(x - xMin)²
  //        (A)  p(x1) = f1 (= a + b*(x1 - xMin)²)
  //        (B)  p(x2) = f2 (= a + b*(x2 - xMin)²)
  //        (C) p'(x1) = g1 (=   2*b*(x1 - xMin) )
  //
  // (B)-(A) => f2-f1 = b*( (x2-xMin)² - (x1-xMin)² )
  //     (C) => b = g1 / (2*(x1-xMin))
  //
  // => (f2-f1)/g1 * 2 * (x1-xMin) = (x2-xMin)² - (x1-xMin)²
  //   = x2² - x1² - 2*(x2-x1)*xMin
  //
  //            x1*(f2-f1)*2 / g1 - (x2²-x1²)       x1*(f2-f1) - (x2-x1)*(x2+x1) * g1/2
  // => xMin = ───────────────────────────────── = ─────────────────────────────────────
  //               (f2-f1)*2 / g1 - (x2 -x1 )*2        (f2-f1) - (x2-x1)         * g1

//  const df = f2-f1,
//       gdx =(x2-x1)*g1;

//  const xMin = (x1*df - gdx*(x2+x1)/2)
//             / (   df - gdx          ); // <- TODO investigate numeric precision

//  const xMin = x1 - 0.5 * (x2-x1)*g1 / ( (f2-f1) / (x2-x1) - g1 );

  const dFdX = (f2-f1) / (x2-x1);

//  const r = - g1 / (dFdX - g1) * 0.5;
//  const xMin = x1*(1-r) + x2*r;

  const w2 = -g1/2,
        w1 = dFdX + w2;
  const xMin = (x1*w1 + x2*w2) / (dFdX - g1);

/*DEBUG*/  const sign = Math.sign(g1) * Math.sign(x1-xMin);
/*DEBUG*/  if( !(0 <= sign) ) throw new Error('Assertion failed.');

  return xMin;
}

export function* _heap_sort( items, isLess = (x,y) => x < y )
{
  if(  items.length%1 !== 0 ) throw new Error('Assertion failed.');
  if(!(items.length   >=  0)) throw new Error('Assertion failed.');
  if(!(isLess instanceof Function)) throw new Error('Assertion failed.');

  let i = 0;
  const len = items.length;

  const less = (i,j) => isLess(items[i], items[j]),
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

  // EXTRACT MINIMA
  while( i < len-1 ) {
    swap(i,  len-1); yield items[i++];
    siftDown(len-1);
  }
}
