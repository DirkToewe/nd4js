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

import {midl} from "../dt/float64_utils";


// REFERENCES
// ----------
// .. [1] https://en.wikipedia.org/wiki/Bisection_method


export function root1d_bisect(F, x1,x2)
{
  x1 *= 1;
  x2 *= 1;
  if( isNaN(x1) ) throw new Error('root1d_bisect(F,x1,x2): x1 must be a number.');
  if( isNaN(x2) ) throw new Error('root1d_bisect(F,x1,x2): x2 must be a number.');

  if( ! (F instanceof Function) ) throw new Error('root1d_bisect(F,x1,x2): F must be a function.');

  let f1 = F(x1),
      f2 = F(x2);

  if( 0 === f1 ) return x1;
  if( 0 === f2 ) return x2;

  if( f1 > f2 ) {
    [f1,f2] = [f2,f1];
    [x1,x2] = [x2,x1];
  }

  if( ! (f1 < 0) ) throw new Error('root1d_bisect(F,x1,x2): F(x1) and F(x2) must have opposite signs.');
  if( ! (f2 > 0) ) throw new Error('root1d_bisect(F,x1,x2): F(x1) and F(x2) must have opposite signs.');
  
  for(;;) {
    const    x = midl(x1, x2),
       f = F(x);
    if(f===0) return x;

    if( x1 === x ) break;

    if( isNaN(f) ) throw new Error('NaN encountered.');
    if( f < 0 ) { x1 = x; f1 = f }
    if( f > 0 ) { x2 = x; f2 = f }
  }

  return Math.abs(f1) <= Math.abs(f2) ? x1 : x2;
}
