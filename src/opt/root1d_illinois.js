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

import {nextDown,
        nextUp} from "../dt/float64_utils";


// REFERENCES
// ----------
// .. [1] https://en.wikipedia.org/wiki/False_position_method
// .. [2] "Improved Algorithms Of Illinois-Type For The Numerical Solution Of Nonlinear Equations"
//         J. A. Ford
//


export function root1d_illinois(F,x1,x2)
{
  x1 *= 1;
  x2 *= 1;
  if( isNaN(x1) ) throw new Error('root1d_falsi(F,x1,x2): x1 must be a number.');
  if( isNaN(x2) ) throw new Error('root1d_falsi(F,x1,x2): x2 must be a number.');

  if( ! (F instanceof Function) ) throw new Error('root1d_falsi(F,x1,x2): F must be a function.');

  let f1 = F(x1),
      f2 = F(x2),
      F1 = f1;

  if( 0 === f1 ) return x1;
  if( 0 === f2 ) return x2;

  if( ! (Math.sign(f1)*f2 < 0) )
    throw new Error('root1d_falsi(F,x1,x2): F(x1) and F(x2) must have opposite signs.');

  for(;;) {
    const l = nextUp  ( Math.min(x1,x2) ),
          r = nextDown( Math.max(x1,x2) );
    if( r < l ) break;

    let x = (x1*f2 - x2*f1) / (f2 - f1);
    if( !(x>=l) ) x = l;
    if( !(x<=r) ) x = r;

    const  f   =   F(x);
    if(0===f) return x;

    if( isNaN(f) ) throw new Error('NaN encountered.');

    if( Math.sign(f2)*f < 0 ) { x1 = x2; x2 = x;
                                f1 = f2; f2 = f;
                                F1 = f1; }
    else {
      // const scale = 0.5; // <- ILLINOIS METHOD
      // let scale = 1 - f/f2; if( scale <= 0 ) scale = 0.5; // <- ANDERSON & BJÖRK

      let   scale = 1 - f/f2/(1-f/f1); // FORD's METHOD 3 (see [2])
      if( !(scale > 0) ) // <- handles NaN
            scale = 0.5;

      f1 *= scale;
      x2 = x;
      f2 = f;
    }
  }

  return Math.abs(F1) <= Math.abs(f2) ? x1 : x2;
}
