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

import {roots1d_polyquad} from '../polyquad'


export function _min1d_interp_gg( x1, x2, g1, g2 )
{
  x1*=1; if( ! isFinite(x1) ) throw new Error('Assertion failed.');
  x2*=1; if( ! isFinite(x2) ) throw new Error('Assertion failed.');
  g1*=1; if( ! isFinite(g1) ) throw new Error('Assertion failed.');
  g2*=1; if( ! isFinite(g2) ) throw new Error('Assertion failed.');

  const dx = x2-x1,
        dg = g2-g1;

  if( ! (Math.sign(dx)*dg > 0) )
    throw new Error('Assertion failed: Bad curvature (i.e. no minimum).');

  //           !
  // g1 + s*dg = 0  =>  s = -g1/dg
  return x1  -  dx/dg * g1;
}


export function _min1d_interp_ffg( x1, x2, f1, f2, g1 )
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

  const             dx = x2-x1, dFdX = (f2-f1) / dx,
    xMin = x1 + 0.5*dx*g1 / (g1-dFdX);

  // TODO: maybe check if xMin is a min. or max.

  return xMin;
}


export function _min1d_interp_ffgg( x1, x2, f1, f2, g1, g2 )
{
  x1*=1; if( ! isFinite(x1) ) throw new Error('Assertion failed.');
  x2*=1; if( ! isFinite(x2) ) throw new Error('Assertion failed.');
  f1*=1; if( ! isFinite(f1) ) throw new Error('Assertion failed.');
  f2*=1; if( ! isFinite(f2) ) throw new Error('Assertion failed.');
  g1*=1; if( ! isFinite(g1) ) throw new Error('Assertion failed.');
  g2*=1; if( ! isFinite(g2) ) throw new Error('Assertion failed.');

  if( x1===x2 ) throw new Error('Assertion failed.');

  f2 -= f1; // shift such that f1 = 0
  x2 -= x1; // shift such that x1 = 0

  // f (x) := a*x³/3 + b*x²/2 + g1*x
  // f'(x)  = a*x²   + b*x    + g1
  //
  //      !
  // f(0) = 0 (check!)
  //
  //       !
  // f'(0) = g1 (check!)
  //
  //       !
  // (A) f(x2) = f2 = a/3*x2*x2*x2 + b/2*x2*x2 + g1*x2
  //
  //        !
  // (B) f'(x2) = g2 = a*x2*x2 + b*x2 + g1
  //           => a*x2*x2 = g2-g1 - b*x2
  //
  // (A+B) => f2 = x2/3*(g2-g1 - b*x2)  +  b/2*x2*x2  +  g1*x2
  //       => b*x2*x2/6 = f2 - x2*(g2 + g1*2)/3
  //       => b = ( 6*f1/x2 - 2*g2 - 4*g1 ) / x2

  // const b = 2*( 3*f2/x2 - g2 - 2*g1 ) / x2,
  //       a = ((g2-g1)/x2 - b) / x2,
  // [z1,z2] = roots1d_polyquad(g1,b,a);
  // return (a < 0 ? z1 : z2) + x1;

  const b = 2*( 3*f2/x2 - g2 - 2*g1 ),
        a = (g2 - g1 - b) / x2,
  [z1,z2] = roots1d_polyquad(g1*x2,b,a);

  return (Math.sign(x2)*a < 0 ? z1 : z2) + x1;
}
