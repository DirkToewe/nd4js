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

import {array, asarray, NDArray} from '../nd_array'
import {TrustRegionSolverLSQ} from './_trust_region_solver'


// TODO implemented double-dogleg method as well as described in
// "An Improved Optimization Method for iSAM2"
// Rana Talaei Shahir and Hamid D. Taghirad Senior Member, IEEE
// Proceeding of the 2nd
// RSI/ISM International Conference on Robotics and Mechatronics
// October 15-17, 2014, Tehran, Iran

// TODO implement dogbox method as well which allows for box-constrained optimization

/* Returns the more positive root of a + bx + cx²
 */
export function roots1d_polyquad( a, b, c )
{
  // https://arxiv.org/abs/1409.8072
  if( isNaN(a) ) throw new Error('roots_quad_poly(a,b,c): a must be number.');
  if( isNaN(b) ) throw new Error('roots_quad_poly(a,b,c): b must be number.');
  if( isNaN(c) ) throw new Error('roots_quad_poly(a,b,c): c must be number.');

  if( 0 === c ) {
    const   x = -a/b;
    return [x,x];
  }

  a /= c;
  b /= c;

  // a + bx + x²
  if( 0 === b ) {
    const x = Math.sqrt(-a);
    return [-x,+x];
  }
  if( 0 === a )
    return 0 <= b ? [-b,0] : [0,-b];

  c  = Math.sqrt(Math.abs(a)) * Math.sign(b);
  b /= 2*c;

  const TOL = Number.EPSILON * 1024 * 64;

  let x1,x2;
  if( 0 < a ) {
    if( b < 1-TOL) throw new Error('Complex roots not yet supported: '+b); // <- TODO: this should be a custom error
    if( b < 1    ) b=1;
    x1 = -b - Math.sqrt((b+1)*(b-1));
    x2 = +1 / x1;
  }
  else {
    x1 = -b - Math.sqrt(b*b + 1);
    x2 = -1 / x1
  }
  x1 *= c;
  x2 *= c;
  return x1 <= x2 ? [x1,x2] : [x2,x1];
}
