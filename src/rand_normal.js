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

//export function rand_normal()
//{
//  // https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
//  const x = Math.random() / (1 - Number.EPSILON) + Number.MIN_VALUE,
//        y = Math.random();
//
//  return Math.sqrt(-2 * Math.log(x)) * Math.cos(Math.PI*2 * y);
//}

export const rand_normal = function(){
  let next = NaN;

  return () => {
    if( ! isNaN(next) ) {
      const  nxt = next;
                   next = NaN;
      return nxt;
    }
    // https://en.wikipedia.org/wiki/Marsaglia_polar_method
    let x,y,r;
    do {
      x = Math.random()*2 - 1;
      y = Math.random()*2 - 1;
      r = x*x + y*y;
    }
    while(r > 1 || r == 0);

    const  z = Math.sqrt(-2 * Math.log(r) / r);
    next = z*x;
    return z*y;
  };
}();

// export const rand_normal = function(){
//   let next = NaN;

//   return () => {
//     if( ! isNaN(next) ) {
//       const  nxt = next;
//                    next = NaN;
//       return nxt;
//     }
//     // https://en.wikipedia.org/wiki/Marsaglia_polar_method
//     let x,y,r;
//     do {
//       x = Math.random()*2 - 1;
//       y = Math.random()*2 - 1;
//       r = Math.hypot(x,y);
//     }
//     while(r > 1 || r == 0);

//     const  z = Math.sqrt( -Math.log(r) ) * 2/r;
//     next = z*x;
//     return z*y;
//   };
// }();
