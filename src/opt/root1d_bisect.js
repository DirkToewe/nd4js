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

export function root1d_bisect(F,x_min,x_max)
{
  let f_min = F(x_min),
      f_max = F(x_max);
  if( f_min > f_max )  [f_min,f_max,x_min,x_max] = [f_max,f_min,x_max,x_min];
  if( f_min > 0 ) { if( f_min <= f_max*Number.EPSILON ) return x_min; else throw new Error('root1d_bisect(F,x_min,x_max): Both F(x_min) and F(x_max) positive.') }
  if( f_max < 0 ) { if( f_max >= f_min*Number.EPSILON ) return x_max; else throw new Error('root1d_bisect(F,x_min,x_max): Both F(x_min) and F(x_max) negaitive.') }
  if( f_min > 0 ) { return x_min; }
  if( f_max < 0 ) { return x_max; }
  
  for(;;) {
    const x = (x_min + x_max) / 2;
    if( x_min == x ||
        x_max == x ) break;
    const f = F(x);
    if( isNaN(f) ) throw new Error('NaN encountered.');
    if( f <= 0 ) { x_min = x; f_min = f }
    if( f >= 0 ) { x_max = x; f_max = f }
  }

  return Math.abs(f_min) <= Math.abs(f_max) ? x_min : x_max;
}
