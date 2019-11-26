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

const φ = Math.sqrt(5/4) + 0.5;


export function opt1d_golden(F,x_min,x_max)
{
  // https://en.wikipedia.org/wiki/Golden_section_search
  if( x_min > x_max )
    throw new Error('opt1d_golden(F,x_min,x_max): x_max must not be less than x_min.');

  for( let a = x_min,
           b = x_max;; )
  {
    const c = b - (b - a) / φ,
          d = a + (b - a) / φ;

    if( c===d ) return (b + a) / 2;

    if( F(c) < F(d) ) b = d;
    else              a = c;
  }
}
