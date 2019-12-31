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


/** Applies a Givens rotation to rows i and j.
 */
export function _giv_rot_rows( W, N, i, j, c, s )
{
  i |= 0;
  j |= 0;
  c = +c;
  s = +s;
  for( let k=N | 0; k-- > 0; ) {
    const W_i = W[i],
          W_j = W[j];
    W[i] = c*W_i + s*W_j;
    W[j] = c*W_j - s*W_i;
    i = i+1 | 0;
    j = j+1 | 0;
  }
}


/** Applies a Givens rotation to columns i and j (to a square matrix).
 */
export function _giv_rot_cols( W, N, i, j, c, s )
{
  i |= 0;
  j |= 0;
  c = +c;
  s = +s;
  for( let k=N | 0; k-- > 0; ) {
    const W_i = W[i],
          W_j = W[j];
    W[i] = c*W_i - s*W_j;
    W[j] = c*W_j + s*W_i;
    i = i+N | 0;
    j = j+N | 0;
  }
}
