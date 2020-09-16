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



export function bitCount( bits )
{
  bits =  bits | 0;
  bits = (bits>>> 1 & 0b01010101010101010101010101010101) + (bits & 0b01010101010101010101010101010101);
  bits = (bits>>> 2 & 0b00110011001100110011001100110011) + (bits & 0b00110011001100110011001100110011);
  bits = (bits>>> 4 & 0b00001111000011110000111100001111) + (bits & 0b00001111000011110000111100001111);
  bits = (bits>>> 8 & 0b00000000111111110000000011111111) + (bits & 0b00000000111111110000000011111111);
  bits = (bits>>>16                                     ) + (bits & 0b00000000000000001111111111111111);
  return  bits;
}
