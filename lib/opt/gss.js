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

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.min1d_gss = min1d_gss;
var γ = 1.5 - Math.sqrt(5 / 4);

function min1d_gss(F, x_min, x_max) {
  // https://en.wikipedia.org/wiki/Golden_section_search
  if (x_min > x_max) throw new Error('opt1d_golden(F,x_min,x_max): x_max must not be less than x_min.');

  for (var a = x_min, b = x_max, c = a + (b - a) * γ, Fc = F(c);;) {
    var d = b - (b - a) * γ;
    if (d === b) return (b + a) / 2;
    var Fd = F(d);

    if (Fc < Fd) {
      b = a;
      a = d;
    } else {
      a = c;
      c = d;
      Fc = Fd;
    }
  }
}