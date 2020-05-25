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
exports.nextUp = nextUp;
exports.nextDown = nextDown;

var _io = require("../io");

var val = Float64Array.of(NaN),
    bits = new Int32Array(val.buffer);
/*DEBUG*/

if (bits.length !== 2) throw new Error('Assertion failed.');

function nextUp(x) {
  if (0 === x) return +Number.MIN_VALUE;
  if (!(x < +Infinity)) // <- handles NaN and Infinity
    return x;
  val[0] = x;
  var i = 1 - _io.IS_LITTLE_ENDIAN,
      j = 1 * _io.IS_LITTLE_ENDIAN;

  if (x > 0) {
    bits[j] += -1 === bits[i];
    bits[i] += 1;
  } else {
    bits[j] -= 0 === bits[i];
    bits[i] -= 1;
  }

  return val[0];
}

;

function nextDown(x) {
  if (0 === x) return -Number.MIN_VALUE;
  if (!(x > -Infinity)) // <- handles NaN and Infinity
    return x;
  val[0] = x;
  var i = 1 - _io.IS_LITTLE_ENDIAN,
      j = 1 * _io.IS_LITTLE_ENDIAN;

  if (x > 0) {
    bits[j] -= 0 === bits[i];
    bits[i] -= 1;
  } else {
    bits[j] += -1 === bits[i];
    bits[i] += 1;
  }

  return val[0];
}