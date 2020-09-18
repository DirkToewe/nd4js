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

var _interopRequireDefault = require("@babel/runtime/helpers/interopRequireDefault");

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.nextUp = nextUp;
exports.nextDown = nextDown;
exports.midl = midl;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _io = require("../io");

var val = Float64Array.of(NaN),
    bits = new Int32Array(val.buffer);
/*DEBUG*/

if (bits.length !== 2) throw new Error('Assertion failed.');

function nextUp(x) {
  x *= 1;
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
  x *= 1;
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

function midl(x, y) {
  x *= 1;
  y *= 1;
  if (isNaN(x)) throw new Error('mid(x,y): x must be number.');
  if (isNaN(y)) throw new Error('mid(x,y): x must be number.');

  if (Math.sign(x) * y < 0) // <- check for opposite signs
    {
      var mid = (x + y) / 2; // <- avoids underflow if e.g x=+MAX_VALUE, y=-MAX_VALUE

      if (!(mid < Math.max(x, y))) throw new Error('Assertion failed.');
      if (!(mid > Math.min(x, y))) throw new Error('Assertion failed.');
      return mid;
    } else {
    var _ref = Math.abs(x) <= Math.abs(y) ? [x, y] : [y, x],
        _ref2 = (0, _slicedToArray2["default"])(_ref, 2),
        a = _ref2[0],
        b = _ref2[1]; // Returns the mid point between two floats (x,y) or x if there is no mid point


    var _mid = a + (b - a) * 0.5; // <- avoids underflow if e.g x=y=MAX_VALUE


    if (!(_mid <= Math.max(x, y))) throw new Error('Assertion failed.');
    if (!(_mid >= Math.min(x, y))) throw new Error('Assertion failed.');
    if (_mid === y) return x;
    return _mid;
  }
}