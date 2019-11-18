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
exports["default"] = exports.math = void 0;

var _complex = require("./dt/complex");

var math = {
  nextUp: function nextUp(x) {
    // FIXME implement this completely
    // https://gist.github.com/Yaffle/4654250
    var EPSILON = Number.EPSILON,
        MIN_VALUE = Number.MIN_VALUE;
    var y = x * (x < 0 ? 1 - EPSILON / 2 : 1 + EPSILON);
    if (!isFinite(x)) throw new Error('Assertion failed!');
    if (!isFinite(y)) throw new Error('Assertion failed!');
    if (y === x) y += MIN_VALUE;
    var b = x + (y - x) / 2;

    if (x < b && b < y) {
      y = b;
    }

    var c = (y + x) / 2;

    if (x < c && c < y) {
      y = c;
    }

    return y === 0 ? -0 : y;
  },
  add: function add(x, y) {
    if (x instanceof _complex.Complex) return x.add(y);
    if (y instanceof _complex.Complex) return y.add(x);
    return x + y;
  },
  sub: function sub(x, y) {
    if (x instanceof _complex.Complex) return x.sub(y);
    if (y instanceof _complex.Complex) return new _complex.Complex(x).sub(y);
    return x - y;
  },
  mul: function mul(x, y) {
    if (x instanceof _complex.Complex) return x.mul(y);
    if (y instanceof _complex.Complex) return y.mul(x);
    return x * y;
  },
  div: function div(x, y) {
    if (x instanceof _complex.Complex) return x.div(y);
    if (y instanceof _complex.Complex) return new _complex.Complex(x).div(y);
    return x / y;
  },
  zero: function zero(dtype) {
    return math.cast(0, dtype);
  },
  one: function one(dtype) {
    return math.cast(1, dtype);
  },
  neg: function neg(x) {
    return x instanceof _complex.Complex ? x.neg() : -x;
  },
  abs: function abs(x) {
    return x instanceof _complex.Complex ? x.abs() : Math.abs(x);
  },
  sqrt: function sqrt(x) {
    return x instanceof _complex.Complex ? x.sqrt() : x >= 0 ? Math.sqrt(x) : new _complex.Complex(x).sqrt();
  },
  exp: function exp(x) {
    return x instanceof _complex.Complex ? x.exp() : Math.exp(x);
  },
  min: function min(x, y) {
    return Math.min(x, y);
  },
  max: function max(x, y) {
    return Math.max(x, y);
  },
  hypot: function hypot(x, y) {
    return Math.hypot(x, y);
  },
  atan2: function atan2(x, y) {
    return Math.atan2(x, y);
  },
  conj: function conj(x) {
    return x instanceof _complex.Complex ? x.conj() : x;
  },
  is_equal: function is_equal(x, y) {
    if (x instanceof _complex.Complex) return x.equals(y);
    if (y instanceof _complex.Complex) return y.equals(x);
    return x == y;
  },
  is_close: function is_close(x, y) {
    var atol = 1e-8,
        rtol = 1e-5,
        tol = atol + rtol * math.max(math.abs(x), math.abs(y));
    return math.abs(math.sub(x, y)) <= tol;
  }
};
exports.math = math;
var _default = math;
exports["default"] = _default;