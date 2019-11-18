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
exports.eye = eye;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _dt = require("../dt");

var _nd_array = require("../nd_array");

function eye() {
  for (var _len = arguments.length, shape = new Array(_len), _key = 0; _key < _len; _key++) {
    shape[_key] = arguments[_key];
  }

  var dtype = shape[0] in _dt.ARRAY_TYPES ? shape.shift() : 'float64';
  if (shape.length < 1) throw new Error('Size parameter missing.');
  if (shape.length == 1) shape.push(shape[shape.length - 1]);
  shape = Int32Array.from(shape, function (s) {
    if (0 !== s % 1) throw new Error("eye(): Invalid shape [".concat(shape, "]."));
    return s;
  });

  var I = new _dt.ARRAY_TYPES[dtype](shape.reduce(function (m, n) {
    return m * n;
  })),
      _shape$slice = shape.slice(-2),
      _shape$slice2 = (0, _slicedToArray2["default"])(_shape$slice, 2),
      M = _shape$slice2[0],
      N = _shape$slice2[1];

  for (var off = I.length; (off -= M * N) >= 0;) {
    for (var i = Math.min(M, N); i-- > 0;) {
      I[off + N * i + i] = 1;
    }
  }

  var result = new _nd_array.NDArray(shape, I);
  return result;
}