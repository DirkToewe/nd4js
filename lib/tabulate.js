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
exports.tabulate = tabulate;

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _nd_array = require("./nd_array");

var _dt = require("./dt");

function tabulate(shape, dtype, idx2val) {
  shape = Int32Array.from(shape, function (s) {
    if (s % 1 !== 0) throw new Error("tabulate(shape, dtype='object', idx2val): Invalid shape [".concat(shape, "]."));
    return s;
  });

  if (null == idx2val) {
    idx2val = dtype;
    dtype = undefined;
  }

  if (null == idx2val) throw new Error("tabulate(shape, dtype='object', idx2val): idx2val missing.");
  if (null == dtype) dtype = 'object';
  (0, _dt._check_dtype)(dtype);
  var multi_idx = new Int32Array(shape.length),
      // <- index in result
  data = new _dt.ARRAY_TYPES[dtype](shape.reduce(function (a, b) {
    return a * b;
  }, 1));
  var flat_idx = 0;

  function write(d) {
    if (d === shape.length) data[flat_idx++] = idx2val.apply(void 0, (0, _toConsumableArray2["default"])(multi_idx));else for (multi_idx[d] = 0; multi_idx[d] < shape[d]; multi_idx[d]++) {
      write(d + 1);
    }
  }

  write(0);
  return new _nd_array.NDArray(shape, data);
}