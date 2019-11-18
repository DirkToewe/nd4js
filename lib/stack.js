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
exports.stack = stack;

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _concat = require("./concat");

function stack(axis, dtype, ndarrays) {
  if (null == ndarrays) {
    if (null == dtype) {
      ndarrays = axis;
      axis = undefined;
    } else {
      ndarrays = dtype;
      dtype = undefined;

      if ('string' === typeof axis) {
        dtype = axis;
        axis = undefined;
      }
    }
  }

  if (!('length' in ndarrays)) ndarrays = (0, _toConsumableArray2["default"])(ndarrays);
  if (null == axis) axis = 0;
  if (0 > axis) axis += ndarrays[0].shape.length + 1;
  if (0 > axis || axis > ndarrays[0].shape.length) throw new Error('Axis out of bounds.');
  ndarrays = ndarrays.map(function (arr) {
    return arr.reshape.apply(arr, (0, _toConsumableArray2["default"])(arr.shape.slice(0, axis)).concat([1], (0, _toConsumableArray2["default"])(arr.shape.slice(axis))));
  });
  return (0, _concat.concat)(axis, dtype, ndarrays);
}