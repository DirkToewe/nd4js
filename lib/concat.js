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
exports.concat = concat;

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _nd_array = require("./nd_array");

var _dt = require("./dt");

function concat(axis, dtype, ndarrays) {
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

  ndarrays = Array.from(ndarrays, function (arr) {
    return (0, _nd_array.asarray)(dtype, arr);
  });
  if (null == axis) axis = 0;
  if (null == dtype) dtype = _dt.super_dtype.apply(void 0, (0, _toConsumableArray2["default"])(ndarrays.map(function (a) {
    return a.dtype;
  })));
  if (0 > axis) axis += ndarrays[0].shape.length;
  if (0 > axis || axis >= ndarrays[0].shape.length) throw new Error('Axis out of bounds.');
  var newShape = Int32Array.from(ndarrays[0].shape);

  var _loop = function _loop(i) {
    var shape = ndarrays[i].shape;
    if (newShape.length != shape.length) throw new Error('All shapes must have the same length.');
    if (!newShape.every(function (len, d) {
      return newShape[d] === shape[d] || axis === d;
    })) throw new Error('Shape along all axes but the concatentation axis must match.');
    newShape[axis] += shape[axis];
  };

  for (var i = ndarrays.length; --i > 0;) {
    _loop(i);
  }

  var rest = newShape.slice(axis + 1).reduce(function (a, b) {
    return a * b;
  }, 1),
      indices = Int32Array.from(ndarrays, function (ndarr) {
    return ndarr.data.length;
  }),
      newData = new _dt.ARRAY_TYPES[dtype](newShape.reduce(function (a, b) {
    return a * b;
  }, 1));
  var newIdx = newData.length;

  function fill(d) {
    if (d === axis) for (var _i = ndarrays.length; _i-- > 0;) {
      for (var j = ndarrays[_i].shape[d] * rest; j-- > 0;) {
        newData[--newIdx] = ndarrays[_i].data[--indices[_i]];
      }
    } else for (var _i2 = newShape[d]; _i2-- > 0;) {
      fill(d + 1);
    }
  }

  fill(0);
  return new _nd_array.NDArray(newShape, newData);
}