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
var _exportNames = {
  ARRAY_TYPES: true,
  eps: true,
  cast_scalar: true,
  _check_dtype: true,
  dtypeof: true,
  super_dtype: true,
  is_subdtype: true
};
exports.eps = eps;
exports.cast_scalar = cast_scalar;
exports._check_dtype = _check_dtype;
exports.dtypeof = dtypeof;
exports.is_subdtype = is_subdtype;
exports.super_dtype = exports.ARRAY_TYPES = void 0;

var _complex = require("./complex");

Object.keys(_complex).forEach(function (key) {
  if (key === "default" || key === "__esModule") return;
  if (Object.prototype.hasOwnProperty.call(_exportNames, key)) return;
  Object.defineProperty(exports, key, {
    enumerable: true,
    get: function get() {
      return _complex[key];
    }
  });
});

var _complex_array = require("./complex_array");

Object.keys(_complex_array).forEach(function (key) {
  if (key === "default" || key === "__esModule") return;
  if (Object.prototype.hasOwnProperty.call(_exportNames, key)) return;
  Object.defineProperty(exports, key, {
    enumerable: true,
    get: function get() {
      return _complex_array[key];
    }
  });
});
var ARRAY_TYPES = {
  'int32': Int32Array,
  'float32': Float32Array,
  'float64': Float64Array,
  'complex128': _complex_array.Complex128Array,
  'object': Array
};
exports.ARRAY_TYPES = ARRAY_TYPES;

function eps(dtype) {
  _check_dtype(dtype);

  switch (dtype) {
    case 'complex128':
      return new _complex.Complex(Number.EPSILON);

    case 'float32':
      return 1.1920928955078125e-7;

    default:
      return Number.EPSILON;
  }
}

function cast_scalar(x, dtype) {
  if (dtype === 'int32') return x & 0xFFFFFFFF;
  if (dtype === 'float32') return Math.fround(x);
  if (dtype === 'complex128') return x instanceof _complex.Complex ? x : new _complex.Complex(x);
  return x * 1;
}

function _check_dtype(dtype) {
  if (!ARRAY_TYPES.hasOwnProperty(dtype)) throw new Error("Invalid dtype '" + dtype + "'. Must be one of {'" + Object.getOwnPropertyNames(ARRAY_TYPES).join("', '") + "'}.");
}

function dtypeof(value) {
  if (value % 1 === 0) {
    if (value <= ~(1 << 31) && value >= 1 << 31) return 'int32';
    return 'object';
  }

  if (value * 1 == value) return 'float64';
  if (value instanceof _complex.Complex) return 'complex128';
  return 'object';
}

var super_dtype = function super_dtype() {
  for (var _len = arguments.length, dtypes = new Array(_len), _key = 0; _key < _len; _key++) {
    dtypes[_key] = arguments[_key];
  }

  return dtypes.reduce(function (dtype1, dtype2) {
    _check_dtype(dtype1);

    _check_dtype(dtype2);

    if (dtype1 === 'object' || dtype2 === 'object') return 'object';
    if (dtype1 === 'complex128' || dtype2 === 'complex128') return 'complex128';
    if (dtype1 === 'float64' || dtype2 === 'float64') return 'float64';
    if (dtype1 === 'float32' || dtype2 === 'float32') return 'float32';
    return 'int32';
  });
};

exports.super_dtype = super_dtype;

function is_subdtype(sub_dtype, sup_dtype) {
  _check_dtype(sub_dtype);

  _check_dtype(sup_dtype);

  var rank = {
    'int32': 0,
    'float32': 1,
    'float64': 2,
    'complex128': 3,
    'object': 4
  };
  return rank[sub_dtype] <= rank[sup_dtype];
}