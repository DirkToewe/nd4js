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
exports.zip_elems = zip_elems;

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _dt = require("./dt");

var _nd_array = require("./nd_array");

function _createForOfIteratorHelper(o) { if (typeof Symbol === "undefined" || o[Symbol.iterator] == null) { if (Array.isArray(o) || (o = _unsupportedIterableToArray(o))) { var i = 0; var F = function F() {}; return { s: F, n: function n() { if (i >= o.length) return { done: true }; return { done: false, value: o[i++] }; }, e: function e(_e) { throw _e; }, f: F }; } throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method."); } var it, normalCompletion = true, didErr = false, err; return { s: function s() { it = o[Symbol.iterator](); }, n: function n() { var step = it.next(); normalCompletion = step.done; return step; }, e: function e(_e2) { didErr = true; err = _e2; }, f: function f() { try { if (!normalCompletion && it["return"] != null) it["return"](); } finally { if (didErr) throw err; } } }; }

function _unsupportedIterableToArray(o, minLen) { if (!o) return; if (typeof o === "string") return _arrayLikeToArray(o, minLen); var n = Object.prototype.toString.call(o).slice(8, -1); if (n === "Object" && o.constructor) n = o.constructor.name; if (n === "Map" || n === "Set") return Array.from(o); if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n)) return _arrayLikeToArray(o, minLen); }

function _arrayLikeToArray(arr, len) { if (len == null || len > arr.length) len = arr.length; for (var i = 0, arr2 = new Array(len); i < len; i++) { arr2[i] = arr[i]; } return arr2; }

function zip_elems(ndarrays, dtype, zip_fn) {
  // PREPROCESS ARGUMENTS
  if (null == zip_fn && dtype instanceof Function) {
    zip_fn = dtype;
    dtype = undefined;
  }

  if (!(ndarrays instanceof Array)) ndarrays = (0, _toConsumableArray2["default"])(ndarrays);
  ndarrays = ndarrays.map(function (a) {
    return (0, _nd_array.asarray)(a);
  });
  if (null == dtype) dtype = 'object';

  if (null == zip_fn) {
    if ('object' !== dtype) throw new Error('If zip_fn is undefined, dtype must be "object" or undefined.');
    var L = ndarrays.length;

    zip_fn = function zip_fn() {
      for (var _len = arguments.length, x = new Array(_len), _key = 0; _key < _len; _key++) {
        x[_key] = arguments[_key];
      }

      return x.slice(0, L);
    };
  }

  (0, _dt._check_dtype)(dtype);
  var ndim = ndarrays.reduce(function (ndim, arr) {
    return Math.max(ndim, arr.ndim);
  }, 0),
      shape = Int32Array.from({
    length: ndim
  }, function () {
    return 1;
  }); // FIND COMMON (BROADCASTED) SHAPE

  var _iterator = _createForOfIteratorHelper(ndarrays),
      _step;

  try {
    for (_iterator.s(); !(_step = _iterator.n()).done;) {
      var arr = _step.value;

      for (var i = ndim, j = arr.ndim; i-- > 0 && j-- > 0;) {
        if (1 === shape[i]) shape[i] = arr.shape[j];else if (shape[i] != arr.shape[j] && arr.shape[j] != 1) throw new Error('Shapes are not broadcast-compatible.');
      }
    } // GENERATE RESULT DATA

  } catch (err) {
    _iterator.e(err);
  } finally {
    _iterator.f();
  }

  var multi_idx = new Int32Array(ndim),
      // <- index in result
  data = new _dt.ARRAY_TYPES[dtype](shape.reduce(function (a, b) {
    return a * b;
  }, 1)),
      values = new Array(ndarrays.length),
      // <- cache of ndarrays[indices]
  indices = new Int32Array(ndarrays.length),
      // <- indices in ndarrays(s)
  strides = new Int32Array(ndarrays.length);
  var flat_idx = 0;

  function write(d) {
    if (d === ndim) {
      strides.fill(1);

      for (var i = ndarrays.length; i-- > 0;) {
        values[i] = ndarrays[i].data[indices[i]++];
      }

      data[flat_idx++] = zip_fn.apply(void 0, values.concat((0, _toConsumableArray2["default"])(multi_idx)));
      return;
    }

    for (multi_idx[d] = 0;;) {
      write(d + 1);
      if (++multi_idx[d] >= shape[d]) break;

      for (var _i = ndarrays.length; _i-- > 0;) {
        if (!(ndarrays[_i].shape[d - ndim + ndarrays[_i].ndim] > 1)) // <- handles undefined (index out of bounds)
          indices[_i] -= strides[_i];
      }
    }

    for (var _i2 = ndarrays.length; _i2-- > 0;) {
      strides[_i2] *= ndarrays[_i2].shape[d - ndim + ndarrays[_i2].ndim] || 1;
    }
  }

  write(0);
  return new _nd_array.NDArray(shape, data);
}