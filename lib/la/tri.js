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
exports.tril = tril;
exports.triu = triu;
exports.tril_solve = tril_solve;
exports.triu_solve = triu_solve;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _dt = require("../dt");

var _nd_array = require("../nd_array");

function tril(m) {
  var k = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0;
  m = (0, _nd_array.asarray)(m);
  if (m.ndim < 2) throw new Error('Input must be at least 2D.');
  return m.mapElems(m.dtype, function (x) {
    for (var _len = arguments.length, indices = new Array(_len > 1 ? _len - 1 : 0), _key = 1; _key < _len; _key++) {
      indices[_key - 1] = arguments[_key];
    }

    var _indices$slice = indices.slice(-2),
        _indices$slice2 = (0, _slicedToArray2["default"])(_indices$slice, 2),
        i = _indices$slice2[0],
        j = _indices$slice2[1];

    return i < j - k ? 0 : x;
  });
}

function triu(m) {
  var k = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0;
  m = (0, _nd_array.asarray)(m);
  if (m.ndim < 2) throw new Error('Input must be at least 2D.');
  return m.mapElems(m.dtype, function (x) {
    for (var _len2 = arguments.length, indices = new Array(_len2 > 1 ? _len2 - 1 : 0), _key2 = 1; _key2 < _len2; _key2++) {
      indices[_key2 - 1] = arguments[_key2];
    }

    var _indices$slice3 = indices.slice(-2),
        _indices$slice4 = (0, _slicedToArray2["default"])(_indices$slice3, 2),
        i = _indices$slice4[0],
        j = _indices$slice4[1];

    return i > j - k ? 0 : x;
  });
}

function tril_solve(L, y) {
  L = (0, _nd_array.asarray)(L);
  if (L.ndim < 2) throw new Error('tril_solve(L,y): L.ndim must be at least 2.');
  y = (0, _nd_array.asarray)(y);
  if (y.ndim < 2) throw new Error('tril_solve(L,y): y.ndim must be at least 2.');

  var _L$shape$slice = L.shape.slice(-2),
      _L$shape$slice2 = (0, _slicedToArray2["default"])(_L$shape$slice, 2),
      N = _L$shape$slice2[0],
      M = _L$shape$slice2[1],
      _y$shape$slice = y.shape.slice(-2),
      _y$shape$slice2 = (0, _slicedToArray2["default"])(_y$shape$slice, 2),
      I = _y$shape$slice2[0],
      J = _y$shape$slice2[1];

  if (N != M) throw new Error('Last two dimensions of L must be quadratic.');
  if (I != M) throw new Error("L and y don't match.");
  var ndim = Math.max(L.ndim, y.ndim),
      shape = Int32Array.from({
    length: ndim
  }, function () {
    return 1;
  });
  shape[ndim - 2] = I;
  shape[ndim - 1] = J; // FIND COMMON (BROADCASTED) SHAPE

  for (var _i = 0, _arr = [L, y]; _i < _arr.length; _i++) {
    var arr = _arr[_i];

    for (var i = ndim - 2, j = arr.ndim - 2; i-- > 0 && j-- > 0;) {
      if (1 === shape[i]) shape[i] = arr.shape[j];else if (shape[i] != arr.shape[j] && arr.shape[j] != 1) throw new Error('Shapes are not broadcast-compatible.');
    }
  } // GENERATE RESULT DATA


  var dtype = [L, y].every(function (a) {
    return a.dtype === 'float32';
  }) ? 'float32' : 'float64',
      x_dat = new _dt.ARRAY_TYPES[dtype](shape.reduce(function (a, b) {
    return a * b;
  })),
      L_dat = L.data,
      y_dat = y.data;
  var L_off = 0,
      L_stride = 1,
      y_off = 0,
      y_stride = 1,
      x_off = 0;

  function solv(d) {
    if (d === ndim - 2) {
      L_stride = N * N;
      y_stride = N * J; // COPYING y

      for (var i = 0; i < y_stride; i++) {
        x_dat[x_off + i] = y_dat[y_off + i];
      } // FORWARD SUBSTITUTION


      for (var _i2 = 0; _i2 < I; _i2++) {
        for (var j = 0; j < J; j++) {
          for (var k = 0; k < _i2; k++) {
            x_dat[x_off + _i2 * J + j] -= L_dat[L_off + N * _i2 + k] * x_dat[x_off + k * J + j];
          }

          x_dat[x_off + _i2 * J + j] /= L_dat[L_off + N * _i2 + _i2];
        }
      }

      L_off += L_stride;
      y_off += y_stride;
      x_off += y_stride;
      return;
    }

    for (var l = shape[d];; l--) {
      solv(d + 1);
      if (l == 1) break;
      if (!(L.shape[d - ndim + L.ndim] > 1)) L_off -= L_stride;
      if (!(y.shape[d - ndim + y.ndim] > 1)) y_off -= y_stride;
    }

    L_stride *= L.shape[d - ndim + L.ndim] || 1;
    y_stride *= y.shape[d - ndim + y.ndim] || 1;
  }

  solv(0);
  return new _nd_array.NDArray(shape, x_dat);
}

function triu_solve(U, y) {
  U = (0, _nd_array.asarray)(U);
  if (U.ndim < 2) throw new Error('triu_solve(U,y): U.ndim must be at least 2.');
  y = (0, _nd_array.asarray)(y);
  if (y.ndim < 2) throw new Error('triu_solve(U,y): y.ndim must be at least 2.');

  var _U$shape$slice = U.shape.slice(-2),
      _U$shape$slice2 = (0, _slicedToArray2["default"])(_U$shape$slice, 2),
      K = _U$shape$slice2[0],
      N = _U$shape$slice2[1],
      _y$shape$slice3 = y.shape.slice(-2),
      _y$shape$slice4 = (0, _slicedToArray2["default"])(_y$shape$slice3, 2),
      I = _y$shape$slice4[0],
      J = _y$shape$slice4[1];

  if (K != N) throw new Error('Last two dimensions of U must be quadratic.');
  if (I != N) throw new Error("U and y don't match.");
  var ndim = Math.max(U.ndim, y.ndim),
      shape = Int32Array.from({
    length: ndim
  }, function () {
    return 1;
  });
  shape[ndim - 2] = I;
  shape[ndim - 1] = J; // FIND COMMON (BROADCASTED) SHAPE

  for (var _i3 = 0, _arr2 = [U, y]; _i3 < _arr2.length; _i3++) {
    var arr = _arr2[_i3];

    for (var i = ndim - 2, j = arr.ndim - 2; i-- > 0 && j-- > 0;) {
      if (1 === shape[i]) shape[i] = arr.shape[j];else if (shape[i] != arr.shape[j] && arr.shape[j] != 1) throw new Error('Shapes are not broadcast-compatible.');
    }
  } // GENERATE RESULT DATA


  var dtype = [U, y].every(function (a) {
    return a.dtype === 'float32';
  }) ? 'float32' : 'float64',
      x_dat = new _dt.ARRAY_TYPES[dtype](shape.reduce(function (a, b) {
    return a * b;
  })),
      U_dat = U.data,
      y_dat = y.data;
  var U_off = 0,
      U_stride = 1,
      y_off = 0,
      y_stride = 1,
      x_off = 0;

  function solv(d) {
    if (d === ndim - 2) {
      U_stride = N * N;
      y_stride = N * J; // COPYING y

      for (var i = 0; i < y_stride; i++) {
        x_dat[x_off + i] = y_dat[y_off + i];
      } // BACKWARD SUBSTITUTION


      for (var _i4 = I; _i4-- > 0;) {
        for (var j = J; j-- > 0;) {
          for (var k = K; --k > _i4;) {
            x_dat[x_off + _i4 * J + j] -= U_dat[U_off + N * _i4 + k] * x_dat[x_off + k * J + j];
          }

          x_dat[x_off + _i4 * J + j] /= U_dat[U_off + N * _i4 + _i4];
        }
      }

      U_off += U_stride;
      y_off += y_stride;
      x_off += y_stride;
      return;
    }

    for (var l = shape[d];; l--) {
      solv(d + 1);
      if (l == 1) break;
      if (!(U.shape[d - ndim + U.ndim] > 1)) U_off -= U_stride;
      if (!(y.shape[d - ndim + y.ndim] > 1)) y_off -= y_stride;
    }

    U_stride *= U.shape[d - ndim + U.ndim] || 1;
    y_stride *= y.shape[d - ndim + y.ndim] || 1;
  }

  solv(0);
  return new _nd_array.NDArray(shape, x_dat);
}