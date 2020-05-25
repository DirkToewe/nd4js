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
exports._cholesky_decomp = _cholesky_decomp;
exports.cholesky_decomp = cholesky_decomp;
exports.cholesky_solve = cholesky_solve;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _dt = require("../dt");

var _kahan_sum = require("../kahan_sum");

var _nd_array = require("../nd_array");

var _tri = require("./tri");

function _cholesky_decomp(M, N, L, L_off) {
  if (!(M <= N)) throw new Error('Assertion failed.');
  var kahan = new _kahan_sum.KahanSum(); // https://de.wikipedia.org/wiki/Cholesky-Zerlegung

  for (var i = 0; i < M; i++) {
    for (var j = 0; j <= i; j++) {
      kahan.set(L[L_off + N * i + j]);

      for (var k = 0; k < j; k++) {
        kahan.add(-L[L_off + N * i + k] * L[L_off + N * j + k]);
      }

      if (i > j) L[L_off + N * i + j] = kahan.sum / L[L_off + N * j + j];else {
        L[L_off + N * i + i] = Math.sqrt(kahan.sum);
        if (isNaN(L[L_off + N * i + i])) throw new Error('Matrix contains NaNs or is (near) singular.');
      }
    }
  }
}

function cholesky_decomp(S) {
  S = (0, _nd_array.asarray)(S);

  var dtype = S.dtype === 'float32' ? 'float32' : 'float64',
      shape = S.shape,
      _shape$slice = shape.slice(-2),
      _shape$slice2 = (0, _slicedToArray2["default"])(_shape$slice, 2),
      N = _shape$slice2[0],
      M = _shape$slice2[1];

  S = S.data;
  var L = new _dt.ARRAY_TYPES[dtype](S.length);
  if (N != M) throw new Error('Last two dimensions must be quadratic.');

  for (var L_off = 0; L_off < L.length; L_off += N * N) {
    for (var i = 0; i < N; i++) {
      for (var j = 0; j <= i; j++) {
        L[L_off + N * i + j] = S[L_off + N * i + j];
      }
    }

    _cholesky_decomp(N, N, L, L_off);
  }

  return new _nd_array.NDArray(shape, L);
}

function cholesky_solve(L, y) {
  L = (0, _nd_array.asarray)(L);
  y = (0, _nd_array.asarray)(y);
  if (L.ndim < 2) throw new Error('L must be at least 2D.');
  if (y.ndim < 2) throw new Error('y must be at least 2D.');

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


  var dtype = L.dtype === 'float32' && y.dtype === 'float32' ? 'float32' : 'float64',
      x_dat = new _dt.ARRAY_TYPES[dtype](shape.reduce(function (a, b) {
    return a * b;
  }, 1)),
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

      for (var _i2 = 0; _i2 < y_stride; _i2++) {
        x_dat[x_off + _i2] = y_dat[y_off + _i2];
      }

      (0, _tri._tril_solve)(I, I, J, L_dat, L_off, x_dat, x_off);
      (0, _tri._tril_t_solve)(I, I, J, L_dat, L_off, x_dat, x_off);
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