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
exports.lu_decomp = lu_decomp;
exports.lu_solve = lu_solve;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _nd_array = require("../nd_array");

var _dt = require("../dt");

var _tri = require("./tri");

function lu_decomp(A) {
  A = (0, _nd_array.asarray)(A);

  var LU = _dt.ARRAY_TYPES[A.dtype === 'float32' ? 'float32' : 'float64'].from(A.data),
      _A$shape$slice = A.shape.slice(-2),
      _A$shape$slice2 = (0, _slicedToArray2["default"])(_A$shape$slice, 2),
      N = _A$shape$slice2[0],
      M = _A$shape$slice2[1],
      P = new Int32Array(LU.length / M);

  if (N != M) throw new Error('Last two dimensions must be quadratic.');

  for (var LU_off = 0, P_off = 0; LU_off < LU.length; LU_off += N * N, P_off += N) {
    for (var i = 0; i < N; i++) {
      P[P_off + i] = i;
    }

    for (var _i = 0; _i < N; _i++) {
      var row_i = LU_off + _i * N; // ROW PIVOTING

      {
        var p = _i;

        for (var j = _i + 1; j < N; j++) {
          if (Math.abs(LU[LU_off + N * j + _i]) > Math.abs(LU[LU_off + N * p + _i])) p = j;
        }

        if (_i != p) {
          var P_p = P[P_off + _i];
          P[P_off + _i] = P[P_off + p];
          P[P_off + p] = P_p; // KEEP TRACK OF ROW SWAPS

          var row_p = LU_off + p * N; // SWAP ROWS

          for (var _j = 0; _j < N; _j++) {
            var tmp = LU[row_i + _j];
            LU[row_i + _j] = LU[row_p + _j];
            LU[row_p + _j] = tmp;
          }
        }
      } // ELIMINATE ELEMENTS BELOW PIVOT

      for (var _j2 = _i + 1; _j2 < N; _j2++) {
        var row_j = LU_off + _j2 * N,
            scale = LU[row_j + _i] / LU[row_i + _i];
        LU[row_j + _i] = scale;

        for (var k = _i + 1; k < N; k++) {
          LU[row_j + k] -= scale * LU[row_i + k];
        }
      }
    }
  }

  return [new _nd_array.NDArray(A.shape, LU), new _nd_array.NDArray(A.shape.slice(0, -1), P)];
}

function lu_solve(LU, P, y) {
  if (undefined == y) {
    y = P;
    var _LU = LU;

    var _LU2 = (0, _slicedToArray2["default"])(_LU, 2);

    LU = _LU2[0];
    P = _LU2[1];
  }

  LU = (0, _nd_array.asarray)(LU);
  if (LU.ndim < 2) throw new Error('LU must be at least 2D.');
  P = (0, _nd_array.asarray)(P);
  if (P.ndim < 1) throw new Error('P must be at least 1D.');
  y = (0, _nd_array.asarray)(y);
  if (y.ndim < 2) throw new Error('y must be at least 2D.');

  var _LU$shape$slice = LU.shape.slice(-2),
      _LU$shape$slice2 = (0, _slicedToArray2["default"])(_LU$shape$slice, 2),
      N = _LU$shape$slice2[0],
      M = _LU$shape$slice2[1],
      _y$shape$slice = y.shape.slice(-2),
      _y$shape$slice2 = (0, _slicedToArray2["default"])(_y$shape$slice, 2),
      I = _y$shape$slice2[0],
      J = _y$shape$slice2[1];

  if (M != N) throw new Error('Last two dimensions of LU must be quadratic.');
  if (M != I) throw new Error("LU and y don't match.");
  if (M != P.shape.slice(-1)) throw new Error("LU and P don't match.");
  var ndim = Math.max(LU.ndim, P.ndim + 1, y.ndim),
      shape = Int32Array.from({
    length: ndim
  }, function () {
    return 1;
  });
  shape[ndim - 2] = I;
  shape[ndim - 1] = J; // FIND COMMON (BROADCASTED) SHAPE

  for (var _i2 = 0, _arr = [LU, y]; _i2 < _arr.length; _i2++) {
    var arr = _arr[_i2];

    for (var _i5 = ndim - 2, _j5 = arr.ndim - 2; _i5-- > 0 && _j5-- > 0;) {
      if (1 === shape[_i5]) shape[_i5] = arr.shape[_j5];else if (shape[_i5] != arr.shape[_j5] && arr.shape[_j5] != 1) throw new Error('LU and y are not broadcast-compatible.');
    }
  }

  for (var i = ndim - 2, j = P.ndim - 1; i-- > 0 && j-- > 0;) {
    if (1 === shape[i]) shape[i] = P.shape[j];else if (shape[i] != P.shape[j] && P.shape[j] != 1) throw new Error('P is not broadcast-compatible.');
  } // GENERATE RESULT DATA


  var x_dat = new _dt.ARRAY_TYPES[[LU, P, y].every(function (a) {
    return a.dtype === 'float32';
  }) ? 'float32' : 'float64'](shape.reduce(function (a, b) {
    return a * b;
  }, 1));
  var LU_dat = LU.data,
      P_dat = P.data,
      y_dat = y.data;
  var LU_off = 0,
      LU_stride = 1,
      P_off = 0,
      P_stride = 1,
      y_off = 0,
      y_stride = 1,
      x_off = 0;

  function solv(d) {
    if (d === ndim - 2) {
      LU_stride = N * N;
      P_stride = N;
      y_stride = N * J; // COPYING PERMUTED y

      for (var _i3 = 0; _i3 < I; _i3++) {
        var row_x = x_off + J * _i3,
            row_y = y_off + J * P_dat[P_off + _i3];

        for (var _j3 = 0; _j3 < J; _j3++) {
          x_dat[row_x + _j3] = y_dat[row_y + _j3];
        }
      } // FORWARD SUBSTITUTION


      for (var _i4 = 0; _i4 < I; _i4++) {
        for (var _j4 = 0; _j4 < J; _j4++) {
          for (var k = 0; k < _i4; k++) {
            x_dat[x_off + _i4 * J + _j4] -= LU_dat[LU_off + N * _i4 + k] * x_dat[x_off + k * J + _j4];
          }
        }
      } // BACKWARD SUBSTITUTION


      (0, _tri._triu_solve)(I, I, J, LU_dat, LU_off, x_dat, x_off);
      LU_off += LU_stride;
      P_off += P_stride;
      y_off += y_stride;
      x_off += y_stride;
      return;
    }

    for (var l = shape[d];; l--) {
      solv(d + 1);
      if (l == 1) break; // RESET ALONG BROADCAST AXES

      if (!(LU.shape[d - ndim + LU.ndim] > 1)) LU_off -= LU_stride;
      if (!(P.shape[d - ndim + P.ndim + 1] > 1)) P_off -= P_stride;
      if (!(y.shape[d - ndim + y.ndim] > 1)) y_off -= y_stride;
    }

    LU_stride *= LU.shape[d - ndim + LU.ndim] || 1;
    P_stride *= P.shape[d - ndim + P.ndim + 1] || 1;
    y_stride *= y.shape[d - ndim + y.ndim] || 1;
  }

  solv(0);
  return new _nd_array.NDArray(shape, x_dat);
}