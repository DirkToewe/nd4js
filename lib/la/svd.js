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
exports.svd_rank = svd_rank;
exports.svd_solve = svd_solve;
exports.svd_lstsq = svd_lstsq;
Object.defineProperty(exports, "svd_decomp", {
  enumerable: true,
  get: function get() {
    return _svd_dc.svd_dc;
  }
});

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _math = _interopRequireDefault(require("../math"));

var _nd_array = require("../nd_array");

var _dt = require("../dt");

var _singular_matrix_solve_error = require("./singular_matrix_solve_error");

var _svd_dc = require("./svd_dc");

function svd_rank(sv) {
  sv = (0, _nd_array.asarray)(sv);

  var N = sv.shape[sv.ndim - 1],
      r_shape = sv.shape.slice(0, -1),
      EPS = _math["default"].sqrt((0, _dt.eps)(sv.dtype));

  sv = sv.data;
  var r = new _dt.ARRAY_TYPES['int32'](sv.length / N);

  for (var off = 0; off < r.length; off++) {
    var T = EPS * _math["default"].abs(sv[N * off]);

    if (r[off] !== 0) throw new Error('Assertion failed.');

    for (; r[off] < N; r[off]++) {
      var sv_r = _math["default"].abs(sv[N * off + r[off]]);

      if (!isFinite(sv_r)) throw new Error('svd_rank(): NaN or Infinity encountered.');
      if (sv_r <= T) break;
    }
  }

  return new _nd_array.NDArray(r_shape, r);
}

function svd_solve(U, sv, V, y) {
  U = (0, _nd_array.asarray)(U);
  sv = (0, _nd_array.asarray)(sv);
  V = (0, _nd_array.asarray)(V);
  if (U.shape[U.ndim - 2] !== V.shape[V.ndim - 1]) throw new Error('rrqr_solve(Q,R,P, y): System not square.');

  var x = svd_lstsq(U, sv, V, y),
      EPS = _math["default"].sqrt((0, _dt.eps)(sv.dtype)),
      N = sv.shape[sv.ndim - 1];

  sv = sv.data;

  for (var sv_off = 0; sv_off < sv.length; sv_off += N) {
    var T = EPS * _math["default"].abs(sv[sv_off]);

    for (var r; r < N; r++) {
      var sv_r = _math["default"].abs(sv[sv_off + r]);

      if (!isFinite(sv_r)) throw new Error('svd_solve(): NaN or Infinity encountered.');
      if (sv_r <= T) throw new _singular_matrix_solve_error.SingularMatrixSolveError(x);
    }
  }

  return x;
}

function svd_lstsq(U, sv, V, y) {
  if (y == undefined) {
    var _Q, _Q2;

    if (V != undefined) throw new Error('svd_lstsq(Q,R,P, y): Either 2 ([Q,R,P], y) or 4 arguments (Q,R,P, y) expected.');
    y = R((_Q = Q, _Q2 = (0, _slicedToArray2["default"])(_Q, 3), U = _Q2[0], sv = _Q2[1], V = _Q2[2], _Q));
  }

  U = (0, _nd_array.asarray)(U);
  if (U.ndim < 2) throw new Error('svd_lstsq(U,sv,V, y): U.ndim must be at least 2.');
  sv = (0, _nd_array.asarray)(sv);
  if (sv.ndim < 1) throw new Error('svd_lstsq(U,sv,V, y): sv.ndim must be at least 1.');
  V = (0, _nd_array.asarray)(V);
  if (V.ndim < 2) throw new Error('svd_lstsq(U,sv,V, y): V.ndim must be at least 2.');
  y = (0, _nd_array.asarray)(y);
  if (y.ndim < 2) throw new Error('svd_lstsq(U,sv,V, y): y.ndim must be at least 2.');

  var _U$shape$slice = U.shape.slice(-2),
      _U$shape$slice2 = (0, _slicedToArray2["default"])(_U$shape$slice, 2),
      N = _U$shape$slice2[0],
      M = _U$shape$slice2[1],
      _V$shape$slice = V.shape.slice(-1),
      _V$shape$slice2 = (0, _slicedToArray2["default"])(_V$shape$slice, 1),
      I = _V$shape$slice2[0],
      _y$shape$slice = y.shape.slice(-1),
      _y$shape$slice2 = (0, _slicedToArray2["default"])(_y$shape$slice, 1),
      J = _y$shape$slice2[0];

  if (N != y.shape[y.ndim - 2]) throw new Error("svd_lstsq(U,sv,V, y): U and y don't match.");
  if (M != sv.shape[sv.ndim - 1]) throw new Error("svd_lstsq(U,sv,V, y): U and sv don't match.");
  if (M != V.shape[V.ndim - 2]) throw new Error("svd_lstsq(U,sv,V, y): V and sv don't match.");
  var ndim = Math.max(U.ndim, sv.ndim + 1, V.ndim, y.ndim),
      shape = Int32Array.from({
    length: ndim
  }, function () {
    return 1;
  });
  shape[ndim - 2] = I;
  shape[ndim - 1] = J; // FIND COMMON (BROADCASTED) SHAPE

  for (var _i = 0, _arr = [U, V, y]; _i < _arr.length; _i++) {
    var arr = _arr[_i];

    for (var _i5 = ndim - 2, _j4 = arr.ndim - 2; _i5-- > 0 && _j4-- > 0;) {
      if (1 === shape[_i5]) shape[_i5] = arr.shape[_j4];else if (shape[_i5] != arr.shape[_j4] && arr.shape[_j4] != 1) throw new Error('svd_lstsq(U,sv,V, y): U,sv,V,y not broadcast-compatible.');
    }
  }

  for (var i = ndim - 2, j = sv.ndim - 1; i-- > 0 && j-- > 0;) {
    if (1 === shape[i]) shape[i] = sv.shape[j];else if (shape[i] != sv.shape[j] && sv.shape[j] != 1) throw new Error('svd_lstsq(U,sv,V, y): U,sv,V,y not broadcast-compatible.');
  } // GENERATE RESULT DATA


  var EPS = _math["default"].sqrt((0, _dt.eps)(sv.dtype)),
      DTypeArray = _dt.ARRAY_TYPES[[U, sv, V, y].every(function (a) {
    return a.dtype === 'float32';
  }) ? 'float32' : 'float64'],
      tmp = new DTypeArray(M * J),
      x_dat = new DTypeArray(shape.reduce(function (a, b) {
    return a * b;
  }, 1)),
      U_dat = U.data,
      sv_dat = sv.data,
      V_dat = V.data,
      y_dat = y.data;

  var U_off = 0,
      U_stride = 1,
      sv_off = 0,
      sv_stride = 1,
      V_off = 0,
      V_stride = 1,
      y_off = 0,
      y_stride = 1,
      x_off = 0;

  function solv(d) {
    if (d === ndim - 2) {
      U_stride = N * M;
      sv_stride = M;
      V_stride = M * I;
      y_stride = N * J;

      var _R = function () {
        var T = EPS * _math["default"].abs(sv_dat[sv_off]);

        for (var r = 0; r < M; r++) {
          var sv_r = _math["default"].abs(sv_dat[sv_off + r]);

          if (!isFinite(sv_r)) throw new Error('svd_solve(): NaN or Infinity encountered.');
          if (sv_r <= T) return r;
        }

        return M;
      }(); // SEE: Gene H. Golub, Charles F. Van Golub
      //     "Matrix Computations", 4th edition
      //      Page 260f, Chap. 5.3.1 (Implications of Full Rank)
      // tmp = U.T @ y


      tmp.fill(0);

      for (var k = 0; k < N; k++) {
        for (var _i2 = 0; _i2 < _R; _i2++) {
          for (var _j = 0; _j < J; _j++) {
            tmp[J * _i2 + _j] += U_dat[U_off + M * k + _i2] * y_dat[y_off + J * k + _j];
          }
        }
      } // tmp \= diag(sv)


      for (var _i3 = 0; _i3 < _R; _i3++) {
        for (var _j2 = 0; _j2 < J; _j2++) {
          tmp[J * _i3 + _j2] /= sv_dat[sv_off + _i3];
        }
      } // x = V.T @ tmp


      for (var _k = 0; _k < _R; _k++) {
        for (var _i4 = 0; _i4 < I; _i4++) {
          for (var _j3 = 0; _j3 < J; _j3++) {
            x_dat[x_off + J * _i4 + _j3] += V_dat[V_off + I * _k + _i4] * tmp[J * _k + _j3];
          }
        }
      }

      U_off += U_stride;
      sv_off += sv_stride;
      V_off += V_stride;
      y_off += y_stride;
      x_off += I * J;
      return;
    }

    for (var l = shape[d];; l--) {
      solv(d + 1);
      if (l == 1) break;
      if (!(U.shape[d - ndim + U.ndim] > 1)) U_off -= U_stride;
      if (!(sv.shape[d - ndim + 1 + sv.ndim] > 1)) sv_off -= sv_stride;
      if (!(V.shape[d - ndim + V.ndim] > 1)) V_off -= V_stride;
      if (!(y.shape[d - ndim + y.ndim] > 1)) y_off -= y_stride;
    }

    U_stride *= U.shape[d - ndim + U.ndim] || 1;
    sv_stride *= sv.shape[d - ndim + 1 + sv.ndim] || 1;
    V_stride *= V.shape[d - ndim + V.ndim] || 1;
    y_stride *= y.shape[d - ndim + y.ndim] || 1;
  }

  solv(0);
  return new _nd_array.NDArray(shape, x_dat);
  throw new Error('Not yet implemented!');
}