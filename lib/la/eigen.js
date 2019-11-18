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
exports.eigen = eigen;
exports.eigenvals = eigenvals;
exports.eigen_balance_pre = eigen_balance_pre;
exports.eigen_balance_post = eigen_balance_post;

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _schur = require("./schur");

var _nd_array = require("../nd_array");

var _dt = require("../dt");

var _zip_elems = require("../zip_elems");

var _math = _interopRequireDefault(require("../math"));

function eigen(A) {
  var _eigen_balance_pre = eigen_balance_pre(A, 2),
      _eigen_balance_pre2 = (0, _slicedToArray2["default"])(_eigen_balance_pre, 2),
      D = _eigen_balance_pre2[0],
      B = _eigen_balance_pre2[1],
      _schur_decomp = (0, _schur.schur_decomp)(B),
      _schur_decomp2 = (0, _slicedToArray2["default"])(_schur_decomp, 2),
      Q = _schur_decomp2[0],
      T = _schur_decomp2[1],
      _schur_eigen = (0, _schur.schur_eigen)(Q, T),
      _schur_eigen2 = (0, _slicedToArray2["default"])(_schur_eigen, 2),
      Λ = _schur_eigen2[0],
      V = _schur_eigen2[1];

  if (!V.dtype.startsWith('complex')) throw new Error('Assertion failed.');
  if (!D.dtype.startsWith('float')) throw new Error('Assertion failed.');
  var N = A.shape[A.ndim - 1],
      norm_sum = new Float64Array(N),
      norm_max = new Float64Array(N),
      V_dat = V.data._array,
      D_dat = D.data; // UNDO BALANCING

  for (var D_off = 0, V_off = 0; V_off < V_dat.length; D_off += N, V_off += 2 * N * N) {
    norm_sum.fill(0);
    norm_max.fill(0); // SCALE ROWS & COMPUTE SCALED COLUMN NORMS

    for (var i = 0; i < N; i++) {
      var D_i = D_dat[D_off + i];

      for (var j = 0; j < 2 * N; j++) {
        // scale row
        var V_ij = Math.abs(V_dat[V_off + 2 * N * i + j] *= D_i); // update norm

        if (V_ij > 0) {
          var k = j >> 1;

          if (V_ij > norm_max[k]) {
            norm_sum[k] *= Math.pow(norm_max[k] / V_ij, 2);
            norm_max[k] = V_ij;
          }

          norm_sum[k] += Math.pow(V_ij / norm_max[k], 2);
        }
      }
    }

    for (var _j = 0; _j < N; _j++) {
      var max = norm_max[_j];
      norm_sum[_j] = isFinite(max) ? Math.sqrt(norm_sum[_j]) * max : max;
    } // NORMALIZE COLUMNS


    for (var _i = 0; _i < N; _i++) {
      for (var _j2 = 0; _j2 < 2 * N; _j2++) {
        V_dat[V_off + 2 * N * _i + _j2] /= norm_sum[_j2 >> 1];
      }
    }
  }

  return [Λ, V];
}

function eigenvals(A) {
  var _eigen_balance_pre3 = eigen_balance_pre(A, 2),
      _eigen_balance_pre4 = (0, _slicedToArray2["default"])(_eigen_balance_pre3, 2),
      D = _eigen_balance_pre4[0],
      B = _eigen_balance_pre4[1],
      _schur_decomp3 = (0, _schur.schur_decomp)(B),
      _schur_decomp4 = (0, _slicedToArray2["default"])(_schur_decomp3, 2),
      Q = _schur_decomp4[0],
      T = _schur_decomp4[1];

  return (0, _schur.schur_eigenvals)(T);
}

function eigen_balance_pre(A, p) {
  // SEE: RODNEY JAMES, JULIEN LANGOU, AND BRADLEY R. LOWERY,
  //      "ON MATRIX BALANCING AND EIGENVECTOR COMPUTATION"
  if (p == null) p = 2;
  if (p > Number.MAX_VALUE) return eigen_balance_pre_inf(A);
  var N = A.shape[A.ndim - 1],
      DTypeArray = _dt.ARRAY_TYPES[A.dtype === 'float32' ? 'float32' : 'float64'],
      A_shape = A.shape,
      D_shape = A_shape.slice(0, -1);
  if (!(p >= 1)) throw new Error("Invalid norm p=".concat(p, ";"));
  if (A.shape[A.ndim - 2] != N) throw new Error('A is not square');
  A = DTypeArray.from(A.data); // <- protection copy

  var D = new DTypeArray(A.length / N),
      TOL = Math.pow(0.95, 1 / p);
  D.fill(1.0);

  for (var A_off = A.length, D_off = D.length; A_off > 0;) {
    A_off -= N * N;
    D_off -= N;

    for (var i, done = false; !done;) {
      for (i = 0, done = true; i < N; i++) {
        var c = 0.0,
            c_max = 0.0,
            r = 0.0,
            r_max = 0.0; // COMPUTE ROW AND COLUMN NORM

        for (var j = 0; j < N; j++) {
          //        if(true)
          if (i !== j) {
            {
              var A_ij = Math.abs(A[A_off + N * i + j]);

              if (A_ij > 0) {
                if (A_ij > r_max) {
                  var _scale = r_max / A_ij;

                  r_max = A_ij;
                  r *= Math.pow(_scale, p);
                }

                var ratio = A_ij / r_max;
                r += Math.pow(ratio, p);
              }
            }
            {
              var A_ji = Math.abs(A[A_off + N * j + i]);

              if (A_ji > 0) {
                if (A_ji > c_max) {
                  var _scale2 = c_max / A_ji;

                  c_max = A_ji;
                  c *= Math.pow(_scale2, p);
                }

                var _ratio = A_ji / c_max;

                c += Math.pow(_ratio, p);
              }
            }
          }
        }

        r = !isFinite(r) ? r : Math.pow(r, 1 / p) * r_max;
        c = !isFinite(c) ? c : Math.pow(c, 1 / p) * c_max;
        if (r * c == 0.0) continue;
        if (!isFinite(r * c)) throw new Error('NaN encountered.');
        var old_norm = c >= r ? Math.pow(1 + Math.pow(r / c, p), 1 / p) * c : Math.pow(1 + Math.pow(c / r, p), 1 / p) * r; // FIND A POWER OF TWO THAT MORE OR LESS BALANCES r AND c

        var scale = 1.0;

        while (r >= c * 2) {
          c *= 2;
          r /= 2;
          scale *= 2;
        }

        while (c >= r * 2) {
          c /= 2;
          r *= 2;
          scale /= 2;
        }

        var new_norm = c >= r ? Math.pow(1 + Math.pow(r / c, p), 1 / p) * c : Math.pow(1 + Math.pow(c / r, p), 1 / p) * r;
        if (new_norm >= TOL * old_norm) continue;
        done = false;
        D[D_off + i] *= scale;

        for (var _j3 = 0; _j3 < N; _j3++) {
          A[A_off + N * i + _j3] /= scale;
          A[A_off + N * _j3 + i] *= scale;
        }
      }
    }
  }

  return [new _nd_array.NDArray(D_shape, D), new _nd_array.NDArray(A_shape, A)];
}

function eigen_balance_pre_inf(A) {
  // SEE: RODNEY JAMES, JULIEN LANGOU, AND BRADLEY R. LOWERY,
  //      "ON MATRIX BALANCING AND EIGENVECTOR COMPUTATION"
  var N = A.shape[A.ndim - 1],
      DTypeArray = _dt.ARRAY_TYPES[A.dtype === 'float32' ? 'float32' : 'float64'],
      A_shape = A.shape,
      D_shape = A_shape.slice(0, -1);
  if (A.shape[A.ndim - 2] != N) throw new Error('A is not square');
  A = DTypeArray.from(A.data);
  var D = new DTypeArray(A.length / N);
  D.fill(1.0);

  for (var A_off = A.length, D_off = D.length; A_off > 0;) {
    A_off -= N * N;
    D_off -= N;

    for (var i, done = false; !done;) {
      for (i = 0, done = true; i < N; i++) {
        var c = 0.0,
            r = 0.0; // COMPUTE ROW AND COLUMN NORM

        for (var j = 0; j < N; j++) {
          if (i !== j) {
            var A_ij = Math.abs(A[A_off + N * i + j]);
            r = Math.max(r, A_ij);
            var A_ji = Math.abs(A[A_off + N * j + i]);
            c = Math.max(c, A_ji);
          }
        }

        if (r * c === 0.0) continue;
        var old_norm = Math.max(c, r);
        if (old_norm == 0.0) continue;
        if (!isFinite(old_norm)) throw new Error('NaN encountered.'); // FIND A POWER OF TWO THAT MORE OR LESS BALANCES r AND c

        var scale = 1.0;

        while (r >= c * 2) {
          c *= 2;
          r /= 2;
          scale *= 2;
        }

        while (c >= r * 2) {
          c /= 2;
          r *= 2;
          scale /= 2;
        }

        var new_norm = Math.max(c, r);
        if (new_norm >= old_norm) continue;
        done = false;
        D[D_off + i] *= scale;

        for (var _j4 = 0; _j4 < N; _j4++) {
          A[A_off + N * i + _j4] /= scale;
          A[A_off + N * _j4 + i] *= scale;
        }
      }
    }
  }

  return [new _nd_array.NDArray(D_shape, D), new _nd_array.NDArray(A_shape, A)];
}

function eigen_balance_post(D, V) {
  if (V.ndim < 2) throw new Error('eigen_balance_post(D,V): V.ndim must be at least 2.');

  var _V$shape$slice = V.shape.slice(-2),
      _V$shape$slice2 = (0, _slicedToArray2["default"])(_V$shape$slice, 2),
      M = _V$shape$slice2[0],
      N = _V$shape$slice2[1];

  if (M !== N) throw new Error('eigen_balance_post(D,V): V must be square.');
  var DTypeArray = _dt.ARRAY_TYPES['float64'],
      V_arr = (0, _zip_elems.zip_elems)([V, D.reshape.apply(D, (0, _toConsumableArray2["default"])(D.shape).concat([1]))], 'complex128', _math["default"].mul);
  V = V_arr.data;
  var norm_sum = new DTypeArray(N),
      norm_max = new DTypeArray(N); // NORMALIZE COLUMNS

  for (var off = 0; off < V.length; off += M * N) {
    norm_sum.fill(0);
    norm_max.fill(0); // COMPUTE COLUMN NORMS

    for (var i = 0; i < N; i++) {
      for (var j = 0; j < N; j++) {
        var V_ij = _math["default"].abs(V[off + N * i + j]);

        if (V_ij > 0) {
          if (V_ij > norm_max[j]) {
            norm_sum[j] *= Math.pow(norm_max[j] / V_ij, 2);
            norm_max[j] = V_ij;
          }

          norm_sum[j] += Math.pow(V_ij / norm_max[j], 2);
        }
      }
    }

    for (var _j5 = 0; _j5 < N; _j5++) {
      var max = norm_max[_j5];
      norm_sum[_j5] = isFinite(max) ? Math.sqrt(norm_sum[_j5]) * max : max;
    } // NORMALIZE COLUMNS


    for (var _i2 = 0; _i2 < N; _i2++) {
      for (var _j6 = 0; _j6 < N; _j6++) {
        V[off + N * _i2 + _j6] = _math["default"].div(V[off + N * _i2 + _j6], norm_sum[_j6]);
      }
    }
  }

  return V_arr;
}