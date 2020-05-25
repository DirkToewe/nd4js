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
exports.qr_decomp_full = qr_decomp_full;
exports.qr_decomp = qr_decomp;
exports.qr_lstsq = qr_lstsq;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _nd_array = require("../nd_array");

var _dt = require("../dt");

var _giv_rot = require("./_giv_rot");

var _transpose_inplace2 = require("./transpose_inplace");

var _tri = require("./tri");

function qr_decomp_full(A) {
  A = (0, _nd_array.asarray)(A);
  if (A.ndim < 2) throw new Error('A must be at least 2D.');

  var DType = A.dtype === 'float32' ? 'float32' : 'float64',
      DTypeArray = _dt.ARRAY_TYPES[DType],
      B = DType === 'float32' ? 64 / 4 : 64 / 8,
      R_shape = A.shape,
      Q_shape = Int32Array.from(R_shape),
      _R_shape$slice = R_shape.slice(-2),
      _R_shape$slice2 = (0, _slicedToArray2["default"])(_R_shape$slice, 2),
      M = _R_shape$slice2[0],
      N = _R_shape$slice2[1],
      R = DTypeArray.from(A.data);

  A = undefined;
  Q_shape[Q_shape.length - 1] = M;
  var Q = new DTypeArray(R.length / N * M); //  Q.fill(0); // <- in case of an object array

  for (var Q_off = 0, R_off = 0; Q_off < Q.length; Q_off += M * M, R_off += M * N) {
    // INIT Q TO IDENTITY MATRIX
    for (var i = 0; i < M; i++) {
      Q[Q_off + M * i + i] = 1.0;
    } // The idea of the blocked loop is that the top B rows


    for (var J = 0; J < N; J += B) {
      // are cached and used to elimate B columns in the remaining
      for (var I = J; I < M; I += B) {
        // rows. This should reduce the number of cache misses.
        for (var _i = I; _i < I + B && _i < M; _i++) {
          // Some quick and dirty benchmarks indicate a ~15% perfomance
          for (var j = J; j < J + B && j < N && j < _i; j++) // improvement for matrix sizes of [300,300] and above.
          {
            // USE GIVENS ROTATION TO ELIMINATE ELEMENT R_ji
            var ij = R_off + N * _i + j,
                R_ij = R[ij];
            if (0 === R_ij) continue;

            var jj = R_off + N * j + j,
                R_jj = R[jj],
                _giv_rot_qr2 = (0, _giv_rot._giv_rot_qr)(R_jj, R_ij),
                _giv_rot_qr3 = (0, _slicedToArray2["default"])(_giv_rot_qr2, 3),
                c = _giv_rot_qr3[0],
                s = _giv_rot_qr3[1],
                norm = _giv_rot_qr3[2];

            R[ij] = 0;
            if (0 === s) continue;
            R[jj] = norm;
            (0, _giv_rot._giv_rot_rows)(R, N - 1 - j, jj + 1, ij + 1, c, s);
            (0, _giv_rot._giv_rot_rows)(Q, 1 + _i, Q_off + M * j, Q_off + M * _i, c, s);
          }
        }
      }
    }

    (0, _transpose_inplace2._transpose_inplace)(M, Q, Q_off);
  }

  return [new _nd_array.NDArray(Q_shape, Q), new _nd_array.NDArray(R_shape, R)];
}

function qr_decomp(A) {
  A = (0, _nd_array.asarray)(A);
  if (A.ndim < 2) throw new Error('qr_decomp(A): A.ndim must be at least 2.');

  var DTypeArray = _dt.ARRAY_TYPES[A.dtype === 'float32' ? 'float32' : 'float64'],
      Q_shape = A.shape,
      R_shape = Int32Array.from(Q_shape),
      _Q_shape$slice = Q_shape.slice(-2),
      _Q_shape$slice2 = (0, _slicedToArray2["default"])(_Q_shape$slice, 2),
      N = _Q_shape$slice2[0],
      M = _Q_shape$slice2[1];

  R_shape[R_shape.length - 2] = M;
  if (N <= M) return qr_decomp_full(A);
  var Q = DTypeArray.from(A.data); // <- we could encourage GC by setting `A = undefined` after this line

  A = undefined;
  var R = new DTypeArray(Q.length / N * M); // <- additional space to temp. store rows of R not contained in the result

  for (var R_off = 0, Q_off = 0; Q_off < Q.length; Q_off += N * M, R_off += M * M) {
    // COMPUTE R (inside of Q)
    for (var i = 1; i < N; i++) {
      var I = Math.min(i, M);

      for (var j = 0; j < I; j++) {
        // USE GIVENS ROTATION TO ELIMINATE ELEMENT R_ji
        var ij = Q_off + M * i + j,
            R_ij = Q[ij];
        if (0 === R_ij) continue;
        var jj = Q_off + M * j + j,
            R_jj = Q[jj];

        var _giv_rot_qr4 = (0, _giv_rot._giv_rot_qr)(R_jj, R_ij),
            _giv_rot_qr5 = (0, _slicedToArray2["default"])(_giv_rot_qr4, 3),
            c = _giv_rot_qr5[0],
            s = _giv_rot_qr5[1],
            norm = _giv_rot_qr5[2];

        if (s !== 0) {
          if (c < 0) {
            c *= -1;
            s *= -1;
            norm *= -1;
          }

          (0, _giv_rot._giv_rot_rows)(Q, M - 1 - j, jj + 1, ij + 1, c, s);
          Q[jj] = norm;
        }

        Q[ij] = s;
      }
    } // MOVE R FROM Q -> R AND INIT Q TO I


    for (var _i2 = 0; _i2 < M; _i2++) {
      for (var _j = _i2; _j < M; _j++) {
        R[R_off + M * _i2 + _j] = Q[Q_off + M * _i2 + _j];
        Q[Q_off + M * _i2 + _j] = +(_i2 === _j);
      }
    } // COMPUTE Q


    for (var _i3 = N; --_i3 > 0;) {
      var _I = Math.min(_i3, M);

      for (var _j2 = _I; _j2-- > 0;) {
        var _s = Q[Q_off + M * _i3 + _j2];
        if (0 === _s) continue;
        Q[Q_off + M * _i3 + _j2] = 0;

        var _c = Math.sqrt((1 - _s) * (1 + _s));

        (0, _giv_rot._giv_rot_rows)(Q, M - _j2, Q_off + M * _i3 + _j2, Q_off + M * _j2 + _j2, _c, _s);
      }
    }
  }

  return [new _nd_array.NDArray(Q_shape, Q), new _nd_array.NDArray(R_shape, R)];
}

function qr_lstsq(Q, R, y) {
  if (undefined == y) {
    y = R;
    var _Q = Q;

    var _Q2 = (0, _slicedToArray2["default"])(_Q, 2);

    Q = _Q2[0];
    R = _Q2[1];
  }

  Q = (0, _nd_array.asarray)(Q);
  if (Q.ndim < 2) throw new Error('qr_lstsq(Q,R,y): Q.ndim must be at least 2.');
  R = (0, _nd_array.asarray)(R);
  if (R.ndim < 2) throw new Error('qr_lstsq(Q,R,y): R.ndim must be at least 2.');
  y = (0, _nd_array.asarray)(y);
  if (y.ndim < 2) throw new Error('qr_lstsq(Q,R,y): y.ndim must be at least 2.'); //  ________________   ______                   ___________
  // |                | |\(MxI)|                 |           |
  // |                | | \ R  |  ___________    |           |
  // |     (NxM)      | |  \   | |   (IxJ)   |   |   (NxJ)   |
  // |       Q        | |   \  | |     X     | = |     Y     |
  // |                | |    \ | |           |   |           |
  // |                | |     \| |___________|   |           |
  // |                | |  0   |                 |           |
  // |________________| |______|                 |___________|

  var _Q$shape$slice = Q.shape.slice(-2),
      _Q$shape$slice2 = (0, _slicedToArray2["default"])(_Q$shape$slice, 2),
      N = _Q$shape$slice2[0],
      M = _Q$shape$slice2[1],
      _R$shape$slice = R.shape.slice(-1),
      _R$shape$slice2 = (0, _slicedToArray2["default"])(_R$shape$slice, 1),
      I = _R$shape$slice2[0],
      _y$shape$slice = y.shape.slice(-1),
      _y$shape$slice2 = (0, _slicedToArray2["default"])(_y$shape$slice, 1),
      J = _y$shape$slice2[0],
      L = Math.min(M, I);

  if (N != y.shape[y.ndim - 2]) throw new Error("qr_lstsq(Q,R,y): Q and y don't match.");
  if (M != R.shape[R.ndim - 2]) throw new Error("qr_lstsq(Q,R,y): Q and R don't match.");
  if (I > N) throw new Error("qr_lstsq(Q,R,y): Under-determined systems not supported. Use rrqr instead.");
  var ndim = Math.max(Q.ndim, R.ndim, y.ndim),
      shape = Int32Array.from({
    length: ndim
  }, function () {
    return 1;
  });
  shape[ndim - 2] = I;
  shape[ndim - 1] = J; // FIND COMMON (BROADCASTED) SHAPE

  for (var _i4 = 0, _arr = [Q, R, y]; _i4 < _arr.length; _i4++) {
    var arr = _arr[_i4];

    for (var i = ndim - 2, j = arr.ndim - 2; i-- > 0 && j-- > 0;) {
      if (1 === shape[i]) shape[i] = arr.shape[j];else if (shape[i] != arr.shape[j] && arr.shape[j] != 1) throw new Error('Q, R, y are not broadcast-compatible.');
    }
  } // GENERATE RESULT DATA


  var DTypeArray = _dt.ARRAY_TYPES[[Q, R, y].every(function (a) {
    return a.dtype === 'float32';
  }) ? 'float32' : 'float64'],
      x_dat = new DTypeArray(shape.reduce(function (a, b) {
    return a * b;
  }, 1)),
      Q_dat = Q.data,
      R_dat = R.data,
      y_dat = y.data;
  var Q_off = 0,
      Q_stride = 1,
      R_off = 0,
      R_stride = 1,
      y_off = 0,
      y_stride = 1,
      x_off = 0;

  function solv(d) {
    if (d === ndim - 2) {
      Q_stride = N * M;
      R_stride = M * I;
      y_stride = N * J; // Q.T @ y

      for (var _i5 = 0; _i5 < L; _i5++) {
        for (var _j3 = 0; _j3 < J; _j3++) {
          for (var k = 0; k < N; k++) {
            x_dat[x_off + _i5 * J + _j3] += Q_dat[Q_off + k * M + _i5] * y_dat[y_off + k * J + _j3];
          }
        }
      }

      (0, _tri._triu_solve)(L, I, J, R_dat, R_off, x_dat, x_off);
      Q_off += Q_stride;
      R_off += R_stride;
      y_off += y_stride;
      x_off += I * J;
      return;
    }

    for (var l = shape[d];; l--) {
      solv(d + 1);
      if (l == 1) break;
      if (!(Q.shape[d - ndim + Q.ndim] > 1)) Q_off -= Q_stride;
      if (!(R.shape[d - ndim + R.ndim] > 1)) R_off -= R_stride;
      if (!(y.shape[d - ndim + y.ndim] > 1)) y_off -= y_stride;
    }

    Q_stride *= Q.shape[d - ndim + Q.ndim] || 1;
    R_stride *= R.shape[d - ndim + R.ndim] || 1;
    y_stride *= y.shape[d - ndim + y.ndim] || 1;
  }

  solv(0);
  return new _nd_array.NDArray(shape, x_dat);
}