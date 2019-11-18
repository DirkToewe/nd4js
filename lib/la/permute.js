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
exports.permute_rows = permute_rows;
exports.permute_cols = permute_cols;
exports.unpermute_rows = unpermute_rows;
exports.unpermute_cols = unpermute_cols;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _nd_array = require("../nd_array");

var _dt = require("../dt");

function permute_rows(A, P) {
  A = (0, _nd_array.asarray)(A);
  if (A.ndim < 2) throw new Error('permute_rows(A,P): A.ndim must be at least 2.');
  P = (0, _nd_array.asarray)(P);
  if (P.ndim < 1) throw new Error('permute_rows(A,P): P.ndim must be at least 1.');
  if (P.dtype !== 'int32') throw new Error('permute_rows(A,P): P.dtype must be "int32".');

  var ndim = Math.max(A.ndim, P.ndim + 1),
      B_shape = Int32Array.from({
    length: ndim
  }, function () {
    return 1;
  }),
      _A$shape$slice = A.shape.slice(-2),
      _A$shape$slice2 = (0, _slicedToArray2["default"])(_A$shape$slice, 2),
      M = _A$shape$slice2[0],
      N = _A$shape$slice2[1];

  if (M !== P.shape[P.ndim - 1]) throw new Error('permute_rows(A,P): A.shape[-2] and P.shape[-1] must match.');

  for (var i = ndim, j = A.ndim; i-- > 0 && j-- > 0;) {
    B_shape[i] = A.shape[j];
  } // FIND COMMON (BROADCASTED) SHAPE


  for (var _i = ndim - 2, _j = P.ndim - 1; _i-- > 0 && _j-- > 0;) {
    if (1 === B_shape[_i]) B_shape[_i] = P.shape[_j];else if (B_shape[_i] != P.shape[_j] && P.shape[_j] != 1) throw new Error('permute_rows(A,P): A and P not broadcast-compatible.');
  } // GENERATE RESULT DATA


  var B = new _dt.ARRAY_TYPES[A.dtype](B_shape.reduce(function (a, b) {
    return a * b;
  })),
      A_shape = A.shape,
      A_ndim = A.ndim,
      P_shape = P.shape,
      P_ndim = P.ndim;
  A = A.data;
  P = P.data;
  var A_off = 0,
      A_stride = -1,
      P_off = 0,
      P_stride = -1,
      B_off = 0;

  function perm(d) {
    if (d === ndim - 2) {
      // Q.T @ y
      for (var _i2 = 0; _i2 < M; _i2++) {
        var k = P[P_off + _i2];
        if (!(k >= 0)) throw new Error('permute_rows(A,P): P.contains invalid indices.');
        if (!(k < M)) throw new Error('permute_rows(A,P): P.contains invalid indices.');

        for (var _j2 = 0; _j2 < N; _j2++) {
          B[B_off + N * _i2 + _j2] = A[A_off + N * k + _j2];
        }
      }

      A_off += A_stride = M * N;
      P_off += P_stride = M;
      B_off += A_stride;
      return;
    }

    for (var l = B_shape[d];; l--) {
      perm(d + 1);
      if (l === 1) break;
      if (!(A_shape[d - ndim + A_ndim] > 1)) A_off -= A_stride;
      if (!(P_shape[d - ndim + 1 + P_ndim] > 1)) P_off -= P_stride;
    }

    A_stride *= A_shape[d - ndim + A_ndim] || 1;
    P_stride *= P_shape[d - ndim + 1 + P_ndim] || 1;
  }

  perm(0);
  return new _nd_array.NDArray(B_shape, B);
}

function permute_cols(A, P) {
  A = (0, _nd_array.asarray)(A);
  if (A.ndim < 2) throw new Error('permute_cols(A,P): A.ndim must be at least 2.');
  P = (0, _nd_array.asarray)(P);
  if (P.ndim < 1) throw new Error('permute_cols(A,P): P.ndim must be at least 1.');
  if (P.dtype !== 'int32') throw new Error('permute_cols(A,P): P.dtype must be "int32".');

  var ndim = Math.max(A.ndim, P.ndim + 1),
      B_shape = Int32Array.from({
    length: ndim
  }, function () {
    return 1;
  }),
      _A$shape$slice3 = A.shape.slice(-2),
      _A$shape$slice4 = (0, _slicedToArray2["default"])(_A$shape$slice3, 2),
      M = _A$shape$slice4[0],
      N = _A$shape$slice4[1];

  if (N !== P.shape[P.ndim - 1]) throw new Error('permute_cols(A,P): A.shape[-1] and P.shape[-1] must match.');

  for (var i = ndim, j = A.ndim; i-- > 0 && j-- > 0;) {
    B_shape[i] = A.shape[j];
  } // FIND COMMON (BROADCASTED) SHAPE


  for (var _i3 = ndim - 2, _j3 = P.ndim - 1; _i3-- > 0 && _j3-- > 0;) {
    if (1 === B_shape[_i3]) B_shape[_i3] = P.shape[_j3];else if (B_shape[_i3] != P.shape[_j3] && P.shape[_j3] != 1) throw new Error('permute_cols(A,P): A and P not broadcast-compatible.');
  } // GENERATE RESULT DATA


  var B = new _dt.ARRAY_TYPES[A.dtype](B_shape.reduce(function (a, b) {
    return a * b;
  })),
      A_shape = A.shape,
      A_ndim = A.ndim,
      P_shape = P.shape,
      P_ndim = P.ndim;
  A = A.data;
  P = P.data;
  var A_off = 0,
      A_stride = -1,
      P_off = 0,
      P_stride = -1,
      B_off = 0;

  function perm(d) {
    if (d === ndim - 2) {
      // Q.T @ y
      for (var _i4 = 0; _i4 < M; _i4++) {
        for (var _j4 = 0; _j4 < N; _j4++) {
          var k = P[P_off + _j4];
          if (!(k >= 0)) throw new Error('permute_cols(A,P): P.contains invalid indices.');
          if (!(k < N)) throw new Error('permute_cols(A,P): P.contains invalid indices.');
          B[B_off + N * _i4 + _j4] = A[A_off + N * _i4 + k];
        }
      }

      A_off += A_stride = M * N;
      P_off += P_stride = N;
      B_off += A_stride;
      return;
    }

    for (var l = B_shape[d];; l--) {
      perm(d + 1);
      if (l === 1) break;
      if (!(A_shape[d - ndim + A_ndim] > 1)) A_off -= A_stride;
      if (!(P_shape[d - ndim + 1 + P_ndim] > 1)) P_off -= P_stride;
    }

    A_stride *= A_shape[d - ndim + A_ndim] || 1;
    P_stride *= P_shape[d - ndim + 1 + P_ndim] || 1;
  }

  perm(0);
  return new _nd_array.NDArray(B_shape, B);
}

function unpermute_rows(A, P) {
  A = (0, _nd_array.asarray)(A);
  if (A.ndim < 2) throw new Error('unpermute_rows(A,P): A.ndim must be at least 2.');
  P = (0, _nd_array.asarray)(P);
  if (P.ndim < 1) throw new Error('unpermute_rows(A,P): P.ndim must be at least 1.');
  if (P.dtype !== 'int32') throw new Error('unpermute_rows(A,P): P.dtype must be "int32".');

  var ndim = Math.max(A.ndim, P.ndim + 1),
      B_shape = Int32Array.from({
    length: ndim
  }, function () {
    return 1;
  }),
      _A$shape$slice5 = A.shape.slice(-2),
      _A$shape$slice6 = (0, _slicedToArray2["default"])(_A$shape$slice5, 2),
      M = _A$shape$slice6[0],
      N = _A$shape$slice6[1];

  if (M !== P.shape[P.ndim - 1]) throw new Error('unpermute_rows(A,P): A.shape[-2] and P.shape[-1] must match.');

  for (var i = ndim, j = A.ndim; i-- > 0 && j-- > 0;) {
    B_shape[i] = A.shape[j];
  } // FIND COMMON (BROADCASTED) SHAPE


  for (var _i5 = ndim - 2, _j5 = P.ndim - 1; _i5-- > 0 && _j5-- > 0;) {
    if (1 === B_shape[_i5]) B_shape[_i5] = P.shape[_j5];else if (B_shape[_i5] != P.shape[_j5] && P.shape[_j5] != 1) throw new Error('unpermute_rows(A,P): A and P not broadcast-compatible.');
  } // GENERATE RESULT DATA


  var B = new _dt.ARRAY_TYPES[A.dtype](B_shape.reduce(function (a, b) {
    return a * b;
  })),
      A_shape = A.shape,
      A_ndim = A.ndim,
      P_shape = P.shape,
      P_ndim = P.ndim;
  A = A.data;
  P = P.data;
  var A_off = 0,
      A_stride = -1,
      P_off = 0,
      P_stride = -1,
      B_off = 0;

  function perm(d) {
    if (d === ndim - 2) {
      // Q.T @ y
      for (var _i6 = 0; _i6 < M; _i6++) {
        var k = P[P_off + _i6];
        if (!(k >= 0)) throw new Error('unpermute_rows(A,P): P.contains invalid indices.');
        if (!(k < M)) throw new Error('unpermute_rows(A,P): P.contains invalid indices.');

        for (var _j6 = 0; _j6 < N; _j6++) {
          B[B_off + N * k + _j6] = A[A_off + N * _i6 + _j6];
        }
      }

      A_off += A_stride = M * N;
      P_off += P_stride = M;
      B_off += A_stride;
      return;
    }

    for (var l = B_shape[d];; l--) {
      perm(d + 1);
      if (l === 1) break;
      if (!(A_shape[d - ndim + A_ndim] > 1)) A_off -= A_stride;
      if (!(P_shape[d - ndim + 1 + P_ndim] > 1)) P_off -= P_stride;
    }

    A_stride *= A_shape[d - ndim + A_ndim] || 1;
    P_stride *= P_shape[d - ndim + 1 + P_ndim] || 1;
  }

  perm(0);
  return new _nd_array.NDArray(B_shape, B);
}

function unpermute_cols(A, P) {
  A = (0, _nd_array.asarray)(A);
  if (A.ndim < 2) throw new Error('unpermute_cols(A,P): A.ndim must be at least 2.');
  P = (0, _nd_array.asarray)(P);
  if (P.ndim < 1) throw new Error('unpermute_cols(A,P): P.ndim must be at least 1.');
  if (P.dtype !== 'int32') throw new Error('unpermute_cols(A,P): P.dtype must be "int32".');

  var ndim = Math.max(A.ndim, P.ndim + 1),
      B_shape = Int32Array.from({
    length: ndim
  }, function () {
    return 1;
  }),
      _A$shape$slice7 = A.shape.slice(-2),
      _A$shape$slice8 = (0, _slicedToArray2["default"])(_A$shape$slice7, 2),
      M = _A$shape$slice8[0],
      N = _A$shape$slice8[1];

  if (N !== P.shape[P.ndim - 1]) throw new Error('unpermute_cols(A,P): A.shape[-1] and P.shape[-1] must match.');

  for (var i = ndim, j = A.ndim; i-- > 0 && j-- > 0;) {
    B_shape[i] = A.shape[j];
  } // FIND COMMON (BROADCASTED) SHAPE


  for (var _i7 = ndim - 2, _j7 = P.ndim - 1; _i7-- > 0 && _j7-- > 0;) {
    if (1 === B_shape[_i7]) B_shape[_i7] = P.shape[_j7];else if (B_shape[_i7] != P.shape[_j7] && P.shape[_j7] != 1) throw new Error('unpermute_cols(A,P): A and P not broadcast-compatible.');
  } // GENERATE RESULT DATA


  var B = new _dt.ARRAY_TYPES[A.dtype](B_shape.reduce(function (a, b) {
    return a * b;
  })),
      A_shape = A.shape,
      A_ndim = A.ndim,
      P_shape = P.shape,
      P_ndim = P.ndim;
  A = A.data;
  P = P.data;
  var A_off = 0,
      A_stride = -1,
      P_off = 0,
      P_stride = -1,
      B_off = 0;

  function perm(d) {
    if (d === ndim - 2) {
      // Q.T @ y
      for (var _i8 = 0; _i8 < M; _i8++) {
        for (var _j8 = 0; _j8 < N; _j8++) {
          var k = P[P_off + _j8];
          if (!(k >= 0)) throw new Error('unpermute_cols(A,P): P.contains invalid indices.');
          if (!(k < N)) throw new Error('unpermute_cols(A,P): P.contains invalid indices.');
          B[B_off + N * _i8 + k] = A[A_off + N * _i8 + _j8];
        }
      }

      A_off += A_stride = M * N;
      P_off += P_stride = N;
      B_off += A_stride;
      return;
    }

    for (var l = B_shape[d];; l--) {
      perm(d + 1);
      if (l === 1) break;
      if (!(A_shape[d - ndim + A_ndim] > 1)) A_off -= A_stride;
      if (!(P_shape[d - ndim + 1 + P_ndim] > 1)) P_off -= P_stride;
    }

    A_stride *= A_shape[d - ndim + A_ndim] || 1;
    P_stride *= P_shape[d - ndim + 1 + P_ndim] || 1;
  }

  perm(0);
  return new _nd_array.NDArray(B_shape, B);
}