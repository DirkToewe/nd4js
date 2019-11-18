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
exports.diag_mat = diag_mat;
exports.diag = diag;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _nd_array = require("../nd_array");

var _dt = require("../dt");

function diag_mat(diag) {
  if (diag.ndim < 1) throw new Error("diag_mat(diag): diag.ndim=".concat(diag.ndim, " is less than 1."));
  diag = (0, _nd_array.asarray)(diag);

  var shape = diag.shape,
      N = shape[shape.length - 1],
      d = diag.data,
      D = _dt.ARRAY_TYPES[diag.dtype].from({
    length: shape.reduce(function (a, b) {
      return a * b;
    }) * N
  }, function () {
    return 0;
  });

  diag = undefined;
  if (shape.length < 1) throw new Error('diag_mat(diag): diag must be at least 1d.');
  if (N <= 0) throw new Error('Assertion Failed!');

  for (var d_off = 0, D_off = 0;;) {
    var d_end = d_off + N;

    while (true) {
      D[D_off] = d[d_off];
      if (++d_off === d_end) break;
      D_off += N + 1;
    }

    if (++D_off == D.length) break;
  }

  return new _nd_array.NDArray(Int32Array.of.apply(Int32Array, (0, _toConsumableArray2["default"])(shape).concat([N])), D);
}

function diag(A) {
  var offset = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0;
  if (A.ndim < 2) throw new Error("diag(A, offset): A.ndim=".concat(A.ndim, " is less than 2."));
  A = (0, _nd_array.asarray)(A);

  var _A$shape$slice = A.shape.slice(-2),
      _A$shape$slice2 = (0, _slicedToArray2["default"])(_A$shape$slice, 2),
      N = _A$shape$slice2[0],
      M = _A$shape$slice2[1],
      shape = A.shape.slice(0, -1),
      DTypeArray = _dt.ARRAY_TYPES[A.dtype];

  A = A.data;
  if (!(offset > -N && offset < +M)) throw new Error("diag(A, offset): offset=".concat(offset, " out of A.shape=[").concat(A.shape, "]."));
  var L = Math.min(N, N + offset, M, M - offset);
  shape[shape.length - 1] = L;
  var d = new DTypeArray(shape.reduce(function (a, b) {
    return a * b;
  })),
      a_off = Math.max(0, offset, -offset * M);

  for (var d_off = 0, A_off = 0; A_off < A.length; d_off += L, A_off += N * M) {
    for (var i = 0; i < L; i++) {
      d[d_off + i] = A[A_off + a_off + M * i + i];
    }
  }

  return new _nd_array.NDArray(shape, d);
}