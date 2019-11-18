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
exports._transpose_inplace = _transpose_inplace;
exports.transpose_inplace = transpose_inplace;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _nd_array = require("../nd_array");

function _transpose_inplace(N, A, A_off) {
  for (var i = N; --i > 0;) {
    for (var j = i; j-- > 0;) {
      var ij = A_off + N * i + j,
          ji = A_off + N * j + i,
          A_ij = A[ij];
      A[ij] = A[ji];
      A[ji] = A_ij;
    }
  }
}

function transpose_inplace(A) {
  var _A$shape$slice = A.shape.slice(-2),
      _A$shape$slice2 = (0, _slicedToArray2["default"])(_A$shape$slice, 2),
      N = _A$shape$slice2[0],
      M = _A$shape$slice2[1];

  if (N != M) throw new Error('In-place transposition is only supported for square matrices.');
  A = A.data;

  for (var A_off = 0; A_off < A.length; A_off += N * N) {
    _transpose_inplace(N, A, A_off);
  }
}