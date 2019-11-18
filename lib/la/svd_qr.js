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
exports._svd_qr_bidiag = _svd_qr_bidiag;
exports._svd_qr = _svd_qr;
exports.svd_qr = svd_qr;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _dt = require("../dt");

var _nd_array = require("../nd_array");

var _transpose_inplace = require("./transpose_inplace");

var _giv_rot = require("./_giv_rot");

var _svd_jac_utils = require("./_svd_jac_utils");

var _norm = require("./norm");

var _bidiag = require("./bidiag");

function _svd_qr_bidiag(N, n, U, U_off, B, B_off, V, V_off) {
  throw new Error('Not yet implemented.');
}

function _svd_qr(M, N, U, U_off, sv, sv_off, V, V_off) {
  throw new Error('Not yet implemented.');
}

function svd_qr(A) {
  // SEE:
  //  - "INTRODUCTION OF DOUBLE DIVIDE AND CONQUER AND THE RECENT PROGRESS", TARO KONDA & YOSHIMASA NAKAMURA
  //  - "A Divide-and-Conquer Approach for Solving Singular Value Decomposition on a Heterogeneous System", Ding Liu & Ruixuan Li & David J. Lilja & Weijun Xiao
  A = (0, _nd_array.asarray)(A);
  if (A.dtype.startsWith('complex')) throw new Error('svd_dc(A): A.dtype must be float.');

  var _A$shape$slice = A.shape.slice(-2),
      _A$shape$slice2 = (0, _slicedToArray2["default"])(_A$shape$slice, 2),
      M = _A$shape$slice2[0],
      N = _A$shape$slice2[1];

  if (M > N) {
    var _svd_dc2 = svd_dc(A.T),
        _svd_dc3 = (0, _slicedToArray2["default"])(_svd_dc2, 3),
        _U = _svd_dc3[0],
        _sv = _svd_dc3[1],
        _V = _svd_dc3[2];

    (0, _transpose_inplace.transpose_inplace)(_U);
    return [_V.T, _sv, _U];
  }

  var V_shape = A.shape,
      U_shape = V_shape.slice(),
      sv_shape = V_shape.slice(0, -1);
  U_shape[U_shape.length - 1] = M;
  var DType = A.dtype === 'float32' ? 'float32' : 'float64',
      DTypeArray = _dt.ARRAY_TYPES[DType];
  A = A.data;
  var len = A.length / (M * N),
      V = A.slice();
  A = undefined;
  var U = new DTypeArray(len * M * M),
      sv = new DTypeArray(len * M),
      F = new DTypeArray(M * (M + 2)
  /*F*/
  + M * 2
  /*B*/
  + (M + 1) * (M + 1)
  /*V1*/
  + (M + 1) * Math.max(1 + M, N)
  /*V2*/
  ),
      I = new Int32Array(M * 3);

  for (var U_off = 0, sv_off = 0, V_off = 0; sv_off < sv.length; U_off += M * M, sv_off += M, V_off += M * N) {
    _svd_dc(M, N, U, U_off, sv, sv_off, V, V_off, I, F);
  }

  return [new _nd_array.NDArray(U_shape, U), new _nd_array.NDArray(sv_shape, sv), new _nd_array.NDArray(V_shape, V)];
}