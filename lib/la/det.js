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
exports.slogdet = exports.det = exports.slogdet_tri = exports.det_tri = void 0;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _nd_array = require("../nd_array");

var _dt = require("../dt");

var _qr = require("./qr");

var det_tri = function det_tri(a) {
  a = (0, _nd_array.asarray)(a);
  if (a.ndim < 2) throw new Error("det_tri(a): a.shape=[".concat(a.shape, "]; a.ndim must be at least 2."));

  var a_shape = a.shape,
      _a_shape$slice = a_shape.slice(-2),
      _a_shape$slice2 = (0, _slicedToArray2["default"])(_a_shape$slice, 2),
      M = _a_shape$slice2[0],
      N = _a_shape$slice2[1];

  if (M !== N) throw new Error("det_tri(a): a must be square matrices.");
  var DTypeArray = _dt.ARRAY_TYPES[a.dtype];
  var A = a.data;
  a = undefined;
  var D = new DTypeArray(A.length / (N * N));

  for (var D_off = 0, A_off = 0; D_off < D.length; D_off++) {
    var i = A_off;
    A_off += N * N;
    var d = 1;

    for (; i < A_off; i += N + 1) {
      d *= A[i];
    }

    D[D_off] = d;
  }

  return new _nd_array.NDArray(a_shape.slice(0, -2), D);
};

exports.det_tri = det_tri;

var slogdet_tri = function slogdet_tri(A) {
  A = (0, _nd_array.asarray)(A);
  if (A.ndim < 2) throw new Error("det_tri(A): A.ndim must be at least 2.");

  var A_shape = A.shape,
      _A_shape$slice = A_shape.slice(-2),
      _A_shape$slice2 = (0, _slicedToArray2["default"])(_A_shape$slice, 2),
      M = _A_shape$slice2[0],
      N = _A_shape$slice2[1];

  if (M !== N) throw new Error("det_tri(A): A must be square matrices.");
  var DTypeArray = _dt.ARRAY_TYPES[A.dtype === 'float32' ? 'float32' : 'float64'];
  A = A.data;
  var S = new DTypeArray(A.length / (N * N)),
      D = new DTypeArray(S.length);

  for (var SD_off = 0, A_off = 0; SD_off < D.length; SD_off++) {
    var i = A_off;
    A_off += N * N;
    var s = 1,
        d = 0;

    for (; i < A_off; i += N + 1) {
      s *= Math.sign(A[i]);
      d += Math.log(Math.abs(A[i]));
    }

    S[SD_off] = s;
    D[SD_off] = d;
  }

  var SD_shape = A_shape.slice(0, -2);
  return [new _nd_array.NDArray(SD_shape, S), new _nd_array.NDArray(SD_shape, D)];
};

exports.slogdet_tri = slogdet_tri;

var det = function det(A) {
  var _qr_decomp = (0, _qr.qr_decomp)(A),
      _qr_decomp2 = (0, _slicedToArray2["default"])(_qr_decomp, 2),
      Q = _qr_decomp2[0],
      R = _qr_decomp2[1];

  return det_tri(R);
};

exports.det = det;

var slogdet = function slogdet(A) {
  var _qr_decomp3 = (0, _qr.qr_decomp)(A),
      _qr_decomp4 = (0, _slicedToArray2["default"])(_qr_decomp3, 2),
      Q = _qr_decomp4[0],
      R = _qr_decomp4[1];

  return slogdet_tri(R);
};

exports.slogdet = slogdet;