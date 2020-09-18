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

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports._hessenberg_decomp = _hessenberg_decomp;
exports.hessenberg_decomp = hessenberg_decomp;

var _nd_array = require("../nd_array");

var _dt = require("../dt");

var _norm = require("./norm");

// TODO: bidiag_solve
function _hessenberg_decomp(N, U, H, off) {
  N |= 0;
  off |= 0; // INIT U TO IDENTITY

  for (var i = N - 1; i-- > 0;) {
    U[off + N * i + i] = 1.0;
  }

  var NORM = new _norm.FrobeniusNorm(),
      lastRow = off + N * (N - 1); // ELIMINATE ELEMENTS LEFT OF (i,i-1) VIA HOUSEHOLDER

  for (var _i = N; --_i > 1;) {
    var rowI = off + N * _i,
        ii = rowI + (_i - 1); // COMPUTE HOUSEHOLDER VECTOR

    NORM.reset();

    for (var j = _i - 1; j-- > 0;) {
      NORM.include(H[rowI + j]);
    }

    if (NORM.max === 0) continue;
    var norm = NORM.resultIncl(H[ii]) * (H[ii] > 0 ? -1 : +1); // <- avoid cancellation error

    NORM.include(H[ii] -= norm);
    var max = NORM.max,
        div = Math.sqrt(NORM.sum);

    for (var _j = _i; _j-- > 0;) {
      H[rowI + _j] = H[rowI + _j] / max * Math.SQRT2 / div;
    } // APPLY HOUSEHOLDER TO RIGHT OF H


    for (var _j2 = _i; _j2-- > 0;) {
      var sum = 0;

      for (var k = _i; k-- > 0;) {
        sum += H[off + N * _j2 + k] * H[rowI + k];
      }

      for (var _k = _i; _k-- > 0;) {
        H[off + N * _j2 + _k] -= H[rowI + _k] * sum;
      }
    } // APPLY HOUSEHOLDER TO LEFT OF H


    U.fill(0.0, lastRow, off + N * N);

    for (var _j3 = _i; _j3-- > 0;) {
      for (var _k2 = N; _k2-- > 0;) {
        U[lastRow + _k2] += H[off + N * _j3 + _k2] * H[rowI + _j3];
      }
    }

    for (var _j4 = _i; _j4-- > 0;) {
      for (var _k3 = N; _k3-- > 0;) {
        H[off + N * _j4 + _k3] -= H[rowI + _j4] * U[lastRow + _k3];
      }
    } // APPLY HOUSEHOLDER TO RIGHT OF U


    for (var _j5 = N - 1; _j5-- > 0;) {
      var _sum = 0;

      for (var _k4 = _i; _k4-- > 0;) {
        _sum += U[off + N * _j5 + _k4] * H[rowI + _k4];
      }

      for (var _k5 = _i; _k5-- > 0;) {
        U[off + N * _j5 + _k5] -= H[rowI + _k5] * _sum;
      }
    } // FINISH ROW i


    H.fill(0.0, rowI, ii);
    H[ii] = norm;
  } // LAST ROW OF U WAS USED AS TEMP MEMORY


  U.fill(0.0, lastRow, off + N * N - 1);
  U[off + N * N - 1] = 1;
}

function hessenberg_decomp(A) {
  A = (0, _nd_array.asarray)(A);
  if (A.ndim < 2) throw new Error('hessenberg_decomp(A): A must at least be 2D.');
  var DType = A.dtype === 'float32' ? 'float32' : 'float64',
      DTypeArray = _dt.ARRAY_TYPES[A.dtype === 'float32' ? 'float32' : 'float64'],
      // <- ensure at least double precision
  shape = A.shape,
      N = shape[shape.length - 1] | 0;
  if (N != shape[shape.length - 2]) throw new Error('hessenberg_decomp(A): A must be square.');
  var H = DTypeArray.from(A.data);
  A = undefined;
  var U = new DTypeArray(H.length);

  for (var off = H.length; (off -= N * N) >= 0;) {
    _hessenberg_decomp(N, U, H, off);
  }

  return [new _nd_array.NDArray(shape, U), new _nd_array.NDArray(shape, H)];
}

;