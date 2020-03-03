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
exports.rand_ortho = rand_ortho;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _dt = require("../dt");

var _nd_array = require("../nd_array");

var _rand_normal = require("../rand_normal");

var _norm = require("./norm");

function rand_ortho(dtype) {
  for (var _len = arguments.length, shape = new Array(_len > 1 ? _len - 1 : 0), _key = 1; _key < _len; _key++) {
    shape[_key - 1] = arguments[_key];
  }

  // REFERENCES:
  //   - https://blogs.sas.com/content/iml/2012/03/28/generating-a-random-orthogonal-matrix.html
  if (!(dtype in _dt.ARRAY_TYPES)) {
    shape.unshift(dtype);
    dtype = 'float64';
  }

  if (!dtype.startsWith('float')) throw new Error("Unsupported dtype: '".concat(dtype, "'."));
  if (shape.length === 1) shape.push(shape[0]);

  for (var i = shape.length; i-- > 0;) {
    if (shape[i] % 1 !== 0 || shape[i] < 1) throw new Error("shape[".concat(i, "] = ").concat(shape[i], " not a positive integer."));
  }

  shape = Int32Array.from(shape);

  var _shape$slice = shape.slice(-2),
      _shape$slice2 = (0, _slicedToArray2["default"])(_shape$slice, 2),
      M = _shape$slice2[0],
      N = _shape$slice2[1],
      L = Math.max(M, N);

  var DTypeArray = _dt.ARRAY_TYPES[dtype],
      U = new DTypeArray(shape.reduce(function (x, y) {
    return x * y;
  })),
      Q = new DTypeArray(L * L),
      x = new DTypeArray(L),
      d = new DTypeArray(L);
  var NORM = new _norm.FrobeniusNorm(); // TODO: this should be computable more efficiently for non-square matrices

  for (var off = 0; off < U.length; off += M * N) {
    // INIT Q TO IDENTITY
    for (var _i = 0; _i < L; _i++) {
      for (var j = 0; j < L; j++) {
        Q[L * _i + j] = _i === j;
      }
    }

    d[0] = Math.random() <= 0.5 ? -1 : +1; // apply series of random householder transformations/reflections

    for (var k = L; k-- > 1;) {
      NORM.reset();

      for (var _i2 = 0;; _i2++) {
        var x_i = x[_i2] = (0, _rand_normal.rand_normal)();
        if (_i2 >= k) break;
        NORM.include(x_i);
      }

      d[k] = x[k] > 0 ? -1 : +1;
      var norm = NORM.resultIncl(x[k]) * d[k];
      if (0 === norm) continue;
      NORM.include(x[k] -= norm);
      var max = NORM.max,
          div = Math.sqrt(NORM.sum);

      for (var _i3 = k; _i3 >= 0; _i3--) {
        x[_i3] = x[_i3] / max / div;
      } // apply householder to right of Q


      for (var _i4 = L; _i4-- > 0;) {
        var sum = 0;

        for (var _j = k; _j >= 0; _j--) {
          sum += Q[L * _i4 + _j] * x[_j];
        }

        sum *= 2;

        for (var _j2 = k; _j2 >= 0; _j2--) {
          Q[L * _i4 + _j2] -= x[_j2] * sum;
        }
      }
    } // flip row signs according to d


    for (var _i5 = 0; _i5 < L; _i5++) {
      for (var _j3 = 0; _j3 < L; _j3++) {
        Q[L * _i5 + _j3] *= d[_i5];
      }
    } // MOVE Q -> U


    for (var _i6 = 0; _i6 < M; _i6++) {
      for (var _j4 = 0; _j4 < N; _j4++) {
        U[off + N * _i6 + _j4] = Q[L * _i6 + _j4];
      }
    }
  }

  return new _nd_array.NDArray(shape, U);
}