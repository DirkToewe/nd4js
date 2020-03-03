'use strict';
/* This file is part of ND4JS.
 *
 * ND4JS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ND4JS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ND4JS. If not, see <http://www.gnu.org/licenses/>.
 */

var _interopRequireDefault = require("@babel/runtime/helpers/interopRequireDefault");

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.fit_lin = fit_lin;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _nd_array = require("../nd_array");

var _lstsq2 = require("../la/lstsq");

function fit_lin(x, y, regularization, funcs) {
  if (funcs == null) {
    funcs = regularization;
    regularization = 0;
  } else if (null == regularization) regularization = 0;else if (!(0 <= regularization && regularization < Infinity)) throw new Error("fit_lin(x,y, regularization, funcs): Invalid regularization: ".concat(regularization, "."));

  funcs = (0, _toConsumableArray2["default"])(funcs);
  x = (0, _nd_array.asarray)('float64', x);
  y = (0, _nd_array.asarray)('float64', y);
  if (y.ndim !== 1) throw new Error('fit_lin(x,y, (regularization,) funcs): y.ndim must be 1.');
  var M = y.shape[0],
      N = funcs.length;
  var x_ndim = x.ndim;
  if (x.ndim !== 1 && x.ndim !== 2) throw new Error('fit_lin(x,y, (regularization,) funcs): x.ndim must be 1 or 2.');
  if (x.shape[0] !== y.shape[0]) throw new Error('fit_lin(x,y, (regularization,) funcs): x.shape[0] and y.shape[0] must equal.');
  if (x.ndim === 1) x = x.data;else {
    var x_shape = x.shape.slice(1),
        _x_shape = (0, _slicedToArray2["default"])(x_shape, 1),
        L = _x_shape[0];

    x = x.data.slice();
    Object.freeze(x.buffer);
    x = Array.from({
      length: M
    }, function (_, i) {
      return Object.freeze(new _nd_array.NDArray(x_shape, x.subarray(L * i, L * (i + 1))));
    });
  }
  var A_shape = Int32Array.of((regularization !== 0) * N + M, N);
  var A = new Float64Array(A_shape[0] * A_shape[1]),
      z = y;

  if (regularization !== 0) {
    regularization = Math.sqrt(regularization);

    for (var i = N; i-- > 0;) {
      A[M * N + N * i + i] = regularization;
    }
  }

  for (var _i = M; _i-- > 0;) {
    for (var j = N; j-- > 0;) {
      A[N * _i + j] = funcs[j](x[_i]);
    }
  }

  if (regularization !== 0) {
    y = y.data;

    var _z = new Float64Array(M + N);

    for (var _i2 = M; _i2-- > 0;) {
      _z[_i2] = y[_i2];
    }

    y = new _nd_array.NDArray(Int32Array.of(M + N, 1), _z);
  } else y = new _nd_array.NDArray(Int32Array.of(M, 1), y.data);

  A = new _nd_array.NDArray(A_shape, A);

  var _lstsq = (0, _lstsq2.lstsq)(A, y),
      coeffs = _lstsq.data;

  var param_lin_func = function param_lin_func(x) {
    var coeffs = param_lin_func.coeffs,
        funcs = param_lin_func.funcs;
    if (coeffs.length !== funcs.length) throw new Error('Assertion failed.');
    var result = 0;

    for (var _i3 = coeffs.length; _i3-- > 0;) {
      result += coeffs[_i3] * funcs[_i3](x);
    }

    return result;
  };

  param_lin_func.coeffs = coeffs;
  param_lin_func.funcs = funcs;
  return param_lin_func;
}