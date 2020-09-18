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
exports.num_grad_forward = exports.num_grad = void 0;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _nd_array = require("../nd_array");

var _dt = require("../dt");

// https://en.wikipedia.org/wiki/Finite_difference#Forward,_backward,_and_central_differences
// https://www.geometrictools.com/Documentation/FiniteDifferences.pdf
var NUM_GRAD_D = Float64Array.of(+2, +1, -1, -2),
    NUM_GRAD_W = Float64Array.of(-1, 8, -8, +1),
    NUM_GRAD_S = 12;

var num_grad = function num_grad(f) {
  var _ref = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : {},
      _ref$h_rel = _ref.h_rel,
      h_rel = _ref$h_rel === void 0 ? undefined : _ref$h_rel,
      _ref$h_abs = _ref.h_abs,
      h_abs = _ref$h_abs === void 0 ? undefined : _ref$h_abs;

  return function (x) {
    x = (0, _nd_array.array)('float', x);
    var dtype = x.dtype === 'float32' ? 'float32' : 'float64',
        DTypeArray = _dt.ARRAY_TYPES[dtype],
        x_shape = x.shape;
    x = x.data;
    var g = null,
        g_shape = null;
    var epsRel = null != h_rel ? h_rel : Math.pow((0, _dt.eps)(dtype), 1 / 3),
        epsAbs = null != h_abs ? h_abs : Math.pow((0, _dt.eps)(dtype), 1 / 3);
    var h = Math.max(Math.abs(x[0]) * epsRel, epsAbs); // first function evaluation to determine output shape

    var y = function () {
      var X = new _nd_array.NDArray(x_shape, DTypeArray.from(x));
      X.data[0] += h * NUM_GRAD_D[0];
      return f(X);
    }();

    if ('number' === typeof y) {
      // SCALAR OUTPUT
      // -------------
      g_shape = x_shape;
      g = new Float64Array(x.length);
      g[0] = y * NUM_GRAD_W[0];

      outer_loop: for (var i = 0, j = 1;; j = 0) {
        for (; j < NUM_GRAD_D.length; j++) {
          var X = new _nd_array.NDArray(x_shape, DTypeArray.from(x));
          X.data[i] += h * NUM_GRAD_D[j];
          y = f(X);
          g[i] += y * NUM_GRAD_W[j];
        }

        g[i] /= h * NUM_GRAD_S;
        if (!(++i < x.length)) break outer_loop;
        h = Math.max(Math.abs(x[i]) * epsRel, epsAbs);
      }
    } else {
      // ND-ARRAY OUTPUT
      // ---------------
      y = (0, _nd_array.array)('float', y);
      var y_shape = y.shape;
      y = y.data;
      var M = y.length,
          N = x.length;
      g_shape = Int32Array.of.apply(Int32Array, (0, _toConsumableArray2["default"])(y_shape).concat((0, _toConsumableArray2["default"])(x_shape)));
      g = new Float64Array(y.length * x.length);

      for (var k = 0; k < M; k++) {
        g[N * k] = y[k] * NUM_GRAD_W[0];
      }

      outer_loop: for (var _i = 0, _j = 1;; _j = 0) {
        for (; _j < NUM_GRAD_D.length; _j++) {
          var _X = new _nd_array.NDArray(x_shape, DTypeArray.from(x));

          _X.data[_i] += h * NUM_GRAD_D[_j];
          y = (0, _nd_array.array)('float', f(_X));
          if (y.ndim !== y_shape.length) throw new Error('num_grad(f)(x): return shape of f must be consistent.');

          for (var _i2 = y_shape.length; _i2-- > 0;) {
            if (y.shape[_i2] !== y_shape[_i2]) throw new Error('num_grad(f)(x): return shape of f must be consistent.');
          }

          y = y.data;

          for (var _k = 0; _k < M; _k++) {
            g[N * _k + _i] += y[_k] * NUM_GRAD_W[_j];
          }
        }

        for (var _k2 = 0; _k2 < M; _k2++) {
          g[N * _k2 + _i] /= h * NUM_GRAD_S;
        }

        if (!(++_i < N)) break outer_loop;
        h = Math.max(Math.abs(x[_i]) * epsRel, epsAbs);
      }
    }

    return new _nd_array.NDArray(g_shape, g);
  };
};

exports.num_grad = num_grad;

var num_grad_forward = function num_grad_forward(f) {
  var _ref2 = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : {},
      _ref2$h_rel = _ref2.h_rel,
      h_rel = _ref2$h_rel === void 0 ? undefined : _ref2$h_rel,
      _ref2$h_abs = _ref2.h_abs,
      h_abs = _ref2$h_abs === void 0 ? undefined : _ref2$h_abs;

  return function (x) {
    x = (0, _nd_array.array)(x);
    var dtype = x.dtype === 'float32' ? 'float32' : 'float64',
        DTypeArray = _dt.ARRAY_TYPES[dtype],
        shape = x.shape;
    x = x.data;
    var g = new DTypeArray(x.length),
        epsRel = null != h_rel ? h_rel : Math.pow((0, _dt.eps)(dtype), 1 / 3),
        epsAbs = null != h_abs ? h_abs : Math.pow((0, _dt.eps)(dtype), 1 / 3);
    var f0 = f(new _nd_array.NDArray(shape, x.slice()));
    g.fill(-1.5 * f0); // http://macs.citadel.edu/chenm/344.dir/temp.dir/lect4_1.pdf

    for (var i = x.length; i-- > 0;) {
      var h = Math.max(Math.abs(x[i]) * epsRel, epsAbs); // TODO maybe a check for h===0 might be sensible

      for (var _i3 = 0, _arr = [[1 * h, +2.0], [2 * h, -0.5]]; _i3 < _arr.length; _i3++) {
        var _arr$_i = (0, _slicedToArray2["default"])(_arr[_i3], 2),
            d = _arr$_i[0],
            w = _arr$_i[1];

        var X = new _nd_array.NDArray(shape, DTypeArray.from(x));
        X.data[i] += d;
        g[i] += f(X) * w;
      }

      g[i] /= h;
    }

    return new _nd_array.NDArray(shape, g);
  };
};

exports.num_grad_forward = num_grad_forward;