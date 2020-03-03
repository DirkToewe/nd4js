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

var _nd_array = require("../nd_array");

var _dt = require("../dt");

var num_grad = function num_grad(f) {
  var _ref = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : {},
      _ref$h_rel = _ref.h_rel,
      h_rel = _ref$h_rel === void 0 ? undefined : _ref$h_rel,
      _ref$h_abs = _ref.h_abs,
      h_abs = _ref$h_abs === void 0 ? undefined : _ref$h_abs;

  return function (x) {
    x = (0, _nd_array.array)(x);
    var dtype = x.dtype === 'float32' ? 'float32' : 'float64',
        DTypeArray = _dt.ARRAY_TYPES[dtype],
        shape = x.shape;
    x = x.data;
    var g = new DTypeArray(x.length),
        epsRel = null != h_rel ? h_rel : Math.pow((0, _dt.eps)(dtype), 1 / 3),
        epsAbs = null != h_abs ? h_abs : Math.pow((0, _dt.eps)(dtype), 1 / 3); // https://en.wikipedia.org/wiki/Finite_difference#Forward,_backward,_and_central_differences
    // https://www.geometrictools.com/Documentation/FiniteDifferences.pdf

    for (var i = x.length; i-- > 0;) {
      var h = Math.max(Math.abs(x[i]) * epsRel, epsAbs); // TODO maybe a check for h===0 might be sensible

      for (var _i = 0, _arr = [[+2 * h, -1], [+1 * h, +8], [-1 * h, -8], [-2 * h, +1]]; _i < _arr.length; _i++) {
        var _arr$_i = (0, _slicedToArray2["default"])(_arr[_i], 2),
            d = _arr$_i[0],
            w = _arr$_i[1];

        var X = new _nd_array.NDArray(shape, DTypeArray.from(x));
        X.data[i] += d;
        g[i] += f(X) * w;
      }

      g[i] /= 12 * h;
    }

    return new _nd_array.NDArray(shape, g);
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

      for (var _i2 = 0, _arr2 = [[1 * h, +2.0], [2 * h, -0.5]]; _i2 < _arr2.length; _i2++) {
        var _arr2$_i = (0, _slicedToArray2["default"])(_arr2[_i2], 2),
            d = _arr2$_i[0],
            w = _arr2$_i[1];

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