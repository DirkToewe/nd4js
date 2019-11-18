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
exports.num_grad = void 0;

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
    x = (0, _nd_array.asarray)(x);
    var dtype = x.dtype === 'float32' ? 'float32' : 'float64',
        DTypeArray = _dt.ARRAY_TYPES[dtype],
        shape = x.shape;
    x = x.data;
    var g = new DTypeArray(x.length),
        epsRel = null != h_rel ? h_rel : Math.pow((0, _dt.eps)(dtype), 1 / 3),
        epsAbs = null != h_abs ? h_abs : Math.pow((0, _dt.eps)(dtype), 1 / 3); // https://en.wikipedia.org/wiki/Finite_difference#Forward,_backward,_and_central_differences
    // https://www.geometrictools.com/Documentation/FiniteDifferences.pdf

    var _loop = function _loop(i) {
      var h = Math.max(Math.abs(x[i]) * epsRel, epsAbs); // TODO maybe a check for h===0 might be sensible

      var _loop2 = function _loop2() {
        var _arr$_i = (0, _slicedToArray2["default"])(_arr[_i], 2),
            d = _arr$_i[0],
            w = _arr$_i[1];

        var X = new _nd_array.NDArray(shape, DTypeArray.from(x, function (xj, j) {
          return xj + d * (i === j);
        }));
        g[i] += w * f(X);
      };

      for (var _i = 0, _arr = [[+2 * h, -1], [+1 * h, +8], [-1 * h, -8], [-2 * h, +1]]; _i < _arr.length; _i++) {
        _loop2();
      }

      g[i] /= 12 * h;
    };

    for (var i = x.length; i-- > 0;) {
      _loop(i);
    }

    return new _nd_array.NDArray(shape, g);
  };
};

exports.num_grad = num_grad;