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
exports.root_newton_gen = root_newton_gen;

var _regenerator = _interopRequireDefault(require("@babel/runtime/regenerator"));

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _nd_array = require("../nd_array");

var _lstsq = require("../la/lstsq");

var _marked = /*#__PURE__*/_regenerator["default"].mark(root_newton_gen);

function root_newton_gen(fJ, x0) {
  var x, _x$shape, N, _fJ$map, _fJ$map2, f, J, X, dx, DX, i;

  return _regenerator["default"].wrap(function root_newton_gen$(_context) {
    while (1) {
      switch (_context.prev = _context.next) {
        case 0:
          x = (0, _nd_array.array)('float64', x0);
          x0 = undefined;

          if (!(x.ndim !== 1)) {
            _context.next = 4;
            break;
          }

          throw new Error('root_newton_gen(fJ, x0): x0.ndim must be 1.');

        case 4:
          _x$shape = (0, _slicedToArray2["default"])(x.shape, 1), N = _x$shape[0];

        case 5:
          _fJ$map = fJ(x).map(function (a) {
            return (0, _nd_array.asarray)(a);
          }), _fJ$map2 = (0, _slicedToArray2["default"])(_fJ$map, 2), f = _fJ$map2[0], J = _fJ$map2[1];
          _context.next = 8;
          return [x, f, J];

        case 8:
          if (!(J.ndim !== 2 || J.shape[0] !== N || J.shape[1] !== N || f.ndim !== 1 || f.shape[0] !== N)) {
            _context.next = 10;
            break;
          }

          throw new Error('root_newton_gen(fJ, x0: float[N]): fJ return type must be [float[N],float[N,N]].');

        case 10:
          X = x.data, dx = (0, _lstsq.lstsq)(J, f.reshape(N, 1)).reshape(N), DX = dx.data;

          for (i = DX.length; i-- > 0;) {
            DX[i] = X[i] - DX[i];
          }

          x = dx;

        case 13:
          _context.next = 5;
          break;

        case 15:
        case "end":
          return _context.stop();
      }
    }
  }, _marked);
}