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
exports.min_lbfgs_gen = min_lbfgs_gen;

var _regenerator = _interopRequireDefault(require("@babel/runtime/regenerator"));

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _strong_wolfe = require("./line_search/strong_wolfe");

var _nd_array = require("../nd_array");

var _line_search_error = require("./line_search/line_search_error");

var _marked = /*#__PURE__*/_regenerator["default"].mark(min_lbfgs_gen);

function dot(u, v) {
  if (u.length !== v.length) throw new Error();
  var result = 0;

  for (var i = u.length; i-- > 0;) {
    result += u[i] * v[i];
  }

  return result;
}

function min_lbfgs_gen(fg, x0) {
  var _ref,
      _ref$historySize,
      historySize,
      _ref$lineSearch,
      lineSearch,
      _ref$negDir,
      negDir0,
      x,
      _fg,
      _fg2,
      f,
      g,
      _x$shape,
      L,
      shape,
      dX,
      dGdX,
      dG,
      negDir,
      step,
      _args = arguments;

  return _regenerator["default"].wrap(function min_lbfgs_gen$(_context) {
    while (1) {
      switch (_context.prev = _context.next) {
        case 0:
          _ref = _args.length > 2 && _args[2] !== undefined ? _args[2] : {}, _ref$historySize = _ref.historySize, historySize = _ref$historySize === void 0 ? 8 : _ref$historySize, _ref$lineSearch = _ref.lineSearch, lineSearch = _ref$lineSearch === void 0 ? (0, _strong_wolfe.strong_wolfe)() : _ref$lineSearch, _ref$negDir = _ref.negDir0, negDir0 = _ref$negDir === void 0 ? function (g) {
            return g;
          } : _ref$negDir;
          x = (0, _nd_array.array)('float64', x0);
          x0 = undefined;
          _fg = fg(x), _fg2 = (0, _slicedToArray2["default"])(_fg, 2), f = _fg2[0], g = _fg2[1];

          if (f.dtype.startsWith('float')) {
            _context.next = 6;
            break;
          }

          throw new Error('min_lbfgs_gen(fg, x0, opt): fg must return [float,float[]].');

        case 6:
          if (g.dtype.startsWith('float')) {
            _context.next = 8;
            break;
          }

          throw new Error('min_lbfgs_gen(fg, x0, opt): fg must return [float,float[]].');

        case 8:
          _x$shape = (0, _slicedToArray2["default"])(x.shape, 1), L = _x$shape[0], shape = x.shape;
          lineSearch = lineSearch(fg);
          dX = [], dGdX = [], dG = [];

          negDir = function negDir(g) {
            //SEE:
            // Jorge Nocedal "Updating Quasi-Newton Matrices with Limited Storage"
            // MATHEMATICS OF COMPUTATION, VOL. 35, NO. 151, JULY 1980, PAGES 773-78
            // https://courses.engr.illinois.edu/ece544na/fa2014/nocedal80.pdf
            g = Float64Array.from(g.data);
            var α = new Float64Array(dX.length);

            for (var i = dGdX.length; i-- > 0;) {
              α[i] = dot(dX[i], g) / dGdX[i];

              for (var j = L; j-- > 0;) {
                g[j] -= α[i] * dG[i][j];
              }
            }

            g = negDir0(new _nd_array.NDArray(shape, g));
            g = Float64Array.from(g.data);

            for (var _i = 0; _i < dGdX.length; _i++) {
              var βi = dot(dG[_i], g) / dGdX[_i];

              for (var _j = L; _j-- > 0;) {
                g[_j] += (α[_i] - βi) * dX[_i][_j];
              }
            }

            return new _nd_array.NDArray(shape, g);
          };

          step = function step() {
            var _lineSearch = lineSearch(x, f, g, negDir(g)),
                _lineSearch2 = (0, _slicedToArray2["default"])(_lineSearch, 3),
                X = _lineSearch2[0],
                F = _lineSearch2[1],
                G = _lineSearch2[2];

            x = x.data;
            g = g.data;
            var dg = Float64Array.from(G.data, function (Gi, i) {
              return Gi - g[i];
            }),
                dx = Float64Array.from(X.data, function (Xi, i) {
              return Xi - x[i];
            });
            dG.push(dg);
            dGdX.push(dot(dg, dx));
            dX.push(dx); // LIMIT THE NUMBER OF MEMOIZED GRADIENTS
            // (while loop in case historySize was changed)

            if (dX.length > historySize) {
              dX.shift();
              dGdX.shift();
              dG.shift();
            }

            if (dX.length > historySize) throw new Error('Assertion#1 failed.');
            if (dX.length !== dGdX.length) throw new Error('Assertion#2 failed.');
            if (dX.length !== dG.length) throw new Error('Assertion#3 failed.');
            x = X;
            f = F;
            g = G;
          };

        case 13:
          if (!(x.ndim !== 1)) {
            _context.next = 15;
            break;
          }

          throw new Error('Assertion#4 failed.');

        case 15:
          if (!isNaN(f * 1)) {
            _context.next = 17;
            break;
          }

          throw new Error('Assertion#5 failed.');

        case 17:
          if (!(g.ndim !== 1)) {
            _context.next = 19;
            break;
          }

          throw new Error('Assertion#6 failed.');

        case 19:
          if (!(x.shape[0] !== L)) {
            _context.next = 21;
            break;
          }

          throw new Error('Assertion#7 failed.');

        case 21:
          if (!(g.shape[0] !== L)) {
            _context.next = 23;
            break;
          }

          throw new Error('Assertion#8 failed.');

        case 23:
          _context.next = 25;
          return [new _nd_array.NDArray(shape, x.data.slice()), f * 1, new _nd_array.NDArray(shape, g.data.slice())];

        case 25:
          _context.prev = 25;
          step();
          _context.next = 37;
          break;

        case 29:
          _context.prev = 29;
          _context.t0 = _context["catch"](25);

          if (!(_context.t0 instanceof _line_search_error.LineSearchNoProgressError)) {
            _context.next = 36;
            break;
          }

          // single attempt to restart
          dX = [], dGdX = [], dG = [];
          step();
          _context.next = 37;
          break;

        case 36:
          throw _context.t0;

        case 37:
          _context.next = 13;
          break;

        case 39:
        case "end":
          return _context.stop();
      }
    }
  }, _marked, null, [[25, 29]]);
}