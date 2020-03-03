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
exports.strong_wolfe = void 0;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _zip_elems = require("../../zip_elems");

var _line_search_error = require("./line_search_error");

var _nd_array = require("../../nd_array");

var strong_wolfe = function strong_wolfe() {
  var _ref = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : {},
      _ref$c = _ref.c1,
      c1 = _ref$c === void 0 ? 0.4 : _ref$c,
      _ref$c2 = _ref.c2,
      c2 = _ref$c2 === void 0 ? 0.8 : _ref$c2,
      _ref$c3 = _ref.c3,
      c3 = _ref$c3 === void 0 ? 1.6 : _ref$c3;

  // SEE:
  //   "Numerical Optimization" 2n Edition,
  //   Jorge Nocedal Stephen J. Wright,
  //   Chapter 3. Line Search Methods, page 60.
  // CHECK 0 < c1 < c2 < 1 < c3
  if (c1 <= 0) throw new Error('strong_wolfe({c1,c2,c3}): c1 must be positive.');
  if (c1 >= c2) throw new Error('strong_wolfe({c1,c2,c3}): c1 must less than c2.');
  if (1 <= c2) throw new Error('strong_wolfe({c1,c2,c3}): c2 must less than 1.');
  if (1 >= c3) throw new Error('strong_wolfe({c1,c2,c3}): c3 must larger than 1.');
  return function (fg) {
    return function (X0, f0, G0, negDir) {
      X0 = (0, _nd_array.asarray)(X0);
      G0 = (0, _nd_array.asarray)(G0);
      if (X0.ndim !== 1) throw new Error('strong_wolfe()(fg)(X0,f0,G0): X0.ndim must be 1.');
      if (G0.ndim !== 1) throw new Error('strong_wolfe()(fg)(X0,f0,G0): G0.ndim must be 1.');
      if (X0.shape[0] !== G0.shape[0]) throw new Error('strong_wolfe()(fg)(X0,f0,G0): X0 and G0 must have the same shape.');
      if (isNaN(f0)) throw new Error('strong_wolfe()(fg)(X0,f0,G0): f0 is NaN.');
      var dtype = [X0, G0].every(function (x) {
        return x.dtype === 'float32';
      }) ? 'float32' : 'float64';

      var projGrad = function projGrad(g) {
        // <- projected gradient
        var pg = 0;

        for (var _negDir$shape = (0, _slicedToArray2["default"])(negDir.shape, 1), i = _negDir$shape[0]; i-- > 0;) {
          pg -= negDir.data[i] * g.data[i];
        }

        return pg;
      };

      var p0 = projGrad(G0); // <- projected gradient

      if (p0 >= 0) throw new Error('strong_wolfe: Initial projected gradient not negative.');
      var αMin = 0,
          α = 1,
          αMax = Infinity,
          fMin = f0; // STEP 1: BRACKETING PHASE
      //   Find a range guaranteed to contain an α satisfying strong Wolfe.

      bracketing: while (true) {
        var X = (0, _zip_elems.zip_elems)([X0, negDir], dtype, function (x, nDir) {
          return x - α * nDir;
        }),
            _fg = fg(X),
            _fg2 = (0, _slicedToArray2["default"])(_fg, 2),
            f = _fg2[0],
            G = _fg2[1],
            p = projGrad(G);

        if (f - f0 > c1 * α * p0 || 0 < αMin && f >= fMin) {
          αMax = α;
          break bracketing;
        }

        if (Math.abs(p) <= -c2 * p0) return [X, f, G];

        if (p >= 0) {
          αMax = αMin;
          αMin = α;
          fMin = f;
          break bracketing;
        }

        if (!(α < αMax)) throw new _line_search_error.LineSearchError('strong_wolfe: Strong Wolfe condition not satisfiable in range.');
        αMin = α;
        α *= c3;
        fMin = f;
      }

      if (αMin === αMax) throw new _line_search_error.LineSearchError('strong_wolfe: bracketing failed.');

      var noProgress = function noProgress() {
        if (αMin === 0) throw new _line_search_error.LineSearchNoProgressError('strongWolfeLineSearch(): no progress.');
        throw new _line_search_error.LineSearchError("strong_wolfe: bisection failed ".concat(JSON.stringify({
          αMin: αMin,
          α: α,
          αMax: αMax
        }), "."));
      }; // STEP 2: BISECTION PHASE
      //   Given a range that is guaranteed to contain a valid
      //   strong Wolfe α values, this method finds such a value.


      while (true) {
        α = (αMin + αMax) / 2; // <- TODO: use quadratic polynomial to find new point

        var _X = (0, _zip_elems.zip_elems)([X0, negDir], dtype, function (x, nDir) {
          return x - α * nDir;
        }),
            _fg3 = fg(_X),
            _fg4 = (0, _slicedToArray2["default"])(_fg3, 2),
            _f = _fg4[0],
            _G = _fg4[1],
            _p = projGrad(_G);

        if (_f - f0 > c1 * α * p0 || _f >= fMin) {
          if (αMax === α) noProgress();
          αMax = α;
        } else {
          if (Math.abs(_p) <= -c2 * p0) return [_X, _f, _G];
          if (_p * (αMax - αMin) >= 0) αMax = αMin;
          if (αMin === α) noProgress();
          αMin = α;
          fMin = _f;
        }
      }
    };
  };
};

exports.strong_wolfe = strong_wolfe;