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

var _nd_array = require("../../nd_array");

var _opt_utils = require("../_opt_utils");

var _line_search_error = require("./line_search_error");

//*DEBUG*/import {num_grad} from '../num_grad'
var strong_wolfe = function strong_wolfe() {
  var _ref = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : {},
      _ref$c = _ref.c1,
      c1 = _ref$c === void 0 ? 0.01 : _ref$c,
      _ref$c2 = _ref.c2,
      c2 = _ref$c2 === void 0 ? 0.9 : _ref$c2,
      _ref$c3 = _ref.c3,
      c3 = _ref$c3 === void 0 ? 2 : _ref$c3,
      _ref$minRed = _ref.minRed,
      minRed = _ref$minRed === void 0 ? 0.2 : _ref$minRed;

  // SEE:
  //   "Numerical Optimization" 2n Edition,
  //   Jorge Nocedal Stephen J. Wright,
  //   Chapter 3. Line Search Methods, page 60.
  c1 *= 1; // (1) Wolfe cond. 1: Sufficient decrease of objective.

  c2 *= 1; // (2) Wolfe cond. 2: Sufficient decrease of projected gradient.

  c3 *= 1; // (3) Growth factor for 1st phase of line search (bracketing).

  minRed *= 1; // (4) Minimum reduction of search range during 2nd phase of line search (zoom).
  // CHECK 0 < c1 < c2 < 1 < c3

  if (!(c1 > 0)) throw new Error('strong_wolfe(opt): opt.c1 must be positive.');
  if (!(c1 < c2)) throw new Error('strong_wolfe(opt): opt.c1 must less than opt.c2.');
  if (!(1 > c2)) throw new Error('strong_wolfe(opt): opt.c2 must less than 1.');
  if (!(1 < c3)) throw new Error('strong_wolfe(opt): opt.c3 must larger than 1.');
  if (!(minRed >= 0.0)) throw new Error('Assertion failed.');
  if (!(minRed <= 0.5)) throw new Error('Assertion failed.');
  var Λ = 1 / minRed,
      λ = Λ - 1;
  return function (fg_raw) {
    return function (X0, f0, G0, negDir) {
      X0 = (0, _nd_array.asarray)(X0);
      G0 = (0, _nd_array.asarray)(G0);
      f0 *= 1;
      if (!(negDir instanceof _nd_array.NDArray)) throw new Error('Assertion failed.');
      var shape = negDir.shape;
      if (X0.ndim !== 1) throw new Error('strong_wolfe()(fg)(X0,f0,G0): X0.ndim must be 1.');
      if (G0.ndim !== 1) throw new Error('strong_wolfe()(fg)(X0,f0,G0): G0.ndim must be 1.');
      if (negDir.ndim !== 1) throw new Error('strong_wolfe()(fg)(X0,f0,G0): negDir.ndim must be 1.');
      if (X0.shape[0] !== shape[0]) throw new Error('strong_wolfe()(fg)(X0,f0,G0): negDir and X0 must have the same shape.');
      if (G0.shape[0] !== shape[0]) throw new Error('strong_wolfe()(fg)(X0,f0,G0): negDir and G0 must have the same shape.');
      if (isNaN(f0)) throw new Error('strong_wolfe()(fg)(X0,f0,G0): f0 is NaN.');
      var dtype = [X0, G0].every(function (x) {
        return x.dtype === 'float32';
      }) ? 'float32' : 'float64';

      var fg = function fg(x) {
        var _fg_raw = fg_raw(x),
            _fg_raw2 = (0, _slicedToArray2["default"])(_fg_raw, 2),
            f = _fg_raw2[0],
            g = _fg_raw2[1];

        f *= 1;
        g = (0, _nd_array.asarray)(g);
        if (isNaN(f)) throw new Error('strong_wolfe()(fg)(X0,f0,G0): fg returned NaN.');
        if (g.ndim !== 1) throw new Error('strong_wolfe()(fg)(X0,f0,G0): fg returned gradient with invalid shape.');
        if (g.shape[0] !== shape[0]) throw new Error('strong_wolfe()(fg)(X0,f0,G0): fg returned gradient with invalid shape.');
        return [f, g];
      };

      negDir = negDir.data; //*DEBUG*/    const projGradTest = num_grad( α => {
      //*DEBUG*/      const X = zip_elems([X0,negDir], dtype, (x,nDir) => x - α*nDir);
      //*DEBUG*/      return fg(X)[0];
      //*DEBUG*/    });

      var projGrad = function projGrad(g) {
        // <- projected gradient
        var pg = 0;

        for (var i = negDir.length; i-- > 0;) {
          pg -= negDir[i] * g.data[i];
        }

        return pg;
      };

      var p0 = projGrad(G0); // <- projected gradient

      if (p0 >= 0) throw new Error('strong_wolfe: Initial projected gradient not negative.');
      var αMin = 0,
          α = 1,
          fMin = f0,
          pMin = p0,
          αMax = Infinity,
          fMax = NaN; // STEP 1: BRACKETING PHASE
      //   Find a range guaranteed to contain an α satisfying strong Wolfe.

      bracketing: for (;;) {
        var X = new _nd_array.NDArray(shape, X0.data.map(function (x, i) {
          return x - α * negDir[i];
        })),
            _fg = fg(X),
            _fg2 = (0, _slicedToArray2["default"])(_fg, 2),
            f = _fg2[0],
            G = _fg2[1],
            p = projGrad(G); //*DEBUG*/      const q = projGradTest(α);
        //*DEBUG*/      if( !(Math.abs(p-q) <= Math.max(Math.abs(p),Math.abs(q),1)*1e-5) ) throw Error(`Assertion failed: ${p} !== ${q}`);


        if (f - f0 > c1 * α * p0 || 0 < αMin && f >= fMin) {
          αMax = α;
          fMax = f;
          break bracketing;
        }

        if (Math.abs(p) <= -c2 * p0) return [X, f, G];

        if (p >= 0) {
          αMax = αMin;
          fMax = fMin;
          αMin = α;
          fMin = f;
          pMin = p;
          break bracketing;
        }

        if (!(α < αMax)) throw new _line_search_error.LineSearchError('strong_wolfe: Strong Wolfe condition not satisfiable in range.');
        αMin = α;
        α *= c3;
        fMin = f;
        pMin = p;
      }

      if (αMin === αMax) throw new _line_search_error.LineSearchError('strong_wolfe: bracketing failed.');

      var noProgress = function noProgress() {
        if (αMin === 0) throw new _line_search_error.LineSearchNoProgressError('strongWolfeLineSearch(): no progress.');
        throw new _line_search_error.LineSearchError("strong_wolfe: bisection failed ".concat(JSON.stringify({
          αMin: αMin,
          α: α,
          αMax: αMax
        }), "."));
      }; // STEP 2: ZOOM PHASE
      //   Given a range that is guaranteed to contain a valid
      //   strong Wolfe α values, this method finds such a value.


      for (var run = 1;; run++) {
        //*DEBUG*/      const q = projGradTest(αMin);
        //*DEBUG*/      if( !(Math.abs(pMin-q) <= Math.max(Math.abs(pMin),Math.abs(q),1)*1e-5) ) throw Error(`Assertion failed: ${pMin} !== ${q}`);
        if (2 === Λ) α = (αMin + αMax) / 2;else {
          α = (0, _opt_utils._min1d_interp_quad)(αMin, αMax, fMin, fMax, pMin);
          var αLo = (αMax + λ * αMin) / Λ,
              αHi = (αMin + λ * αMax) / Λ; //*DEBUG*/        if( !(Math.min(αMin,αMax) <= α) ) throw new Error('Assertion failed.');
          //*DEBUG*/        if( !(Math.max(αMin,αMax) >= α) ) throw new Error('Assertion failed.');

          if (!(αLo <= α)) α = αLo; // < handles NaN
          else if (!(αHi >= α)) α = αHi; // < handles NaN
        }

        var _X = new _nd_array.NDArray(shape, X0.data.map(function (x, i) {
          return x - α * negDir[i];
        })),
            _fg3 = fg(_X),
            _fg4 = (0, _slicedToArray2["default"])(_fg3, 2),
            _f = _fg4[0],
            _G = _fg4[1],
            _p = projGrad(_G);

        if (_f - f0 > c1 * α * p0 || _f >= fMin) {
          if (αMax === α) noProgress();
          αMax = α;
          fMax = _f;
        } else {
          if (Math.abs(_p) <= -c2 * p0) return [_X, _f, _G];

          if (_p * (αMax - αMin) >= 0) {
            αMax = αMin;
            fMax = fMin;
          }

          if (αMin === α) noProgress();
          αMin = α;
          fMin = _f;
          pMin = _p;
        }
      }
    };
  };
};

exports.strong_wolfe = strong_wolfe;