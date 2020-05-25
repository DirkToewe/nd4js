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
exports.albaali_fletcher = void 0;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _nd_array = require("../../nd_array");

var _float64_utils = require("../../dt/float64_utils");

var _line_search_error = require("./line_search_error");

var _line_search_utils = require("./_line_search_utils");

// References
// ----------
// .. [1] "Numerical Optimization", 2nd Edition,
//         Jorge Nocedal & Stephen J. Wright,
//         Chapter 3. Line Search Methods, page 60f
// .. [2] "An Efficient Line Search for Nonlinear Least Squares",
//         M. Al-Baali & R. Fletcher
//         Journal Of Optimization Theory and Applications: Vol. 48, No. 3, MARCH 1986
// .. [3] "Line Search Algorithms with Guaranteed Sufficient Decrease"
//         Jorge J. Moré and David J. Thuente
//         ACM Transactions on Mathematical Software,
//         Vol 20, No. 3, Septermber 1994, Pages 286-307
// Implementation of Scheme S2 from [2].
var albaali_fletcher = function albaali_fletcher() {
  var opt = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : {};
  var _opt = opt,
      _opt$fRed = _opt.fRed,
      fRed = _opt$fRed === void 0 ? 0.1 : _opt$fRed,
      _opt$gRed = _opt.gRed,
      gRed = _opt$gRed === void 0 ? 0.9 : _opt$gRed,
      _opt$grow = _opt.grow,
      grow = _opt$grow === void 0 ? Math.PI / 3 : _opt$grow,
      _opt$shrinkLeast = _opt.shrinkLeast,
      shrinkLeast = _opt$shrinkLeast === void 0 ? 0.2 : _opt$shrinkLeast;
  opt = {
    fRed: fRed,
    gRed: gRed,
    grow: grow,
    shrinkLeast: shrinkLeast
  };
  fRed *= 1; // (1) Wolfe cond. 1: Sufficient decrease of objective.

  gRed *= 1; // (2) Wolfe cond. 2: Sufficient decrease of projected gradient.

  grow *= 1; // (3) Growth factor for 1st phase of line search (bracketing).

  shrinkLeast *= 1; // (4) Minimum reduction of search range during 2nd phase of line search (zoom).
  // CHECK 0 < fRed < gRed < 1 < c3

  if (!(fRed > 0)) throw new Error('albaali_fletcher(opt): opt.fRed must be positive.');
  if (!(fRed < gRed)) throw new Error('albaali_fletcher(opt): opt.fRed must less than opt.gRed.');
  if (!(1 > gRed)) throw new Error('albaali_fletcher(opt): opt.gRed must less than 1.');
  if (!(1 < grow)) throw new Error('albaali_fletcher(opt): opt.grow must larger than 1.');
  if (!(fRed < 0.5)) console.warn('albaali_fletcher(opt): opt.fRed should be less than 0.5 to work properly with (quasi-)Newton methods.');
  if (!(shrinkLeast >= 0.0)) throw new Error('albaali_fletcher(opt): opt.shrinkLeast must be non-negative.');
  if (!(shrinkLeast <= 0.5)) throw new Error('albaali_fletcher(opt): opt.shrinkLeast must not be greater than 0.5.');
  var Λ = 1 / shrinkLeast,
      λ = Λ - 1;

  var line_search = function line_search(fg_raw) {
    return function (X0, f0, G0, negDir) {
      var αMin = arguments.length > 4 && arguments[4] !== undefined ? arguments[4] : 0;
      var α0 = arguments.length > 5 && arguments[5] !== undefined ? arguments[5] : undefined;
      var αMax = arguments.length > 6 && arguments[6] !== undefined ? arguments[6] : Infinity;
      X0 = (0, _nd_array.array)('float64', X0);
      G0 = (0, _nd_array.array)('float64', G0);
      f0 *= 1;
      if (!(negDir instanceof _nd_array.NDArray)) throw new Error('Assertion failed.');

      var shape = negDir.shape,
          _shape = (0, _slicedToArray2["default"])(shape, 1),
          L = _shape[0];

      if (X0.ndim !== 1) throw new Error('albaali_fletcher()(fg)(X0,f0,G0): X0.ndim must be 1.');
      if (G0.ndim !== 1) throw new Error('albaali_fletcher()(fg)(X0,f0,G0): G0.ndim must be 1.');
      if (negDir.ndim !== 1) throw new Error('albaali_fletcher()(fg)(X0,f0,G0): negDir.ndim must be 1.');
      if (X0.shape[0] !== L) throw new Error('albaali_fletcher()(fg)(X0,f0,G0): negDir and X0 must have the same shape.');
      if (G0.shape[0] !== L) throw new Error('albaali_fletcher()(fg)(X0,f0,G0): negDir and G0 must have the same shape.');
      if (isNaN(f0)) throw new Error('albaali_fletcher()(fg)(X0,f0,G0): f0 is NaN.');
      if (αMin !== 0) throw new Error('albaali_fletcher()(fg)(X0,f0,G0, αMin,α0,αMax): αMin !== 0 not (yet) supported.');
      negDir = negDir.data;
      if (null == α0) α0 = Math.min(1, αMax / 2);
      if (!(αMin <= α0)) throw new Error('albaali_fletcher()(fg)(X0,f0,G0, αMin,α0,αMax): αMin !== 0 not (yet) supported.');
      if (!(αMax >= α0)) throw new Error('albaali_fletcher()(fg)(X0,f0,G0, αMin,α0,αMax): αMax must not be less than αMin.');
      if (0 === αMax) throw new _line_search_error.LineSearchNoProgressError();

      var projGrad = function projGrad(g) {
        return g.data.reduce(function (p, gi, i) {
          return p - negDir[i] * gi;
        }, 0);
      }; // wrap some checks aroung fg_raw


      var xfgp = function xfgp(α) {
        var X = new _nd_array.NDArray(shape, X0.data.map(function (x, i) {
          return x - α * negDir[i];
        }));

        var _fg_raw = fg_raw(new _nd_array.NDArray(shape, X.data.slice())),
            _fg_raw2 = (0, _slicedToArray2["default"])(_fg_raw, 2),
            f = _fg_raw2[0],
            G = _fg_raw2[1];

        f *= 1;
        G = (0, _nd_array.array)('float64', G);
        if (isNaN(f)) throw new Error('more_thuente_u123()(fg)(X0,f0,G0): fg returned NaN.');
        if (G.ndim !== 1) throw new Error('more_thuente_u123()(fg)(X0,f0,G0): fg returned gradient with invalid shape.');
        if (G.shape[0] !== L) throw new Error('more_thuente_u123()(fg)(X0,f0,G0): fg returned gradient with invalid shape.');
        var p = projGrad(G);
        return [X, f, G, p];
      };

      var p0 = projGrad(G0); // <- projected gradient

      if (0 === p0) throw new _line_search_error.LineSearchNoProgressError();
      if (0 < p0) throw new Error('albaali_fletcher: Initial projected gradient not negative.');
      var αLo = 0,
          αHi = Infinity,
          α = α0,
          fLo = f0,
          fHi = NaN,
          pLo = p0,
          pHi = NaN; // PHASE 1: BRACKETING PHASE
      // ------------------------
      //   Find a range guaranteed to contain an α satisfying strong Wolfe.

      bracketing: for (;;) {
        if (!isFinite(α)) throw new Error('Assertion failed.');

        var _xfgp = xfgp(α),
            _xfgp2 = (0, _slicedToArray2["default"])(_xfgp, 4),
            _X = _xfgp2[0],
            _f = _xfgp2[1],
            _G = _xfgp2[2],
            p = _xfgp2[3];

        if (_f - f0 > fRed * α * p0 || 0 < αLo && _f >= fLo) {
          αHi = α;
          fHi = _f;
          pHi = p;
          break bracketing;
        }

        if (Math.abs(p) <= -gRed * p0) return [_X, _f, _G];

        if (p >= 0) {
          αHi = αLo;
          αLo = α;
          fHi = fLo;
          fLo = _f;
          pHi = pLo;
          pLo = p;
          break bracketing;
        }

        if (!(α < αMax)) throw new _line_search_error.LineSearchBoundReachedError(_X, _f, _G);
        if (!(α < αHi)) throw new _line_search_error.LineSearchError('albaali_fletcher: Strong Wolfe condition not satisfiable in range.');
        αLo = α;
        α = Math.min(Math.max((0, _float64_utils.nextUp)(α), α * grow), αMax);
        fLo = _f;
        pLo = p;
      }

      if (αLo === αHi) throw new _line_search_error.LineSearchError('albaali_fletcher: bracketing failed.');

      var noProgress = function noProgress() {
        if (0 === αLo) // <- TODO FIXME
          throw new _line_search_error.LineSearchNoProgressError('strongWolfeLineSearch(): no progress.');
        throw new _line_search_error.LineSearchBisectionError(X, f, G);
      }; // PHASE 2: ZOOM PHASE
      // ------------------
      //   Given a range that is guaranteed to contain a valid
      //   strong Wolfe α values, this method finds such a value.


      for (var run = 1;; run++) {
        // SELECT A NEW TRIAL VALUE FOR α
        {
          var αLil = Math.min(αLo, αHi),
              αBig = Math.max(αLo, αHi);
          var αLst = Math.max((0, _float64_utils.nextUp)(αLil), (αBig + λ * αLil) / Λ),
              // <- Least
          αMst = Math.min((0, _float64_utils.nextDown)(αBig), (αLil + λ * αBig) / Λ); // <- Most

          if (2 === Λ || !(αLst < αMst) || !isFinite(fHi) || !isFinite(pHi)) α = αLo + (αHi - αLo) / 2;else {
            // TODO: the minimum will satisfy gRed but might not satisfy fRed.
            //       It might be smarter to use the interpolation to find
            //       a point that potentially satisfies both.
            α = (0, _line_search_utils._min1d_interp_ffg)(αLo, αHi, fLo, fHi, pLo);
            if (!(αLst <= α)) α = αLst; // < handles NaN
            else if (!(αMst >= α)) α = αMst; // < handles NaN
          }
          /*DEBUG*/

          if (!(αLil <= α)) throw new Error('Assertion failed.');
          /*DEBUG*/

          if (!(αBig >= α)) throw new Error('Assertion failed.');
        }

        var _xfgp3 = xfgp(α),
            _xfgp4 = (0, _slicedToArray2["default"])(_xfgp3, 4),
            _X2 = _xfgp4[0],
            _f2 = _xfgp4[1],
            _G2 = _xfgp4[2],
            _p = _xfgp4[3];

        if (_f2 - f0 > fRed * α * p0 || _f2 >= fLo) {
          if (αHi === α) noProgress();
          αHi = α;
          fHi = _f2;
          pHi = _p;
        } else {
          if (Math.abs(_p) <= -gRed * p0) return [_X2, _f2, _G2];

          if (Math.sign(αHi - αLo) * _p >= 0) {
            αHi = αLo;
            fHi = fLo;
            pHi = pLo;
          }

          if (αLo === α) noProgress();
          αLo = α;
          fLo = _f2;
          pLo = _p;
        }
      }
    };
  };

  Object.defineProperty(line_search, 'name', {
    value: "albaali_fletcher(".concat(JSON.stringify(opt), ")"),
    writable: false
  });
  Object.assign(line_search, opt);
  Object.freeze(line_search);
  return line_search;
};

exports.albaali_fletcher = albaali_fletcher;