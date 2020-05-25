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
exports.more_thuente_u123 = void 0;

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
// .. [2] "Line Search Algorithms with Guaranteed Sufficient Decrease"
//         Jorge J. Moré and David J. Thuente
//         ACM Transactions on Mathematical Software,
//         Vol 20, No. 3, Septermber 1994, Pages 286-307
// TODO: study dcsrch.f and dcstep.f from MINPACK2 for possible improvements
//       https://github.com/scipy/scipy/tree/master/scipy/optimize/minpack2
var DEFAULT_OPTIONS = Object.freeze({
  fRed: 1e-2,
  gRed: 0.9,
  growMin: Math.PI / 3,
  growMax: Math.E - 1.5,
  shrinkLeast: 0.1
}); // implements the pseudocode lines (U1,U2,U3) from [2].

var more_thuente_u123 = function more_thuente_u123() {
  var opt = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : {};

  for (var _i = 0, _Object$keys = Object.keys(opt); _i < _Object$keys.length; _i++) {
    var key = _Object$keys[_i];
    if (!(key in DEFAULT_OPTIONS)) console.warn("more_thuente_u123(opt): unknown parameter \"".concat(key, "\" in opt."));
  }

  var _opt = opt = Object.assign({}, DEFAULT_OPTIONS, opt),
      fRed = _opt.fRed,
      gRed = _opt.gRed,
      growMin = _opt.growMin,
      growMax = _opt.growMax,
      shrinkLeast = _opt.shrinkLeast;

  fRed *= 1; // (1) Wolfe cond. 1: Sufficient decrease of objective.

  gRed *= 1; // (2) Wolfe cond. 2: Sufficient decrease of projected gradient.

  growMin *= 1; // (3) Lower growth factor bound for 1st phase of line search (bracketing).

  growMax *= 1; // (4) Upper growth factor bound for 1st phase of line search (bracketing).

  shrinkLeast *= 1; // (5) Minimum reduction of search range during 2nd phase of line search (zoom).
  // CHECK 0 < fRed < gRed < 1 < c3

  if (!(fRed > 0)) throw new Error('more_thuente_u123(opt): opt.fRed must be positive.');
  if (!(fRed < gRed)) throw new Error('more_thuente_u123(opt): opt.fRed must less than opt.gRed.');
  if (!(1 > gRed)) throw new Error('more_thuente_u123(opt): opt.gRed must less than 1.');
  if (!(1 < growMin)) throw new Error('more_thuente_u123(opt): opt.growMin must larger than 1.');
  if (!(growMax >= growMin)) throw new Error('more_thuente_u123(opt): opt.growMax must not be less than opt.growMin.');
  if (!(fRed < 0.5)) console.warn('more_thuente_u123(opt): opt.fRed should be less than 0.5 to work properly with (quasi-)Newton methods.');
  if (!(shrinkLeast >= 0.0)) throw new Error('more_thuente_u123(opt): opt.shrinkLeast must be non-negative.');
  if (!(shrinkLeast <= 0.5)) throw new Error('more_thuente_u123(opt): opt.shrinkLeast must not be greater than 0.5.');

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

      if (X0.ndim !== 1) throw new Error('more_thuente_u123()(fg)(X0,f0,G0): X0.ndim must be 1.');
      if (G0.ndim !== 1) throw new Error('more_thuente_u123()(fg)(X0,f0,G0): G0.ndim must be 1.');
      if (negDir.ndim !== 1) throw new Error('more_thuente_u123()(fg)(X0,f0,G0): negDir.ndim must be 1.');
      if (X0.shape[0] !== L) throw new Error('more_thuente_u123()(fg)(X0,f0,G0): negDir and X0 must have the same shape.');
      if (G0.shape[0] !== L) throw new Error('more_thuente_u123()(fg)(X0,f0,G0): negDir and G0 must have the same shape.');
      if (αMin !== 0) throw new Error('more_thuente_u123()(fg)(X0,f0,G0, αMin,α0,αMax): αMin !== 0 not (yet) supported.');
      if (isNaN(f0)) throw new Error('more_thuente_u123()(fg)(X0,f0,G0): f0 is NaN.');
      if (X0.data.some(isNaN)) throw new Error('more_thuente_u123()(fg)(X0,f0,G0): X0 contains NaN.');
      if (G0.data.some(isNaN)) throw new Error('more_thuente_u123()(fg)(X0,f0,G0): G0 contains NaN.');
      negDir = negDir.data;
      if (negDir.some(isNaN)) throw new Error('Assertion failed.'); // TODO: Instead of a specified αMin, we could compute
      //       the smallest αMin that still results in an
      //       X different from X0, i.e. an αMin that
      //       does not result in complete underflow.

      if (null == α0) α0 = Math.min(1, αMax / 2);
      if (!(αMin <= α0)) throw new Error('more_thuente_u123()(fg)(X0,f0,G0, αMin,α0,αMax): αMin must not be greater than α0.');
      if (!(αMax >= α0)) throw new Error('more_thuente_u123()(fg)(X0,f0,G0, αMin,α0,αMax): αMax must not be less than α0.');
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
      if (0 < p0) throw new Error('more_thuente_u123: Initial projected gradient not negative.');
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
            X = _xfgp2[0],
            f = _xfgp2[1],
            G = _xfgp2[2],
            p = _xfgp2[3]; // check convergence


        if (f - f0 <= fRed * α * p0 && Math.abs(p) <= -gRed * p0) return [X, f, G];

        if (f - fLo > fRed * (α - αLo) * p0) {
          αHi = α;
          fHi = f;
          pHi = p;
          break bracketing;
        }

        if (p > 0) {
          αHi = αLo;
          αLo = α;
          fHi = fLo;
          fLo = f;
          pHi = pLo;
          pLo = p;
          break bracketing;
        }

        if (!(α < αMax)) throw new _line_search_error.LineSearchBoundReachedError(X, f, G);
        if (!(α < Number.MAX_VALUE)) throw new Error('more_thuente_u123: Infinity reached.');
        var αTrial = void 0;
        if (pLo < p) αTrial = (0, _line_search_utils._min1d_interp_gg)(αLo, α, pLo, p);else {
          αTrial = α * growMin;
        }
        αTrial = Math.min(αTrial, α * growMax);
        αTrial = Math.max(αTrial, α * growMin);
        if (!(α < αTrial)) αTrial = (0, _float64_utils.nextUp)(α); // <- handles NaN

        if (!(αMax >= αTrial)) αTrial = αMax; // <- handles NaN

        αLo = α;
        α = αTrial;
        fLo = f;
        pLo = p;
      }

      if (αLo === αHi) throw new _line_search_error.LineSearchError('more_thuente_u123: bracketing failed.');

      var noProgress = function noProgress(X, f, G) {
        if (αLo === 0) // <- TODO FIXME
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
          var αLst = Math.max((0, _float64_utils.nextUp)(αLil), shrinkLeast * αBig + (1 - shrinkLeast) * αLil),
              // <- Least
          αMst = Math.min((0, _float64_utils.nextDown)(αBig), shrinkLeast * αLil + (1 - shrinkLeast) * αBig); // <- Most

          var αMid = αLo + (αHi - αLo) / 2;
          if (shrinkLeast === 0.5 || !(αLst < αMst) || !isFinite(fHi) || !isFinite(pHi)) α = αMid;else {
            // TODO: the minimum will satisfy gRed but might not satisfy fRed.
            //       It might be smarter to use the interpolation to find
            //       a point that potentially satisfies both.
            if (fLo < fHi) {
              // See [2], section "Trial Value Selection", Case 1
              var αc = (0, _line_search_utils._min1d_interp_ffgg)(αLo, αHi, fLo, fHi, pLo, pHi),
                  αq = (0, _line_search_utils._min1d_interp_ffg)(αLo, αHi, fLo, fHi, pLo);
              α = Math.abs(αc - αLo) < Math.abs(αq - αLo) ? αc : (αc + αq) / 2; // <- TODO: Looking for a minimum might not be smart, we should consider trying to find a strong wolfe range instead
            } else if (Math.sign(pLo) * pHi < 0) {
              // See [2], section "Trial Value Selection", Case 1
              var _c = (0, _line_search_utils._min1d_interp_ffgg)(αLo, αHi, fLo, fHi, pLo, pHi),
                  αs = (0, _line_search_utils._min1d_interp_gg)(αLo, αHi, pLo, pHi);

              α = Math.abs(αs - αHi) <= Math.abs(_c - αHi) ? _c : αs; // <- TODO: Looking for a minimum might not be smart, we should consider trying to find a strong wolfe range instead
            } else {
              // α = αMid;
              α = (0, _line_search_utils._min1d_interp_ffg)(αLo, αHi, fLo, fHi, pLo);
            }

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
            _X = _xfgp4[0],
            _f = _xfgp4[1],
            _G = _xfgp4[2],
            _p = _xfgp4[3]; // check convergence


        if (_f - f0 <= fRed * α * p0 && Math.abs(_p) <= -gRed * p0) return [_X, _f, _G];

        if (_f - fLo > fRed * (α - αLo) * p0) {
          if (αHi === α) noProgress(_X, _f, _G);
          αHi = α;
          fHi = _f;
          pHi = _p;
        } else {
          if (Math.sign(αHi - αLo) * _p > 0) {
            αHi = αLo;
            fHi = fLo;
            pHi = pLo;
          }

          if (αLo === α) noProgress(_X, _f, _G);
          αLo = α;
          fLo = _f;
          pLo = _p;
        }
      }
    };
  };

  Object.defineProperty(line_search, 'name', {
    value: "more_thuente_u123(".concat(JSON.stringify(opt), ")"),
    writable: false
  });
  Object.assign(line_search, opt);
  Object.freeze(line_search);
  return line_search;
};

exports.more_thuente_u123 = more_thuente_u123;