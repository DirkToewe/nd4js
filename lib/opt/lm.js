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
exports.fit_lm_gen = fit_lm_gen;
exports.odr_lm_gen = exports.tls_lm_gen = exports.lsq_lm_gen = void 0;

var _regenerator = _interopRequireDefault(require("@babel/runtime/regenerator"));

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _nd_array = require("../nd_array");

var _float64_utils = require("../dt/float64_utils");

var _norm = require("../la/norm");

var _optimization_error = require("./optimization_error");

var _trust_region_solver_lsq = require("./_trust_region_solver_lsq");

var _trust_region_solver_tls = require("./_trust_region_solver_tls");

var _marked = /*#__PURE__*/_regenerator["default"].mark(_lm),
    _marked2 = /*#__PURE__*/_regenerator["default"].mark(fit_lm_gen);

// References
// ----------
// .. [1] "THE LEVENBERG-MARQUARDT ALGORITHM: IMPLEMENTATION AND THEORY"
//         Jorge J. Moré.

/** Solves a nonlinear least-squares problem using a
 *  trust-region variant Levenberg-Marquard algorithm
 *  as described in:
 *
 *    "THE LEVENBERG-MARQUARDT ALGORITHM: IMPLEMENTATION AND THEORY"
 *    by Jorge J. Moré.
 */
var lsq_lm_gen = function lsq_lm_gen(fJ, x0, opt) {
  return _lm(new _trust_region_solver_lsq.TrustRegionSolverLSQ(fJ, x0), opt);
};

exports.lsq_lm_gen = lsq_lm_gen;

function _lm(solver) {
  var _ref,
      _ref$r,
      r0,
      _ref$rMin,
      rMin,
      _ref$rMax,
      rMax,
      _ref$rNewton,
      rNewton,
      _ref$rTol,
      rTol,
      _ref$lmLower,
      lmLower,
      _ref$shrinkLower,
      shrinkLower,
      _ref$shrinkUpper,
      shrinkUpper,
      _ref$grow,
      grow,
      _ref$expectGainMin,
      expectGainMin,
      _ref$expectGainMax,
      expectGainMax,
      _ref$stuckLimit,
      stuckLimit,
      N,
      D,
      dX,
      G,
      NORM,
      R,
      loss,
      stuckometer,
      _solver$computeNewton,
      _solver$computeNewton2,
      r,
      dr,
      gnInRadius,
      i,
      d,
      g,
      λMin,
      λMax,
      λ,
      nIter,
      _solver$computeNewton3,
      _solver$computeNewton4,
      _solver$considerMove,
      _solver$considerMove2,
      loss_predict,
      loss_consider,
      rScale,
      predict,
      actual,
      grad,
      _i,
      shrink,
      _r,
      rs,
      _args = arguments;

  return _regenerator["default"].wrap(function _lm$(_context) {
    while (1) {
      switch (_context.prev = _context.next) {
        case 0:
          _ref = _args.length > 1 && _args[1] !== undefined ? _args[1] : {}, _ref$r = _ref.r0, r0 = _ref$r === void 0 ? 1.1 : _ref$r, _ref$rMin = _ref.rMin, rMin = _ref$rMin === void 0 ? 0 : _ref$rMin, _ref$rMax = _ref.rMax, rMax = _ref$rMax === void 0 ? Infinity : _ref$rMax, _ref$rNewton = _ref.rNewton, rNewton = _ref$rNewton === void 0 ? 2 : _ref$rNewton, _ref$rTol = _ref.rTol, rTol = _ref$rTol === void 0 ? 0.05 : _ref$rTol, _ref$lmLower = _ref.lmLower, lmLower = _ref$lmLower === void 0 ? 0.001 : _ref$lmLower, _ref$shrinkLower = _ref.shrinkLower, shrinkLower = _ref$shrinkLower === void 0 ? 0.05 : _ref$shrinkLower, _ref$shrinkUpper = _ref.shrinkUpper, shrinkUpper = _ref$shrinkUpper === void 0 ? 0.95 : _ref$shrinkUpper, _ref$grow = _ref.grow, grow = _ref$grow === void 0 ? 1.4142135623730951 : _ref$grow, _ref$expectGainMin = _ref.expectGainMin, expectGainMin = _ref$expectGainMin === void 0 ? 0.25 : _ref$expectGainMin, _ref$expectGainMax = _ref.expectGainMax, expectGainMax = _ref$expectGainMax === void 0 ? 0.75 : _ref$expectGainMax, _ref$stuckLimit = _ref.stuckLimit, stuckLimit = _ref$stuckLimit === void 0 ? 64 : _ref$stuckLimit;

          if (0 <= lmLower) {
            _context.next = 3;
            break;
          }

          throw new Error('lsq_lm_gen(fJ, x0, opt): opt.lmLower must be within [0,1).');

        case 3:
          if (lmLower < 1) {
            _context.next = 5;
            break;
          }

          throw new Error('lsq_lm_gen(fJ, x0, opt): opt.lmLower must be within [0,1).');

        case 5:
          if (1 > rTol) {
            _context.next = 7;
            break;
          }

          throw new Error('lsq_lm_gen(fJ, x0, opt): opt.rTol must be less than 1.');

        case 7:
          if (0 < rTol) {
            _context.next = 9;
            break;
          }

          throw new Error('lsq_lm_gen(fJ, x0, opt): opt.rTol must be positive number.');

        case 9:
          if (0 < r0) {
            _context.next = 11;
            break;
          }

          throw new Error('lsq_lm_gen(fJ, x0, opt): opt.r0 must be positive number.');

        case 11:
          if (0 <= rMin) {
            _context.next = 13;
            break;
          }

          throw new Error('lsq_lm_gen(fJ, x0, opt): opt.rMin must be non-negative.');

        case 13:
          if (rMin <= rMax) {
            _context.next = 15;
            break;
          }

          throw new Error('lsq_lm_gen(fJ, x0, opt): opt.rMin must not be greater than opt.rMax.');

        case 15:
          if (rMin <= r0) {
            _context.next = 17;
            break;
          }

          throw new Error('lsq_lm_gen(fJ, x0, opt): opt.r0 must not be less than opt.rMin.');

        case 17:
          if (rMax >= r0) {
            _context.next = 19;
            break;
          }

          throw new Error('lsq_lm_gen(fJ, x0, opt): opt.r0 must not be greater than opt.rMax.');

        case 19:
          if (1 / (1 + rTol) < rNewton) {
            _context.next = 21;
            break;
          }

          throw new Error('lsq_lm_gen(fJ, x0, opt): opt.rGN must be greater than 1/(1+rTol).');

        case 21:
          if (0 < shrinkLower) {
            _context.next = 23;
            break;
          }

          throw new Error('lsq_lm_gen(fJ, x0, opt): opt.shrinkLower must be positive number.');

        case 23:
          if (shrinkLower <= shrinkUpper) {
            _context.next = 25;
            break;
          }

          throw new Error('lsq_lm_gen(fJ, x0, opt): opt.shrinkLower must not be greater than opt.shrinkUpper.');

        case 25:
          if (shrinkUpper < 1) {
            _context.next = 27;
            break;
          }

          throw new Error('lsq_lm_gen(fJ, x0, opt): opt.shrinkUpper must be less than 1.');

        case 27:
          if (1 < grow) {
            _context.next = 29;
            break;
          }

          throw new Error('lsq_lm_gen(fJ, x0, opt): opt.grow must be number greater than 1.');

        case 29:
          if (0 < expectGainMin) {
            _context.next = 31;
            break;
          }

          throw new Error('lsq_lm_gen(fJ, x0, opt): opt.expectGainMin must be positive number.');

        case 31:
          if (expectGainMin < expectGainMax) {
            _context.next = 33;
            break;
          }

          throw new Error('lsq_lm_gen(fJ, x0, opt): opt.expectGainMin must be less than expectGainMax.');

        case 33:
          if (!(expectGainMax < 1)) console.warn('lsq_lm_gen(fJ, x0, opt): opt.expectGainMax should be less than 1.');

          if (!(0 !== stuckLimit % 1)) {
            _context.next = 36;
            break;
          }

          throw new Error('lsq_lm_gen(fJ, x0, opt): opt.stuckLimit must be an integer.');

        case 36:
          if (0 <= stuckLimit) {
            _context.next = 38;
            break;
          }

          throw new Error('lsq_lm_gen(fJ, x0, opt): opt.stuckLimit must not be negative.');

        case 38:
          N = solver.N, D = solver.D, dX = solver.regularized_dX, G = solver.G0, NORM = new _norm.FrobeniusNorm(); //  let R = r0; // <- TODO maybe there is a better starting value for R0 like r0*solver.scaledNorm(G) order something like that

          R = -r0 * solver.cauchyTravel() * solver.scaledNorm(G), loss = solver.loss, stuckometer = 0;

        case 40:
          if (!(0 === stuckometer)) {
            _context.next = 43;
            break;
          }

          _context.next = 43;
          return solver.report();

        case 43:
          _solver$computeNewton = solver.computeNewtonRegularized(0), _solver$computeNewton2 = (0, _slicedToArray2["default"])(_solver$computeNewton, 2), r = _solver$computeNewton2[0], dr = _solver$computeNewton2[1];

          if (isFinite(r)) {
            _context.next = 46;
            break;
          }

          throw new Error('Assertion failed.');

        case 46:
          gnInRadius = r <= R * (1 + rTol);

          if (gnInRadius) {
            _context.next = 85;
            break;
          }

          NORM.reset();

          for (i = N; i-- > 0;) {
            d = D[i];
            g = G[i];
            if (0 !== d) g /= d;
            NORM.include(g);
          }

          r -= R;

          if (isFinite(r)) {
            _context.next = 53;
            break;
          }

          throw new Error('Assertion failed.');

        case 53:
          if (0 < r) {
            _context.next = 55;
            break;
          }

          throw new Error('Assertion failed.');

        case 55:
          if (0 > dr) {
            _context.next = 57;
            break;
          }

          throw new Error('Assertion failed.');

        case 57:
          λMin = -r / dr, λMax = NORM.result / R, λ = 0;
          /*DEBUG*/

          nIter = -1;

        case 59:
          if (!(64 < ++nIter)) {
            _context.next = 61;
            break;
          }

          throw new Error('Assertion failed.');

        case 61:
          if (λMin >= 0) {
            _context.next = 63;
            break;
          }

          throw new Error('Assertion failed.');

        case 63:
          if (0 <= λMax) {
            _context.next = 65;
            break;
          }

          throw new Error('Assertion failed.');

        case 65:
          if (λMin < λMax) {
            _context.next = 67;
            break;
          }

          throw new Error('Assertion failed.');

        case 67:
          // Algorithm (5.5) (a)
          if (!(λMin < λ && λ < λMax)) λ = Math.max(lmLower * λMax, Math.sqrt(λMin * λMax));
          /*DEBUG*/

          if (λMin < λ && λ < λMax) {
            _context.next = 70;
            break;
          }

          throw new Error('Assertion failed.');

        case 70:
          _solver$computeNewton3 = solver.computeNewtonRegularized(λ);
          _solver$computeNewton4 = (0, _slicedToArray2["default"])(_solver$computeNewton3, 2);
          r = _solver$computeNewton4[0];
          dr = _solver$computeNewton4[1];

          if (isFinite(r)) {
            _context.next = 76;
            break;
          }

          throw new Error('Assertion failed.');

        case 76:
          if (!(Math.abs(R - r) <= R * rTol)) {
            _context.next = 78;
            break;
          }

          return _context.abrupt("break", 83);

        case 78:
          if (r < R) λMax = λ;
          λMin = Math.max(λMin, λ + (R - r) / dr); // Algorithm (5.5) (c)

          λ += (R - r) / dr * (r / R);

        case 81:
          _context.next = 59;
          break;

        case 83:
          if (Math.abs(solver.scaledNorm(dX) - R) <= R * rTol) {
            _context.next = 85;
            break;
          }

          throw new Error('Assertion failed.');

        case 85:
          if (solver.scaledNorm(dX) <= R * (1 + rTol)) {
            _context.next = 87;
            break;
          }

          throw new Error('Assertion failed.');

        case 87:
          _solver$considerMove = solver.considerMove(dX), _solver$considerMove2 = (0, _slicedToArray2["default"])(_solver$considerMove, 2), loss_predict = _solver$considerMove2[0], loss_consider = _solver$considerMove2[1]; //*DEBUG*/    if( !(loss_predict <= loss+1e-6) ) throw new Error('Assertion failed.');

          rScale = function rScale() {
            return 1;
          }; // const rScale = () => solver.scaledNorm(G); // <- TODO: examine scaling alternatives for the radius limits


          predict = loss - loss_predict, actual = loss - loss_consider;

          if (!(actual < predict * expectGainMin)) {
            _context.next = 103;
            break;
          }

          // PREDICTION BAD -> REDUCE TRUST RADIUS
          // (by fitting a quadratic polynomial through err_now, grad_now and err_next)
          grad = 0;

          for (_i = N; _i-- > 0;) {
            grad += dX[_i] * G[_i];
          }

          grad *= 2;
          shrink = grad / (2 * grad + loss - loss_consider);
          if (!(shrink >= shrinkLower)) shrink = shrinkLower;
          if (!(shrink <= shrinkUpper)) shrink = shrinkUpper;
          _r = Math.max(rMin * rScale(), R * shrink);

          if (!(_r === R && !(loss > loss_consider))) {
            _context.next = 100;
            break;
          }

          throw new _optimization_error.OptimizationNoProgressError();

        case 100:
          R = _r;
          _context.next = 104;
          break;

        case 103:
          if (gnInRadius) {
            // See [1] Algorithm (7.1) (d). If Gauss-Newton point is inside the trust radius and
            // prediction is not bad, set the radius to rGN times the distance to the Gauss-Newton
            // point. This also avoids that the trust radius grows uncontrollably while the
            // Gauss-Newton point is inside the trust region for a few iterations.
            rs = rScale();
            R = Math.max(rs * rMin, Math.min(rs * rMax, rNewton * r));
          } else if (actual > predict * expectGainMax || actual === 0) {
            // PREDICTION GREAT -> INCREASE TRUST RADIUS
            R = Math.min(rMax * rScale(), Math.max((0, _float64_utils.nextUp)(R), R * grow));
          }

        case 104:
          if (!(loss > loss_consider)) {
            _context.next = 110;
            break;
          }

          loss = loss_consider;
          solver.makeConsideredMove();
          stuckometer = 0;
          _context.next = 112;
          break;

        case 110:
          if (!(++stuckometer > stuckLimit)) {
            _context.next = 112;
            break;
          }

          throw new _optimization_error.OptimizationNoProgressError('Too many unsuccessfull iterations.');

        case 112:
          _context.next = 40;
          break;

        case 114:
        case "end":
          return _context.stop();
      }
    }
  }, _marked);
}

var tls_lm_gen = function tls_lm_gen(fjj, p0, dx0, opt) {
  return _lm(new _trust_region_solver_tls.TrustRegionSolverTLS(fjj, p0, dx0), opt);
};

exports.tls_lm_gen = tls_lm_gen;
var odr_lm_gen = (0, _trust_region_solver_tls.odr_gen)(_lm);
exports.odr_lm_gen = odr_lm_gen;

function fit_lm_gen(x, //: float[nSamples,nDim]
y, //: float[nSamples]
fg, //: (params: float[nParam]) => (x: float[nDim]) => [f: float, g: float[nParam]],
p0, //: float[nParam]
opt) {
  var FloatArray, _p0$shape, P, _x$shape, M, N, x_shape, R, J, RJ, fJ;

  return _regenerator["default"].wrap(function fit_lm_gen$(_context2) {
    while (1) {
      switch (_context2.prev = _context2.next) {
        case 0:
          if (fg instanceof Function) {
            _context2.next = 2;
            break;
          }

          throw new Error('fit_lm_gen(x,y, fg, p0): fg must be a function.');

        case 2:
          x = (0, _nd_array.asarray)(x); // <- TODO allow non-ndarray x ?

          y = (0, _nd_array.asarray)('float64', y);
          p0 = (0, _nd_array.asarray)('float64', p0);
          FloatArray = Float64Array;

          if (!(x.ndim !== 2)) {
            _context2.next = 8;
            break;
          }

          throw new Error('fit_lm_gen(x,y, fg, p0): x.ndim must be 2.');

        case 8:
          if (!(y.ndim !== 1)) {
            _context2.next = 10;
            break;
          }

          throw new Error('fit_lm_gen(x,y, fg, p0): y.ndim must be 1.');

        case 10:
          if (!(p0.ndim !== 1)) {
            _context2.next = 12;
            break;
          }

          throw new Error('fit_lm_gen(x,y, fg, p0): p0.ndim must be 1.');

        case 12:
          _p0$shape = (0, _slicedToArray2["default"])(p0.shape, 1), P = _p0$shape[0], _x$shape = (0, _slicedToArray2["default"])(x.shape, 2), M = _x$shape[0], N = _x$shape[1], x_shape = Int32Array.of(N);

          if (!(M != y.shape[0])) {
            _context2.next = 15;
            break;
          }

          throw new Error('fit_lm_gen(x,y, fg, p0): x.shape[0] must be equal to y.shape[0].');

        case 15:
          x = x.data.slice(); // <- TODO: TypedArray.subarray could be used here

          y = y.data; // if x is frozen, no protection copies are necessary while passing it to fg

          Object.freeze(x.buffer);
          x = Array.from({
            length: M
          }, function (_, i) {
            return Object.freeze(new _nd_array.NDArray(x_shape, x.subarray(N * i, N * i + 1)));
          });
          Object.freeze(x);
          R = new FloatArray(M), J = new FloatArray(M * P), RJ = [new _nd_array.NDArray(Int32Array.of(M), R), new _nd_array.NDArray(Int32Array.of(M, P), J)];

          fJ = function fJ(p) {
            var fgp = fg(p);
            if (!(fgp instanceof Function)) throw new Error();

            for (var i = M; i-- > 0;) {
              var _fgp = fgp(x[i]),
                  _fgp2 = (0, _slicedToArray2["default"])(_fgp, 2),
                  f = _fgp2[0],
                  g = _fgp2[1];

              f = (0, _nd_array.asarray)(f);
              g = (0, _nd_array.asarray)(g);
              if (f.ndim !== 0) throw new Error('fit_lm_gen(x,y, fg, p0): fg must have signature float[nParams] => float[nDim] => [float, float[nParams]].');
              if (g.ndim !== 1) throw new Error('fit_lm_gen(x,y, fg, p0): fg must have signature float[nParams] => float[nDim] => [float, float[nParams]].');
              if (g.shape[0] !== P) throw new Error('fit_lm_gen(x,y, fg, p0): fg must have signature float[nParams] => float[nDim] => [float, float[nParams]].');
              f = f.data[0];
              g = g.data;
              if (isNaN(f)) throw new Error('fit_lm_gen(x,y, fg, p0): NaN encountered.');
              if (g.some(isNaN)) throw new Error('fit_lm_gen(x,y, fg, p0): NaN encountered.');
              R[i] = f - y[i];

              for (var j = P; j-- > 0;) {
                J[P * i + j] = g[j];
              }
            }

            return RJ;
          };

          return _context2.delegateYield(lsq_lm_gen(fJ, p0, opt), "t0", 23);

        case 23:
        case "end":
          return _context2.stop();
      }
    }
  }, _marked2);
}