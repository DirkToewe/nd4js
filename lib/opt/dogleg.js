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
exports.min_dogleg_gen = min_dogleg_gen;
exports.fit_dogleg_gen = fit_dogleg_gen;
exports.odr_dogleg_gen = exports.tls_dogleg_gen = exports.lsq_dogleg_gen = void 0;

var _regenerator = _interopRequireDefault(require("@babel/runtime/regenerator"));

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _nd_array = require("../nd_array");

var _float64_utils = require("../dt/float64_utils");

var _optimization_error = require("./optimization_error");

var _polyquad = require("./polyquad");

var _trust_region_solver_lbfgs = require("./_trust_region_solver_lbfgs");

var _trust_region_solver_lsq = require("./_trust_region_solver_lsq");

var _trust_region_solver_tls = require("./_trust_region_solver_tls");

var _marked = /*#__PURE__*/_regenerator["default"].mark(_dogleg),
    _marked2 = /*#__PURE__*/_regenerator["default"].mark(fit_dogleg_gen);

function min_dogleg_gen(fg, x0) {
  var opt = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : {};
  var _opt$updateTol = opt.updateTol,
      updateTol = _opt$updateTol === void 0 ? 1e-14 : _opt$updateTol,
      _opt$historySize = opt.historySize,
      historySize = _opt$historySize === void 0 ? 8 : _opt$historySize,
      _opt$scaleInit = opt.scaleInit,
      scaleInit = _opt$scaleInit === void 0 ? 1e-2 : _opt$scaleInit;
  return _dogleg(new _trust_region_solver_lbfgs.TrustRegionSolverLBFGS(fg, x0, {
    updateTol: updateTol,
    historySize: historySize,
    scaleInit: scaleInit
  }), opt);
} // References
// ----------
// .. [1] "An Improved Optimization Method for iSAM2"
//         Rana Talaei Shahir and Hamid D. Taghirad Senior Member, IEEE
//         Proceeding of the 2nd RSI/ISM International Conference on Robotics and Mechatronics
//         October 15-17, 2014, Tehran, Iran
// .. [2] "THE LEVENBERG-MARQUARDT ALGORITHM: IMPLEMENTATION AND THEORY"
//         Jorge J. Moré.
// TODO implemented double-dogleg method as well as described in [1]
// TODO implement dogbox method as well which allows for box-constrained optimization

/** Solves a nonlinear least-squares problem using the dogleg method.
 */


var lsq_dogleg_gen = function lsq_dogleg_gen(fJ, x0, opt) {
  return _dogleg(new _trust_region_solver_lsq.TrustRegionSolverLSQ(fJ, x0), opt);
};

exports.lsq_dogleg_gen = lsq_dogleg_gen;

function _dogleg(solver) {
  var _ref,
      _ref$r,
      r0,
      _ref$rMin,
      rMin,
      _ref$rMax,
      rMax,
      _ref$rNewton,
      rNewton,
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
      newton_dX,
      G,
      dX,
      gNorm,
      cp,
      R,
      loss,
      stuckometer,
      t,
      i,
      gnInRadius,
      a,
      b,
      c,
      _i,
      d,
      u,
      v,
      _t,
      _i2,
      _v,
      _solver$considerMove,
      _solver$considerMove2,
      loss_predict,
      loss_consider,
      rScale,
      predict,
      actual,
      grad,
      _i3,
      shrink,
      r,
      distanceToNewton,
      rs,
      _args = arguments;

  return _regenerator["default"].wrap(function _dogleg$(_context) {
    while (1) {
      switch (_context.prev = _context.next) {
        case 0:
          _ref = _args.length > 1 && _args[1] !== undefined ? _args[1] : {}, _ref$r = _ref.r0, r0 = _ref$r === void 0 ? 1.1 : _ref$r, _ref$rMin = _ref.rMin, rMin = _ref$rMin === void 0 ? 0 : _ref$rMin, _ref$rMax = _ref.rMax, rMax = _ref$rMax === void 0 ? Infinity : _ref$rMax, _ref$rNewton = _ref.rNewton, rNewton = _ref$rNewton === void 0 ? 2 : _ref$rNewton, _ref$shrinkLower = _ref.shrinkLower, shrinkLower = _ref$shrinkLower === void 0 ? 0.05 : _ref$shrinkLower, _ref$shrinkUpper = _ref.shrinkUpper, shrinkUpper = _ref$shrinkUpper === void 0 ? 0.95 : _ref$shrinkUpper, _ref$grow = _ref.grow, grow = _ref$grow === void 0 ? 1.4142135623730951 : _ref$grow, _ref$expectGainMin = _ref.expectGainMin, expectGainMin = _ref$expectGainMin === void 0 ? 0.25 : _ref$expectGainMin, _ref$expectGainMax = _ref.expectGainMax, expectGainMax = _ref$expectGainMax === void 0 ? 0.75 : _ref$expectGainMax, _ref$stuckLimit = _ref.stuckLimit, stuckLimit = _ref$stuckLimit === void 0 ? 64 : _ref$stuckLimit;

          if (0 < r0) {
            _context.next = 3;
            break;
          }

          throw new Error('dogleg_gen(fJ, x0, opt): opt.r0 must be positive number.');

        case 3:
          if (0 <= rMin) {
            _context.next = 5;
            break;
          }

          throw new Error('dogleg_gen(fJ, x0, opt): opt.rMin must be non-negative.');

        case 5:
          if (rMin <= rMax) {
            _context.next = 7;
            break;
          }

          throw new Error('dogleg_gen(fJ, x0, opt): opt.rMin must not be greater than opt.rMax.');

        case 7:
          if (rMin <= r0) {
            _context.next = 9;
            break;
          }

          throw new Error('dogleg_gen(fJ, x0, opt): opt.r0 must not be less than opt.rMin.');

        case 9:
          if (rMax >= r0) {
            _context.next = 11;
            break;
          }

          throw new Error('dogleg_gen(fJ, x0, opt): opt.r0 must not be greater than opt.rMax.');

        case 11:
          if (1 < rNewton) {
            _context.next = 13;
            break;
          }

          throw new Error('dogleg_gen(fJ, x0, opt): opt.rGN must be greater than 1.');

        case 13:
          if (0 < shrinkLower) {
            _context.next = 15;
            break;
          }

          throw new Error('dogleg_gen(fJ, x0, opt): opt.shrinkLower must be positive number.');

        case 15:
          if (shrinkLower <= shrinkUpper) {
            _context.next = 17;
            break;
          }

          throw new Error('dogleg_gen(fJ, x0, opt): opt.shrinkLower must not be greater than opt.shrinkUpper.');

        case 17:
          if (shrinkUpper < 1) {
            _context.next = 19;
            break;
          }

          throw new Error('dogleg_gen(fJ, x0, opt): opt.shrinkUpper must be less than 1.');

        case 19:
          if (1 < grow) {
            _context.next = 21;
            break;
          }

          throw new Error('dogleg_gen(fJ, x0, opt): opt.grow must be number greater than 1.');

        case 21:
          if (0 < expectGainMin) {
            _context.next = 23;
            break;
          }

          throw new Error('dogleg_gen(fJ, x0, opt): opt.expectGainMin must be positive number.');

        case 23:
          if (expectGainMin < expectGainMax) {
            _context.next = 25;
            break;
          }

          throw new Error('dogleg_gen(fJ, x0, opt): opt.expectGainMin must be less than expectGainMax.');

        case 25:
          if (!(expectGainMax < 1)) console.warn('dogleg_gen(fJ, x0, opt): opt.expectGainMax should be less than 1.');

          if (!(0 !== stuckLimit % 1)) {
            _context.next = 28;
            break;
          }

          throw new Error('dogleg_gen(fJ, x0, opt): opt.stuckLimit must be an integer.');

        case 28:
          if (0 <= stuckLimit) {
            _context.next = 30;
            break;
          }

          throw new Error('dogleg_gen(fJ, x0, opt): opt.stuckLimit must not be negative.');

        case 30:
          N = solver.N, D = solver.D, newton_dX = solver.newton_dX, G = solver.G0, dX = new Float64Array(N);
          gNorm = solver.scaledNorm(G), cp = solver.cauchyTravel(), R = r0 * -gNorm * cp, loss = solver.loss, stuckometer = 0; // <- indicates how stuck we are, i.e. counts the number of consecutive operations without progress

        case 32:
          if (cp <= 0) {
            _context.next = 34;
            break;
          }

          throw new Error('Assertion failed.');

        case 34:
          if (!(0 === stuckometer)) {
            _context.next = 37;
            break;
          }

          _context.next = 37;
          return solver.report();

        case 37:
          t = Math.max(-R / gNorm, cp); // <- don't go beyond trust region radius

          for (i = N; i-- > 0;) {
            dX[i] = t * G[i];
          } // <- FIXME: Shouldn't the scaling be used here?!?!


          ;
          gnInRadius = false;

          if (!(-gNorm * cp < R)) {
            _context.next = 52;
            break;
          }

          solver.computeNewton(); // <- computes Gauss-Newton point
          // TODO In [2], a method is described to adjust trust radius when global min is
          //      inside trust region. Maybe it could be used in the Dogleg method as well.
          // let's find the travel t from cauchy point (cp) to gauss-newton point (gp) which intersects
          // with the trust region boundary. t=0 means cauchy point t=1 mean gauss-newton point.
          // The travel is the solution to following the quadratic equation:
          //   R² = a + 2bt + ct²

          a = 0, b = 0, c = 0;

          for (_i = N; _i-- > 0;) {
            d = D[_i], u = d * dX[_i], v = d * (newton_dX[_i] - dX[_i]);
            a += u * u;
            b += u * v;
            c += v * v;
          }

          a -= R * R;
          b *= 2;
          _t = (0, _polyquad.roots1d_polyquad)(a, b, c)[1];
          /*DEBUG*/

          if (0 <= _t) {
            _context.next = 50;
            break;
          }

          throw new Error('Assertion failed.');

        case 50:
          if (_t > 1) {
            // <- don't go beyond gauss-newton point (i.e. when it lies inside rust region)
            _t = 1;
            gnInRadius = true;
          }

          for (_i2 = N; _i2-- > 0;) {
            _v = newton_dX[_i2] - dX[_i2];
            dX[_i2] += _v * _t;
          }

        case 52:
          if (gnInRadius) {
            _context.next = 57;
            break;
          }

          if (Math.abs(solver.scaledNorm(dX) - R) <= R * 1e-5 + 1e-8) {
            _context.next = 55;
            break;
          }

          throw new Error('Assertion failed: ' + JSON.stringify({
            cp: cp,
            err: Math.abs(solver.scaledNorm(dX) - R),
            D: (0, _toConsumableArray2["default"])(D),
            G: (0, _toConsumableArray2["default"])(G),
            dX: (0, _toConsumableArray2["default"])(dX)
          }));

        case 55:
          _context.next = 59;
          break;

        case 57:
          if (solver.scaledNorm(dX) <= R * (1 + 1e-5)) {
            _context.next = 59;
            break;
          }

          throw new Error('Assertion failed.');

        case 59:
          // <- check that solution is in ellipsoid
          _solver$considerMove = solver.considerMove(dX), _solver$considerMove2 = (0, _slicedToArray2["default"])(_solver$considerMove, 2), loss_predict = _solver$considerMove2[0], loss_consider = _solver$considerMove2[1]; //*DEBUG*/    if( !(loss_predict <= loss+1e-12) ) throw new Error('Assertion failed: ' + JSON.stringify({loss, loss_consider, loss_predict, gnInRadius}));

          rScale = function rScale() {
            return 1;
          }; // const rScale = () => gNorm; // <- TODO: examine scaling alternatives


          predict = loss - loss_predict, actual = loss - loss_consider;

          if (!(actual < predict * expectGainMin)) {
            _context.next = 74;
            break;
          }

          // PREDICTION BAD -> REDUCE TRUST RADIUS
          // (by fitting a quadratic polynomial through err_now, grad_now and err_next)
          grad = 0;

          for (_i3 = N; _i3-- > 0;) {
            grad += dX[_i3] * G[_i3];
          }

          shrink = grad / (2 * grad + loss - loss_consider);
          if (!(shrink >= shrinkLower)) shrink = shrinkLower;
          if (!(shrink <= shrinkUpper)) shrink = shrinkUpper;
          r = Math.max(rMin * rScale(), R * shrink);

          if (!(r === R && !(loss > loss_consider))) {
            _context.next = 71;
            break;
          }

          throw new _optimization_error.OptimizationNoProgressError();

        case 71:
          R = r;
          _context.next = 75;
          break;

        case 74:
          if (gnInRadius) {
            // See [2] Algorithm (7.1) (d). If Gauss-Newton point is inside the trust radius and
            // prediction is not bad, set the radius to rGN times the distance to the Gauss-Newton
            // point. This also avoids that the trust radius grows uncontrollably while the
            // Gauss-Newton point is inside the trust region for a few iterations.
            distanceToNewton = solver.scaledNorm(dX), rs = rScale();
            R = Math.max(rMin * rs, Math.min(rMax * rs, rNewton * distanceToNewton));
          } else if (actual > predict * expectGainMax || actual === 0) {
            // PREDICTION GREAT -> INCREASE TRUST RADIUS
            R = Math.min(rMax * rScale(), Math.max((0, _float64_utils.nextUp)(R), R * grow));
          }

        case 75:
          // ONLY ACCEPT NEW X IF BETTER
          if (loss > loss_consider) {
            loss = loss_consider;
            solver.makeConsideredMove();
            stuckometer = 0;
          } else if (++stuckometer > stuckLimit) solver.wiggle();

          gNorm = solver.scaledNorm(G);
          cp = solver.cauchyTravel();

        case 78:
          _context.next = 32;
          break;

        case 80:
        case "end":
          return _context.stop();
      }
    }
  }, _marked);
}

var tls_dogleg_gen = function tls_dogleg_gen(fjj, p0, dx0, opt) {
  return _dogleg(new _trust_region_solver_tls.TrustRegionSolverTLS(fjj, p0, dx0), opt);
};

exports.tls_dogleg_gen = tls_dogleg_gen;
var odr_dogleg_gen = (0, _trust_region_solver_tls.odr_gen)(_dogleg);
exports.odr_dogleg_gen = odr_dogleg_gen;

function fit_dogleg_gen(x, //: float[nSamples,nDim]
y, //: float[nSamples]
fg, //: (params: float[nParam]) => (x: float[nDim]) => [f: float, g: float[nParam]],
p0, //: float[nParam]
opt) {
  var FloatArray, _p0$shape, P, _x$shape, M, N, x_shape, R, J, RJ, fJ;

  return _regenerator["default"].wrap(function fit_dogleg_gen$(_context2) {
    while (1) {
      switch (_context2.prev = _context2.next) {
        case 0:
          if (fg instanceof Function) {
            _context2.next = 2;
            break;
          }

          throw new Error('fit_dogleg_gen(x,y, fg, p0): fg must be a function.');

        case 2:
          x = (0, _nd_array.array)(x); // <- TODO allow non-ndarray x ?

          y = (0, _nd_array.array)('float64', y);
          p0 = (0, _nd_array.array)('float64', p0);
          FloatArray = Float64Array;

          if (!(x.ndim !== 2)) {
            _context2.next = 8;
            break;
          }

          throw new Error('fit_dogleg_gen(x,y, fg, p0): x.ndim must be 2.');

        case 8:
          if (!(y.ndim !== 1)) {
            _context2.next = 10;
            break;
          }

          throw new Error('fit_dogleg_gen(x,y, fg, p0): y.ndim must be 1.');

        case 10:
          if (!(p0.ndim !== 1)) {
            _context2.next = 12;
            break;
          }

          throw new Error('fit_dogleg_gen(x,y, fg, p0): p0.ndim must be 1.');

        case 12:
          _p0$shape = (0, _slicedToArray2["default"])(p0.shape, 1), P = _p0$shape[0], _x$shape = (0, _slicedToArray2["default"])(x.shape, 2), M = _x$shape[0], N = _x$shape[1], x_shape = Int32Array.of(N);

          if (!(M != y.shape[0])) {
            _context2.next = 15;
            break;
          }

          throw new Error('fit_dogleg_gen(x,y, fg, p0): x.shape[0] must be equal to y.shape[0].');

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
              if (f.ndim !== 0) throw new Error('fit_dogleg_gen(x,y, fg, p0): fg must have signature float[nParams] => float[nDim] => [float, float[nParams]].');
              if (g.ndim !== 1) throw new Error('fit_dogleg_gen(x,y, fg, p0): fg must have signature float[nParams] => float[nDim] => [float, float[nParams]].');
              if (g.shape[0] !== P) throw new Error('fit_dogleg_gen(x,y, fg, p0): fg must have signature float[nParams] => float[nDim] => [float, float[nParams]].');
              f = f.data[0];
              g = g.data;
              if (isNaN(f)) throw new Error('fit_dogleg_gen(x,y, fg, p0): NaN encountered.');
              if (g.some(isNaN)) throw new Error('fit_dogleg_gen(x,y, fg, p0): NaN encountered.');
              R[i] = f - y[i];

              for (var j = P; j-- > 0;) {
                J[P * i + j] = g[j];
              }
            }

            return RJ;
          };

          return _context2.delegateYield(lsq_dogleg_gen(fJ, p0, opt), "t0", 23);

        case 23:
        case "end":
          return _context2.stop();
      }
    }
  }, _marked2);
}