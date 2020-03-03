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
exports.lsq_lm_gen = lsq_lm_gen;
exports.fit_lm_gen = fit_lm_gen;

var _regenerator = _interopRequireDefault(require("@babel/runtime/regenerator"));

var _construct2 = _interopRequireDefault(require("@babel/runtime/helpers/construct"));

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _nd_array = require("../nd_array");

var _norm = require("../la/norm");

var _trust_region_solver = require("./_trust_region_solver");

var _marked = /*#__PURE__*/_regenerator["default"].mark(lsq_lm_gen),
    _marked2 = /*#__PURE__*/_regenerator["default"].mark(fit_lm_gen);

/** Solves a nonlinear least-squares problem using a
 *  trust-region variant Levenberg-Marquard algorithm
 *  as described in:
 *
 *    "THE LEVENBERG-MARQUARDT ALGORITHM: IMPLEMENTATION AND THEORY"
 *    by Jorge J. Moré.
 */
function lsq_lm_gen(fJ, //: (x0: float[nInputs]) => [f: float[nSamples], j: float[nSamples,nInputs]] - THE OPTIMIZATION FUNCTION AND ITS JACOBIAN
x0) {
  var _ref,
      _ref$r,
      r0,
      _ref$rMin,
      rMin,
      _ref$rMax,
      rMax,
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
      R,
      X,
      W,
      NORM,
      X_shape,
      mse_shape,
      _fJ,
      _fJ2,
      ff,
      JJ,
      solver,
      M,
      N,
      D,
      dX,
      G,
      f,
      J,
      err_now,
      i,
      s,
      _i,
      d,
      g,
      _solver$computeMinReg,
      _solver$computeMinReg2,
      r,
      dr,
      λMin,
      λMax,
      λ,
      nIter,
      _solver$computeMinReg3,
      _solver$computeMinReg4,
      err_predict,
      _i2,
      _s,
      j,
      _i3,
      _fJ3,
      _fJ4,
      _f,
      _J,
      err_next,
      _i4,
      _s2,
      predict,
      actual,
      grad_now,
      _i5,
      shrink,
      _ref2,
      _args = arguments;

  return _regenerator["default"].wrap(function lsq_lm_gen$(_context) {
    while (1) {
      switch (_context.prev = _context.next) {
        case 0:
          _ref = _args.length > 2 && _args[2] !== undefined ? _args[2] : {}, _ref$r = _ref.r0, r0 = _ref$r === void 0 ? 1e-2 : _ref$r, _ref$rMin = _ref.rMin, rMin = _ref$rMin === void 0 ? 1e-4 : _ref$rMin, _ref$rMax = _ref.rMax, rMax = _ref$rMax === void 0 ? Infinity : _ref$rMax, _ref$rTol = _ref.rTol, rTol = _ref$rTol === void 0 ? 0.05 : _ref$rTol, _ref$lmLower = _ref.lmLower, lmLower = _ref$lmLower === void 0 ? 0.001 : _ref$lmLower, _ref$shrinkLower = _ref.shrinkLower, shrinkLower = _ref$shrinkLower === void 0 ? 0.05 : _ref$shrinkLower, _ref$shrinkUpper = _ref.shrinkUpper, shrinkUpper = _ref$shrinkUpper === void 0 ? 0.95 : _ref$shrinkUpper, _ref$grow = _ref.grow, grow = _ref$grow === void 0 ? 1.5 : _ref$grow, _ref$expectGainMin = _ref.expectGainMin, expectGainMin = _ref$expectGainMin === void 0 ? 0.25 : _ref$expectGainMin, _ref$expectGainMax = _ref.expectGainMax, expectGainMax = _ref$expectGainMax === void 0 ? 0.75 : _ref$expectGainMax;
          x0 = (0, _nd_array.array)('float64', x0);

          if (!(x0.ndim !== 1)) {
            _context.next = 4;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0): x0.ndim must be 1.');

        case 4:
          if (fJ instanceof Function) {
            _context.next = 6;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0): fJ must be Function.');

        case 6:
          if (0 <= lmLower) {
            _context.next = 8;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0, opt): opt.lmLower must be within [0,1).');

        case 8:
          if (lmLower < 1) {
            _context.next = 10;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0, opt): opt.lmLower must be within [0,1).');

        case 10:
          if (0 < r0) {
            _context.next = 12;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0, opt): opt.r0 must be positive number.');

        case 12:
          if (0 < rMin) {
            _context.next = 14;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0, opt): opt.rMin must be positive number.');

        case 14:
          if (rMin <= rMax) {
            _context.next = 16;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0, opt): opt.rMin must not be greater than opt.rMax.');

        case 16:
          if (rMin <= r0) {
            _context.next = 18;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0, opt): opt.r0 must not be less than opt.rMin.');

        case 18:
          if (rMax >= r0) {
            _context.next = 20;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0, opt): opt.r0 must not be greater than opt.rMax.');

        case 20:
          if (0 < shrinkLower) {
            _context.next = 22;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0, opt): opt.shrinkLower must be positive number.');

        case 22:
          if (shrinkLower <= shrinkUpper) {
            _context.next = 24;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0, opt): opt.shrinkLower must not be greater than opt.shrinkUpper.');

        case 24:
          if (shrinkUpper < 1) {
            _context.next = 26;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0, opt): opt.shrinkUpper must be less than 1.');

        case 26:
          if (1 < grow) {
            _context.next = 28;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0, opt): opt.grow must be number greater than 1.');

        case 28:
          if (0 < expectGainMin) {
            _context.next = 30;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0, opt): opt.expectGainMin must be positive number.');

        case 30:
          if (expectGainMin < expectGainMax) {
            _context.next = 32;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0, opt): opt.expectGainMin must be less than expectGainMax.');

        case 32:
          if (!(expectGainMax < 1)) console.warn('lsq_dogleg_gen(fJ, x0, opt): opt.expectGainMax should be less than 1.');
          R = r0, X = x0.data, W = new Float64Array(X.length);
          NORM = new _norm.FrobeniusNorm();
          X_shape = x0.shape, mse_shape = new Int32Array(0);
          _fJ = fJ(new _nd_array.NDArray(X_shape, X.slice())), _fJ2 = (0, _slicedToArray2["default"])(_fJ, 2), ff = _fJ2[0], JJ = _fJ2[1];
          ff = (0, _nd_array.array)(ff);
          JJ = (0, _nd_array.array)(JJ);

          if (!(ff.ndim !== 1)) {
            _context.next = 41;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0): fJ(x) returned [f,J] with f.ndim !== 1.');

        case 41:
          if (!(JJ.ndim !== 2)) {
            _context.next = 43;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0): fJ(x) returned [f,J] with J.ndim !== 2.');

        case 43:
          solver = (0, _construct2["default"])(_trust_region_solver.TrustRegionSolverLSQ, (0, _toConsumableArray2["default"])(JJ.shape)), M = solver.M, N = solver.N, D = solver.D, dX = solver.X, G = solver.G;

          if (!(x0.shape[0] !== N)) {
            _context.next = 46;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0): fJ(x) returned [f,J] where J.shape[1] did not match x0.shape[0].');

        case 46:
          f = ff.data, J = JJ.data;
          err_now = 0;

          for (i = M; i-- > 0;) {
            s = f[i];
            err_now += s * s;
          }

        case 49:
          solver.update(ff, JJ);
          _context.next = 52;
          return [new _nd_array.NDArray(X_shape, X.slice()), new _nd_array.NDArray(mse_shape, Float64Array.of(err_now / M)), new _nd_array.NDArray(X_shape, G.map(function (g) {
            return g * 2 / M;
          })), new _nd_array.NDArray(ff.shape, ff.data.slice()), new _nd_array.NDArray(JJ.shape, JJ.data.slice())];

        case 52:
          solver.computeMinGlobal();

          if (solver.scaledNorm(dX) <= R * (1 + rTol)) {
            _context.next = 92;
            break;
          }

          NORM.reset();

          for (_i = N; _i-- > 0;) {
            d = D[_i];
            g = G[_i];
            if (0 !== d) g /= d;
            NORM.include(g);
          } // GLOBAL MIN OUTSIDE TRUST REGION -> FIND LAMBDA


          _solver$computeMinReg = solver.computeMinRegularized(0), _solver$computeMinReg2 = (0, _slicedToArray2["default"])(_solver$computeMinReg, 2), r = _solver$computeMinReg2[0], dr = _solver$computeMinReg2[1];
          r -= R;

          if (isFinite(r)) {
            _context.next = 60;
            break;
          }

          throw new Error('Assertion failed.');

        case 60:
          if (0 < r) {
            _context.next = 62;
            break;
          }

          throw new Error('Assertion failed.');

        case 62:
          if (0 > dr) {
            _context.next = 64;
            break;
          }

          throw new Error('Assertion failed.');

        case 64:
          λMin = -r / dr, λMax = NORM.result / R, λ = 0;
          /*DEBUG*/

          nIter = -1;

        case 66:
          if (!(8 < ++nIter)) {
            _context.next = 68;
            break;
          }

          throw new Error('Assertion failed.');

        case 68:
          if (λMin >= 0) {
            _context.next = 70;
            break;
          }

          throw new Error('Assertion failed.');

        case 70:
          if (0 <= λMax) {
            _context.next = 72;
            break;
          }

          throw new Error('Assertion failed.');

        case 72:
          if (λMin < λMax) {
            _context.next = 74;
            break;
          }

          throw new Error('Assertion failed.');

        case 74:
          // Algorithm (5.5) (a)
          if (!(λMin < λ && λ < λMax)) λ = Math.max(lmLower * λMax, Math.sqrt(λMin * λMax));
          /*DEBUG*/

          if (λMin < λ && λ < λMax) {
            _context.next = 77;
            break;
          }

          throw new Error('Assertion failed.');

        case 77:
          _solver$computeMinReg3 = solver.computeMinRegularized(λ);
          _solver$computeMinReg4 = (0, _slicedToArray2["default"])(_solver$computeMinReg3, 2);
          r = _solver$computeMinReg4[0];
          dr = _solver$computeMinReg4[1];

          if (isFinite(r)) {
            _context.next = 83;
            break;
          }

          throw new Error('Assertion failed.');

        case 83:
          if (!(Math.abs(R - r) <= R * rTol)) {
            _context.next = 85;
            break;
          }

          return _context.abrupt("break", 90);

        case 85:
          if (r < R) λMax = λ;
          λMin = Math.max(λMin, λ + (R - r) / dr); // Algorithm (5.5) (c)

          λ += (R - r) / dr * (r / R);

        case 88:
          _context.next = 66;
          break;

        case 90:
          if (Math.abs(solver.scaledNorm(dX) - R) <= R * rTol) {
            _context.next = 92;
            break;
          }

          throw new Error('Assertion failed.');

        case 92:
          if (solver.scaledNorm(dX) <= R * (1 + rTol)) {
            _context.next = 94;
            break;
          }

          throw new Error('Assertion failed.');

        case 94:
          err_predict = 0;

          for (_i2 = M; _i2-- > 0;) {
            _s = 0;

            for (j = N; j-- > 0;) {
              _s += J[N * _i2 + j] * dX[j];
            }

            _s += f[_i2];
            err_predict += _s * _s;
          }
          /*DEBUG*/


          if (err_predict <= err_now) {
            _context.next = 98;
            break;
          }

          throw new Error('Assertion failed.');

        case 98:
          for (_i3 = N; _i3-- > 0;) {
            W[_i3] = dX[_i3] + X[_i3];
          }

          _fJ3 = fJ(new _nd_array.NDArray(X_shape, W.slice())), _fJ4 = (0, _slicedToArray2["default"])(_fJ3, 2), _f = _fJ4[0], _J = _fJ4[1];
          _f = (0, _nd_array.array)(_f);
          _J = (0, _nd_array.array)(_J);

          if (!(_f.ndim !== 1)) {
            _context.next = 104;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0): fJ(x) returned [f,J] with f.ndim !== 1.');

        case 104:
          if (!(_J.ndim !== 2)) {
            _context.next = 106;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0): fJ(x) returned [f,J] with J.ndim !== 2.');

        case 106:
          if (!(M !== _f.shape[0])) {
            _context.next = 108;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0): fJ must return [f,J] with consistent shapes');

        case 108:
          if (!(M !== _J.shape[0])) {
            _context.next = 110;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0): fJ must return [f,J] with consistent shapes');

        case 110:
          if (!(N !== _J.shape[1])) {
            _context.next = 112;
            break;
          }

          throw new Error('lsq_dogleg_gen(fJ, x0): fJ must return [f,J] with consistent shapes');

        case 112:
          f = _f.data;
          err_next = 0;

          for (_i4 = M; _i4-- > 0;) {
            _s2 = f[_i4];
            err_next += _s2 * _s2;
          }

          predict = err_now - err_predict, actual = err_now - err_next;
          if (actual > predict * expectGainMax) R = Math.min(rMax, R * grow);else if (actual < predict * expectGainMin) {
            grad_now = 0;

            for (_i5 = N; _i5-- > 0;) {
              grad_now += dX[_i5] * G[_i5];
            }

            grad_now *= 2; //*DEBUG*/      const func = s => {
            //*DEBUG*/        const Z = new Float64Array(N);
            //*DEBUG*/        for( let i=N; i-- > 0; )
            //*DEBUG*/          Z[i] = X[i] + dX[i] * s;
            //*DEBUG*/        let [{data: f}] = fJ( new NDArray(X_shape,Z) );
            //*DEBUG*/        return f.reduce((x,y) => x + y*y, 0);
            //*DEBUG*/      };
            //*DEBUG*/      const grad = num_grad(func);
            //*DEBUG*/
            //*DEBUG*/      if( !(Math.abs(func(0) - err_now ) <= 1e-6) ) throw new Error('Assertion failed.');
            //*DEBUG*/      if( !(Math.abs(func(1) - err_next) <= 1e-6) ) throw new Error('Assertion failed.');
            //*DEBUG*/      if( !(Math.abs(grad(0) - grad_now) <= 1e-6) ) throw new Error('Assertion failed.');
            //*DEBUG*/
            //*DEBUG*/      // minimize quadratic polynomial through {err_now, grad_now, err_next} to compute shrink
            //*DEBUG*/      const a = err_now,
            //*DEBUG*/            b = grad_now,
            //*DEBUG*/            c = err_next - err_now - grad_now;
            //*DEBUG*/
            //*DEBUG*/      const poly = s => a + b*s + c*s*s,
            //*DEBUG*/            goly = s =>     b   + c*s*2;
            //*DEBUG*/
            //*DEBUG*/      if( !(Math.abs(goly(0) -grad_now ) <= 1e-6) ) throw new Error('Assertion failed.');
            //*DEBUG*/      if( !(Math.abs(poly(0) - err_now ) <= 1e-6) ) throw new Error('Assertion failed.');
            //*DEBUG*/      if( !(Math.abs(poly(1) - err_next) <= 1e-6) ) throw new Error('Assertion failed.');

            shrink = 0.5 * grad_now / (grad_now + err_now - err_next); //*DEBUG*/      if( !(Math.abs(goly(shrink)) <= 1e-6) ) {
            //*DEBUG*/        console.log({a,b,err_now,err_next,grad_now})
            //*DEBUG*/        throw new Error('Assertion failed: ' + Math.abs(grad(shrink)));
            //*DEBUG*/      }

            if (!(shrinkLower <= shrink)) shrink = shrinkLower;
            if (!(shrinkUpper >= shrink)) shrink = shrinkUpper;
            R = Math.max(rMin, R * shrink);
          } // ONLY ACCEPT NEW X IF BETTER

          if (err_now > err_next) {
            err_now = err_next;
            _ref2 = [X, W];
            W = _ref2[0];
            X = _ref2[1];
            ff = _f;
            JJ = _J;
          }

          f = ff.data;
          J = JJ.data;

        case 120:
          _context.next = 49;
          break;

        case 122:
        case "end":
          return _context.stop();
      }
    }
  }, _marked);
}

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