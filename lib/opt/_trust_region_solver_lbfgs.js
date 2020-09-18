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
exports.TrustRegionSolverLBFGS = void 0;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _classCallCheck2 = _interopRequireDefault(require("@babel/runtime/helpers/classCallCheck"));

var _createClass2 = _interopRequireDefault(require("@babel/runtime/helpers/createClass"));

var _nd_array = require("../nd_array");

var _norm = require("../la/norm");

var _lbfgsb_solver = require("./_lbfgsb_solver");

var _optimization_error = require("./optimization_error");

var REPORT_STATE_READY = 1,
    REPORT_STATE_NA = 2,
    REPORT_STATE_CONSIDER = 3;
/** Computes the dot product of two arrays.
 */

var dot = function dot(u, v) {
  if (u.length % 1 !== 0) throw new Error('Assertion failed.');
  if (v.length % 1 !== 0) throw new Error('Assertion failed.');
  if (u.length !== v.length) throw new Error('Assertion failed.');
  var uv = 0;

  for (var i = u.length; i-- > 0;) {
    uv += u[i] * v[i];
  }

  return uv;
};

var TrustRegionSolverLBFGS = /*#__PURE__*/function () {
  function TrustRegionSolverLBFGS(fg, x0, _ref) {
    var updateTol = _ref.updateTol,
        historySize = _ref.historySize,
        scaleInit = _ref.scaleInit;
    (0, _classCallCheck2["default"])(this, TrustRegionSolverLBFGS);
    if (!(fg instanceof Function)) throw new Error('new TrustRegionSolverLBFGS(fg, x0, updateTol): fg must be Function.');
    scaleInit *= 1;
    if (!(0 < scaleInit)) throw new Error('TrustRegionSolverLBFGS(fg,x0,opt): opt.scaleInit must be greater than 0.');
    if (!isFinite(scaleInit)) throw new Error('TrustRegionSolverLBFGS(fg,x0,opt): opt.scaleInit must be greater than 0.');
    this.fg = fg;
    this.report_x = x0 = (0, _nd_array.array)('float64', x0);
    if (x0.ndim !== 1) throw new Error('new TrustRegionSolverLBFGS(fg, x0, updateTol): x0.ndim must be 1.');

    var _fg = fg(new _nd_array.NDArray(x0.shape, x0.data.slice())),
        _fg2 = (0, _slicedToArray2["default"])(_fg, 2),
        f = _fg2[0],
        g = _fg2[1];

    this.report_loss = f *= 1;
    this.report_loss_grad = g = (0, _nd_array.array)('float64', g);
    if (isNaN(f)) throw new Error('new TrustRegionSolverLBFGS(fg, x0, updateTol): fg must have signature float[N] => [float,float[N]].');
    if (g.ndim !== 1) throw new Error('new TrustRegionSolverLBFGS(fg, x0, updateTol): fg must have signature float[N] => [float,float[N]].');

    var _x0$shape = (0, _slicedToArray2["default"])(x0.shape, 1),
        N = _x0$shape[0];

    if (g.shape[0] !== N) throw new Error('new TrustRegionSolverLBFGS(fg, x0, updateTol): fg must have signature float[N] => [float,float[N]].');
    this.M = historySize;
    this.N = N;
    this.scaleInit = scaleInit = 1 / scaleInit;
    this.lbfgs = new _lbfgsb_solver.LBFGSB_Solver(historySize, N);
    this.lbfgs.scale = scaleInit;
    this.D = new Float64Array(N);
    this.D.fill(Math.sqrt(scaleInit)); // <- TODO: improve/rethink ... this scaling is awful

    this.X0 = Float64Array.from(x0);
    this.G0 = Float64Array.from(g); // the last X,G that was used as LBFGS update

    this.last_X = Float64Array.from(x0);
    this.last_G = Float64Array.from(g);
    this.dg = new Float64Array(N);
    this.loss = f;
    this.bv = new Float64Array(2 * historySize);
    this.newton_dX = new Float64Array(N); // this.regularized_dX = new Float64Array(N);

    this.updateTol = updateTol;
    this._report_state = REPORT_STATE_READY;
    Object.seal(this);
  }

  (0, _createClass2["default"])(TrustRegionSolverLBFGS, [{
    key: "wiggle",
    value: function wiggle() {
      var lbfgs = this.lbfgs,
          m = lbfgs.m;
      if (0 === m) throw new _optimization_error.OptimizationNoProgressError('Too many unsuccessfull iterations.');
      if (!(0 < m)) throw new Error('Assertion failed.');
      lbfgs.forget(lbfgs.m + 1 >>> 1); // lbfgs.scale = scaleInit;
      // D.fill( Math.sqrt(scaleInit) );
    }
  }, {
    key: "report",
    value: function report() {
      if (this._report_state !== REPORT_STATE_READY) throw new Error('TrustRegionSolverLSQ::report: can only be called once after each makeConsideredMove() but not directly after considerMove(dX).');
      this._report_state = REPORT_STATE_NA;
      var result = [this.report_x, this.report_loss, this.report_loss_grad];
      this.report_x = this.report_loss_grad = null;
      this.report_loss = NaN;
      return result;
    }
  }, {
    key: "considerMove",
    value: function considerMove(dX) {
      var N = this.N,
          D = this.D,
          updateTol = this.updateTol,
          lbfgs = this.lbfgs,
          loss = this.loss,
          bv = this.bv,
          X0 = this.X0,
          G0 = this.G0,
          dg = this.dg,
          last_X = this.last_X,
          last_G = this.last_G,
          scaleSharpness = this.scaleSharpness;
      this._report_state = REPORT_STATE_CONSIDER;
      if (N !== dX.length) throw new Error('Assertion failed.');
      var X_shape = Int32Array.of(N),
          X = new Float64Array(N); // EVALUATE FUNCTION AT CONSIDERED X
      // compute considered x

      for (var i = N; i-- > 0;) {
        X[i] = X0[i] + dX[i];
      } // evaluate function at considered x


      var _this$fg = this.fg(new _nd_array.NDArray(X_shape, X.slice())),
          _this$fg2 = (0, _slicedToArray2["default"])(_this$fg, 2),
          f = _this$fg2[0],
          g = _this$fg2[1];

      f *= 1;
      g = (0, _nd_array.array)('float64', g);
      if (isNaN(f)) throw new Error('Assertion failed.');
      if (g.ndim !== 1) throw new Error('Assertion failed.');
      if (g.shape[0] !== N) throw new Error('Assertion failed.');
      this.report_x = new _nd_array.NDArray(X_shape, X);
      this.report_loss = f;
      this.report_loss_grad = g;
      g = g.data; // compute actual deltaX

      for (var _i = N; _i-- > 0;) {
        X[_i] -= X0[_i];
      } // COMPUTE PREDICTION


      var predict_loss = loss;

      for (var _i2 = N; _i2-- > 0;) {
        predict_loss += G0[_i2] * X[_i2];
      }

      lbfgs.compute_bv(X, bv);
      predict_loss += 0.5 * lbfgs.compute_ubbv(bv, dot(X, X), bv); // UPDATE LBFGSB MODEL

      ;
      {
        var dx = X;

        for (var _i3 = N; _i3-- > 0;) {
          dg[_i3] = g[_i3] - last_G[_i3];
        }

        for (var _i4 = N; _i4-- > 0;) {
          dx[_i4] = X0[_i4] + dX[_i4];
          dx[_i4] -= last_X[_i4];
        }

        var mx = dg.reduce(function (mx, x) {
          return Math.max(mx, Math.abs(x));
        }, 0); // <- avoid underflow

        if (0 < mx) {
          var gg = 0,
              gx = 0;

          for (var _i5 = N; _i5-- > 0;) {
            var xi = dx[_i5] / mx,
                gi = dg[_i5] / mx;
            gg += gi * gi;
            gx += xi * gi;
          }

          if (1 <= gg) // <- checks NaN
            {
              if (!isFinite(gx)) throw new Error('Assertion failed.');

              if (updateTol * gg < gx) {
                lbfgs.update(dx, dg);
                lbfgs.scale = gg / gx; // lbfgs.scale = Math.max(lbfgs.scale, gg/gx);
                // UPDATE SCALING

                for (var _i6 = N; _i6-- > 0;) {
                  lbfgs.compute_be(_i6, bv);
                  var B_ii = lbfgs.compute_ubbv(bv, 1, bv);
                  if (!(0 < B_ii)) throw new Error('Assertion failed: LBFGS model not positive definite.');
                  D[_i6] = Math.max(D[_i6], Math.sqrt(B_ii)); // <- TODO: improve/rethink ... this scaling is awful
                }

                for (var _i7 = N; _i7-- > 0;) {
                  last_X[_i7] = X0[_i7] + dX[_i7];
                }

                for (var _i8 = N; _i8-- > 0;) {
                  last_G[_i8] = g[_i8];
                }
              }
            }
        }
      } // compute considered x (again)

      for (var _i9 = N; _i9-- > 0;) {
        X[_i9] = X0[_i9] + dX[_i9];
      }

      return [predict_loss, this.report_loss];
    }
  }, {
    key: "makeConsideredMove",
    value: function makeConsideredMove() {
      if (this._report_state !== REPORT_STATE_CONSIDER) throw new Error('Assertion failed.');
      this._report_state = REPORT_STATE_READY; // discard consideration

      this.loss = this.report_loss;
      var N = this.N,
          X0 = this.X0,
          G0 = this.G0;
      var x0 = this.report_x.data,
          g0 = this.report_loss_grad.data;

      for (var i = N; i-- > 0;) {
        X0[i] = x0[i];
      }

      for (var _i10 = N; _i10-- > 0;) {
        G0[_i10] = g0[_i10];
      }
    }
  }, {
    key: "scaledNorm",
    value: function scaledNorm(X) {
      var N = this.N,
          D = this.D;
      if (X.length !== N) throw new Error('Assertion failed.');
      var norm = new _norm.FrobeniusNorm();

      for (var i = N; i-- > 0;) {
        norm.include(D[i] * X[i]);
      }

      return norm.result;
    }
  }, {
    key: "cauchyTravel",
    value: function cauchyTravel() {
      var N = this.N,
          lbfgs = this.lbfgs,
          G0 = this.G0,
          bv = this.bv; // polynomial along the gradient const1 + 2ax + bx² = const2 + (a/b + x)² * b

      var a = 0;

      for (var i = N; i-- > 0;) {
        a += G0[i] * G0[i];
      }

      lbfgs.compute_bv(G0, bv);
      var b = lbfgs.compute_ubbv(bv, dot(G0, G0), bv);
      if (isNaN(b)) throw new Error('Assertion failed: ' + b);
      return 0 === a ? 0 : 0 === b ? -Infinity : -a / b;
    }
  }, {
    key: "computeNewton",
    value: function computeNewton() {
      var N = this.N,
          lbfgs = this.lbfgs,
          G0 = this.G0,
          newton_dX = this.newton_dX;
      lbfgs.compute_Hv(G0, newton_dX);
      if (newton_dX.some(isNaN)) throw new Error('Assertion failed.');

      for (var i = N; i-- > 0;) {
        newton_dX[i] *= -1;
      }
    }
  }]);
  return TrustRegionSolverLBFGS;
}();

exports.TrustRegionSolverLBFGS = TrustRegionSolverLBFGS;