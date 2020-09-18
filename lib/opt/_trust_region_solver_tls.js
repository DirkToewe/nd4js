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
exports.odr_gen = odr_gen;
exports.TrustRegionSolverTLS = void 0;

var _classCallCheck2 = _interopRequireDefault(require("@babel/runtime/helpers/classCallCheck"));

var _createClass2 = _interopRequireDefault(require("@babel/runtime/helpers/createClass"));

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _nd_array = require("../nd_array");

var _giv_rot = require("../la/_giv_rot");

var _norm = require("../la/norm");

var _qr = require("../la/qr");

var _rrqr = require("../la/rrqr");

var _tri = require("../la/tri");

var _optimization_error = require("./optimization_error");

function odr_gen(trust_region) {
  var fit_odr_gen = function fit_odr_gen(x, y, fgg, p0) {
    var opt = arguments.length > 4 && arguments[4] !== undefined ? arguments[4] : {};
    if (!(fgg instanceof Function)) throw new Error("".concat(NAME, "(x,y, fgg,p0): fgg must be function."));
    x = (0, _nd_array.array)('float64', x);
    y = (0, _nd_array.array)('float64', y);
    p0 = (0, _nd_array.asarray)('float64', p0);
    var dx0 = 'dx0' in opt ? (0, _nd_array.array)('float64', opt.dx0) : x.mapElems('float64', function () {
      return 0;
    });
    var NAME = fit_odr_gen.name;
    if (p0.ndim !== 1) throw new Error("".concat(NAME, "(x,y, fgg,p0): p0.ndim must be 1."));
    if (x.ndim !== 1 && x.ndim !== 2) throw new Error("".concat(NAME, "(x,y, fgg,p0): x.ndim must be 1 or 2."));
    if (y.ndim !== 1 && y.ndim !== 2) throw new Error("".concat(NAME, "(x,y, fgg,p0): y.ndim must be 1 or 2."));
    if (x.ndim !== dx0.ndim) throw new Error("".concat(NAME, "(x,y, fgg,p0, opt): opt.dx0 and x must have same ndim."));

    var _x$shape = (0, _slicedToArray2["default"])(x.shape, 2),
        MX = _x$shape[0],
        _x$shape$ = _x$shape[1],
        NX = _x$shape$ === void 0 ? 1 : _x$shape$,
        _y$shape = (0, _slicedToArray2["default"])(y.shape, 2),
        MY = _y$shape[0],
        _y$shape$ = _y$shape[1],
        NY = _y$shape$ === void 0 ? 1 : _y$shape$,
        _p0$shape = (0, _slicedToArray2["default"])(p0.shape, 1),
        NP = _p0$shape[0];

    if (x.ndim == 2 && NX !== dx0.shape[1]) throw new Error("".concat(NAME, "(x,y, fgg,p0, opt): opt.dx0 and x must have same shape."));
    if (MX !== dx0.shape[0]) throw new Error("".concat(NAME, "(x,y, fgg,p0, opt): opt.dx0 and x must have same shape."));
    if (MX !== MY) throw new Error("".concat(NAME, "(x,y, fgg,p0, opt): x.shape[0] must equal y.shape[0]."));
    var x_ndim = x.ndim,
        y_ndim = y.ndim,
        xi_shape = x.shape.slice(1);
    var dy = new Float64Array(MX * NY),
        dy_dp = new Float64Array(MX * NY * NP),
        dy_dx = new Float64Array(MX * NY * NX),
        result_dy = new _nd_array.NDArray(y.shape, dy),
        result_dy_dp = new _nd_array.NDArray(Int32Array.of.apply(Int32Array, (0, _toConsumableArray2["default"])(y.shape).concat([NP])), dy_dp),
        result_dy_dx = new _nd_array.NDArray(Int32Array.of.apply(Int32Array, (0, _toConsumableArray2["default"])(y.shape).concat((0, _toConsumableArray2["default"])(x.shape.subarray(1)))), dy_dx);
    x = x.data;
    y = y.data;

    var fjj = function fjj(p, dx) {
      if (!p.dtype.startsWith('float')) throw new Error('Assertion failed.');
      if (!dx.dtype.startsWith('float')) throw new Error('Assertion failed.');
      if (p.ndim !== 1) throw new Error('Assertion failed.');
      if (p.shape[0] !== NP) throw new Error('Assertion failed.');
      if (x_ndim !== dx.ndim) throw new Error('Assertion failed.');
      if (x_ndim !== 1 && dx.shape[1] !== NX) throw new Error('Assertion failed.');
      if (dx.shape[0] !== MX) throw new Error('Assertion failed.');
      var fgg_p = fgg(p);
      dx = dx.data;

      for (var i = 0; i < MX; i++) {
        var xi = dx.slice(NX * i, NX * (i + 1));

        for (var j = NX; j-- > 0;) {
          xi[j] += x[NX * i + j];
        }

        var _fgg_p = fgg_p(new _nd_array.NDArray(xi_shape, xi)),
            _fgg_p2 = (0, _slicedToArray2["default"])(_fgg_p, 3),
            dyi = _fgg_p2[0],
            dyi_dp = _fgg_p2[1],
            dyi_dx = _fgg_p2[2];

        dyi = (0, _nd_array.asarray)('float', dyi);
        dyi_dp = (0, _nd_array.asarray)('float', dyi_dp);
        dyi_dx = (0, _nd_array.asarray)('float', dyi_dx);
        if (dyi.ndim !== y_ndim - 1) throw new Error('Assertion failed.');
        if (dyi_dp.ndim !== y_ndim) throw new Error('Assertion failed.');
        if (dyi_dx.ndim !== y_ndim - 2 + x_ndim) throw new Error('Assertion failed.');

        if (y_ndim !== 1) {
          if (dyi.shape[0] !== NY) throw new Error('Assertion failed.');
          if (dyi_dp.shape[0] !== NY) throw new Error('Assertion failed.');
          if (dyi_dx.shape[0] !== NY) throw new Error('Assertion failed.');
        }

        if (x_ndim !== 1 && dyi_dx.shape[y_ndim - 1] !== NX) throw new Error('Assertion failed.');
        if (dyi_dp.shape[y_ndim - 1] !== NP) throw new Error('Assertion failed.');
        dyi = dyi.data;
        dyi_dp = dyi_dp.data;
        dyi_dx = dyi_dx.data;

        for (var _j = NY; _j-- > 0;) {
          dy[NY * i + _j] = dyi[_j] - y[NY * i + _j];
        }

        for (var _j2 = NP * NY; _j2-- > 0;) {
          dy_dp[NP * NY * i + _j2] = dyi_dp[_j2];
        }

        for (var _j3 = NX * NY; _j3-- > 0;) {
          dy_dx[NX * NY * i + _j3] = dyi_dx[_j3];
        }
      } // we know that the TrustRegionSolverTLS performs protection copies so we can reuse memory


      return [result_dy, result_dy_dp, result_dy_dx];
    };

    return trust_region(new TrustRegionSolverTLS(fjj, p0, dx0), opt);
  };

  Object.defineProperty(fit_odr_gen, 'name', {
    value: "odr".concat(trust_region.name, "_gen")
  });
  return fit_odr_gen;
}

var REPORT_STATE_READY = 1,
    REPORT_STATE_NA = 2,
    REPORT_STATE_CONSIDER = 3;

var TrustRegionSolverTLS = /*#__PURE__*/function () {
  function TrustRegionSolverTLS(fgg, p0, dx0) {
    (0, _classCallCheck2["default"])(this, TrustRegionSolverTLS);
    if (!(fgg instanceof Function)) throw new Error('Assertion failed.');
    p0 = (0, _nd_array.array)('float64', p0);
    dx0 = (0, _nd_array.array)('float64', dx0);
    if (p0.ndim !== 1) throw new Error('Assertion failed.');
    if (dx0.ndim !== 1 && dx0.ndim !== 2) throw new Error('Assertion failed.');

    var _dx0$shape = (0, _slicedToArray2["default"])(dx0.shape, 2),
        MX = _dx0$shape[0],
        _dx0$shape$ = _dx0$shape[1],
        NX = _dx0$shape$ === void 0 ? 1 : _dx0$shape$,
        _p0$shape2 = (0, _slicedToArray2["default"])(p0.shape, 1),
        NP = _p0$shape2[0];

    var _fgg = fgg(new _nd_array.NDArray(p0.shape, p0.data.slice()), new _nd_array.NDArray(dx0.shape, dx0.data.slice())),
        _fgg2 = (0, _slicedToArray2["default"])(_fgg, 3),
        dy = _fgg2[0],
        dy_dp = _fgg2[1],
        dy_dx = _fgg2[2];

    dy = (0, _nd_array.array)('float64', dy);
    dy_dp = (0, _nd_array.array)('float64', dy_dp);
    dy_dx = (0, _nd_array.array)('float64', dy_dx);
    if (dy.ndim !== 1 && dy.ndim !== 2) throw new Error('Assertion failed.');
    if (dy_dp.ndim !== dy.ndim + 1) throw new Error('Assertion failed.');
    if (dy_dx.ndim !== dy.ndim + dx0.ndim - 1) throw new Error('Assertion failed.');

    var _dy$shape = (0, _slicedToArray2["default"])(dy.shape, 2),
        MY = _dy$shape[0],
        _dy$shape$ = _dy$shape[1],
        NY = _dy$shape$ === void 0 ? 1 : _dy$shape$;

    if (MX !== MY) throw new Error('Assertion failed.');
    if (dy_dp.shape[0] !== MX) throw new Error('Assertion failed.');
    if (dy_dx.shape[0] !== MX) throw new Error('Assertion failed.');

    if (dy.ndim === 2) {
      if (dy_dp.shape[1] !== NY) throw new Error('Assertion failed.');
      if (dy_dx.shape[1] !== NY) throw new Error('Assertion failed.');
    }

    if (dx0.ndim === 2 && dy_dx.shape[dy.ndim] !== NX) throw new Error('Assertion failed.');
    if (dy_dp.shape[dy.ndim] !== NP) throw new Error('Assertion failed.');

    var M = MX * NX + MX * NY,
        N = MX * NX + NP,
        L = Math.min(NX, NY),
        K = Math.min(MX * NX + NP + 1, MX * NY),
        // <- +1 as temp. memory for QR-decomp.
    J11 = new Float64Array(MX * NX),
        J21 = new Float64Array(MX * NY * NX),
        J22 = new Float64Array(MX * NY * NP),
        D = new Float64Array(N),
        X0 = new Float64Array(N),
        F0 = new Float64Array(M),
        G0 = new Float64Array(N),
        tmp = new Float64Array(L * (L + 3) >>> 1),
        // If NX < NY, we can use the prepare() step to reduce the work computeNewtonRegularized(λ) which should speed up Levenberg-Marquardt
    prepared_J21 = NY <= NX ? J21 : new Float64Array(MX * NX * NX + NX),
        // <- +NX as temp. memory for QR decomp.
    prepared_J22 = NY <= NX ? J22 : new Float64Array(K * NP),
        prepared_QF = NY <= NX ? F0.subarray(MX * NX) : new Float64Array(K),
        // Working memory and result of computeNewton()
    newton_R11 = new Float64Array(MX * NX),
        // <- after computeNewton(), contains diagonal of R11
    newton_R21 = new Float64Array(MX * L * NX),
        // <- after computeNewton(), contains sines of givens rotations used to compute off diagonal entries of R11 and R12
    newton_R22 = new Float64Array(K * NP),
        // <- after computeNewton(), contains R22 which is dense
    newton_P = new Int32Array(NP),
        newton_dX = new Float64Array(N),
        // Working memory and result of computeNewtonRegularized(λ)
    regularized_R11 = new Float64Array(MX * NX),
        regularized_R21 = new Float64Array(MX * L * NX),
        regularized_R22 = new Float64Array(Math.max(K, NP + 1) * NP),
        regularized_P = new Int32Array(NP),
        regularized_dX = new Float64Array(N),
        QF = new Float64Array(Math.max(MX * NX + K, N + 1)),
        norm = new Float64Array(2 * NP),
        _consider_J11 = J11.slice(),
        _consider_J21 = dy_dx.data,
        _consider_J22 = dy_dp.data;

    _consider_J11.fill(1.0);

    Object.assign(this, {
      MX: MX,
      NX: NX,
      NY: NY,
      NP: NP,
      M: M,
      N: N,
      loss: 0.0,
      rank: -1,
      _report_state: REPORT_STATE_NA,
      fgg: fgg,
      report_p: p0,
      report_dx: dx0,
      report_loss: NaN,
      report_dloss_dp: null,
      report_dloss_ddx: null,
      report_dy: dy,
      p_shape: p0.shape,
      x_shape: dx0.shape,
      y_shape: dy.shape,
      QF: QF,
      J11: J11,
      J21: J21,
      J22: J22,
      tmp: tmp,
      prepared_QF: prepared_QF,
      prepared_J21: prepared_J21,
      prepared_J22: prepared_J22,
      prepared: false,
      newton_R11: newton_R11,
      newton_R21: newton_R21,
      newton_R22: newton_R22,
      newton_P: newton_P,
      newton_dX: newton_dX,
      regularized_R11: regularized_R11,
      regularized_R21: regularized_R21,
      regularized_R22: regularized_R22,
      regularized_P: regularized_P,
      regularized_dX: regularized_dX,
      _consider_J11: _consider_J11,
      _consider_J21: _consider_J21,
      _consider_J22: _consider_J22,
      D: D,
      norm: norm,
      X0: X0,
      F0: F0,
      G0: G0
    });
    Object.seal(this);

    this._considerMove_computeLoss();

    this.makeConsideredMove();
  }

  (0, _createClass2["default"])(TrustRegionSolverTLS, [{
    key: "_considerMove_computeLoss",
    value: function _considerMove_computeLoss() {
      var M = this.M,
          MX = this.MX,
          NX = this.NX,
          NY = this.NY,
          NP = this.NP,
          x_shape = this.x_shape,
          p_shape = this.p_shape,
          J11 = this._consider_J11,
          J21 = this._consider_J21,
          J22 = this._consider_J22;
      this._report_state = REPORT_STATE_CONSIDER;
      var report_dloss_dp = new Float64Array(NP),
          report_dloss_ddx = new Float64Array(MX * NX),
          report_dx = this.report_dx.data,
          report_dy = this.report_dy.data; // COMPUTE LOSS GRADIENT w.r.t. P

      for (var i = MX * NY; i-- > 0;) {
        for (var j = NP; j-- > 0;) {
          report_dloss_dp[j] += report_dy[i] * J22[NP * i + j] / M * 2;
        }
      } // COMPUTE LOSS GRADIENT w.r.t. ΔX


      for (var _i = MX; _i-- > 0;) {
        for (var _j4 = NX; _j4-- > 0;) {
          report_dloss_ddx[NX * _i + _j4] = J11[NX * _i + _j4] * report_dx[NX * _i + _j4] / M * 2;
        }
      }

      for (var _i2 = MX; _i2-- > 0;) {
        for (var _j5 = NY; _j5-- > 0;) {
          for (var k = NX; k-- > 0;) {
            report_dloss_ddx[NX * _i2 + k] += J21[NX * (NY * _i2 + _j5) + k] * report_dy[NY * _i2 + _j5] / M * 2;
          }
        }
      } // COMPUTE LOSS (mean squared error)


      var report_loss = 0.0;

      for (var _i3 = MX * NX; _i3-- > 0;) {
        var s = report_dx[_i3];
        report_loss += s * s / M;
      }

      for (var _i4 = MX * NY; _i4-- > 0;) {
        var _s = report_dy[_i4];
        report_loss += _s * _s / M;
      }

      this.report_dloss_dp = new _nd_array.NDArray(p_shape, report_dloss_dp);
      this.report_dloss_ddx = new _nd_array.NDArray(x_shape, report_dloss_ddx);
      this.report_loss = report_loss;
    }
  }, {
    key: "considerMove",
    value: function considerMove(dX) {
      var M = this.M,
          N = this.N,
          MX = this.MX,
          NX = this.NX,
          NY = this.NY,
          NP = this.NP,
          x_shape = this.x_shape,
          y_shape = this.y_shape,
          p_shape = this.p_shape,
          _consider_J11 = this._consider_J11,
          _consider_J21 = this._consider_J21,
          _consider_J22 = this._consider_J22,
          J11 = this.J11,
          J21 = this.J21,
          J22 = this.J22,
          F0 = this.F0,
          X0 = this.X0,
          fgg = this.fgg;
      if (dX.length !== N) throw new Error('Assertion failed.');
      var report_p = new Float64Array(NP),
          report_dx = new Float64Array(MX * NX);
      this.report_p = new _nd_array.NDArray(p_shape, report_p);
      this.report_dx = new _nd_array.NDArray(x_shape, report_dx);

      for (var i = NP; i-- > 0;) {
        var I = MX * NX + i;
        report_p[i] = X0[I] + dX[I];
      }

      for (var _i5 = MX * NX; _i5-- > 0;) {
        report_dx[_i5] = X0[_i5] + dX[_i5];
      }

      var _fgg3 = fgg(new _nd_array.NDArray(p_shape, report_p.slice()), new _nd_array.NDArray(x_shape, report_dx.slice())),
          _fgg4 = (0, _slicedToArray2["default"])(_fgg3, 3),
          dy = _fgg4[0],
          dy_dp = _fgg4[1],
          dy_dx = _fgg4[2];

      dy = (0, _nd_array.array)('float64', dy);
      dy_dp = (0, _nd_array.asarray)('float64', dy_dp);
      dy_dx = (0, _nd_array.asarray)('float64', dy_dx);
      if (dy.ndim !== y_shape.length) throw new Error('Assertion failed.');
      if (dy_dp.ndim !== y_shape.length + 1) throw new Error('Assertion failed.');
      if (dy_dx.ndim !== y_shape.length + x_shape.length - 1) throw new Error('Assertion failed.');
      if (dy.shape[0] !== MX) throw new Error('Assertion failed.');
      if (dy_dp.shape[0] !== MX) throw new Error('Assertion failed.');
      if (dy_dx.shape[0] !== MX) throw new Error('Assertion failed.');
      if (dy_dp.shape[y_shape.length] !== NP) throw new Error('Assertion failed.');
      if (x_shape.length === 2 && dy_dx.shape[y_shape.length] !== NX) throw new Error('Assertion failed.');

      if (y_shape.length === 2) {
        if (dy.shape[1] !== NY) throw new Error('Assertion failed.');
        if (dy_dp.shape[1] !== NY) throw new Error('Assertion failed.');
        if (dy_dx.shape[1] !== NY) throw new Error('Assertion failed.');
      }

      this.report_dy = dy;
      dy_dx = dy_dx.data;
      dy_dp = dy_dp.data;

      for (var _i6 = MX * NY * NX; _i6-- > 0;) {
        _consider_J21[_i6] = dy_dx[_i6];
      }

      for (var _i7 = MX * NY * NP; _i7-- > 0;) {
        _consider_J22[_i7] = dy_dp[_i7];
      }

      _consider_J11.fill(1.0);

      this._considerMove_computeLoss();

      var predict_loss = 0.0;

      for (var _i8 = MX * NX; _i8-- > 0;) {
        var f = F0[_i8] + J11[_i8] * (report_dx[_i8] - X0[_i8]);
        predict_loss += f * f / M;
      }

      for (var _i9 = MX; _i9-- > 0;) {
        for (var j = NY; j-- > 0;) {
          var _f = F0[MX * NX + NY * _i9 + j];

          for (var k = NX; k-- > 0;) {
            _f += J21[NX * (NY * _i9 + j) + k] * (report_dx[NX * _i9 + k] - X0[NX * _i9 + k]);
          }

          for (var _k = NP; _k-- > 0;) {
            _f += J22[NP * (NY * _i9 + j) + _k] * (report_p[_k] - X0[MX * NX + _k]);
          }

          predict_loss += _f * _f / M;
        }
      }

      return [predict_loss, this.report_loss];
    }
  }, {
    key: "makeConsideredMove",
    value: function makeConsideredMove() {
      if (this._report_state !== REPORT_STATE_CONSIDER) throw new Error('Assertion failed.');
      this._report_state = REPORT_STATE_READY;
      this.loss = this.report_loss;
      this.rank = -1;
      this.prepared = false; // swap in consideration

      var MX = this.MX,
          NX = this.NX,
          NY = this.NY,
          NP = this.NP,
          _consider_J11 = this._consider_J11,
          J11 = this.J11,
          _consider_J21 = this._consider_J21,
          J21 = this.J21,
          _consider_J22 = this._consider_J22,
          J22 = this.J22,
          D = this.D,
          X0 = this.X0,
          F0 = this.F0,
          G0 = this.G0;

      for (var i = MX * NX; i-- > 0;) {
        J11[i] = _consider_J11[i];
      }

      for (var _i10 = MX * NY * NX; _i10-- > 0;) {
        J21[_i10] = _consider_J21[_i10];
      }

      for (var _i11 = MX * NY * NP; _i11-- > 0;) {
        J22[_i11] = _consider_J22[_i11];
      }

      ;
      {
        var p = this.report_p.data,
            dx = this.report_dx.data,
            dy = this.report_dy.data;

        for (var _i12 = NP; _i12-- > 0;) {
          X0[MX * NX + _i12] = p[_i12];
        }

        for (var _i13 = MX * NX; _i13-- > 0;) {
          X0[_i13] = dx[_i13];
        }

        for (var _i14 = MX * NY; _i14-- > 0;) {
          F0[MX * NX + _i14] = dy[_i14];
        }

        for (var _i15 = MX * NX; _i15-- > 0;) {
          F0[_i15] = X0[_i15];
        }
      } // COMPUTE GRADIENT OF (HALF) SQUARED ERROR (MIND THE HALF :P)

      for (var _i16 = MX * NX; _i16-- > 0;) {
        G0[_i16] = J11[_i16] * F0[_i16];
      }

      for (var _i17 = MX; _i17-- > 0;) {
        for (var j = NY; j-- > 0;) {
          for (var k = NX; k-- > 0;) {
            G0[NX * _i17 + k] += J21[NX * (NY * _i17 + j) + k] * F0[MX * NX + NY * _i17 + j];
          }
        }
      }

      G0.fill(0.0, MX * NX, MX * NX + NP);

      for (var _i18 = MX * NY; _i18-- > 0;) {
        for (var _j6 = NP; _j6-- > 0;) {
          G0[MX * NX + _j6] += J22[NP * _i18 + _j6] * F0[MX * NX + _i18];
        }
      }

      var norm = new _norm.FrobeniusNorm();

      for (var _i19 = MX; _i19-- > 0;) {
        for (var _k2 = NX; _k2-- > 0;) {
          norm.reset();
          norm.include(J11[NX * _i19 + _k2]);

          for (var _j7 = NY; _j7-- > 0;) {
            norm.include(J21[NX * (NY * _i19 + _j7) + _k2]);
          }

          var I = NX * _i19 + _k2;
          D[I] = Math.max(D[I], norm.result);
        }
      }

      for (var _j8 = NP; _j8-- > 0;) {
        norm.reset();

        for (var _i20 = MX * NY; _i20-- > 0;) {
          norm.include(J22[NP * _i20 + _j8]);
        }

        var J = MX * NX + _j8;
        D[J] = Math.max(D[J], norm.result);
      }
    }
  }, {
    key: "cauchyTravel",
    value: function cauchyTravel() {
      var N = this.N,
          MX = this.MX,
          NX = this.NX,
          NY = this.NY,
          NP = this.NP,
          J11 = this.J11,
          J21 = this.J21,
          J22 = this.J22,
          G = this.G0; // polynomial along the gradient const1 + 2ax + bx² = const2 + (a/b + x)² * b

      var a = 0,
          b = 0;

      for (var i = N; i-- > 0;) {
        a += G[i] * G[i];
      }

      for (var _i21 = MX; _i21-- > 0;) {
        for (var j = NY; j-- > 0;) {
          var Jg = 0;

          for (var k = NP; k-- > 0;) {
            Jg += J22[NP * (NY * _i21 + j) + k] * G[MX * NX + k];
          }

          for (var _k3 = NX; _k3-- > 0;) {
            Jg += J21[NX * (NY * _i21 + j) + _k3] * G[NX * _i21 + _k3];
          }

          b += Jg * Jg;
        }
      }

      for (var _i22 = MX * NX; _i22-- > 0;) {
        var _Jg = J11[_i22] * G[_i22];

        b += _Jg * _Jg;
      }

      return 0 === a ? 0 : 0 === b ? -Infinity : -a / b;
    }
  }, {
    key: "report",
    value: function report() {
      if (this._report_state !== REPORT_STATE_READY) throw new Error('TrustRegionSolverLSQ::report: can only be called once after each makeConsideredMove() but not directly after considerMove(dX).');
      this._report_state = REPORT_STATE_NA;
      var result = [this.report_p, this.report_dx, this.report_loss, this.report_dloss_dp, this.report_dloss_ddx, this.report_dy];
      this.report_p = this.report_dx = this.report_dy = this.report_dloss_dp = this.report_dloss_ddx = null;
      this.report_loss = NaN;
      return result;
    }
  }, {
    key: "wiggle",
    value: function wiggle() {
      throw new _optimization_error.OptimizationNoProgressError('Too many unsuccessfull iterations.');
    }
  }, {
    key: "__DEBUG_J",
    value: function __DEBUG_J(i, j) // <- meant for debugging only!
    {
      if (i % 1 !== 0) throw new Error('Assertion failed.');
      if (j % 1 !== 0) throw new Error('Assertion failed.');
      i |= 0;
      j |= 0;
      if (!(0 <= i)) throw new Error('Assertion failed.');
      if (!(0 <= j)) throw new Error('Assertion failed.');
      var M = this.M,
          N = this.N,
          MX = this.MX,
          NX = this.NX,
          NY = this.NY,
          NP = this.NP,
          J11 = this.J11,
          J21 = this.J21,
          J22 = this.J22;
      if (!(i < M)) throw new Error('Assertion failed.');
      if (!(j < N)) throw new Error('Assertion failed.');
      if (i < MX * NX) return i === j ? J11[i] : 0;
      i -= MX * NX;

      if (j < MX * NX) {
        var I = i / NY | 0,
            J = j / NX | 0;
        if (I !== J) return 0;
        i = i % NY;
        j = j % NX;
        return J21[NX * (NY * I + i) + j];
      }

      j -= MX * NX;
      return J22[NP * i + j];
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
    } // The (ODR) trust-region solver decomposes J into:
    //
    // J = Q·R·V·D
    //
    // Q: float[M,M], orthogonal
    // R: float[M,N], upper triangular
    // V: float[N,N], orthogonal
    // D: float[N,N], diagonal scaling matrix
    //
    // Where Q is only available implicitly via QF which is the product (Q·F).
    // This method returns [R,V,D] for debugging purposes, where D is returned
    // as row vector.

  }, {
    key: "_qr_decomp",
    value: function _qr_decomp(R11, R21, R22, P, QF) {
      var MX = this.MX,
          NX = this.NX,
          NY = this.NY,
          NP = this.NP,
          norm = this.norm,
          Q = this.tmp,
          J21 = this.prepared_J21,
          J22 = this.prepared_J22;
      var L = Math.min(NX, NY),
          K = Math.min(MX * NY, MX * NX + NP);
      if (!(P instanceof Int32Array)) throw new Error('Assertion failed.');
      if (P.length !== NP) throw new Error('Assertion failed.');
      if (R11.length !== MX * NX) throw new Error('Assertion failed.');
      if (R21.length !== MX * L * NX) throw new Error('Assertion failed.');
      if (R22.length !== (NP + 1) * NP && R22.length !== Math.min(MX * NY, MX * NX + NP + 1) * NP) throw new Error('Assertion failed: ' + JSON.stringify({
        len: R22.length,
        MX: MX,
        NX: NX,
        NP: NP
      })); //
      // STEP 2.1: ELIMINATE R21 USING GIVENS ROTATIONS
      //           O( MX*(NX+NP) ) operations
      //

      if (1 === L) {
        // THIS BRANCH IS PURELY FOR PERFORMANCE REASONS
        //   - JIT compilers seem to be unable to optimize for L===1
        //   - With this branch takes half as long as without it for L===1
        //   - TODO: remove once JIT compilers become better
        for (var i = 0; i < MX; i++) {
          var cc = 1;

          for (var j = 0; j < NX; j++) {
            var ij = NX * i + j,
                r1 = R11[ij],
                r2 = J21[ij] * cc,
                _giv_rot_qr2 = (0, _giv_rot._giv_rot_qr)(r1, r2),
                _giv_rot_qr3 = (0, _slicedToArray2["default"])(_giv_rot_qr2, 3),
                c = _giv_rot_qr3[0],
                s = _giv_rot_qr3[1],
                nrm = _giv_rot_qr3[2];

            if (0 === nrm) throw new Error('Assertion failed: Sparse part of J must not be singular.');
            R21[ij] = cc * s;
            if (s === 0) continue;
            cc *= c;
            R11[ij] = nrm;
            (0, _giv_rot._giv_rot_rows)(QF, 1, ij, MX * NX + i, c, s);
          }

          for (var _j9 = 0; _j9 < NP; _j9++) {
            R22[NP * i + _j9] = cc * J22[NP * i + _j9];
          }
        }
      } else {
        for (var _i23 = 0; _i23 < MX; _i23++) // <- for each block i
        {
          // Q keeps track of rotations in block
          Q.fill(0.0,
          /*start=*/
          L,
          /*end=*/
          L * (L + 3) >>> 1); // init Q to:
          //     ┏                  ┓
          //     ┃ 0                ┃
          //     ┃ 1                ┃
          //     ┃    1             ┃
          // Q = ┃       .          ┃
          //     ┃          .       ┃
          //     ┃             .    ┃
          //     ┃                1 ┃
          //     ┗                  ┛
          // (Q[1:] is stored sparsly as the upper off-diagonal entries will always be 0)

          for (var _i24 = L; _i24 > 0; _i24--) {
            Q[L - 1 + (_i24 * (_i24 + 1) >>> 1)] = 1;
          }

          for (var k = 0; k < NX; k++) // <- for each column in block i
          {
            Q.fill(0.0, 0, L);
            var _r = R11[NX * _i23 + k];

            for (var _j10 = 0; _j10 < L; _j10++) // <- for each entry j in column k of block i
            {
              var Q_off = L + (_j10 * (_j10 + 1) >>> 1);
              var _r2 = 0.0;

              for (var l = -1; l++ < _j10;) {
                _r2 += Q[Q_off + l] * J21[NX * (L * _i23 + l) + k];
              }

              if (0 === _r2) continue;

              var _giv_rot_qr4 = (0, _giv_rot._giv_rot_qr)(_r, _r2),
                  _giv_rot_qr5 = (0, _slicedToArray2["default"])(_giv_rot_qr4, 3),
                  _c = _giv_rot_qr5[0],
                  _s2 = _giv_rot_qr5[1],
                  _nrm = _giv_rot_qr5[2];

              _r = _nrm;
              if (0 === _s2) continue;
              (0, _giv_rot._giv_rot_rows)(Q, _j10 + 1, 0, Q_off, _c, _s2);
              (0, _giv_rot._giv_rot_rows)(QF, 1, NX * _i23 + k, MX * NX + L * _i23 + _j10, _c, _s2);
            }

            R11[NX * _i23 + k] = _r; // write finished row of Q to R21

            for (var _j11 = 0; _j11 < L; _j11++) {
              R21[NX * (L * _i23 + _j11) + k] = Q[_j11];
            }
          } // apply Q to J22


          R22.fill(0.0, NP * L * _i23, NP * L * (_i23 + 1));

          for (var _j12 = 0; _j12 < L; _j12++) {
            var _Q_off = L + (_j12 * (_j12 + 1) >>> 1);

            for (var _k4 = 0; _k4 <= _j12; _k4++) {
              for (var _l = 0; _l < NP; _l++) {
                R22[NP * (L * _i23 + _j12) + _l] += Q[_Q_off + _k4] * J22[NP * (L * _i23 + _k4) + _l];
              }
            }
          }
        }
      } // copy remaining rows of R22


      for (var _i25 = MX * L * NP; _i25 < K * NP; _i25++) {
        R22[_i25] = J22[_i25];
      } //
      // STEP 2.3: RRQR-DECOMPOSE R22
      //           O( min(MX,NP)*NP*MX ) operations
      //


      for (var _i26 = NP; _i26-- > 0;) {
        P[_i26] = _i26;
      }

      (0, _rrqr._rrqr_decomp_inplace)(K, NP, 1, R22, 0, QF, MX * NX, P, 0, norm);
      var rnk = (0, _rrqr._rrqr_rank)(K, NP, R22, 0, norm);
      R22.fill(0.0, NP * rnk, NP * Math.min(K, NP)); // <- zero out the rank-deficient rows after RRQR

      return rnk;
    }
  }, {
    key: "_qr_solve",
    value: function _qr_solve(R11, R21, R22, P, rnk, X) {
      var MX = this.MX,
          NX = this.NX,
          NY = this.NY,
          NP = this.NP,
          J21 = this.prepared_J21,
          J22 = this.prepared_J22,
          D = this.D,
          Jx = this.tmp,
          tmp = this.norm;
      var L = Math.min(NX, NY);
      if (!(P instanceof Int32Array)) throw new Error('Assertion failed.');
      if (P.length !== NP) throw new Error('Assertion failed.');
      if (R11.length !== MX * NX) throw new Error('Assertion failed.');
      if (R21.length !== MX * L * NX) throw new Error('Assertion failed.');
      if (!(R22.length >= MX * L * NP)) throw new Error('Assertion failed.');
      if (0 !== rnk % 1) throw new Error('Assertion failed.');
      if (!(0 <= rnk)) throw new Error('Assertion failed.');
      if (!(NP >= rnk)) throw new Error('Assertion failed.');
      rnk |= 0; //
      // STEP 3.1: BACKWARDS SUBSITUTION OF R2-PART
      //           O( rnk² ) operations
      //

      (0, _tri._triu_solve)(rnk, NP, 1, R22, 0, X, MX * NX);

      if (rnk != NP) {
        //
        // STEP 3.2: APPLY GIVENS ROTATIONS TO X (UNDO 2.4)
        //           O( rnk * (NP-rnk) ) operations
        //
        for (var i = 0; i < rnk; i++) {
          for (var j = rnk; j < NP; j++) {
            var s = R22[NP * i + j];
            if (s === 0) continue;
            var c = Math.sqrt(1 - s * s);
            (0, _giv_rot._giv_rot_rows)(X, 1, MX * NX + j, MX * NX + i, c, s);
          }
        }
      } //
      // STEP 3.3: APPLY COLUMN PERMUTATIONS P TO X (UNDO 2.3)
      //           O( NP ) operations
      //


      for (var _i27 = NP; _i27-- > 0;) {
        tmp[P[_i27]] = X[MX * NX + _i27];
      }

      for (var _i28 = NP; _i28-- > 0;) {
        X[MX * NX + _i28] = tmp[_i28];
      } // factor out scaling


      if (rnk !== NP) {
        for (var _i29 = NP; _i29-- > 0;) {
          var d = D[MX * NX + _i29];
          if (0 !== d) X[MX * NX + _i29] /= d;
        }
      } //
      // STEP 3.4: MOVE SOLVED R22-PART TO THE RIGHT
      //           O( MX*(NX+NP) ) operations
      // -------------------------------------------
      // As a result of steps 2.1 and 2.2, each row in R12 can be described as
      // a scaled row of J22, i.e.:
      //
      // R12[MX*j+i,:] = s[i,j] * J22[i,:]
      //
      // Where:
      //
      // s[i,j] = R21[i,j] * c[i,j-1]
      // s[i,0] = R21[i,0]
      // c[i,j] = sqrt( 1 - (s[i,j])² )
      //


      for (var _i30 = MX; _i30-- > 0;) {
        for (var _j13 = L; _j13-- > 0;) {
          var xj = 0.0;

          for (var k = NP; k-- > 0;) {
            xj += J22[NP * (L * _i30 + _j13) + k] * X[MX * NX + k];
          }

          for (var _k5 = NX; _k5-- > 0;) {
            X[NX * _i30 + _k5] -= xj * R21[NX * (L * _i30 + _j13) + _k5];
          }
        }
      } // TODO: The following can be done more efficient as the off-diagonal rows of R11 are scaled rows of J21
      //
      // STEP 3.5: BACKWARD SUBSTITUTION OF R11
      //           O( MX*NX ) operations


      for (var _i31 = MX; _i31-- > 0;) // <- for each block bottom to top
      {
        Jx.fill(0.0, 0, L);

        for (var _k6 = NX; _k6-- > 0;) // <- for each diagonal entry in block bottom to top
        {
          var ik = NX * _i31 + _k6;

          for (var _j14 = L; _j14-- > 0;) {
            X[ik] -= R21[NX * (L * _i31 + _j14) + _k6] * Jx[_j14];
          } // <- move partial solution to the right side


          X[ik] /= R11[ik];

          for (var _j15 = L; _j15-- > 0;) {
            Jx[_j15] += J21[NX * (L * _i31 + _j15) + _k6] * X[ik];
          } // <- Jx accumulates the part to be moved to the right side

        }
      }
    }
  }, {
    key: "_rt_solve",
    value: function _rt_solve(R11, R21, R22, P, rnk, X) {
      var MX = this.MX,
          NX = this.NX,
          NY = this.NY,
          NP = this.NP,
          J21 = this.prepared_J21,
          J22 = this.prepared_J22,
          D = this.D,
          sx = this.tmp,
          tmp = this.norm;
      var L = Math.min(NX, NY);
      if (!(P instanceof Int32Array)) throw new Error('Assertion failed.');
      if (P.length !== NP) throw new Error('Assertion failed.');
      if (R11.length !== MX * NX) throw new Error('Assertion failed.');
      if (R21.length !== MX * L * NX) throw new Error('Assertion failed.');
      if (!(R22.length >= MX * NP)) throw new Error('Assertion failed.');
      if (!(X.length >= MX * NX + NP)) throw new Error('Assertion failed.');
      if (0 !== rnk % 1) throw new Error('Assertion failed.');
      if (!(0 <= rnk)) throw new Error('Assertion failed.');
      if (!(NP >= rnk)) throw new Error('Assertion failed.');
      rnk |= 0; //
      // STEP 4.1: SOLVE SPARSE PART
      //           O( MX*(NX+NP) ) operations
      //

      tmp.fill(0.0, 0, NP); // <- accumulates sparse part of solution to be moved to right side (because we have to apply givens rotations to it before moving it to the right side)

      for (var i = 0; i < MX; i++) {
        sx.fill(0.0, 0, L); // forward substitute block

        for (var k = 0; k < NX; k++) {
          var ik = NX * i + k;

          for (var j = 0; j < L; j++) {
            X[ik] -= J21[NX * (L * i + j) + k] * sx[j];
          }

          X[ik] /= R11[ik];

          for (var _j16 = 0; _j16 < L; _j16++) {
            sx[_j16] += R21[NX * (L * i + _j16) + k] * X[ik];
          }
        } // move solved block to right hand side


        for (var _j17 = 0; _j17 < L; _j17++) {
          for (var _k7 = 0; _k7 < NP; _k7++) {
            tmp[_k7] += sx[_j17] * J22[NP * (L * i + _j17) + P[_k7]];
          }
        }
      }

      if (rnk < NP) {
        for (var _i32 = NP; _i32-- > 0;) {
          var d = D[MX * NX + P[_i32]];
          if (0 !== d) tmp[_i32] /= d;
        }

        for (var _i33 = rnk; _i33-- > 0;) {
          for (var _j18 = NP; _j18-- > rnk;) {
            var s = R22[NP * _i33 + _j18];
            if (s === 0) continue;
            var c = Math.sqrt(1 - s * s);
            (0, _giv_rot._giv_rot_rows)(tmp, 1, _i33, _j18, c, s);
          }
        }
      }

      for (var _j19 = rnk; _j19-- > 0;) {
        X[MX * NX + _j19] -= tmp[_j19];
      }

      (0, _tri._triu_t_solve)(rnk, NP, 1, R22, 0, X, MX * NX);
    }
  }, {
    key: "prepare",
    value: function prepare() {
      var MX = this.MX,
          NX = this.NX,
          NY = this.NY,
          NP = this.NP,
          J21 = this.prepared_J21,
          J22 = this.prepared_J22,
          QF = this.prepared_QF,
          raw_J21 = this.J21,
          raw_J22 = this.J22,
          raw_F = this.F0;
      this.prepared = this.prepared || NX >= NY; // <- FIXME: NX >= NY

      if (this.prepared) return;
      this.prepared = true; // IF BLOCKS IN J21 ARE PORTRAIT-SHAPED, WE CAN PRECOMPUTE SOME WORK TO MAKE computeNewtonRegularized(λ) CHEAPER
      // ┏           ╷   ┓      ┏           ╷   ┓
      // ┃ ╲         ┊   ┃      ┃ ╲         ┊   ┃
      // ┃   ╲       ┊   ┃      ┃   ╲       ┊   ┃
      // ┃     ╲     ┊ 0 ┃      ┃     ╲     ┊ 0 ┃
      // ┃       ╲   ┊   ┃      ┃       ╲   ┊   ┃
      // ┃         ╲ ┊   ┃      ┃         ╲ ┊   ┃
      // ┃┄┄┄┄┄┄┄┄┄┄┄┼┄┄┄┃      ┃┄┄┄┄┄┄┄┄┄┄┄┼┄┄┄┃
      // ┃ ██        ┊██ ┃      ┃ ▜█        ┊██ ┃
      // ┃ ██        ┊██ ┃  =>  ┃  ▜        ┊██ ┃
      // ┃ ██        ┊██ ┃      ┃   ▜█      ┊██ ┃
      // ┃   ██      ┊██ ┃      ┃    ▜      ┊██ ┃
      // ┃   ██      ┊██ ┃      ┃     .     ┊██ ┃
      // ┃   ██      ┊██ ┃      ┃      .    ┊██ ┃
      // ┃     .     ┊██ ┃      ┃       .   ┊██ ┃
      // ┃      .    ┊██ ┃      ┃         ▜█┊██ ┃
      // ┃       .   ┊██ ┃      ┃          ▜┊██ ┃
      // ┃         ██┊██ ┃      ┃           ┊▜█ ┃
      // ┃         ██┊██ ┃      ┃           ┊ ▜ ┃
      // ┃         ██┊██ ┃      ┃           ┊   ┃
      // ┗           ╵   ┛      ┗           ╵   ┛

      if (J21.length !== MX * NX * NX + NX) throw new Error('Assertion failed.');
      if (J22.length !== Math.min(MX * NY, MX * NX + NP + 1) * NP) throw new Error('Assertion failed.');
      if (QF.length !== Math.min(MX * NY, MX * NX + NP + 1)) throw new Error('Assertion failed.');

      for (var l = -1, i = 0; i < MX; i++) {
        var J21_off = NX * NX * i,
            J22_off = NP * NX * i,
            QF_off = NX * i; // copy upper square region

        for (var j = 0; j < NX * NX; j++) {
          J21[J21_off + j] = raw_J21[NX * NY * i + j];
        }

        for (var _j20 = 0; _j20 < NX * NP; _j20++) {
          J22[J22_off + _j20] = raw_J22[NP * NY * i + _j20];
        }

        for (var _j21 = 0; _j21 < NX; _j21++) {
          QF[QF_off + _j21] = raw_F[MX * NX + NY * i + _j21];
        } // QR decomp. square region


        for (var _j22 = 1; _j22 < NX; _j22++) {
          for (var k = 0; k < _j22; k++) {
            var jk = J21_off + NX * _j22 + k,
                R_jk = J21[jk];
            if (R_jk === 0) continue;

            var kk = J21_off + NX * k + k,
                R_kk = J21[kk],
                _giv_rot_qr6 = (0, _giv_rot._giv_rot_qr)(R_kk, R_jk),
                _giv_rot_qr7 = (0, _slicedToArray2["default"])(_giv_rot_qr6, 3),
                c = _giv_rot_qr7[0],
                s = _giv_rot_qr7[1],
                nrm = _giv_rot_qr7[2];

            J21[jk] = 0;
            if (0 === s) continue;
            J21[kk] = nrm;
            (0, _giv_rot._giv_rot_rows)(J21, NX - 1 - k, kk + 1, jk + 1, c, s);
            (0, _giv_rot._giv_rot_rows)(J22, NP, J22_off + k * NP, J22_off + _j22 * NP, c, s);
            (0, _giv_rot._giv_rot_rows)(QF, 1, QF_off + k, QF_off + _j22, c, s);
          }
        } // QR decomp. remaining rows


        for (var _j23 = NX; _j23 < NY; _j23++) {
          if (++l === NP) // <- more than NP rows at bottom of J22 -> start QR decomp.
            (0, _qr._qr_decomp_inplace)(NP, NP, 1, J22, MX * NX * NP, QF, MX * NX);
          l = Math.min(l, NP); // copy row

          for (var _k8 = 0; _k8 < NX; _k8++) {
            J21[J21_off + NX * NX + _k8] = raw_J21[NX * NY * i + NX * _j23 + _k8];
          }

          for (var _k9 = 0; _k9 < NP; _k9++) {
            J22[(MX * NX + l) * NP + _k9] = raw_J22[NP * NY * i + NP * _j23 + _k9];
          }

          QF[MX * NX + l] = raw_F[MX * NX + NY * i + _j23]; // eliminate entries in J21

          for (var _k10 = 0; _k10 < NX; _k10++) {
            var _jk = J21_off + NX * NX + _k10,
                _R_jk = J21[_jk];

            if (_R_jk === 0) continue;

            var _kk = J21_off + NX * _k10 + _k10,
                _R_kk = J21[_kk],
                _giv_rot_qr8 = (0, _giv_rot._giv_rot_qr)(_R_kk, _R_jk),
                _giv_rot_qr9 = (0, _slicedToArray2["default"])(_giv_rot_qr8, 3),
                _c2 = _giv_rot_qr9[0],
                _s3 = _giv_rot_qr9[1],
                _nrm2 = _giv_rot_qr9[2];

            J21[_jk] = 0;
            if (0 === _s3) continue;
            J21[_kk] = _nrm2;
            (0, _giv_rot._giv_rot_rows)(J21, NX - 1 - _k10, _kk + 1, _jk + 1, _c2, _s3);
            (0, _giv_rot._giv_rot_rows)(J22, NP, J22_off + NP * _k10, (MX * NX + l) * NP, _c2, _s3);
            (0, _giv_rot._giv_rot_rows)(QF, 1, QF_off + _k10, MX * NX + l, _c2, _s3);
          }

          if (!(l <= NP)) throw new Error('Assertion failed.'); // more than NP trailing rows in J22 -> eliminate further rows using QR decomp.

          if (l === NP) for (var _k11 = 0; _k11 < NP; _k11++) {
            var lk = MX * NX * NP + NP * l + _k11,
                R_lk = J22[lk];
            if (R_lk === 0) continue;

            var _kk2 = MX * NX * NP + NP * _k11 + _k11,
                _R_kk2 = J22[_kk2],
                _giv_rot_qr10 = (0, _giv_rot._giv_rot_qr)(_R_kk2, R_lk),
                _giv_rot_qr11 = (0, _slicedToArray2["default"])(_giv_rot_qr10, 3),
                _c3 = _giv_rot_qr11[0],
                _s4 = _giv_rot_qr11[1],
                _nrm3 = _giv_rot_qr11[2];

            J22[lk] = 0;
            if (0 === _s4) continue;
            J22[_kk2] = _nrm3;
            (0, _giv_rot._giv_rot_rows)(J22, NP - 1 - _k11, _kk2 + 1, lk + 1, _c3, _s4);
            (0, _giv_rot._giv_rot_rows)(QF, 1, MX * NX + _k11, MX * NX + l, _c3, _s4);
          }
        }
      }
    } // The Jacobian of the orthogonal least squares problem is sparse with the following structure:
    //
    //     ┏                   ╷    ┓   ┏                 ╷     ┓   
    //     ┃  ╲                ┊    ┃   ┃                 ┊     ┃
    //     ┃    ╲              ┊    ┃   ┃                 ┊     ┃
    //     ┃      ╲            ┊    ┃   ┃                 ┊     ┃
    //     ┃        ╲          ┊ 0  ┃   ┃       J11       ┊     ┃
    // J = ┃          ╲        ┊    ┃ = ┃                 ┊     ┃
    //     ┃            ╲      ┊    ┃   ┃                 ┊     ┃
    //     ┃              ╲    ┊    ┃   ┃                 ┊     ┃
    //     ┃                ╲  ┊    ┃   ┃                 ┊     ┃
    //     ┃┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┃   ┃┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┃
    //     ┃ ████              ┊ ██ ┃   ┃                 ┊     ┃
    //     ┃     ████          ┊ ██ ┃   ┃       J21       ┊ J22 ┃
    //     ┃          ...      ┊ ██ ┃   ┃                 ┊     ┃
    //     ┃              ████ ┊ ██ ┃   ┃                 ┊     ┃
    //     ┗                   ╵    ┛   ┗                 ╵     ┛
    //
    //  J  : float[MX*(NX+NY),MX*NX+NP];
    //  J11: float[MX*NX    ,MX*NX]; J11 is a Diagonal Matrix; Diag. represent the weights on Δ; Assumed to NOT be rank-deficient;
    //  J21: float[MX*NX    ,   NP]; J21 is a diagonal block matrix (possibly non-square).
    //  J22: float[MX*NY    ,   NP]; J22 is a dense matrix
    //
    // If there are no rank deficiencies in the upper left MX*NX columns, J can be sparsely QR decomposed as follows to solve ODR problem:
    // 
    //       ┏              ╷     ┓     ┏                 ╷     ┓
    //       ┃ ▜▓▓          ┊ ▓▓▓ ┃     ┃                 ┊     ┃
    //       ┃  ▜▓          ┊ ▓▓▓ ┃     ┃                 ┊     ┃
    //       ┃   ▜          ┊ ▓▓▓ ┃     ┃                 ┊     ┃
    //       ┃    ▜▓▓       ┊ ▓▓▓ ┃     ┃                 ┊     ┃
    //       ┃     ▜▓       ┊ ▓▓▓ ┃     ┃                 ┊     ┃
    //       ┃      ▜       ┊ ▓▓▓ ┃     ┃       R11       ┊ R21 ┃
    //       ┃       .      ┊ ▓▓▓ ┃     ┃                 ┊     ┃
    // J = Q·┃        .     ┊ ▓▓▓ ┃ = Q·┃                 ┊     ┃
    //       ┃         .    ┊ ▓▓▓ ┃     ┃                 ┊     ┃
    //       ┃          ▜▓▓ ┊ ▓▓▓ ┃     ┃                 ┊     ┃
    //       ┃           ▜▓ ┊ ▓▓▓ ┃     ┃                 ┊     ┃
    //       ┃            ▜ ┊ ▓▓▓ ┃     ┃                 ┊     ┃
    //       ┃┄┄┄┄┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┃     ┃┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┃
    //       ┃              ┊ ▜██ ┃     ┃                 ┊     ┃
    //       ┃       0      ┊  ▜█ ┃     ┃                 ┊     ┃
    //       ┃              ┊   ▜ ┃     ┃        0        ┊ R22 ┃
    //       ┗              ╵     ┛     ┗                 ╵     ┛
    //
    // R11: float[MX*NX               , MX*NX]; Block diagonal matrix where each block is of size [NX,NX] and upper-triangular. Off-diagonal entries can be computed "on-demand".
    // R21: float[MX*NX               ,    NP]; Dense matrix. Entries can be computed implicitly.
    // R22: float[min(MX*NY, MX*NX+NP),    NP]; Dense upper triangular matrix.
    //

  }, {
    key: "computeNewton",
    value: function computeNewton() {
      var MX = this.MX,
          NX = this.NX,
          NY = this.NY,
          NP = this.NP,
          QF = this.QF,
          J11 = this.J11,
          F1 = this.prepared_QF,
          R11 = this.newton_R11,
          R21 = this.newton_R21,
          R22 = this.newton_R22,
          P = this.newton_P,
          X = this.newton_dX,
          F0 = this.F0,
          D = this.D;
      if (this.rank >= 0) return;
      if (this.rank !== -1) throw new Error('Assertion failed.');
      this.prepare();
      var K = Math.min(MX * NY, MX * NX + NP); //
      // STEP 1: MEMORY INITIALIZATION
      //
      // for R11, only the diagonal is stored in memory explicitly. Off-diagonal entries are computed "on-demand"

      for (var i = MX * NX; i-- > 0;) {
        R11[i] = J11[i];
      }

      for (var _i34 = K; _i34-- > 0;) {
        QF[MX * NX + _i34] = F1[_i34];
      }

      for (var _i35 = MX * NX; _i35-- > 0;) {
        QF[_i35] = F0[_i35];
      }

      var rnk = this._qr_decomp(R11, R21, R22, P, QF);

      this.rank = rnk + MX * NX;

      for (var _i36 = MX * NX + rnk; _i36-- > 0;) {
        X[_i36] = -QF[_i36];
      }

      if (rnk !== NP) {
        //
        // STEP 2.4: ELIMINATE RANK-DEFICIENT COLUMS
        //           O( rnk² * (NP-rnk) ) operations
        X.fill(0.0, MX * NX + rnk, MX * NX + NP); // factor in scaling into R22

        for (var j = NP; j-- > 0;) {
          var d = D[MX * NX + P[j]];
          if (0 !== d) for (var _i37 = rnk; _i37-- > 0;) {
            R22[NP * _i37 + j] /= d;
          }
        } // eliminate lower part of linear dependent columns of R22


        for (var _i38 = rnk; _i38-- > 0;) {
          var ii = NP * _i38 + _i38;

          for (var _j24 = NP; _j24-- > rnk;) {
            var ij = NP * _i38 + _j24,
                R_ij = R22[ij];
            if (0 === R_ij) continue;
            var R_ii = R22[ii]; // compute Givens rot.

            var _giv_rot_qr12 = (0, _giv_rot._giv_rot_qr)(R_ii, R_ij),
                _giv_rot_qr13 = (0, _slicedToArray2["default"])(_giv_rot_qr12, 3),
                c = _giv_rot_qr13[0],
                s = _giv_rot_qr13[1],
                nrm = _giv_rot_qr13[2];

            if (s !== 0) {
              if (c < 0) {
                c *= -1;
                s *= -1;
                nrm *= -1;
              } // apply Givens rot.


              for (var k = _i38; k-- > 0;) {
                var R_ki = R22[NP * k + _i38],
                    R_kj = R22[NP * k + _j24];
                R22[NP * k + _i38] = R_kj * s + c * R_ki;
                R22[NP * k + _j24] = R_kj * c - s * R_ki;
              }

              R22[ii] = nrm;
            }

            R22[ij] = s;
          }
        }
      }

      this._qr_solve(R11, R21, R22, P, rnk, X);
    }
  }, {
    key: "computeNewtonRegularized",
    value: function computeNewtonRegularized(λ) {
      λ *= 1;
      if (!(λ >= 0)) throw new Error('Assertion failed.');
      var N = this.N,
          MX = this.MX,
          NX = this.NX,
          NY = this.NY,
          NP = this.NP,
          R11 = this.regularized_R11,
          J11 = this.J11,
          R21 = this.regularized_R21,
          R22 = this.regularized_R22,
          P = this.regularized_P,
          X = this.regularized_dX,
          D = this.D,
          QF = this.QF,
          F0 = this.F0,
          F1 = this.prepared_QF;
      var Y = QF; // <- alias to reuse memory

      var K = Math.min(MX * NY, MX * NX + NP);
      if (R22.length !== Math.max(Math.min(MX * NY, MX * NX + NP + 1), NP + 1) * NP) throw new Error('Assertion failed.');
      if (QF.length !== Math.max(Math.min(MX * NY, MX * NX + NP + 1), NP + 1) + MX * NX) throw new Error('Assertion failed.');

      if (0 === λ) {
        this.computeNewton();
        var _R = this.newton_R11,
            _R2 = this.newton_R21,
            _R3 = this.newton_R22,
            _P = this.newton_P,
            newton_dX = this.newton_dX,
            rank = this.rank;

        for (var i = N; i-- > 0;) {
          X[i] = newton_dX[i];
        }

        var _r3 = this.scaledNorm(X);

        if (0 === _r3) return [0, 0];

        var _rnk = rank - MX * NX;

        if (rank < N) {
          for (var _i39 = NP; _i39-- > 0;) {
            var j = MX * NX + _P[_i39];
            Y[MX * NX + _i39] = X[j] * D[j];
          }

          for (var _i40 = _rnk; _i40-- > 0;) {
            for (var _j25 = NP; _j25-- > _rnk;) {
              var s = _R3[NP * _i40 + _j25];
              if (s === 0) continue;
              var c = Math.sqrt(1 - s * s);
              (0, _giv_rot._giv_rot_rows)(Y, 1, MX * NX + _i40, MX * NX + _j25, c, s);
            }
          }
        } else {
          for (var _i41 = NP; _i41-- > 0;) {
            var _j26 = MX * NX + _P[_i41];

            Y[MX * NX + _i41] = X[_j26] * D[_j26] * D[_j26];
          }
        }

        for (var _i42 = MX * NX; _i42-- > 0;) {
          Y[_i42] = X[_i42] * D[_i42] * D[_i42];
        }

        this._rt_solve(_R, _R2, _R3, _P, _rnk, Y);

        var _dr = 0;

        for (var _i43 = rank; _i43-- > 0;) {
          var d = Y[_i43];
          _dr += d * d;
        }

        _dr /= -_r3;
        return [_r3, _dr];
      }

      this.prepare(); //
      // STEP 4: MEMORY INITIALIZATION
      //
      // for R11, only the diagonal is stored in memory explicitly. Off-diagonal entries are computed "on-demand"

      for (var _i44 = MX * NX; _i44-- > 0;) {
        R11[_i44] = J11[_i44];
      }

      for (var _i45 = K; _i45-- > 0;) {
        QF[MX * NX + _i45] = F1[_i45];
      }

      for (var _i46 = MX * NX; _i46-- > 0;) {
        QF[_i46] = F0[_i46];
      }

      var λSqrt = Math.sqrt(λ); // eliminate the upper part of regularization before _qr_decomp
      // O(MX*NX) operations

      for (var _i47 = MX * NX; _i47-- > 0;) {
        var Dλ = D[_i47] * λSqrt;
        if (!(Dλ > 0)) throw new Error('Assertion failed.');

        var _giv_rot_qr14 = (0, _giv_rot._giv_rot_qr)(R11[_i47], Dλ),
            _giv_rot_qr15 = (0, _slicedToArray2["default"])(_giv_rot_qr14, 3),
            _c4 = _giv_rot_qr15[0],
            _s5 = _giv_rot_qr15[1],
            nrm = _giv_rot_qr15[2];

        R11[_i47] = nrm;
        QF[_i47] *= _c4;
      }

      var rnk = this._qr_decomp(R11, R21, R22, P, QF),
          R22_end = NP * NP;

      R22.fill(0.0, K * NP, (NP + 1) * NP); // <- zero out entries from previous calls to computeNewtonRegularized()

      for (var _i48 = NP; _i48-- > 0;) {
        var _D = D[MX * NX + P[_i48]];
        if (0 === _D) _D = 1;else _D *= λSqrt;
        if (!(_D > 0)) throw new Error('Assertion failed.');

        if (rnk <= _i48) {
          // fill up the rank-deficient rows with regularization
          R22[NP * _i48 + _i48] = _D;
          QF[MX * NX + _i48] = 0;
        } else {
          // eliminate remaining regularization entries (using Givens QR)
          R22[R22_end + _i48] = _D;
          QF[N] = 0;

          for (var _j27 = _i48; _j27 < NP; _j27++) {
            var jj = NP * _j27 + _j27,
                ij = R22_end + _j27,
                R_jj = R22[jj],
                R_ij = R22[ij],
                _giv_rot_qr16 = (0, _giv_rot._giv_rot_qr)(R_jj, R_ij),
                _giv_rot_qr17 = (0, _slicedToArray2["default"])(_giv_rot_qr16, 3),
                _c5 = _giv_rot_qr17[0],
                _s6 = _giv_rot_qr17[1],
                _nrm4 = _giv_rot_qr17[2];

            R22[ij] = 0;
            if (_s6 === 0) continue;
            R22[jj] = _nrm4;
            (0, _giv_rot._giv_rot_rows)(R22, NP - _j27 - 1, jj + 1, ij + 1, _c5, _s6);
            (0, _giv_rot._giv_rot_rows)(QF, 1, MX * NX + _j27, N, _c5, _s6);
          }
        }
      }

      for (var _i49 = N; _i49-- > 0;) {
        X[_i49] = -QF[_i49];
      }

      this._qr_solve(R11, R21, R22, P, NP, X);

      var r = this.scaledNorm(X);
      if (0 === r) return [0, 0];

      for (var _i50 = NP; _i50-- > 0;) {
        var _j28 = MX * NX + P[_i50];

        Y[MX * NX + _i50] = X[_j28] * D[_j28] * D[_j28];
      }

      for (var _i51 = MX * NX; _i51-- > 0;) {
        Y[_i51] = X[_i51] * D[_i51] * D[_i51];
      }

      this._rt_solve(R11, R21, R22, P, NP, Y);

      var dr = 0;

      for (var _i52 = N; _i52-- > 0;) {
        var _d = Y[_i52];
        dr += _d * _d;
      }

      dr /= -r;
      return [r, dr];
    }
  }, {
    key: "__DEBUG_RVD",
    get: function get() // <- meant for debugging only!
    {
      var M = this.M,
          N = this.N,
          MX = this.MX,
          NX = this.NX,
          NY = this.NY,
          NP = this.NP,
          R11 = this.newton_R11,
          R21 = this.newton_R21,
          R22 = this.newton_R22,
          P = this.newton_P,
          J21 = this.prepared_J21,
          J22 = this.prepared_J22;
      if (!(0 <= this.rank)) throw new Error('Assertion failed.');
      var L = Math.min(NX, NY),
          rnk = this.rank - MX * NX,
          R = new Float64Array(M * N),
          V = new Float64Array(N * N),
          D = new Float64Array(1 * N);
      D.fill(1.0);

      for (var i = 0; i < MX * NX; i++) {
        V[N * i + i] = 1;
      }

      for (var _i53 = 0; _i53 < NP; _i53++) {
        V[N * (MX * NX + _i53) + (MX * NX + P[_i53])] = 1;
      } // diag(R11)


      for (var _i54 = 0; _i54 < MX * NX; _i54++) {
        R[N * _i54 + _i54] = R11[_i54];
      } // triu(R11,+1)


      for (var _i55 = 0; _i55 < MX; _i55++) {
        for (var j = 0; j < NX; j++) {
          for (var k = 1 + j; k < NX; k++) {
            for (var l = 0; l < L; l++) {
              R[N * (NX * _i55 + j) + (NX * _i55 + k)] += R21[NX * (L * _i55 + l) + j] * J21[NX * (L * _i55 + l) + k];
            }
          }
        }
      } // R12


      for (var _i56 = 0; _i56 < MX; _i56++) {
        for (var _j29 = 0; _j29 < NX; _j29++) {
          for (var _k12 = 0; _k12 < NP; _k12++) {
            for (var _l2 = 0; _l2 < L; _l2++) {
              R[N * (NX * _i56 + _j29) + (NX * MX + _k12)] += R21[NX * (L * _i56 + _l2) + _j29] * J22[NP * (L * _i56 + _l2) + P[_k12]];
            }
          }
        }
      } // R22


      for (var _i57 = 0; _i57 < rnk; _i57++) {
        for (var _j30 = 0; _j30 < rnk; _j30++) {
          R[N * (MX * NX + _i57) + (MX * NX + _j30)] = R22[NP * _i57 + _j30];
        }
      }

      if (rnk !== NP) {
        for (var _i58 = MX * NX + NP; _i58-- > MX * NX;) {
          D[_i58] = this.D[_i58] || 1;
        } // APPLY SCALING TO R12 


        for (var _j31 = NP; _j31-- > 0;) {
          var D_j = this.D[MX * NX + P[_j31]];
          if (0 !== D_j) for (var _i59 = MX * NX; _i59-- > 0;) {
            R[N * _i59 + (MX * NX + _j31)] /= D_j;
          }
        } // APPLY GIVENS ROTATIONS TO R12


        for (var _i60 = rnk; _i60-- > 0;) {
          for (var _j32 = NP; _j32-- > rnk;) {
            var s = R22[NP * _i60 + _j32];
            if (s === 0) continue;
            var c = Math.sqrt(1 - s * s);

            for (var _k13 = MX * NX; _k13-- > 0;) {
              var ki = N * _k13 + (MX * NX + _i60),
                  kj = N * _k13 + (MX * NX + _j32),
                  R_ki = R[ki],
                  R_kj = R[kj];
              R[ki] = R_kj * s + c * R_ki;
              R[kj] = R_kj * c - s * R_ki;
            }
          }
        }

        for (var _i61 = MX * NX; _i61-- > 0;) {
          for (var _j33 = MX * NX + NP; _j33-- > MX * NX + rnk;) {
            R[N * _i61 + _j33] = 0;
          }
        } // APPLY GIVENS ROTATIONS TO V


        for (var _i62 = rnk; _i62-- > 0;) {
          for (var _j34 = NP; _j34-- > rnk;) {
            var _s7 = R22[NP * _i62 + _j34];
            if (0 === _s7) continue;

            var _c6 = Math.sqrt(1 - _s7 * _s7);

            (0, _giv_rot._giv_rot_rows)(V, N, N * (MX * NX + _i62), N * (MX * NX + _j34), _c6, _s7);
          }
        }
      }

      return [new _nd_array.NDArray(Int32Array.of(M, N), R), new _nd_array.NDArray(Int32Array.of(N, N), V), new _nd_array.NDArray(Int32Array.of(1, N), D)];
    }
  }]);
  return TrustRegionSolverTLS;
}();

exports.TrustRegionSolverTLS = TrustRegionSolverTLS;