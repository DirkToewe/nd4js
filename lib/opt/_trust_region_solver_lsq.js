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
exports.TrustRegionSolverLSQ = void 0;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _classCallCheck2 = _interopRequireDefault(require("@babel/runtime/helpers/classCallCheck"));

var _createClass2 = _interopRequireDefault(require("@babel/runtime/helpers/createClass"));

var _nd_array = require("../nd_array");

var _giv_rot = require("../la/_giv_rot");

var _norm2 = require("../la/norm");

var _rrqr = require("../la/rrqr");

var _tri = require("../la/tri");

var _urv = require("../la/urv");

var _optimization_error = require("./optimization_error");

var REPORT_STATE_READY = 1,
    REPORT_STATE_NA = 2,
    REPORT_STATE_CONSIDER = 3;
/** 
 */

var TrustRegionSolverLSQ = /*#__PURE__*/function () {
  function TrustRegionSolverLSQ(fJ, x0) {
    (0, _classCallCheck2["default"])(this, TrustRegionSolverLSQ);
    // EVAL fJ(x0)
    if (!(fJ instanceof Function)) throw new Error('new TrustRegionSolverLSQ(fJ, x0): fJ must be Function.');
    this.fJ = fJ;
    this.report_x = x0 = (0, _nd_array.array)('float64', x0);
    if (x0.ndim !== 1) throw new Error('new TrustRegionSolverLSQ(fJ, x0): x0.ndim must be 1.');

    var _fJ = fJ(new _nd_array.NDArray(x0.shape, x0.data.slice())),
        _fJ2 = (0, _slicedToArray2["default"])(_fJ, 2),
        f = _fJ2[0],
        J = _fJ2[1];

    this.report_f = f = (0, _nd_array.array)('float64', f);
    this.report_J = J = (0, _nd_array.array)('float64', J);
    if (f.ndim !== 1) throw new Error('new TrustRegionSolverLSQ(fJ, x0): fJ(x) returned [f,J] with f.ndim !== 1.');
    if (J.ndim !== 2) throw new Error('new TrustRegionSolverLSQ(fJ, x0): fJ(x) returned [f,J] with J.ndim !== 2.');

    var _J$shape = (0, _slicedToArray2["default"])(J.shape, 2),
        M = _J$shape[0],
        N = _J$shape[1];

    if (f.shape[0] !== M) throw new Error('new TrustRegionSolverLSQ(fJ, x0): fJ(x) returned [f,J] where J.shape[0] did not match f.shape[0].');
    if (x0.shape[0] !== N) throw new Error('new TrustRegionSolverLSQ(fJ, x0): fJ(x) returned [f,J] where J.shape[1] did not match x0.shape[0].');
    this.loss = this.report_loss = NaN;
    this.report_loss_grad = null; // ALLOCATE MEMORY

    this.M = M;
    this.N = N;
    var L = Math.min(M, N),
        K = Math.max(M, N + 1);
    this.D = new Float64Array(N), // <- trust region scaling diagonal
    // STORED FUNCTION INPUT AND OUTPUT AT CURRENT TRUST REGION CENTER
    this.X0 = new Float64Array(N * 1); // <- center of the current trust radius (point where fJ was last evaluated)

    this.F0 = new Float64Array(M * 1);
    this.J0 = new Float64Array(M * N);
    this.G0 = new Float64Array(N); // USED BY computeNewton() and computeNewtonRegularized

    this.tmp = new Float64Array(2 * N); // <- used to compute the column norms of J for the diagonal scaling

    this.TMP = new Float64Array(K * N); // <- stores R of regularized RRQR
    // USED IN computeNewton()

    this.newton_P = new Int32Array(N); // <- column permutation by RRQR decomp.

    this.newton_Q = new Int32Array(N); // <- column permutation by URV  decomp.

    this.newton_R = new Float64Array(L * L);
    this.newton_V = new Float64Array(N * N); // <- stores complete orthogonal decomp V TODO: In unlikely case (M < N), economic SRRQR could be used instead of full

    this.newton_UF = new Float64Array(M * 1);
    this.newton_dX = new Float64Array(N * 1); // <- stores gauss-newton point relative to X0
    // USED IN computeRegularized(λ)

    this.regularized_R0 = new Float64Array(L * N); // <- stores R of RRQR (backed up from Newton)

    this.regularized_QF = new Float64Array(K * 1);
    this.regularized_dX = new Float64Array(N * 1);
    this.rank = -1;
    this._report_state = REPORT_STATE_CONSIDER;
    Object.seal(this); // INIT DATA

    this._considerMove_computeLoss();

    this.makeConsideredMove();
  }

  (0, _createClass2["default"])(TrustRegionSolverLSQ, [{
    key: "wiggle",
    value: function wiggle() {
      throw new _optimization_error.OptimizationNoProgressError('Too many unsuccessfull iterations.');
    }
  }, {
    key: "_considerMove_computeLoss",
    value: function _considerMove_computeLoss() {
      var M = this.M,
          N = this.N,
          F = this.report_f.data,
          J = this.report_J.data;
      var mse = 0.0;

      for (var i = M; i-- > 0;) {
        var s = F[i];
        mse += s * s / M;
      }

      var mse_grad = new Float64Array(N);

      for (var _i = M; _i-- > 0;) {
        for (var j = N; j-- > 0;) {
          mse_grad[j] += J[N * _i + j] * F[_i] / M * 2;
        }
      }

      this.report_loss = mse;
      this.report_loss_grad = new _nd_array.NDArray(Int32Array.of(N), mse_grad); // <- TODO: computation of loss_grad could probably be delayed until makeConsideredMove() is called
    }
  }, {
    key: "scaledNorm",
    value: function scaledNorm(X) {
      var N = this.N,
          D = this.D;
      if (X.length !== N) throw new Error('Assertion failed.');
      var norm = new _norm2.FrobeniusNorm();

      for (var i = N; i-- > 0;) {
        norm.include(D[i] * X[i]);
      }

      return norm.result;
    }
  }, {
    key: "cauchyTravel",
    value: function cauchyTravel() {
      var M = this.M,
          N = this.N,
          J = this.J0,
          G = this.G0; // polynomial along the gradient const1 + 2ax + bx² = const2 + (a/b + x)² * b

      var a = 0,
          b = 0;

      for (var i = N; i-- > 0;) {
        a += G[i] * G[i];
      }

      for (var _i2 = M; _i2-- > 0;) {
        var Jg = 0;

        for (var j = N; j-- > 0;) {
          Jg += J[N * _i2 + j] * G[j];
        }

        b += Jg * Jg;
      }

      return 0 === a ? 0 : 0 === b ? -Infinity : -a / b;
    }
  }, {
    key: "report",
    value: function report() {
      if (this._report_state !== REPORT_STATE_READY) throw new Error('TrustRegionSolverLSQ::report: can only be called once after each makeConsideredMove() but not directly after considerMove(dX).');
      this._report_state = REPORT_STATE_NA;
      var result = [this.report_x, this.report_loss, this.report_loss_grad, this.report_f, this.report_J];
      this.report_x = this.report_f = this.report_J = this.report_loss_grad = null;
      this.report_loss = NaN;
      return result;
    }
  }, {
    key: "considerMove",
    value: function considerMove(dX) {
      var M = this.M,
          N = this.N,
          X0 = this.X0,
          F0 = this.F0,
          J0 = this.J0;
      this._report_state = REPORT_STATE_CONSIDER;
      if (N !== dX.length) throw new Error('Assertion failed.');
      var X_shape = Int32Array.of(N),
          X = new Float64Array(N);

      for (var i = N; i-- > 0;) {
        X[i] = X0[i] + dX[i];
      }

      var _this$fJ = this.fJ(new _nd_array.NDArray(X_shape, X.slice())),
          _this$fJ2 = (0, _slicedToArray2["default"])(_this$fJ, 2),
          f = _this$fJ2[0],
          J = _this$fJ2[1];

      f = (0, _nd_array.array)('float64', f);
      J = (0, _nd_array.array)('float64', J);
      if (f.ndim !== 1) throw new Error('Assertion failed.');
      if (J.ndim !== 2) throw new Error('Assertion failed.');
      if (f.shape[0] !== M) throw new Error('Assertion failed.');
      if (J.shape[0] !== M) throw new Error('Assertion failed.');
      if (J.shape[1] !== N) throw new Error('Assertion failed.');
      this.report_x = new _nd_array.NDArray(X_shape, X);
      this.report_f = f;
      this.report_J = J;

      this._considerMove_computeLoss();

      var predict_loss = 0.0;

      for (var _i3 = M; _i3-- > 0;) {
        var s = 0;

        for (var j = N; j-- > 0;) {
          s += J0[N * _i3 + j] * (X[j] - X0[j]);
        }

        s += F0[_i3];
        predict_loss += s * s / M;
      }

      return [predict_loss, this.report_loss];
    }
  }, {
    key: "makeConsideredMove",
    value: function makeConsideredMove() {
      if (this._report_state !== REPORT_STATE_CONSIDER) throw new Error('Assertion failed.');
      this._report_state = REPORT_STATE_READY;
      this.rank = -1;
      this.loss = this.report_loss;
      var M = this.M,
          N = this.N,
          D = this.D,
          norm = this.tmp,
          X0 = this.X0,
          F0 = this.F0,
          J0 = this.J0,
          G0 = this.G0;
      ;
      {
        var f = this.report_f.data,
            J = this.report_J.data,
            x0 = this.report_x.data; // COPY (x,f,J)

        for (var i = N; i-- > 0;) {
          X0[i] = x0[i];
        }

        for (var _i4 = M; _i4-- > 0;) {
          F0[_i4] = f[_i4];
        }

        for (var _i5 = M * N; _i5-- > 0;) {
          J0[_i5] = J[_i5];
        }
      } // COMPUTE GRADIENT OF (HALF) SQUARED ERROR (MIND THE HALF :P)

      G0.fill(0.0);

      for (var _i6 = M; _i6-- > 0;) {
        for (var j = N; j-- > 0;) {
          G0[j] += J0[N * _i6 + j] * F0[_i6];
        }
      } // UPDATE SCALING
      // compute column norms of J


      if (norm.length !== 2 * N) throw new Error('Assertion failed.');
      norm.fill(0);

      for (var _i7 = 0; _i7 < M; _i7++) {
        (0, _rrqr._norm_update)(norm, J0, N * _i7, 0);
      } // use column norms of J to update D


      for (var _j = 0; _j < N; _j++) {
        D[_j] = Math.max(D[_j], (0, _rrqr._norm)(norm, _j));
      }
    }
  }, {
    key: "computeNewton",
    value: function computeNewton() {
      if (0 <= this.rank) return;
      var M = this.M,
          N = this.N,
          D = this.D,
          J0 = this.J0,
          F0 = this.F0,
          Y = this.tmp,
          T = this.TMP,
          P = this.newton_P,
          Q = this.newton_Q,
          V = this.newton_V,
          X = this.newton_dX,
          UF = this.newton_UF;
      var L = Math.min(M, N); // INIT

      for (var i = N; i-- > 0;) {
        P[i] = i;
      }

      for (var _i8 = M * N; _i8-- > 0;) {
        T[_i8] = J0[_i8];
      }

      for (var _i9 = M; _i9-- > 0;) {
        UF[_i9] = F0[_i9];
      } // SOLVE GLOBAL MINIMUM


      (0, _rrqr._rrqr_decomp_inplace)(M, N, 1, T, 0, UF, 0, P, 0, Y); // backup the rrqr decomposition in R (for finding regularized minima)

      ;
      {
        var R0 = this.regularized_R0;

        for (var _i10 = L * N; _i10-- > 0;) {
          R0[_i10] = T[_i10];
        }
      }
      var rank = this.rank = (0, _rrqr._rrqr_rank)(M, N, T, 0, Y);

      if (rank < N) {
        // RANK DEFICIENT CASE -> USE COMPLETE ORTHOGONAL DECOMPOSITION
        // factor in scaling
        for (var j = N; j-- > 0;) {
          var d = D[P[j]];
          if (0 !== d) for (var _i11 = rank; _i11-- > 0;) {
            T[N * _i11 + j] /= d;
          }
        } // complete orthogonal decompositon


        for (var _i12 = N; _i12-- > 0;) {
          Q[_i12] = P[_i12];
        }

        V.fill(0.0);
        (0, _urv._urv_decomp_full)(rank, N, T, 0, V, 0, Q, 0); // compute least squares solution

        for (var _i13 = rank; _i13-- > 0;) {
          Y[_i13] = -UF[_i13];
        }

        ;
        {
          var R = this.newton_R; // backup complete orthogonal T to R

          for (var _i14 = rank; _i14-- > 0;) {
            for (var _j2 = rank; _j2-- > 0;) {
              R[L * _i14 + _j2] = T[N * _i14 + _j2];
            }
          }

          (0, _tri._triu_solve)(rank, L, 1, R, 0, Y, 0); // Y = S \ Y
        }
        X.fill(0.0);

        for (var k = rank; k-- > 0;) {
          for (var _i15 = N; _i15-- > 0;) {
            X[_i15] += V[N * k + _i15] * Y[k];
          }
        } // X = V.T @ Y
        // factor out scaling


        for (var _i16 = N; _i16-- > 0;) {
          var _d = D[_i16];
          if (0 !== _d) X[_i16] /= _d;
        }
      } else {
        // FULL RANK CASE -> SCALING NOT RELEVANT TO GLOBAL SOLUTION
        for (var _i17 = N; _i17-- > 0;) {
          Y[_i17] = -UF[_i17];
        }

        (0, _tri._triu_solve)(N, N, 1, T, 0, Y, 0);

        for (var _i18 = N; _i18-- > 0;) {
          X[P[_i18]] = Y[_i18];
        }
      }
    } // Regularized least squares (RLS) solver. Solves the follwing equation system:
    //
    //   (JᵀJ + λDᵀD)x = Jᵀf
    //
    // The equation above can be solved as the following least squares problem:
    //
    //   min ║Ax - z║₂
    //    x
    //
    // where
    //       ┏              ┓        ┏   ┓
    //       ┃              ┃        ┃   ┃
    //       ┃              ┃        ┃   ┃
    //       ┃       J      ┃        ┃-f ┃
    //       ┃              ┃        ┃   ┃
    //       ┃              ┃        ┃   ┃
    //   A = ┃┄┄┄┄┄┄┄┄┄┄┄┄┄┄┃    z = ┃┄┄┄┃
    //       ┃√λ D[0]       ┃        ┃ 0 ┃
    //       ┃   ⋱         ┃        ┃ ⋮ ┃
    //       ┃    √λ D[N-1] ┃        ┃ 0 ┃
    //       ┗              ┛        ┗   ┛
    //
    // On top of that we can utilize the fact that J has already
    // been RRQR decomposed in computeNewton().

  }, {
    key: "computeNewtonRegularized",
    value: function computeNewtonRegularized(λ) {
      this.computeNewton();
      var M = this.M,
          N = this.N,
          D = this.D,
          rnk = this.rank,
          Y = this.tmp,
          T = this.TMP,
          P = this.newton_P,
          QF0 = this.newton_UF,
          newton_R = this.newton_R,
          V = this.newton_V,
          newton_dX = this.newton_dX,
          R0 = this.regularized_R0,
          X = this.regularized_dX,
          QF = this.regularized_QF;
      var L = Math.min(M, N);

      if (0 === λ) {
        for (var i = N; i-- > 0;) {
          X[i] = newton_dX[i];
        }

        if (rnk < N) {
          // UNDER-DETERMINED/RANK-DEFICIENT CASE -> USE COMPLETE ORTHOGONAL DECOMP. FROM computeNewton()
          // COMPUTE DISTANCE
          var _r = this.scaledNorm(X); // <- the scaled length


          if (0 === _r) return [0, 0]; // COMPUTE DISTANCE DERIVATIVE dr/dλ, see:
          //  "THE LEVENBERG-MARQUARDT ALGORITHM: IMPLEMENTATION AND THEORY"
          //   by Jorge J. Moré
          //   Chapter 5 "The Levenberg-Marquardt Parameter", Equation (5.8)

          Y.fill(0.0, 0, rnk);

          for (var _i19 = rnk; _i19-- > 0;) {
            for (var j = N; j-- > 0;) {
              Y[_i19] += V[N * _i19 + j] * D[j] * X[j];
            }
          }

          (0, _tri._triu_t_solve)(rnk, L, 1, newton_R, 0, Y, 0);
          var _dr = 0;

          for (var _i20 = rnk; _i20-- > 0;) {
            var d = Y[_i20];
            _dr += d * d;
          }

          _dr /= -_r;
          return [_r, _dr];
        }
      }

      if (!(0 <= λ)) throw new Error('TrustRegionSolverLSQ.prototype.computeMinRegularized(λ): λ must be non-negative.'); // START AT RRQR FROM computeNewton() (by copying R0 -> T)

      for (var _i21 = rnk * N; _i21-- > 0;) {
        T[_i21] = R0[_i21];
      } // TODO scale T to norm(T,'fro')=1 to avoid underflow ?


      if (0 !== λ) {
        for (var _i22 = rnk; _i22-- > 0;) {
          QF[_i22] = -QF0[_i22];
        }

        T.fill(0.0, rnk * N, (N + 1) * N);
        QF.fill(0.0, rnk, N + 1);
        var λSqrt = Math.sqrt(λ);

        for (var _j3 = N; _j3-- > 0;) {
          var Dλ = D[P[_j3]];
          if (0 === Dλ) Dλ = 1;else Dλ *= λSqrt;
          if (!(Dλ > 0)) throw new Error('TrustRegionSolverLSQ.prototype.computeMinRegularized(λ): λ too small (caused underflow).');
          if (rnk <= _j3) // fill up the rank-deficient rows with regularization
            T[N * _j3 + _j3] = Dλ;else {
            // eliminate remaining regularization entries (using Givens QR)
            T[N * N + _j3] = Dλ;
            QF[N] = 0;

            for (var _i23 = _j3; _i23 < N; _i23++) {
              var ii = N * _i23 + _i23,
                  Ni = N * N + _i23,
                  T_Ni = T[Ni];
              if (0 === T_Ni) continue;

              var T_ii = T[ii],
                  _giv_rot_qr2 = (0, _giv_rot._giv_rot_qr)(T_ii, T_Ni),
                  _giv_rot_qr3 = (0, _slicedToArray2["default"])(_giv_rot_qr2, 3),
                  c = _giv_rot_qr3[0],
                  s = _giv_rot_qr3[1],
                  norm = _giv_rot_qr3[2];

              T[Ni] = 0;
              if (0 === s) continue;
              T[ii] = norm;
              (0, _giv_rot._giv_rot_rows)(T, N - 1 - _i23, ii + 1, Ni + 1, c, s);
              (0, _giv_rot._giv_rot_rows)(QF, 1, _i23, N, c, s);
            }
          }
        }

        for (var _i24 = N; _i24-- > 0;) {
          Y[_i24] = QF[_i24];
        }

        (0, _tri._triu_solve)(N, N, 1, T, 0, Y, 0);

        for (var _i25 = N; _i25-- > 0;) {
          X[P[_i25]] = Y[_i25];
        }
      } else if (rnk !== N) throw new Error('Assertion failed.'); // COMPUTE DISTANCE


      var r = this.scaledNorm(X); // <- the scaled length

      if (0 === r) return [0, 0]; // COMPUTE DISTANCE DERIVATIVE dr/dλ, see:
      //  "THE LEVENBERG-MARQUARDT ALGORITHM: IMPLEMENTATION AND THEORY"
      //   by Jorge J. Moré
      //   Chapter 5 "The Levenberg-Marquardt Parameter", Equation (5.8)

      for (var _i26 = N; _i26-- > 0;) {
        var _j4 = P[_i26];
        Y[_i26] = X[_j4] * D[_j4] * D[_j4];
      }

      (0, _tri._triu_t_solve)(N, N, 1, T, 0, Y, 0);
      var dr = 0;

      for (var _i27 = N; _i27-- > 0;) {
        var _d2 = Y[_i27];
        dr += _d2 * _d2;
      }

      dr /= -r;
      return [r, dr];
    }
  }]);
  return TrustRegionSolverLSQ;
}();

exports.TrustRegionSolverLSQ = TrustRegionSolverLSQ;