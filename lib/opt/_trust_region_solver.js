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

var TrustRegionSolverLSQ = /*#__PURE__*/function () {
  function TrustRegionSolverLSQ(M, N) {
    (0, _classCallCheck2["default"])(this, TrustRegionSolverLSQ);
    if (M % 1 !== 0) throw new Error('Assertion failed.');
    if (N % 1 !== 0) throw new Error('Assertion failed.');
    if (!(0 < M)) throw new Error('Assertion failed.');
    if (!(0 < N)) throw new Error('Assertion failed.');
    this.M = M | 0;
    this.N = N | 0;
    var L = Math.min(M, N),
        K = Math.max(L + N, M);
    this.P = new Int32Array(N); // <- column permutation by RRQR

    this.Q = new Int32Array(N); // <- column permutation by RRQR

    this.D = new Float64Array(N), // <- trust region scaling diagonal
    this.R = new Float64Array(L * N); // <- stores result of RRQR

    this.S = new Float64Array(L * L); // <- stores complete orthogonal decomp R (in rank-deficient case)

    this.T = new Float64Array(K * N); // <- working memory for decompositions (RRQR/URV)

    this.V = new Float64Array(N * N); // <- stores complete orthogonal decomp V TODO: In unlikely case (M < N), economic SRRQR could be used instead of full

    this.F = new Float64Array(M * 1);
    this.Y = new Float64Array(K * 1); // <- (L+N)*1 should be enough ... right?

    this.X = new Float64Array(N * 1);
    this.norm = new Float64Array(2 * N); // <- couldn't Y be used instead?

    this.G = new Float64Array(N);
    this.rank = -1;
    this._state = 0;
    Object.seal(this);
  }

  (0, _createClass2["default"])(TrustRegionSolverLSQ, [{
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
    key: "cauchyPointTravel",
    value: function cauchyPointTravel() {
      if (this._state !== 1) throw new Error('TrustRegionSolverLSQ.prototype.cauchyPointTravel(): can only be called directly after update(f,J) call.');
      var M = this.M,
          N = this.N,
          X = this.X,
          F = this.F,
          J = this.T,
          G = this.G; // polynomial along the gradient const1 + 2ax + bx² = const2 + (a/b + x)² * b

      var a = 0,
          b = 0;

      for (var i = N; i-- > 0;) {
        a += G[i] * G[i];
      }

      for (var _i = M; _i-- > 0;) {
        var Jg = 0;

        for (var j = N; j-- > 0;) {
          Jg += J[N * _i + j] * G[j];
        }

        b += Jg * Jg;
      }

      return 0 === a ? 0 : 0 === b ? -Infinity : -a / b;
    } // METHODS NEED TO BE CALLED IN THE FOLLOWING ORDER:
    //   1) update(f,J)
    //   2) computeMinGlobal()
    //   3) computeMinRegularized() (can be called multiple times)
    // computeMinGlobal() and computeMinRegularized() store their results
    // in X overwriting the previous result.

  }, {
    key: "update",
    value: function update(f, J) {
      this._state = 1;
      var M = this.M,
          N = this.N,
          D = this.D,
          P = this.P,
          Q = this.Q,
          T = this.T,
          F = this.F,
          G = this.G,
          norm = this.norm; // ASSERTIONS

      if (!(f instanceof _nd_array.NDArray)) throw new Error('Assertion failed.');
      if (!(J instanceof _nd_array.NDArray)) throw new Error('Assertion failed.');
      if (f.ndim !== 1) throw new Error('Assertion failed.');
      if (J.ndim !== 2) throw new Error('Assertion failed.');
      if (f.shape[0] !== M) throw new Error('Assertion failed.');
      if (J.shape[0] !== M) throw new Error('Assertion failed.');
      if (J.shape[1] !== N) throw new Error('Assertion failed.');
      f = f.data;
      J = J.data; // COPY f -> F

      for (var i = M * 1; i-- > 0;) {
        F[i] = f[i];
      } // COPY J -> T


      for (var _i2 = M * N; _i2-- > 0;) {
        T[_i2] = J[_i2];
      } // INIT P


      for (var _i3 = N; _i3-- > 0;) {
        P[_i3] = _i3;
      } // UPDATE SCALING
      // compute column norms of J


      norm.fill(0);

      for (var _i4 = 0; _i4 < M; _i4++) {
        (0, _rrqr._norm_update)(norm, J, N * _i4, 0);
      } // use column norms of J to update D


      for (var j = 0; j < N; j++) {
        D[j] = Math.max(D[j], (0, _rrqr._norm)(norm, j));
      } // COMPUTE GRADIENT OF THE (HALF) SQUARED ERROR (MIND THE HALF :P)


      G.fill(0.0);

      for (var _i5 = M; _i5-- > 0;) {
        for (var _j = N; _j-- > 0;) {
          G[_j] += J[N * _i5 + _j] * F[_i5];
        }
      }
    }
  }, {
    key: "computeMinGlobal",
    value: function computeMinGlobal() {
      if (this._state !== 1) throw new Error('TrustRegionSolverLSQ.prototype.computeMinGlobal(): can only be called once after every update(f,J) call.');
      this._state = 2;
      var M = this.M,
          N = this.N,
          P = this.P,
          Q = this.Q,
          D = this.D,
          R = this.R,
          S = this.S,
          T = this.T,
          V = this.V,
          X = this.X,
          F = this.F,
          Y = this.Y,
          norm = this.norm;
      var L = Math.min(M, N); // SOLVE GLOBAL MINIMUM

      (0, _rrqr._rrqr_decomp_inplace)(M, N, 1, T, 0, F, 0, P, 0, norm); // backup the rrqr decomposition in R (for finding regularized minima)

      for (var i = L * N; i-- > 0;) {
        R[i] = T[i];
      }

      var rank = this.rank = (0, _rrqr._rrqr_rank)(M, N, T, 0,
      /*tmp=*/
      norm);

      if (rank < N) {
        // RANK DEFICIENT CASE -> USE COMPLETE ORTHOGONAL DECOMPOSITION
        // factor in scaling
        for (var _i6 = L; _i6-- > 0;) {
          for (var j = N; j-- > 0;) {
            var d = D[P[j]];
            if (0 !== d) T[N * _i6 + j] /= d;
          }
        } // complete orthogonal decompositon


        for (var _i7 = N; _i7-- > 0;) {
          Q[_i7] = P[_i7];
        }

        V.fill(0.0);
        (0, _urv._urv_decomp_full)(rank, N, T, 0, V, 0, Q, 0); // backup complete orthogonal T to S

        for (var _i8 = rank; _i8-- > 0;) {
          for (var _j2 = rank; _j2-- > 0;) {
            S[L * _i8 + _j2] = T[N * _i8 + _j2];
          }
        } // compute least squares solution


        for (var _i9 = rank; _i9-- > 0;) {
          Y[_i9] = -F[_i9];
        }

        (0, _tri._triu_solve)(rank, L, 1, S, 0, Y, 0); // Y = S \ Y

        X.fill(0.0);

        for (var k = rank; k-- > 0;) {
          for (var _i10 = N; _i10-- > 0;) {
            X[_i10] += V[N * k + _i10] * Y[k];
          }
        } // X = V.T @ Y
        // factor out scaling


        for (var _i11 = N; _i11-- > 0;) {
          var _d = D[_i11];
          if (0 !== _d) X[_i11] /= _d;
        }
      } else {
        // FULL RANK CASE -> SCALING NOT RELEVANT TO GLOBAL SOLUTION
        for (var _i12 = N; _i12-- > 0;) {
          Y[_i12] = -F[_i12];
        }

        (0, _tri._triu_solve)(N, N, 1, T, 0, Y, 0);

        for (var _i13 = N; _i13-- > 0;) {
          X[P[_i13]] = Y[_i13];
        }
      }
    } // Regularized least squares (RLS) solver. Solves the follwing equation system:
    //
    //   (JᵀJ + λI)x = Jᵀf
    //
    // The equation above can be solved as the following least squares problem:
    //
    //   min ║Ax - z║₂
    //    x
    //
    // where
    //       ┏       ┓        ┏   ┓
    //       ┃       ┃        ┃   ┃
    //       ┃       ┃        ┃   ┃
    //       ┃   J   ┃        ┃-f ┃
    //       ┃       ┃        ┃   ┃
    //       ┃       ┃        ┃   ┃
    //   A = ┃┄┄┄┄┄┄┄┃    z = ┃┄┄┄┃
    //       ┃√λ     ┃        ┃ 0 ┃
    //       ┃   ⋱   ┃        ┃ ⋮ ┃
    //       ┃    √λ ┃        ┃ 0 ┃
    //       ┗       ┛        ┗   ┛
    //
    // On top of that we can utilize the fact that J has already
    // been RRQR decomposed in computeMinGlobal().

  }, {
    key: "computeMinRegularized",
    value: function computeMinRegularized(λ) {
      if (this._state !== 2 && this._state !== 3) throw new Error('TrustRegionSolverLSQ.prototype.computeMinRegularized(λ): can only be (repeatedly) called after computeMinGlobal().');
      if (0 !== λ) this._state = 3;
      var M = this.M,
          N = this.N,
          P = this.P,
          D = this.D,
          R = this.R,
          S = this.S,
          T = this.T,
          V = this.V,
          X = this.X,
          F = this.F,
          Y = this.Y,
          rnk = this.rank;
      var L = Math.min(M, N);

      if (0 === λ && rnk < N) {
        // UNDER-DETERMINED/RANK-DEFICIENT CASE -> USE COMPLETE ORTHOGONAL DECOMPOSTION
        if (this._state !== 2) {
          this._state = 2; // GLOBAL MIN. NEEDS TO BE RECOMPUTED
          // compute least squares solution

          for (var i = rnk; i-- > 0;) {
            Y[i] = -F[i];
          }

          (0, _tri._triu_solve)(rnk, L, 1, S, 0, Y, 0); // Y = S \ Y

          X.fill(0.0);

          for (var k = rnk; k-- > 0;) {
            for (var _i14 = N; _i14-- > 0;) {
              X[_i14] += V[N * k + _i14] * Y[k];
            }
          } // X = V.T @ Y
          // factor out scaling


          for (var _i15 = N; _i15-- > 0;) {
            var d = D[_i15];
            if (0 !== d) X[_i15] /= d;
          }
        } // COMPUTE DISTANCE


        var _r = this.scaledNorm(X); // <- the scaled length


        if (0 === _r) return [0, 0]; // COMPUTE DISTANCE DERIVATIVE dr/dλ, see:
        //  "THE LEVENBERG-MARQUARDT ALGORITHM: IMPLEMENTATION AND THEORY"
        //   by Jorge J. Moré
        //   Chapter 5 "The Levenberg-Marquardt Parameter", Equation (5.8)

        Y.fill(0.0, 0, rnk);

        for (var _i16 = rnk; _i16-- > 0;) {
          for (var j = N; j-- > 0;) {
            Y[_i16] += X[j] * D[j] * V[N * _i16 + j];
          }
        } // FORWARD SUBSTITUTION


        for (var _i17 = 0; _i17 < rnk; _i17++) {
          Y[_i17] /= S[L * _i17 + _i17];

          for (var _j3 = _i17; ++_j3 < rnk;) {
            Y[_j3] -= S[L * _i17 + _j3] * Y[_i17];
          }
        }

        var _dr = 0;

        for (var _i18 = rnk; _i18-- > 0;) {
          var _d2 = Y[_i18];
          _dr += _d2 * _d2;
        }

        _dr /= -_r;
        return [_r, _dr];
      }

      if (!(0 <= λ)) throw new Error('TrustRegionSolverLSQ.prototype.computeMinRegularized(λ): λ must be positive.');

      for (var _i19 = rnk; _i19-- > 0;) {
        Y[_i19] = -F[_i19];
      }

      Y.fill(0.0, rnk, rnk + N); // COPY R->T TO REUSE THE RRQR OF J

      for (var _i20 = rnk * N; _i20-- > 0;) {
        T[_i20] = R[_i20];
      }

      T.fill(0.0, rnk * N, (rnk + N) * N); // TODO scale T to norm(T,'fro')=1 to avoid underflow ?

      if (0 !== λ) {
        var λSqrt = Math.sqrt(λ);

        for (var _i21 = N; _i21-- > 0;) {
          var Dλ = D[P[_i21]];
          if (Dλ === 0) Dλ = 1;else Dλ *= λSqrt;
          if (!(0 < Dλ)) throw new Error('TrustRegionSolverLSQ.prototype.computeMinRegularized(λ): λ too small (caused underflow).');

          var _j4 = (_i21 - rnk + N) % N;

          T[rnk * N + N * _j4 + _i21] = Dλ;
        } // COMPLETE QR DECOMPOSITION


        for (var _i22 = 0; _i22 < N; _i22++) {
          var J = Math.min(_i22 + 1, rnk) + N;

          for (var _j5 = N; _j5 < J; _j5++) {
            // USE GIVENS ROTATION TO ELIMINATE ELEMENT R_ji
            var ji = N * _j5 + _i22,
                T_ji = T[ji];
            if (0 === T_ji) continue;

            var ii = N * _i22 + _i22,
                T_ii = T[ii],
                _giv_rot_qr2 = (0, _giv_rot._giv_rot_qr)(T_ii, T_ji),
                _giv_rot_qr3 = (0, _slicedToArray2["default"])(_giv_rot_qr2, 3),
                c = _giv_rot_qr3[0],
                s = _giv_rot_qr3[1],
                norm = _giv_rot_qr3[2];

            T[ji] = 0;
            if (0 === s) continue;
            T[ii] = norm;
            (0, _giv_rot._giv_rot_rows)(T, N - 1 - _i22, ii + 1, ji + 1, c, s);
            (0, _giv_rot._giv_rot_rows)(Y, 1, _i22, _j5, c, s);
          }
        }
      } else if (rnk !== N) throw new Error('Assertion failed.');

      if (0 !== λ || this._state !== 2) // <- AVOID RECOMPUTING GLOBAL MIN.
        {
          if (0 === λ) this._state = 2;
          (0, _tri._triu_solve)(N, N, 1, T, 0, Y, 0);

          for (var _i23 = N; _i23-- > 0;) {
            X[P[_i23]] = Y[_i23];
          }
        } // COMPUTE DISTANCE


      var r = this.scaledNorm(X); // <- the scaled length

      if (0 === r) return [0, 0]; // COMPUTE DISTANCE DERIVATIVE dr/dλ, see:
      //  "THE LEVENBERG-MARQUARDT ALGORITHM: IMPLEMENTATION AND THEORY"
      //   by Jorge J. Moré
      //   Chapter 5 "The Levenberg-Marquardt Parameter", Equation (5.8)

      for (var _i24 = N; _i24-- > 0;) {
        var _j6 = P[_i24];
        Y[_i24] = X[_j6] * D[_j6] * D[_j6];
      } // FORWARD SUBSTITUTION


      for (var _i25 = 0; _i25 < N; _i25++) {
        Y[_i25] /= T[N * _i25 + _i25];

        for (var _j7 = _i25; ++_j7 < N;) {
          Y[_j7] -= T[N * _i25 + _j7] * Y[_i25];
        }
      }

      var dr = 0;

      for (var _i26 = N; _i26-- > 0;) {
        var _d3 = Y[_i26];
        dr += _d3 * _d3;
      }

      dr /= -r;
      return [r, dr];
    }
  }]);
  return TrustRegionSolverLSQ;
}();

exports.TrustRegionSolverLSQ = TrustRegionSolverLSQ;