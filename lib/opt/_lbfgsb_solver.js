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
exports.LBFGSB_Solver = void 0;

var _classCallCheck2 = _interopRequireDefault(require("@babel/runtime/helpers/classCallCheck"));

var _createClass2 = _interopRequireDefault(require("@babel/runtime/helpers/createClass"));

var _heap_sort_gen = require("../arrays/heap_sort_gen");

var _float64_utils = require("../dt/float64_utils");

var _cholesky = require("../la/cholesky");

var _ldl = require("../la/ldl");

var _pldlp = require("../la/pldlp");

var _tri = require("../la/tri");

function _createForOfIteratorHelper(o, allowArrayLike) { var it; if (typeof Symbol === "undefined" || o[Symbol.iterator] == null) { if (Array.isArray(o) || (it = _unsupportedIterableToArray(o)) || allowArrayLike && o && typeof o.length === "number") { if (it) o = it; var i = 0; var F = function F() {}; return { s: F, n: function n() { if (i >= o.length) return { done: true }; return { done: false, value: o[i++] }; }, e: function e(_e) { throw _e; }, f: F }; } throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method."); } var normalCompletion = true, didErr = false, err; return { s: function s() { it = o[Symbol.iterator](); }, n: function n() { var step = it.next(); normalCompletion = step.done; return step; }, e: function e(_e2) { didErr = true; err = _e2; }, f: function f() { try { if (!normalCompletion && it["return"] != null) it["return"](); } finally { if (didErr) throw err; } } }; }

function _unsupportedIterableToArray(o, minLen) { if (!o) return; if (typeof o === "string") return _arrayLikeToArray(o, minLen); var n = Object.prototype.toString.call(o).slice(8, -1); if (n === "Object" && o.constructor) n = o.constructor.name; if (n === "Map" || n === "Set") return Array.from(o); if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n)) return _arrayLikeToArray(o, minLen); }

function _arrayLikeToArray(arr, len) { if (len == null || len > arr.length) len = arr.length; for (var i = 0, arr2 = new Array(len); i < len; i++) { arr2[i] = arr[i]; } return arr2; }

var USE_LDL = false; // References
// ----------
// .. [1] "REPRESENTATIONS OF QUASI-NEWTON MATRICES AND THEIR USE IN LIMITED MEMORY METHODS"
//         by Richard H. Byrd, Jorge Nocedal and Robert B. Schnabel
//         Technical Report NAM-03, June 1992; revised January 21, 1996
//         http://users.iems.northwestern.edu/~nocedal/PDFfiles/representations.pdf
//
// .. [2] "A LIMITED MEMORY ALGORITHM FOR BOUND CONSTRAINED OPTIMIZATION"
//         by Richard H. Byrd, Peihuang Lu, Jorge Nocedal and Ciyou Zhu
//         Technical Report NAM-08; Revised: May 1994
//         http://users.iems.northwestern.edu/~nocedal/PDFfiles/limited.pdf
// Compares to [2], the following symbols have been renamed in the source:
//
//   Y -> dG
//   S -> dX
//
// dG is the delta in gradients between iterations. dX is the delta in X.
// The names dG and dX therefore seemed to be more intuitive.

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

var LBFGSB_Solver = /*#__PURE__*/function () {
  function LBFGSB_Solver(M, N) {
    (0, _classCallCheck2["default"])(this, LBFGSB_Solver);
    if (M % 1 !== 0) throw new Error('Assertion failed.');
    if (N % 1 !== 0) throw new Error('Assertion failed.');
    if (!(0 < M)) throw new Error('Assertion failed.');
    if (!(0 < N)) throw new Error('Assertion failed.');
    this.m = 0;
    this.M = M | 0;
    this.N = N | 0;
    this._scale = 1;
    this.dX = new Float64Array(M * N);
    this.dG = new Float64Array(M * N);
    this.dXdX = new Float64Array(M * M);
    this.dXdG = new Float64Array(M * M);
    this.dGdG = new Float64Array(M * M);
    this.D = new Float64Array(M);
    this._J = new Float64Array(M * M);
    this._J_invalid = true;
    this.travels = new Float64Array(N);
    this.Bdx = new Float64Array(2 * M);
    this.Bg = new Float64Array(2 * M);
    this.Bei = new Float64Array(2 * M);
    this.M_WAAW = new Float64Array(4 * M * M);
    if (!USE_LDL) this.P = new Int32Array(2 * M);
  }

  (0, _createClass2["default"])(LBFGSB_Solver, [{
    key: "forget",
    value: function forget(k) {
      // TODO:
      //   Moving all rows/columns up by just to forget old entries not efficient.
      //   A better way would be to use the matrix rows/columns like a circular buffer
      //   restarting in the rop row each time the bottom is reached. That will
      //   however make the implementation of all the member methods quite a bit
      //   more complicated.
      var m = this.m,
          M = this.M,
          N = this.N,
          dX = this.dX,
          dXdX = this.dXdX,
          dG = this.dG,
          dXdG = this.dXdG,
          dGdG = this.dGdG,
          D = this.D;
      if (m > M) throw new Error('Assertion failed.');
      if (k % 1 !== 0) throw new Error('Assertion failed.');
      if (!(k > 0)) throw new Error('Assertion failed.');
      if (!(k <= m)) throw new Error('Assertion failed.'); // MAKE SPACE IN LAST ROW BY OVERWRITING 1ST ROW

      for (var _i = 0, _arr = [dX, dG]; _i < _arr.length; _i++) {
        var A = _arr[_i];

        for (var i = k; i < m; i++) {
          for (var j = 0; j < N; j++) {
            A[N * (i - k) + j] = A[N * i + j];
          }
        }
      }

      for (var _i2 = 0, _arr2 = [dXdX, dXdG, dGdG]; _i2 < _arr2.length; _i2++) {
        var _A = _arr2[_i2];

        for (var _i3 = k; _i3 < m; _i3++) {
          for (var _j = k; _j < m; _j++) {
            _A[M * (_i3 - k) + (_j - k)] = _A[M * _i3 + _j];
          }
        }
      }

      for (var _i4 = k; _i4 < m; _i4++) {
        D[_i4 - k] = D[_i4];
      }

      this.m = m - k;
      this._J_invalid = true;
    }
  }, {
    key: "update",
    value: function update(dx, dg) {
      var M = this.M,
          N = this.N,
          dX = this.dX,
          dXdX = this.dXdX,
          dG = this.dG,
          dXdG = this.dXdG,
          dGdG = this.dGdG,
          D = this.D;
      if (N !== dx.length) throw new Error('Assertion failed.');
      if (N !== dg.length) throw new Error('Assertion failed.');
      if (this.m >= M) this.forget(1);
      this.m++;
      var m = this.m,
          i = m - 1;
      if (m > M) throw new Error('Assertion failed.');

      for (var j = 0; j < N; j++) {
        dX[N * i + j] = dx[j];
      }

      for (var _j2 = 0; _j2 < N; _j2++) {
        dG[N * i + _j2] = dg[_j2];
      } // UPDATE dX @ dX.T


      for (var _j3 = 0; _j3 < m; _j3++) {
        var dxdx = 0;

        for (var k = 0; k < N; k++) {
          dxdx += dX[N * i + k] * dX[N * _j3 + k];
        }

        dXdX[M * i + _j3] = dXdX[M * _j3 + i] = dxdx;
      } // UPDATE dX @ dG.T


      for (var _j4 = 0; _j4 < m; _j4++) {
        dXdG[M * i + _j4] = 0;

        for (var _k = 0; _k < N; _k++) {
          dXdG[M * i + _j4] += dX[N * i + _k] * dG[N * _j4 + _k];
        }

        dXdG[M * _j4 + i] = 0;

        for (var _k2 = 0; _k2 < N; _k2++) {
          dXdG[M * _j4 + i] += dX[N * _j4 + _k2] * dG[N * i + _k2];
        }
      } // UPDATE sqrt( diag(dX @ dG.T) )


      D[i] = Math.sqrt(dXdG[M * i + i]); // UPDATE dG @ dG.T

      for (var _j5 = 0; _j5 < m; _j5++) {
        var dgdg = 0;

        for (var _k3 = 0; _k3 < N; _k3++) {
          dgdg += dG[N * i + _k3] * dG[N * _j5 + _k3];
        }

        dGdG[M * i + _j5] = dGdG[M * _j5 + i] = dgdg;
      } // INVALIDATE J


      this._J_invalid = true; // <- set J to be recomputed on demand
    }
  }, {
    key: "_compute_Bv",
    value: function _compute_Bv(bv, Bv) {
      var m = this.m,
          M = this.M,
          N = this.N,
          dG = this.dG,
          dXdG = this.dXdG,
          dX = this.dX,
          J = this.J,
          D = this.D,
          scale = this.scale;
      if (bv.length !== M * 2) throw new Error('Assertion failed.');
      if (Bv.length !== N) throw new Error('Assertion failed.');

      for (var i = 2 * m; i-- > m;) {
        bv[i] *= -1;
      }

      (0, _tri._tril_t_solve)(m, M, 1, J, 0, bv, m);

      for (var _i5 = 0; _i5 < m; _i5++) {
        bv[_i5] /= D[_i5];
      }

      for (var j = 1; j < m; j++) {
        for (var _i6 = 0; _i6 < j; _i6++) {
          bv[_i6] += dXdG[M * j + _i6] * bv[m + j] / dXdG[M * _i6 + _i6];
        }
      }

      for (var _i7 = N; _i7-- > 0;) {
        for (var _j6 = m; _j6-- > 0;) {
          Bv[_i7] += dX[N * _j6 + _i7] * bv[m + _j6] * scale;
        }
      }

      for (var _i8 = N; _i8-- > 0;) {
        for (var _j7 = m; _j7-- > 0;) {
          Bv[_i8] += dG[N * _j7 + _i8] * bv[_j7];
        }
      }
    }
  }, {
    key: "compute_Bv",
    value: function compute_Bv(v, Bv) {
      var N = this.N,
          scale = this.scale,
          bv = this.Bg;
      if (v.length !== N) throw new Error('Assertion failed.');
      this.compute_bv(v, bv);

      for (var i = N; i-- > 0;) {
        Bv[i] = v[i] * scale;
      }

      this._compute_Bv(bv, Bv);
    }
    /** This method "sort of" computes half an inner product over B.
     *  To compute an inner product u.T @ B @ v, You can use:
     *
     *    compute_bv(u, ub);
     *    compute_bv(v, vb);
     *    compute_ubbv(ub, dot(u,v), bv)
     */

  }, {
    key: "compute_bv",
    value: function compute_bv(v, bv) {
      //               ┏                 ┓
      //               ┃        ┊        ┃
      //               ┃        ┊        ┃ ┏                 ┓ -1  ┏                 ┓   ┏                 ┓ -1  ┏                         ┓
      //               ┃        ┊        ┃ ┃        ┊        ┃     ┃        ┊        ┃   ┃        ┊        ┃     ┃                         ┃
      //               ┃        ┊        ┃ ┃   √D   ┊-√D \ Lᵀ┃     ┃   -I   ┊    0   ┃   ┃   √D   ┊   0    ┃     ┃           dG            ┃
      //               ┃        ┊        ┃ ┃        ┊        ┃     ┃        ┊        ┃   ┃        ┊        ┃     ┃                         ┃
      // uᵀ⋅WMW⋅v = uᵀ⋅┃   dGᵀ  ┊ B0⋅dXᵀ ┃⋅┃┈┈┈┈┈┈┈┈┼┈┈┈┈┈┈┈┈┃  ⋅  ┃┈┈┈┈┈┈┈┈┼┈┈┈┈┈┈┈┈┃ ⋅ ┃┈┈┈┈┈┈┈┈┼┈┈┈┈┈┈┈┈┃  ⋅  ┃┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┃⋅v
      //               ┃        ┊        ┃ ┃        ┊        ┃     ┃        ┊        ┃   ┃        ┊        ┃     ┃                         ┃
      //               ┃        ┊        ┃ ┃    0   ┊   Jᵀ   ┃     ┃    0   ┊    I   ┃   ┃-L / √D ┊   J    ┃     ┃          dX⋅B0          ┃
      //               ┃        ┊        ┃ ┃        ┊        ┃     ┃        ┊        ┃   ┃        ┊        ┃     ┃                         ┃
      //               ┃        ┊        ┃ ┗                 ┛     ┗                 ┛   ┗                 ┛     ┗                         ┛
      //               ┃        ┊        ┃
      //               ┗                 ┛
      //                 ┏                 ┓                  ┏                 ┓ -1  ┏                      ┓
      //                 ┃        ┊        ┃                  ┃        ┊        ┃     ┃                      ┃
      //                 ┃   -I   ┊    0   ┃                  ┃   √D   ┊   0    ┃     ┃         dG           ┃
      //                 ┃        ┊        ┃                  ┃        ┊        ┃     ┃                      ┃
      //        =: uᵀ⋅bᵀ⋅┃┈┈┈┈┈┈┈┈┼┈┈┈┈┈┈┈┈┃⋅b⋅v    =>    b = ┃┈┈┈┈┈┈┈┈┼┈┈┈┈┈┈┈┈┃  ⋅  ┃┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┃
      //                 ┃        ┊        ┃                  ┃        ┊        ┃     ┃                      ┃
      //                 ┃    0   ┊    I   ┃                  ┃-L / √D ┊   J    ┃     ┃        dX⋅B0         ┃
      //                 ┃        ┊        ┃                  ┃        ┊        ┃     ┃                      ┃
      //                 ┗                 ┛                  ┗                 ┛     ┗                      ┛
      //
      // 
      // B = scale⋅I - WMW =>  scale⋅uᵀ⋅v - uᵀ⋅WMW⋅v
      //
      // Where:
      //   L = tril(dX⋅dGᵀ, -1)
      //   D = diag(dX⋅dGᵀ)
      //   I = eye(m)
      //   J = cholesky( dX⋅B0⋅dXᵀ + (L/D)⋅Lᵀ )
      //  B0 = I⋅scale
      var m = this.m,
          M = this.M,
          N = this.N,
          dG = this.dG,
          dXdG = this.dXdG,
          dX = this.dX,
          J = this.J,
          D = this.D,
          scale = this.scale;
      if (v.length !== N) throw new Error('Assertion failed.');
      if (bv.length !== 2 * M) throw new Error('Assertion failed.');
      bv.fill(0.0, 0, 2 * m);

      for (var i = 0; i < m; i++) {
        for (var j = 0; j < N; j++) {
          bv[i] += dG[N * i + j] * v[j];
        }
      }

      for (var _i9 = 1; _i9 < m; _i9++) {
        for (var _j8 = 0; _j8 < _i9; _j8++) {
          bv[m + _i9] += dXdG[M * _i9 + _j8] * bv[_j8] / dXdG[M * _j8 + _j8];
        }
      }

      for (var _i10 = 0; _i10 < m; _i10++) {
        bv[_i10] /= D[_i10];
      }

      for (var _i11 = 0; _i11 < m; _i11++) {
        for (var _j9 = 0; _j9 < N; _j9++) {
          bv[m + _i11] += dX[N * _i11 + _j9] * v[_j9] * scale;
        }
      }

      (0, _tri._tril_solve)(m, M, 1, J, 0, bv, m);
    }
  }, {
    key: "compute_be",
    value: function compute_be(j, bej) {
      var m = this.m,
          M = this.M,
          N = this.N,
          dG = this.dG,
          dXdG = this.dXdG,
          dX = this.dX,
          J = this.J,
          D = this.D,
          scale = this.scale;
      if (j % 1 !== 0) throw new Error('Assertion failed.');
      j |= 0;
      if (!(j >= 0)) throw new Error('Assertion failed.');
      if (!(j < N)) throw new Error('Assertion failed.');
      if (bej.length !== 2 * M) throw new Error('Assertion failed.');
      bej.fill(0.0, 0, 2 * m);

      for (var i = 0; i < m; i++) {
        bej[i] += dG[N * i + j];
      }

      for (var _i12 = 1; _i12 < m; _i12++) {
        for (var _j10 = 0; _j10 < _i12; _j10++) {
          bej[m + _i12] += dXdG[M * _i12 + _j10] * bej[_j10] / dXdG[M * _j10 + _j10];
        }
      }

      for (var _i13 = 0; _i13 < m; _i13++) {
        bej[_i13] /= D[_i13];
      }

      for (var _i14 = 0; _i14 < m; _i14++) {
        bej[m + _i14] += dX[N * _i14 + j] * scale;
      }

      (0, _tri._tril_solve)(m, M, 1, J, 0, bej, m);
    }
  }, {
    key: "compute_ubbv",
    value: function compute_ubbv(hu, uv, hv) {
      var m = this.m,
          M = this.M,
          scale = this.scale;
      uv *= scale;
      if (isNaN(uv)) throw new Error('Assertion failed.');
      if (hu.length !== 2 * M) throw new Error('Assertion failed.');
      if (hv.length !== 2 * M) throw new Error('Assertion failed.');
      var result = 0;

      for (var i = 2 * m; i-- > m;) {
        result -= hu[i] * hv[i];
      }

      for (var _i15 = m; _i15-- > 0;) {
        result += hu[_i15] * hv[_i15];
      }

      return result + uv;
    }
  }, {
    key: "compute_cauchyGeneralized",
    value: function compute_cauchyGeneralized(G, X, bounds, indices) {
      var complete = arguments.length > 4 && arguments[4] !== undefined ? arguments[4] : true;
      var m = this.m,
          N = this.N,
          Bg = this.Bg,
          Bdx = this.Bdx,
          Bei = this.Bei,
          scale = this.scale,
          travels = this.travels;
      if (G.length !== N) throw new Error('Assertion failed.');
      if (X.length !== N) throw new Error('Assertion failed.'); //    if(     dX.length !==  N) throw new Error('Assertion failed.');

      if (bounds.length !== 2 * N) throw new Error('Assertion failed.');
      if (indices.length !== N) throw new Error('Assertion failed.');
      /*DEBUG*/

      if (!indices.slice().sort(function (i, j) {
        return i - j;
      }).every(function (i, I) {
        return i === I;
      }))
        /*DEBUG*/
        throw new Error('Assertion failed.');

      for (var i = N; i-- > 0;) {
        if (!(X[i] >= bounds[2 * i + 0])) throw new Error('Assertion failed.');
        if (!(X[i] <= bounds[2 * i + 1])) throw new Error('Assertion failed.');
      }

      for (var _i16 = N; _i16-- > 0;) {
        if (0 === G[_i16]) travels[_i16] = Infinity;else {
          var end = bounds[2 * _i16 + (G[_i16] < 0)];
          travels[_i16] = (X[_i16] - end) / G[_i16];
          if (!(0 <= travels[_i16])) throw new Error('Assertion failed: ' + travels[_i16]);
          var n = 0;

          while (Math.sign(G[_i16]) * (X[_i16] - travels[_i16] * G[_i16] - end) < 0) {
            if (++n > 2) throw new Error('Assertion failed.');
            travels[_i16] = (0, _float64_utils.nextDown)(travels[_i16]);
          }

          if (G[_i16] < 0 && X[_i16] - travels[_i16] * G[_i16] > end) throw new Error('Assertion failed.');
          if (G[_i16] > 0 && X[_i16] - travels[_i16] * G[_i16] < end) throw new Error('Assertion failed.');
        }
      }

      var t = 0,
          n_fix = 0;
      this.compute_bv(G, Bg);
      Bdx.fill(0.0, 0, 2 * m);
      var fp = dot(G, G),
          fpp = this.compute_ubbv(Bg, fp, Bg),
          df = 0; // Make sure zero gradient entries are the very last
      // TODO: entries with travels[i]===0 could be skipped (moved to the back) before heap_sort_gen()

      var order = (0, _heap_sort_gen.heap_sort_gen)(indices, function (i, j) {
        var ti = travels[i],
            tj = travels[j];
        if (ti < tj) return -1;
        if (ti > tj) return +1;
        return Math.abs(G[j]) - Math.abs(G[i]); // move entries with G[i] === 0 to the very back
      });

      var _iterator = _createForOfIteratorHelper(order),
          _step;

      try {
        loop: for (_iterator.s(); !(_step = _iterator.n()).done;) {
          var _i17 = _step.value;

          if (0 === G[_i17]) {
            var _iterator3 = _createForOfIteratorHelper(indices.subarray(n_fix)),
                _step3;

            try {
              for (_iterator3.s(); !(_step3 = _iterator3.n()).done;) {
                var k = _step3.value;
                if (!(Number.MAX_VALUE < travels[k])) throw new Error('Assertion failed.');
                if (0 !== G[k]) throw new Error('Assertion failed.');
              }
            } catch (err) {
              _iterator3.e(err);
            } finally {
              _iterator3.f();
            }

            break loop;
          }

          var _end = bounds[2 * _i17 + (G[_i17] < 0)],
              dt = travels[_i17] - t;
          if (!(0 <= dt)) throw new Error('Assertion failed.');

          if (0 === dt) {
            // <- TODO: entries with travels[i]===0 could be skipped before this loop and _heap_sort()
            X[_i17] = _end;
            t = travels[_i17];
            ++n_fix;
            continue loop;
          } //*DEBUG*/      if( !(0 <= fpp) ) throw new Error('Assertion failed.'); // <- B is supposed to be positive definite


          var cp = Math.min(dt, fp / fpp);
          if (!complete && !(cp >= dt)) // <- handles NaN
            break loop;

          if (0 < cp) {
            for (var j = 2 * m; j-- > 0;) {
              Bdx[j] -= cp * Bg[j];
            }

            t += cp;
          }

          if (complete && !(cp >= dt)) // <- handles NaN
            break loop;
          t = travels[_i17];
          this.compute_be(_i17, Bei);
          var eBx = this.compute_ubbv(Bei, -t * G[_i17], Bdx),
              eBe = this.compute_ubbv(Bei, 1, Bei),
              eBg = this.compute_ubbv(Bei, G[_i17], Bg);
          fp -= G[_i17] * (G[_i17] + eBx) + dt * fpp;
          fpp += G[_i17] * (G[_i17] * eBe - eBg * 2);

          for (var _j11 = 2 * m; _j11-- > 0;) {
            Bg[_j11] -= G[_i17] * Bei[_j11];
          }

          df += G[_i17] * G[_i17] * t * (t * scale / 2 - 1); // df = g.T @ dx + 0.5 * dx.T @ B @ dx = g.T @ dx + 0.5 * scale * dx.T @ dx + compute_ubbv(B @ dx, B @ dx)

          X[_i17] = _end;
          G[_i17] *= 1 - t * scale;
          ++n_fix;
        }
      } catch (err) {
        _iterator.e(err);
      } finally {
        _iterator.f();
      }

      var _iterator2 = _createForOfIteratorHelper(indices.subarray(n_fix)),
          _step2;

      try {
        for (_iterator2.s(); !(_step2 = _iterator2.n()).done;) {
          var _i18 = _step2.value;
          df += G[_i18] * G[_i18] * t * (t * scale / 2 - 1); // df = g.T @ dx + 0.5 * dx.T @ B @ dx = g.T @ dx + 0.5 * scale * dx.T @ dx + compute_ubbv(B @ dx, B @ dx)

          X[_i18] -= t * G[_i18];
          G[_i18] *= 1 - t * scale;
          if (!(X[_i18] >= bounds[2 * _i18 + 0])) throw new Error('Assertion failed.');
          if (!(X[_i18] <= bounds[2 * _i18 + 1])) throw new Error('Assertion failed.');
        }
      } catch (err) {
        _iterator2.e(err);
      } finally {
        _iterator2.f();
      }

      df += 0.5 * this.compute_ubbv(Bdx, 0, Bdx); // compute new gradient (G = G0 + B @ dx)

      this._compute_Bv(Bdx, G); // <- keep in mind this overwrites both Bdx and g


      return [n_fix, df];
    }
  }, {
    key: "compute_subspace_Hv",
    value: function compute_subspace_Hv(v, Hv, indices, n_fix) {
      // For more information, see [2].
      // The goal is to solve the following subspace_problem:
      //        _
      // Hv = Z⋅B⋅Zᵀ⋅v
      //       
      // Where v and Hv are the input and output parameters of
      // of compute_subspace_Hv (see above).
      // _
      // B is the subspace Hessian (approximation):
      // _
      // B = Zᵀ⋅B⋅Z
      //                _
      // The inverse of B can be written as:
      // _                     _           _       _
      // B⁻¹ = I/scale + Zᵀ⋅W⋅(M⁻¹ + scale⋅Wᵀ⋅A⋅Aᵀ⋅W)⋅Wᵀ⋅Z
      //
      //       ┏                                  ┓
      //       ┃                ┊                 ┃
      //       ┃                ┊                 ┃
      //       ┃-(D+dGᵀdG/scale)┊       -Rᵀ       ┃
      //       ┃                ┊                 ┃
      //       ┃                ┊                 ┃
      // _     ┃                ┊                 ┃
      // M⁻¹ = ┃┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┼┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┃
      //       ┃                ┊                 ┃
      //       ┃                ┊                 ┃
      //       ┃                ┊                 ┃
      //       ┃       -R       ┊        0        ┃
      //       ┃                ┊                 ┃
      //       ┃                ┊                 ┃
      //       ┗                                  ┛
      //
      // R = triu(dXᵀdG)
      // D = diag(dXᵀdG)
      //
      var m = this.m,
          M = this.M,
          N = this.N,
          dX = this.dX,
          dXdG = this.dXdG,
          dG = this.dG,
          dGdG = this.dGdG,
          M_WAAW = this.M_WAAW,
          u = this.Bei,
          scale = this.scale;
      if (!(indices instanceof Int32Array)) throw new Error('Assertion failed.');
      if (v.length !== N) throw new Error('Assertion failed.');
      if (Hv.length !== N) throw new Error('Assertion failed.');
      if (indices.length !== N) throw new Error('Assertion failed.');
      if (n_fix % 1 !== 0) throw new Error('Assertion failed.');
      if (!(n_fix >= 0)) throw new Error('Assertion failed.');
      if (!(n_fix <= N)) throw new Error('Assertion failed.');
      /*DEBUG*/

      if (!indices.slice().sort(function (i, j) {
        return i - j;
      }).every(function (i, I) {
        return i === I;
      }))
        /*DEBUG*/
        throw new Error('Assertion failed.'); // /*DEBUG*/    for( let i=0; i < n_fix; i++ )
      // /*DEBUG*/      if( v[indices[i]] !== 0 ) throw new Error('Assertion failed.');
      // indices of the inactive bound constraints

      var free_indices = indices.subarray(n_fix); //
      // RIGHT PART
      //

      for (var i = 0; i < m; i++) {
        var sum = 0;

        var _iterator4 = _createForOfIteratorHelper(free_indices),
            _step4;

        try {
          for (_iterator4.s(); !(_step4 = _iterator4.n()).done;) {
            var j = _step4.value;
            sum += dG[N * i + j] * v[j];
          }
        } catch (err) {
          _iterator4.e(err);
        } finally {
          _iterator4.f();
        }

        u[i] = sum / scale;
      }

      for (var _i19 = 0; _i19 < m; _i19++) {
        var _sum = 0;

        var _iterator5 = _createForOfIteratorHelper(free_indices),
            _step5;

        try {
          for (_iterator5.s(); !(_step5 = _iterator5.n()).done;) {
            var _j12 = _step5.value;
            _sum += dX[N * _i19 + _j12] * v[_j12];
          }
        } catch (err) {
          _iterator5.e(err);
        } finally {
          _iterator5.f();
        }

        u[m + _i19] = _sum;
      } //              _
      // MIDDLE PART (M PART)
      //


      if (0 < n_fix) {
        // indices of the active bound constraints
        var fix_indices = indices.subarray(0, n_fix); // UPPER LEFT BLOCK

        for (var _i20 = 0; _i20 < m; _i20++) {
          for (var _j13 = 0; _j13 <= _i20; _j13++) {
            var _sum2 = 0;

            var _iterator6 = _createForOfIteratorHelper(fix_indices),
                _step6;

            try {
              for (_iterator6.s(); !(_step6 = _iterator6.n()).done;) {
                var k = _step6.value;
                _sum2 += dG[N * _i20 + k] * dG[N * _j13 + k];
              }
            } catch (err) {
              _iterator6.e(err);
            } finally {
              _iterator6.f();
            }

            M_WAAW[2 * m * _i20 + _j13] = _sum2;
          }
        } // LOWER LEFT BLOCK


        for (var _i21 = 0; _i21 < m; _i21++) {
          for (var _j14 = 0; _j14 < m; _j14++) {
            var _sum3 = 0;

            var _iterator7 = _createForOfIteratorHelper(fix_indices),
                _step7;

            try {
              for (_iterator7.s(); !(_step7 = _iterator7.n()).done;) {
                var _k4 = _step7.value;
                _sum3 += dX[N * _i21 + _k4] * dG[N * _j14 + _k4];
              }
            } catch (err) {
              _iterator7.e(err);
            } finally {
              _iterator7.f();
            }

            M_WAAW[2 * m * (m + _i21) + _j14] = _sum3;
          }
        } // LOWER RIGHT BLOCK


        for (var _i22 = 0; _i22 < m; _i22++) {
          for (var _j15 = 0; _j15 <= _i22; _j15++) {
            var _sum4 = 0;

            var _iterator8 = _createForOfIteratorHelper(fix_indices),
                _step8;

            try {
              for (_iterator8.s(); !(_step8 = _iterator8.n()).done;) {
                var _k5 = _step8.value;
                _sum4 += dX[N * _i22 + _k5] * dX[N * _j15 + _k5];
              }
            } catch (err) {
              _iterator8.e(err);
            } finally {
              _iterator8.f();
            }

            M_WAAW[2 * m * (m + _i22) + (m + _j15)] = _sum4 * scale;
          }
        }
      } else {
        for (var _i23 = 0; _i23 < m * 2; _i23++) {
          for (var _j16 = 0; _j16 <= _i23; _j16++) {
            M_WAAW[2 * m * _i23 + _j16] = 0;
          }
        }
      } // UPPER LEFT BLOCK


      for (var _i24 = 0; _i24 < m; _i24++) {
        for (var _j17 = 0; _j17 <= _i24; _j17++) {
          M_WAAW[2 * m * _i24 + _j17] -= dGdG[M * _i24 + _j17];
          M_WAAW[2 * m * _i24 + _j17] /= scale;
        }
      }

      for (var _i25 = 0; _i25 < m; _i25++) {
        M_WAAW[2 * m * _i25 + _i25] -= dXdG[M * _i25 + _i25];
      } // LOWER LEFT BLOCK


      for (var _i26 = 0; _i26 < m; _i26++) {
        for (var _j18 = _i26; _j18 < m; _j18++) {
          M_WAAW[2 * m * (m + _i26) + _j18] -= dXdG[M * _i26 + _j18];
        }
      } // FIXME: It remains to be seen whether or not LDL is accurate enough.
      //        If it isn't, we have to switch to PLDLP (Bunch-Kaufman).


      if (0 < m) {
        if (USE_LDL) {
          (0, _ldl._ldl_decomp)(2 * m, 2 * m, M_WAAW, 0);
          (0, _ldl._ldl_solve)(2 * m, 2 * m, 1, M_WAAW, 0, u, 0);
        } else {
          var tmp = this.Bdx,
              P = this.P;
          (0, _pldlp._pldlp_decomp)(2 * m, 2 * m, M_WAAW, 0, P, 0);
          (0, _pldlp._pldlp_solve)(2 * m, 2 * m, 1, M_WAAW, 0, P, 0, u, 0, tmp);
        }
      } //
      // LEFT PART
      //


      Hv.fill(0.0, 0, N);

      for (var _j19 = 0; _j19 < m; _j19++) {
        var _iterator9 = _createForOfIteratorHelper(free_indices),
            _step9;

        try {
          for (_iterator9.s(); !(_step9 = _iterator9.n()).done;) {
            var _i27 = _step9.value;
            Hv[_i27] += dG[N * _j19 + _i27] * u[_j19];
          }
        } catch (err) {
          _iterator9.e(err);
        } finally {
          _iterator9.f();
        }
      }

      var _iterator10 = _createForOfIteratorHelper(free_indices),
          _step10;

      try {
        for (_iterator10.s(); !(_step10 = _iterator10.n()).done;) {
          var _i29 = _step10.value;
          Hv[_i29] /= scale;
        }
      } catch (err) {
        _iterator10.e(err);
      } finally {
        _iterator10.f();
      }

      for (var _j20 = 0; _j20 < m; _j20++) {
        var _iterator11 = _createForOfIteratorHelper(free_indices),
            _step11;

        try {
          for (_iterator11.s(); !(_step11 = _iterator11.n()).done;) {
            var _i28 = _step11.value;
            Hv[_i28] += dX[N * _j20 + _i28] * u[m + _j20];
          }
        } catch (err) {
          _iterator11.e(err);
        } finally {
          _iterator11.f();
        }
      }

      var _iterator12 = _createForOfIteratorHelper(free_indices),
          _step12;

      try {
        for (_iterator12.s(); !(_step12 = _iterator12.n()).done;) {
          var _i30 = _step12.value;
          Hv[_i30] += v[_i30] / scale;
        }
      } catch (err) {
        _iterator12.e(err);
      } finally {
        _iterator12.f();
      }
    }
  }, {
    key: "compute_Hv",
    value: function compute_Hv(v, Hv) {
      var m = this.m,
          M = this.M,
          N = this.N,
          dX = this.dX,
          dXdG = this.dXdG,
          dG = this.dG,
          dGdG = this.dGdG,
          M_WAAW = this.M_WAAW,
          u = this.Bei,
          scale = this.scale;
      if (v.length !== N) throw new Error('Assertion failed.');
      if (Hv.length !== N) throw new Error('Assertion failed.'); //
      // RIGHT PART
      //

      for (var i = 0; i < m; i++) {
        var sum = 0;

        for (var j = 0; j < N; j++) {
          sum += dG[N * i + j] * v[j];
        }

        u[i] = sum / scale;
      }

      for (var _i31 = 0; _i31 < m; _i31++) {
        var _sum5 = 0;

        for (var _j21 = 0; _j21 < N; _j21++) {
          _sum5 += dX[N * _i31 + _j21] * v[_j21];
        }

        u[m + _i31] = _sum5;
      } //              _
      // MIDDLE PART (M PART)
      //


      for (var _i32 = 0; _i32 < m * 2; _i32++) {
        for (var _j22 = 0; _j22 <= _i32; _j22++) {
          M_WAAW[2 * m * _i32 + _j22] = 0;
        }
      } // UPPER LEFT BLOCK


      for (var _i33 = 0; _i33 < m; _i33++) {
        for (var _j23 = 0; _j23 <= _i33; _j23++) {
          M_WAAW[2 * m * _i33 + _j23] -= dGdG[M * _i33 + _j23];
          M_WAAW[2 * m * _i33 + _j23] /= scale;
        }
      }

      for (var _i34 = 0; _i34 < m; _i34++) {
        M_WAAW[2 * m * _i34 + _i34] -= dXdG[M * _i34 + _i34];
      } // LOWER LEFT BLOCK


      for (var _i35 = 0; _i35 < m; _i35++) {
        for (var _j24 = _i35; _j24 < m; _j24++) {
          M_WAAW[2 * m * (m + _i35) + _j24] -= dXdG[M * _i35 + _j24];
        }
      } // FIXME: It remains to be seen whether or not LDL is accurate enough.
      //        If it isn't, we have to switch to PLDLP (Bunch-Kaufman).


      if (0 < m) {
        if (USE_LDL) {
          (0, _ldl._ldl_decomp)(2 * m, 2 * m, M_WAAW, 0);
          (0, _ldl._ldl_solve)(2 * m, 2 * m, 1, M_WAAW, 0, u, 0);
        } else {
          var tmp = this.Bdx,
              P = this.P;
          (0, _pldlp._pldlp_decomp)(2 * m, 2 * m, M_WAAW, 0, P, 0);
          (0, _pldlp._pldlp_solve)(2 * m, 2 * m, 1, M_WAAW, 0, P, 0, u, 0, tmp);
        }
      } //
      // LEFT PART
      //


      Hv.fill(0.0, 0, N);

      for (var _j25 = 0; _j25 < m; _j25++) {
        for (var _i36 = 0; _i36 < N; _i36++) {
          Hv[_i36] += dG[N * _j25 + _i36] * u[_j25];
        }
      }

      for (var _i37 = 0; _i37 < N; _i37++) {
        Hv[_i37] /= scale;
      }

      for (var _j26 = 0; _j26 < m; _j26++) {
        for (var _i38 = 0; _i38 < N; _i38++) {
          Hv[_i38] += dX[N * _j26 + _i38] * u[m + _j26];
        }
      }

      for (var _i39 = 0; _i39 < N; _i39++) {
        Hv[_i39] += v[_i39] / scale;
      }
    }
  }, {
    key: "J",
    get: function get() {
      var m = this.m,
          M = this.M,
          N = this.N,
          dXdX = this.dXdX,
          dXdG = this.dXdG,
          scale = this.scale,
          J = this._J;

      if (this._J_invalid) {
        this._J_invalid = false;

        for (var i = 0; i < m; i++) {
          for (var j = 0; j <= i; j++) {
            J[M * i + j] = scale * dXdX[M * i + j];
          }
        }

        for (var _i40 = 0; _i40 < m; _i40++) {
          for (var _j27 = 0; _j27 <= _i40; _j27++) {
            for (var k = 0; k < _j27; k++) {
              J[M * _i40 + _j27] += dXdG[M * _i40 + k] / dXdG[M * k + k] * dXdG[M * _j27 + k];
            }
          }
        }

        (0, _cholesky._cholesky_decomp)(m, M, J, 0);
      }

      return J;
    }
  }, {
    key: "scale",
    get: function get() {
      return this._scale;
    },
    set: function set(s) {
      s *= 1;
      if (isNaN(s)) throw new Error('Assertion failed.');
      if (!(0 < s)) throw new Error('Assertion failed.');

      if (this._scale !== s) {
        this._scale = s;
        this._J_invalid = true;
      }
    }
  }]);
  return LBFGSB_Solver;
}();

exports.LBFGSB_Solver = LBFGSB_Solver;