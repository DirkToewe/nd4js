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
exports.LBFGS_Solver = void 0;

var _classCallCheck2 = _interopRequireDefault(require("@babel/runtime/helpers/classCallCheck"));

var _createClass2 = _interopRequireDefault(require("@babel/runtime/helpers/createClass"));

var _line_search_error = require("./line_search/line_search_error");

var LBFGS_Solver = /*#__PURE__*/function () {
  function LBFGS_Solver(M, N) {
    (0, _classCallCheck2["default"])(this, LBFGS_Solver);
    this.dXdG = new Float64Array(M);
    this.dX = new Float64Array(M * N);
    this.dG = new Float64Array(M * N);
    this.tmp = new Float64Array(M);
    Object.assign(this, {
      M: M,
      N: N
    });
    this.m = 0;
    this.off = 0;
    Object.seal(this);
  }

  (0, _createClass2["default"])(LBFGS_Solver, [{
    key: "update",
    value: function update(dx, dg) {
      var M = this.M,
          N = this.N,
          dG = this.dG,
          dXdG = this.dXdG,
          dX = this.dX;
      if (dx.length !== N) throw new Error('Assertion failed.');
      if (dg.length !== N) throw new Error('Assertion failed.');
      var dxdg = 0;

      for (var j = N; j-- > 0;) {
        dxdg += dx[j] * dg[j];
      }

      if (!(0 < dxdg)) throw new _line_search_error.LineSearchNoProgressError();
      if (this.m === M) this.off = (this.off + 1) % M;else this.m++;
      var i = (this.off + this.m - 1) % M;
      dXdG[i] = dxdg;

      for (var _j = N; _j-- > 0;) {
        dX[N * i + _j] = dx[_j];
      }

      for (var _j2 = N; _j2-- > 0;) {
        dG[N * i + _j2] = dg[_j2];
      }
    }
  }, {
    key: "forget",
    value: function forget(k) {
      if (0 !== k % 1) throw new Error('Assertion failed.');
      if (!(k >= 0)) throw new Error('Assertion failed.');
      if (!(k <= this.m)) throw new Error('Assertion failed.');
      this.off = (this.off + k) % this.M;
      this.m -= k;
    }
  }, {
    key: "compute_Hv_phase1",
    value: function compute_Hv_phase1(v) {
      var M = this.M,
          N = this.N,
          dG = this.dG,
          m = this.m,
          dXdG = this.dXdG,
          dX = this.dX,
          off = this.off,
          α = this.tmp;
      if (v.length !== N) throw new Error("Assertion failed.");

      for (var I = m; I-- > 0;) {
        var i = off + I;
        i -= (i >= M) * M;
        var αi = 0;

        for (var j = N; j-- > 0;) {
          αi += dX[N * i + j] * v[j];
        }

        α[i] = αi /= dXdG[i];

        for (var _j3 = N; _j3-- > 0;) {
          v[_j3] -= αi * dG[N * i + _j3];
        }
      }
    }
  }, {
    key: "compute_Hv_phase2",
    value: function compute_Hv_phase2(v) {
      var M = this.M,
          N = this.N,
          dG = this.dG,
          m = this.m,
          dXdG = this.dXdG,
          dX = this.dX,
          off = this.off,
          α = this.tmp;
      if (v.length !== N) throw new Error("Assertion failed.");

      for (var I = 0; I < m; I++) {
        var i = off + I;
        i -= (i >= M) * M;
        var βi = 0;

        for (var j = N; j-- > 0;) {
          βi += dG[N * i + j] * v[j];
        }

        βi /= dXdG[i];

        for (var _j4 = N; _j4-- > 0;) {
          v[_j4] += (α[i] - βi) * dX[N * i + _j4];
        }
      }
    }
  }]);
  return LBFGS_Solver;
}();

exports.LBFGS_Solver = LBFGS_Solver;