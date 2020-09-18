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
 * along with ND.JS. If not, see <http://www.gnu.org/licenses/>.
 */

var _interopRequireDefault = require("@babel/runtime/helpers/interopRequireDefault");

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.Rosenbrock = void 0;

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _classCallCheck2 = _interopRequireDefault(require("@babel/runtime/helpers/classCallCheck"));

var _createClass2 = _interopRequireDefault(require("@babel/runtime/helpers/createClass"));

var _inherits2 = _interopRequireDefault(require("@babel/runtime/helpers/inherits"));

var _possibleConstructorReturn2 = _interopRequireDefault(require("@babel/runtime/helpers/possibleConstructorReturn"));

var _getPrototypeOf2 = _interopRequireDefault(require("@babel/runtime/helpers/getPrototypeOf"));

var _wrapNativeSuper2 = _interopRequireDefault(require("@babel/runtime/helpers/wrapNativeSuper"));

var _nd_array = require("../../nd_array");

function _createSuper(Derived) { var hasNativeReflectConstruct = _isNativeReflectConstruct(); return function _createSuperInternal() { var Super = (0, _getPrototypeOf2["default"])(Derived), result; if (hasNativeReflectConstruct) { var NewTarget = (0, _getPrototypeOf2["default"])(this).constructor; result = Reflect.construct(Super, arguments, NewTarget); } else { result = Super.apply(this, arguments); } return (0, _possibleConstructorReturn2["default"])(this, result); }; }

function _isNativeReflectConstruct() { if (typeof Reflect === "undefined" || !Reflect.construct) return false; if (Reflect.construct.sham) return false; if (typeof Proxy === "function") return true; try { Date.prototype.toString.call(Reflect.construct(Date, [], function () {})); return true; } catch (e) { return false; } }

var Rosenbrock = /*#__PURE__*/function (_Function) {
  (0, _inherits2["default"])(Rosenbrock, _Function);

  var _super = _createSuper(Rosenbrock);

  function Rosenbrock(n) {
    var _this;

    (0, _classCallCheck2["default"])(this, Rosenbrock);
    _this = _super.call(this);
    if (0 !== n % 1) throw new Error("new Rosenbrock(n): n not a valid integer.");
    if (!(1 < n)) throw new Error("new Rosenbrock(n): n must be greater than 1.");

    var rosenbrock_fn = function rosenbrock_fn(x) {
      x = (0, _nd_array.asarray)(x);
      var N = rosenbrock_fn.nIn;
      if (x.ndim < 1) throw new Error("".concat(rosenbrock_fn.name, "(x): x.ndim must be at least 1."));
      if (x.shape[x.ndim - 1] !== N) throw new Error("".concat(rosenbrock_fn.name, "(x): x.shape[-1] must be ").concat(N, "."));
      var F_shape = x.shape.slice(0, -1),
          F = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length / N);
      x = x.data;

      for (var x_off = x.length, F_off = F.length; F_off-- > 0;) {
        x_off -= N;

        for (var i = N - 1; i-- > 0;) {
          var xj = x[x_off + i + 1],
              xi = x[x_off + i],
              u = xj - xi * xi,
              v = 1 - xi;
          F[F_off] += 100 * u * u + v * v;
        }
      }

      return new _nd_array.NDArray(F_shape, F);
    };

    if (n == 2) rosenbrock_fn.minima = [[1, 1]];else if (n == 3) rosenbrock_fn.minima = [[1, 1, 1]];
    if (n <= 7) rosenbrock_fn.roots = [Array.from({
      length: n
    }, function () {
      return 1;
    })];
    if (n <= 7) rosenbrock_fn.minima_global = [Array.from({
      length: n
    }, function () {
      return 1;
    })];
    rosenbrock_fn.nIn = n;
    Object.defineProperty(rosenbrock_fn, 'name', {
      value: "rosenbrock".concat(n, "d"),
      writable: false
    });
    Object.setPrototypeOf(rosenbrock_fn, Rosenbrock.prototype);
    return (0, _possibleConstructorReturn2["default"])(_this, Object.freeze(rosenbrock_fn));
  }

  (0, _createClass2["default"])(Rosenbrock, [{
    key: "grad",
    value: function grad(x) {
      x = (0, _nd_array.asarray)(x);
      var N = this.nIn;
      if (x.ndim < 1) throw new Error("".concat(this.name, ".grad(x): x.ndim must be at least 1."));
      if (x.shape[x.ndim - 1] !== N) throw new Error("".concat(this.name, ".grad(x): x.shape[-1] must be ").concat(N, "."));
      var G_shape = x.shape,
          G = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length);
      x = x.data;

      for (var off = x.length; off >= 0;) {
        off -= N;

        for (var i = N; i-- > 0;) {
          var xh = x[off + i - 1],
              xi = x[off + i],
              xj = x[off + i + 1];
          if (i < N - 1) G[off + i] = 400 * (xi * xi - xj) * xi - 2 * (1 - xi);
          if (0 < i) G[off + i] += 200 * (xi - xh * xh);
        }
      }

      return new _nd_array.NDArray(G_shape, G);
    }
  }, {
    key: "hess",
    value: function hess(x) {
      x = (0, _nd_array.asarray)(x);
      var N = this.nIn;
      if (x.ndim < 1) throw new Error("".concat(this.name, ".hess(x): x.ndim must be at least 1."));
      if (x.shape[x.ndim - 1] !== N) throw new Error("".concat(this.name, ".hess(x): x.shape[-1] must be ").concat(N, "."));
      var H_shape = Int32Array.of.apply(Int32Array, (0, _toConsumableArray2["default"])(x.shape).concat([N])),
          H = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length * N);
      x = x.data;

      for (var H_off = H.length, x_off = x.length; (x_off -= N) >= 0;) {
        H_off -= N * N;

        for (var i = N; i-- > 0;) {
          if (0 < i) {
            H[H_off + N * i + i] = 200;
            H[H_off + N * (i - 1) + i] = H[H_off + N * i + (i - 1)] = -400 * x[x_off + i - 1];
          }

          if (i < N - 1) {
            var xi = x[x_off + i];
            H[H_off + N * i + i] -= 400 * x[x_off + i + 1] - 1200 * xi * xi - 2;
          }
        }
      }

      return new _nd_array.NDArray(H_shape, H);
    }
  }, {
    key: "lsq",
    value: function lsq(x) {
      x = (0, _nd_array.asarray)(x);
      var N = this.nIn,
          M = this.nOut;
      if (x.ndim < 1) throw new Error("".concat(this.name, ".lsq(x): x.ndim must be at least 1."));
      if (x.shape[x.ndim - 1] !== N) throw new Error("".concat(this.name, ".lsq(x): x.shape[-1] must be ").concat(N, "."));
      var F_shape = x.shape.slice(),
          F = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length / N * M);
      F_shape[F_shape.length - 1] = M;
      x = x.data;

      for (var x_off = x.length, F_off = F.length; (F_off -= M) >= 0;) {
        x_off -= N;

        for (var i = N - 1; i-- > 0;) {
          var xj = x[x_off + i + 1],
              xi = x[x_off + i],
              u = xj - xi * xi,
              v = 1 - xi;
          F[F_off + 2 * i + 0] = u * 10;
          F[F_off + 2 * i + 1] = v;
        }
      }

      return new _nd_array.NDArray(F_shape, F);
    }
  }, {
    key: "lsq_jac",
    value: function lsq_jac(x) {
      x = (0, _nd_array.asarray)(x);
      var N = this.nIn,
          M = this.nOut;
      if (x.ndim < 1) throw new Error("".concat(this.name, ".lsq_jac(x): x.ndim must be at least 1."));
      if (x.shape[x.ndim - 1] !== N) throw new Error("".concat(this.name, ".lsq_jac(x): x.shape[-1] must be ").concat(N, "."));
      var J_shape = Int32Array.of.apply(Int32Array, (0, _toConsumableArray2["default"])(x.shape.slice(0, -1)).concat([M, N])),
          J = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length * M);
      x = x.data;

      for (var x_off = x.length, J_off = J.length; (J_off -= M * N) >= 0;) {
        x_off -= N;

        for (var i = N - 1; i-- > 0;) {
          var j = i + 1,
              xj = x[x_off + j],
              xi = x[x_off + i];
          J[J_off + N * (2 * i + 0) + i] += -20 * xi;
          J[J_off + N * (2 * i + 0) + j] += 10;
          J[J_off + N * (2 * i + 1) + i] -= 1;
        }
      }

      return new _nd_array.NDArray(J_shape, J);
    }
  }, {
    key: "nOut",
    get: function get() {
      return this.nIn * 2 - 2;
    }
  }]);
  return Rosenbrock;
}( /*#__PURE__*/(0, _wrapNativeSuper2["default"])(Function));

exports.Rosenbrock = Rosenbrock;