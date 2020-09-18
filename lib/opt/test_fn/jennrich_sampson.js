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
exports.JennrichSampson = void 0;

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

// REFERENCES
// ----------
// .. [1] "Testing Unconstrained Optimization Software"
//         Jorge J. Moré, Burton S. Garbow, Kenneth E. Hillstrom
//         ACM Transactions on Mathematical Software, Vol 7, No. 1, March 1982, pp. 17-41
// .. [2] "Application of stepwise regression to nonlinear estimation"
//         R. I. Jennrich, P. F. Sampson
//         Technometrics 10 (1968), pp. 64-72
var JennrichSampson = /*#__PURE__*/function (_Function) {
  (0, _inherits2["default"])(JennrichSampson, _Function);

  var _super = _createSuper(JennrichSampson);

  function JennrichSampson(m) {
    var _this;

    (0, _classCallCheck2["default"])(this, JennrichSampson);
    _this = _super.call(this);
    if (0 !== m % 1) throw new Error("new JennrichSampson(m): m not a valid integer.");
    if (!(1 < m)) throw new Error("new JennrichSampson(m): m must be greater than 1.");

    var jennrich_sampson_fn = function jennrich_sampson_fn(x) {
      x = (0, _nd_array.asarray)(x);
      if (x.ndim < 1) throw new Error("".concat(jennrich_sampson_fn.name, "(x): x.ndim must be at least 1."));
      if (x.shape[x.ndim - 1] !== 2) throw new Error("".concat(jennrich_sampson_fn.name, "(x): x.shape[-1] must be 2."));
      var F_shape = x.shape.slice(0, -1),
          F = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length / 2);
      x = x.data;
      var M = jennrich_sampson_fn.nOut;

      for (var x_off = x.length, F_off = F.length; F_off-- > 0;) {
        x_off -= 2;
        var x1 = x[x_off + 0],
            x2 = x[x_off + 1];
        var sum = 0;

        for (var i = M; i-- > 0;) {
          var fi = 4 + 2 * i - Math.exp((i + 1) * x1) - Math.exp((i + 1) * x2);
          sum += fi * fi;
        }

        F[F_off] = sum;
      }

      return new _nd_array.NDArray(F_shape, F);
    };

    if (m === 10) {
      jennrich_sampson_fn.minima_global = jennrich_sampson_fn.minima = [[0.25782521321500883, 0.25782521321500883]];
    }

    jennrich_sampson_fn.roots = [];
    jennrich_sampson_fn.nOut = m;
    Object.defineProperty(jennrich_sampson_fn, 'name', {
      value: "jennrich_sampson".concat(m, "d"),
      writable: false
    });
    Object.setPrototypeOf(jennrich_sampson_fn, JennrichSampson.prototype);
    return (0, _possibleConstructorReturn2["default"])(_this, Object.freeze(jennrich_sampson_fn));
  }

  (0, _createClass2["default"])(JennrichSampson, [{
    key: "grad",
    value: function grad(x) {
      x = (0, _nd_array.asarray)(x);
      if (x.ndim < 1) throw new Error("".concat(this.name, ".grad(x): x.ndim must be at least 1."));
      if (x.shape[x.ndim - 1] !== 2) throw new Error("".concat(this.name, ".grad(x): x.shape[-1] must be 2."));
      var G_shape = x.shape,
          G = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length);
      x = x.data;
      var M = this.nOut;

      for (var off = x.length; off >= 0;) {
        off -= 2;
        var x1 = x[off + 0],
            x2 = x[off + 1];
        var sum = 0;

        for (var i = M; i-- > 0;) {
          var fi = 4 + 2 * i - Math.exp((i + 1) * x1) - Math.exp((i + 1) * x2);
          G[off + 0] -= 2 * fi * (i + 1) * Math.exp((i + 1) * x1);
          G[off + 1] -= 2 * fi * (i + 1) * Math.exp((i + 1) * x2);
        }
      }

      return new _nd_array.NDArray(G_shape, G);
    }
  }, {
    key: "hess",
    value: function hess(x) {
      x = (0, _nd_array.asarray)(x);
      if (x.ndim < 1) throw new Error("".concat(this.name, ".hess(x): x.ndim must be at least 1."));
      if (x.shape[x.ndim - 1] !== 2) throw new Error("".concat(this.name, ".hess(x): x.shape[-1] must be 2."));
      var H_shape = Int32Array.of.apply(Int32Array, (0, _toConsumableArray2["default"])(x.shape).concat([2])),
          H = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length * 2);
      x = x.data;
      var M = this.nOut;

      for (var H_off = H.length, x_off = x.length; (x_off -= 2) >= 0;) {
        H_off -= 2 * 2;
        var x1 = x[x_off + 0],
            x2 = x[x_off + 1];
        var sum = 0;

        for (var i = M; i-- > 0;) {
          var fi = 4 + 2 * i - Math.exp((i + 1) * x1) - Math.exp((i + 1) * x2);
          H[H_off + 0] += 2 * (i + 1) * (i + 1) * (Math.exp((i + 1) * x1 * 2) - fi * Math.exp((i + 1) * x1));
          H[H_off + 3] += 2 * (i + 1) * (i + 1) * (Math.exp((i + 1) * x2 * 2) - fi * Math.exp((i + 1) * x2));
          H[H_off + 1] += 2 * (i + 1) * (i + 1) * Math.exp((i + 1) * (x1 + x2));
        }

        H[H_off + 2] = H[H_off + 1];
      }

      return new _nd_array.NDArray(H_shape, H);
    }
  }, {
    key: "lsq",
    value: function lsq(x) {
      x = (0, _nd_array.asarray)(x);
      if (x.ndim < 1) throw new Error("".concat(this.name, ".lsq(x): x.ndim must be at least 1."));
      if (x.shape[x.ndim - 1] !== 2) throw new Error("".concat(this.name, ".lsq(x): x.shape[-1] must be 2."));
      var M = this.nOut;
      var F_shape = x.shape.slice(),
          F = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length / 2 * M);
      F_shape[F_shape.length - 1] = M;
      x = x.data;

      for (var x_off = x.length, F_off = F.length; (F_off -= M) >= 0;) {
        x_off -= 2;
        var x1 = x[x_off + 0],
            x2 = x[x_off + 1];

        for (var i = M; i-- > 0;) {
          F[F_off + i] = 4 + 2 * i - Math.exp((i + 1) * x1) - Math.exp((i + 1) * x2);
        }
      }

      return new _nd_array.NDArray(F_shape, F);
    }
  }, {
    key: "lsq_jac",
    value: function lsq_jac(x) {
      x = (0, _nd_array.asarray)(x);
      var M = this.nOut;
      if (x.ndim < 1) throw new Error("".concat(this.name, ".lsq_jac(x): x.ndim must be at least 1."));
      if (x.shape[x.ndim - 1] !== 2) throw new Error("".concat(this.name, ".lsq_jac(x): x.shape[-1] must be 2."));
      var J_shape = Int32Array.of.apply(Int32Array, (0, _toConsumableArray2["default"])(x.shape.slice(0, -1)).concat([M, 2])),
          J = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length * M);
      x = x.data;

      for (var x_off = x.length, J_off = J.length; (J_off -= M * 2) >= 0;) {
        x_off -= 2;
        var x1 = x[x_off + 0],
            x2 = x[x_off + 1];

        for (var i = M; i-- > 0;) {
          J[J_off + 2 * i + 0] = -(i + 1) * Math.exp((i + 1) * x1);
          J[J_off + 2 * i + 1] = -(i + 1) * Math.exp((i + 1) * x2);
        }
      }

      return new _nd_array.NDArray(J_shape, J);
    }
  }, {
    key: "nIn",
    get: function get() {
      return 2;
    }
  }]);
  return JennrichSampson;
}( /*#__PURE__*/(0, _wrapNativeSuper2["default"])(Function));

exports.JennrichSampson = JennrichSampson;