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
exports.Rastrigin = void 0;

var _defineProperty2 = _interopRequireDefault(require("@babel/runtime/helpers/defineProperty"));

var _regenerator = _interopRequireDefault(require("@babel/runtime/regenerator"));

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _classCallCheck2 = _interopRequireDefault(require("@babel/runtime/helpers/classCallCheck"));

var _createClass2 = _interopRequireDefault(require("@babel/runtime/helpers/createClass"));

var _inherits2 = _interopRequireDefault(require("@babel/runtime/helpers/inherits"));

var _possibleConstructorReturn2 = _interopRequireDefault(require("@babel/runtime/helpers/possibleConstructorReturn"));

var _getPrototypeOf2 = _interopRequireDefault(require("@babel/runtime/helpers/getPrototypeOf"));

var _wrapNativeSuper2 = _interopRequireDefault(require("@babel/runtime/helpers/wrapNativeSuper"));

var _iter = require("../../iter");

var _nd_array = require("../../nd_array");

var _root1d_bisect = require("../root1d_bisect");

function _createSuper(Derived) { var hasNativeReflectConstruct = _isNativeReflectConstruct(); return function () { var Super = (0, _getPrototypeOf2["default"])(Derived), result; if (hasNativeReflectConstruct) { var NewTarget = (0, _getPrototypeOf2["default"])(this).constructor; result = Reflect.construct(Super, arguments, NewTarget); } else { result = Super.apply(this, arguments); } return (0, _possibleConstructorReturn2["default"])(this, result); }; }

function _isNativeReflectConstruct() { if (typeof Reflect === "undefined" || !Reflect.construct) return false; if (Reflect.construct.sham) return false; if (typeof Proxy === "function") return true; try { Date.prototype.toString.call(Reflect.construct(Date, [], function () {})); return true; } catch (e) { return false; } }

// REFERENCES
// ----------
// .. [1] https://en.wikipedia.org/wiki/Rastrigin_function
var NONGLOBAL_MINIMA = function () {
  var g = function g(x) {
    return x + Math.PI * 10 * Math.sin(Math.PI * 2 * x);
  };

  var minima = [];

  for (var x = 1; x < 32; x++) {
    if (!(g(x - 0.25) < 0)) throw new Error('Assertion failed.');
    if (!(g(x + 0.25) > 0)) throw new Error('Assertion failed.');
    minima.push((0, _root1d_bisect.root1d_bisect)(g, x - 0.25, x + 0.25));
  }

  Object.freeze(minima);
  return minima;
}();

var Rastrigin = /*#__PURE__*/function (_Function) {
  (0, _inherits2["default"])(Rastrigin, _Function);

  var _super = _createSuper(Rastrigin);

  function Rastrigin(n) {
    var _this;

    (0, _classCallCheck2["default"])(this, Rastrigin);
    _this = _super.call(this);
    if (0 !== n % 1) throw new Error("new Rastrigin(n): n not a valid integer.");
    if (!(0 < n)) throw new Error("new Rastrigin(n): n must be positive.");

    var rastrigin_fn = function rastrigin_fn(x) {
      x = (0, _nd_array.asarray)(x);
      var N = rastrigin_fn.nIn;
      if (x.ndim < 1) throw new Error("".concat(rastrigin_fn.name, "(x): x.ndim must be at least 1."));
      if (x.shape[x.ndim - 1] !== N) throw new Error("".concat(rastrigin_fn.name, "(x): x.shape[-1] must be ").concat(N, "."));
      var F_shape = x.shape.slice(0, -1),
          F = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length / N);
      x = x.data;

      for (var x_off = x.length, F_off = F.length; F_off-- > 0;) {
        x_off -= N;

        for (var i = N; i-- > 0;) {
          var xi = x[x_off + i];
          F[F_off] += xi * xi + 10 * (1 - Math.cos(Math.PI * 2 * xi));
        }
      }

      return new _nd_array.NDArray(F_shape, F);
    };

    rastrigin_fn.nIn = n;
    rastrigin_fn.roots = rastrigin_fn.minima_global = [Array.from({
      length: n
    }, function () {
      return 0;
    })];
    rastrigin_fn.minima = (0, _defineProperty2["default"])({}, Symbol.iterator, /*#__PURE__*/_regenerator["default"].mark(function _callee() {
      var inner, outer, layer, min, ranges, i;
      return _regenerator["default"].wrap(function _callee$(_context) {
        while (1) {
          switch (_context.prev = _context.next) {
            case 0:
              _context.next = 2;
              return Array.from({
                length: n
              }, function () {
                return 0;
              });

            case 2:
              outer = [0]; // grow the minima in a hypercube from the inside out

              layer = 0;

            case 4:
              if (!(layer < NONGLOBAL_MINIMA.length)) {
                _context.next = 20;
                break;
              }

              min = NONGLOBAL_MINIMA[layer];
              inner = outer;
              outer = [-min].concat((0, _toConsumableArray2["default"])(inner), [+min]);
              ranges = Array.from({
                length: n
              }, function () {
                return outer;
              });
              i = 0;

            case 10:
              if (!(i < n)) {
                _context.next = 17;
                break;
              }

              ranges[i] = [-min, +min];
              return _context.delegateYield(_iter.cartesian_prod.apply(void 0, (0, _toConsumableArray2["default"])(ranges)), "t0", 13);

            case 13:
              ranges[i] = inner;

            case 14:
              i++;
              _context.next = 10;
              break;

            case 17:
              layer++;
              _context.next = 4;
              break;

            case 20:
            case "end":
              return _context.stop();
          }
        }
      }, _callee);
    }));
    Object.defineProperty(rastrigin_fn, 'name', {
      value: "rastrigin".concat(n, "d"),
      writable: false
    });
    Object.setPrototypeOf(rastrigin_fn, Rastrigin.prototype);
    return (0, _possibleConstructorReturn2["default"])(_this, Object.freeze(rastrigin_fn));
  }

  (0, _createClass2["default"])(Rastrigin, [{
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
          var xi = x[off + i];
          G[off + i] = 2 * xi + Math.PI * 20 * Math.sin(Math.PI * 2 * xi);
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
          var xi = x[x_off + i];
          H[H_off + N * i + i] += 2 + Math.PI * Math.PI * 40 * Math.cos(Math.PI * 2 * xi);
        }
      }

      return new _nd_array.NDArray(H_shape, H);
    }
  }, {
    key: "lsq",
    value: function lsq(x) {
      x = (0, _nd_array.asarray)(x);
      var N = this.nIn;
      if (x.ndim < 1) throw new Error("".concat(this.name, ".lsq(x): x.ndim must be at least 1."));
      if (x.shape[x.ndim - 1] !== N) throw new Error("".concat(this.name, ".lsq(x): x.shape[-1] must be ").concat(N, "."));
      var F_shape = x.shape,
          F = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length);
      x = x.data;

      for (var x_off = x.length, F_off = F.length; (F_off -= N) >= 0;) {
        x_off -= N;

        for (var i = N; i-- > 0;) {
          var xi = x[x_off + i];
          F[F_off + i] = Math.sqrt(xi * xi + 10 * (1 - Math.cos(Math.PI * 2 * xi)));
        }
      }

      return new _nd_array.NDArray(F_shape, F);
    }
  }, {
    key: "lsq_jac",
    value: function lsq_jac(x) {
      x = (0, _nd_array.asarray)(x);
      var N = this.nIn;
      if (x.ndim < 1) throw new Error("".concat(this.name, ".lsq_jac(x): x.ndim must be at least 1."));
      if (x.shape[x.ndim - 1] !== N) throw new Error("".concat(this.name, ".lsq_jac(x): x.shape[-1] must be ").concat(N, "."));
      var J_shape = Int32Array.of.apply(Int32Array, (0, _toConsumableArray2["default"])(x.shape).concat([N])),
          J = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length * N);
      x = x.data;

      for (var x_off = x.length, J_off = J.length; (J_off -= N * N) >= 0;) {
        x_off -= N;

        for (var i = N; i-- > 0;) {
          var xi = x[x_off + i];
          var denom = Math.sqrt(xi * xi + 10 * (1 - Math.cos(Math.PI * 2 * xi)));
          if (0 !== denom) J[J_off + N * i + i] = (xi + Math.PI * 10 * Math.sin(Math.PI * 2 * xi)) / denom;
        }
      }

      return new _nd_array.NDArray(J_shape, J);
    }
  }, {
    key: "nOut",
    get: function get() {
      return this.nIn;
    }
  }]);
  return Rastrigin;
}( /*#__PURE__*/(0, _wrapNativeSuper2["default"])(Function));

exports.Rastrigin = Rastrigin;