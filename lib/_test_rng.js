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
exports.TestRNG = void 0;

var _classCallCheck2 = _interopRequireDefault(require("@babel/runtime/helpers/classCallCheck"));

var _createClass2 = _interopRequireDefault(require("@babel/runtime/helpers/createClass"));

var _inherits2 = _interopRequireDefault(require("@babel/runtime/helpers/inherits"));

var _possibleConstructorReturn2 = _interopRequireDefault(require("@babel/runtime/helpers/possibleConstructorReturn"));

var _getPrototypeOf2 = _interopRequireDefault(require("@babel/runtime/helpers/getPrototypeOf"));

var _dt = require("./dt");

var _tabulate = require("./tabulate");

var _matmul = require("./la/matmul");

var _alea_rng = require("./rand/alea_rng");

function _createSuper(Derived) { var hasNativeReflectConstruct = _isNativeReflectConstruct(); return function _createSuperInternal() { var Super = (0, _getPrototypeOf2["default"])(Derived), result; if (hasNativeReflectConstruct) { var NewTarget = (0, _getPrototypeOf2["default"])(this).constructor; result = Reflect.construct(Super, arguments, NewTarget); } else { result = Super.apply(this, arguments); } return (0, _possibleConstructorReturn2["default"])(this, result); }; }

function _isNativeReflectConstruct() { if (typeof Reflect === "undefined" || !Reflect.construct) return false; if (Reflect.construct.sham) return false; if (typeof Proxy === "function") return true; try { Date.prototype.toString.call(Reflect.construct(Date, [], function () {})); return true; } catch (e) { return false; } }

var TestRNG = /*#__PURE__*/function (_AleaRNG) {
  (0, _inherits2["default"])(TestRNG, _AleaRNG);

  var _super = _createSuper(TestRNG);

  function TestRNG() {
    (0, _classCallCheck2["default"])(this, TestRNG);
    return _super.apply(this, arguments);
  }

  (0, _createClass2["default"])(TestRNG, [{
    key: "rankDef",
    value: function rankDef() {
      var _this = this;

      for (var _len = arguments.length, shape = new Array(_len), _key = 0; _key < _len; _key++) {
        shape[_key] = arguments[_key];
      }

      if (!(shape.length >= 2)) throw new Error('TestRNG::rankDef(...shape): shape.length must be at least 2.');
      var dtype = shape[0] in _dt.ARRAY_TYPES ? shape.shift() : 'float64';
      if (!shape.every(function (s) {
        return s % 1 === 0;
      })) throw new Error('TestRNG::rankDef(...shape): shape must be all integers.');
      if (!shape.every(function (s) {
        return s > 0;
      })) throw new Error('TestRNG::rankDef(...shape): shape must be all positive integers.');
      var N = shape.pop(),
          M = shape.pop(),
          L = Math.min(M, N); // use random, mostly rank-deficient SVD to generate test matrix A

      var ranks = (0, _tabulate.tabulate)(shape, 'int32', function () {
        return _this["int"](0, L + 1);
      }),
          // <- ranks
      u = this.ortho.apply(this, [dtype].concat(shape, [M, L])),
          v = this.ortho.apply(this, [dtype].concat(shape, [L, N])),
          V = v.data,
          R = ranks.data; // V = S @ V, where S is a batch of diagonal scaling matrices with given `ranks`

      for (var i = R.length; i-- > 0;) {
        for (var j = L; j-- > 0;) {
          var scale = R[i] <= j ? 0 : this.uniform(1e-4, 1e+4);

          for (var k = N; k-- > 0;) {
            V[N * (L * i + j) + k] *= scale;
          }
        }
      }

      var a = (0, _matmul.matmul2)(u, v); // <- A = U @ S @ V

      Object.freeze(ranks.data.buffer);
      Object.freeze(a.data.buffer);
      return [a, ranks];
    }
  }]);
  return TestRNG;
}(_alea_rng.AleaRNG);

exports.TestRNG = TestRNG;