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
exports.compare = exports.Comparator = void 0;

var _classCallCheck2 = _interopRequireDefault(require("@babel/runtime/helpers/classCallCheck"));

var _inherits2 = _interopRequireDefault(require("@babel/runtime/helpers/inherits"));

var _possibleConstructorReturn2 = _interopRequireDefault(require("@babel/runtime/helpers/possibleConstructorReturn"));

var _getPrototypeOf2 = _interopRequireDefault(require("@babel/runtime/helpers/getPrototypeOf"));

var _wrapNativeSuper2 = _interopRequireDefault(require("@babel/runtime/helpers/wrapNativeSuper"));

function _createSuper(Derived) { var hasNativeReflectConstruct = _isNativeReflectConstruct(); return function () { var Super = (0, _getPrototypeOf2["default"])(Derived), result; if (hasNativeReflectConstruct) { var NewTarget = (0, _getPrototypeOf2["default"])(this).constructor; result = Reflect.construct(Super, arguments, NewTarget); } else { result = Super.apply(this, arguments); } return (0, _possibleConstructorReturn2["default"])(this, result); }; }

function _isNativeReflectConstruct() { if (typeof Reflect === "undefined" || !Reflect.construct) return false; if (Reflect.construct.sham) return false; if (typeof Proxy === "function") return true; try { Date.prototype.toString.call(Reflect.construct(Date, [], function () {})); return true; } catch (e) { return false; } }

var Comparator = /*#__PURE__*/function (_Function) {
  (0, _inherits2["default"])(Comparator, _Function);

  var _super = _createSuper(Comparator);

  function Comparator(compare_fn) {
    var _this;

    (0, _classCallCheck2["default"])(this, Comparator);

    var self = function self(x, y) {
      var len = x.length;
      if (len !== y.length) throw new Error('arrays.Comparator(x,y): x and y must have same length.');

      for (var i = 0; i < len; i++) {
        var c = compare_fn(x[i], y[i]);
        if (0 !== c) return c;
      }

      return 0;
    };

    self.prototype = Comparator.prototype;
    return (0, _possibleConstructorReturn2["default"])(_this, Object.freeze(self));
  }

  return Comparator;
}( /*#__PURE__*/(0, _wrapNativeSuper2["default"])(Function));

exports.Comparator = Comparator;
var compare = new Comparator(function (x, y) {
  return (x > y) - (x < y);
});
exports.compare = compare;