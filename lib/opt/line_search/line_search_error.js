'use strict';
/* This file is part of ND.JS.
 *
 * ND.JS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ND.JS is distributed in the hope that it will be useful,
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
exports.LineSearchBoundReachedError = exports.LineSearchBisectionError = exports.LineSearchNoProgressError = exports.LineSearchError = void 0;

var _assertThisInitialized2 = _interopRequireDefault(require("@babel/runtime/helpers/assertThisInitialized"));

var _classCallCheck2 = _interopRequireDefault(require("@babel/runtime/helpers/classCallCheck"));

var _inherits2 = _interopRequireDefault(require("@babel/runtime/helpers/inherits"));

var _possibleConstructorReturn2 = _interopRequireDefault(require("@babel/runtime/helpers/possibleConstructorReturn"));

var _getPrototypeOf2 = _interopRequireDefault(require("@babel/runtime/helpers/getPrototypeOf"));

var _wrapNativeSuper2 = _interopRequireDefault(require("@babel/runtime/helpers/wrapNativeSuper"));

function _createSuper(Derived) { var hasNativeReflectConstruct = _isNativeReflectConstruct(); return function _createSuperInternal() { var Super = (0, _getPrototypeOf2["default"])(Derived), result; if (hasNativeReflectConstruct) { var NewTarget = (0, _getPrototypeOf2["default"])(this).constructor; result = Reflect.construct(Super, arguments, NewTarget); } else { result = Super.apply(this, arguments); } return (0, _possibleConstructorReturn2["default"])(this, result); }; }

function _isNativeReflectConstruct() { if (typeof Reflect === "undefined" || !Reflect.construct) return false; if (Reflect.construct.sham) return false; if (typeof Proxy === "function") return true; try { Date.prototype.toString.call(Reflect.construct(Date, [], function () {})); return true; } catch (e) { return false; } }

var LineSearchError = /*#__PURE__*/function (_Error) {
  (0, _inherits2["default"])(LineSearchError, _Error);

  var _super = _createSuper(LineSearchError);

  function LineSearchError() {
    (0, _classCallCheck2["default"])(this, LineSearchError);
    return _super.apply(this, arguments);
  }

  return LineSearchError;
}( /*#__PURE__*/(0, _wrapNativeSuper2["default"])(Error));

exports.LineSearchError = LineSearchError;

var LineSearchNoProgressError = /*#__PURE__*/function (_LineSearchError) {
  (0, _inherits2["default"])(LineSearchNoProgressError, _LineSearchError);

  var _super2 = _createSuper(LineSearchNoProgressError);

  function LineSearchNoProgressError() {
    (0, _classCallCheck2["default"])(this, LineSearchNoProgressError);
    return _super2.apply(this, arguments);
  }

  return LineSearchNoProgressError;
}(LineSearchError);

exports.LineSearchNoProgressError = LineSearchNoProgressError;

var LineSearchBisectionError = /*#__PURE__*/function (_LineSearchError2) {
  (0, _inherits2["default"])(LineSearchBisectionError, _LineSearchError2);

  var _super3 = _createSuper(LineSearchBisectionError);

  function LineSearchBisectionError(x, f, g) {
    var _this;

    (0, _classCallCheck2["default"])(this, LineSearchBisectionError);
    _this = _super3.call(this);
    Object.assign((0, _assertThisInitialized2["default"])(_this), {
      x: x,
      f: f,
      g: g
    });
    return _this;
  }

  return LineSearchBisectionError;
}(LineSearchError);

exports.LineSearchBisectionError = LineSearchBisectionError;

var LineSearchBoundReachedError = /*#__PURE__*/function (_LineSearchError3) {
  (0, _inherits2["default"])(LineSearchBoundReachedError, _LineSearchError3);

  var _super4 = _createSuper(LineSearchBoundReachedError);

  function LineSearchBoundReachedError(x, f, g) {
    var _this2;

    (0, _classCallCheck2["default"])(this, LineSearchBoundReachedError);
    _this2 = _super4.call(this);
    Object.assign((0, _assertThisInitialized2["default"])(_this2), {
      x: x,
      f: f,
      g: g
    });
    return _this2;
  }

  return LineSearchBoundReachedError;
}(LineSearchError);

exports.LineSearchBoundReachedError = LineSearchBoundReachedError;