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
exports.OptimizationNoProgressError = exports.OptimizationError = void 0;

var _classCallCheck2 = _interopRequireDefault(require("@babel/runtime/helpers/classCallCheck"));

var _inherits2 = _interopRequireDefault(require("@babel/runtime/helpers/inherits"));

var _possibleConstructorReturn2 = _interopRequireDefault(require("@babel/runtime/helpers/possibleConstructorReturn"));

var _getPrototypeOf2 = _interopRequireDefault(require("@babel/runtime/helpers/getPrototypeOf"));

var _wrapNativeSuper2 = _interopRequireDefault(require("@babel/runtime/helpers/wrapNativeSuper"));

function _createSuper(Derived) { var hasNativeReflectConstruct = _isNativeReflectConstruct(); return function _createSuperInternal() { var Super = (0, _getPrototypeOf2["default"])(Derived), result; if (hasNativeReflectConstruct) { var NewTarget = (0, _getPrototypeOf2["default"])(this).constructor; result = Reflect.construct(Super, arguments, NewTarget); } else { result = Super.apply(this, arguments); } return (0, _possibleConstructorReturn2["default"])(this, result); }; }

function _isNativeReflectConstruct() { if (typeof Reflect === "undefined" || !Reflect.construct) return false; if (Reflect.construct.sham) return false; if (typeof Proxy === "function") return true; try { Date.prototype.toString.call(Reflect.construct(Date, [], function () {})); return true; } catch (e) { return false; } }

var OptimizationError = /*#__PURE__*/function (_Error) {
  (0, _inherits2["default"])(OptimizationError, _Error);

  var _super = _createSuper(OptimizationError);

  function OptimizationError() {
    (0, _classCallCheck2["default"])(this, OptimizationError);
    return _super.apply(this, arguments);
  }

  return OptimizationError;
}( /*#__PURE__*/(0, _wrapNativeSuper2["default"])(Error));

exports.OptimizationError = OptimizationError;

var OptimizationNoProgressError = /*#__PURE__*/function (_OptimizationError) {
  (0, _inherits2["default"])(OptimizationNoProgressError, _OptimizationError);

  var _super2 = _createSuper(OptimizationNoProgressError);

  function OptimizationNoProgressError() {
    (0, _classCallCheck2["default"])(this, OptimizationNoProgressError);
    return _super2.apply(this, arguments);
  }

  return OptimizationNoProgressError;
}(OptimizationError);

exports.OptimizationNoProgressError = OptimizationNoProgressError;