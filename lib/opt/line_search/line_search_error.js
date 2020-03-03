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
exports.LineSearchNoProgressError = exports.LineSearchError = void 0;

var _classCallCheck2 = _interopRequireDefault(require("@babel/runtime/helpers/classCallCheck"));

var _possibleConstructorReturn2 = _interopRequireDefault(require("@babel/runtime/helpers/possibleConstructorReturn"));

var _getPrototypeOf2 = _interopRequireDefault(require("@babel/runtime/helpers/getPrototypeOf"));

var _inherits2 = _interopRequireDefault(require("@babel/runtime/helpers/inherits"));

var _wrapNativeSuper2 = _interopRequireDefault(require("@babel/runtime/helpers/wrapNativeSuper"));

var LineSearchError = /*#__PURE__*/function (_Error) {
  (0, _inherits2["default"])(LineSearchError, _Error);

  function LineSearchError() {
    (0, _classCallCheck2["default"])(this, LineSearchError);
    return (0, _possibleConstructorReturn2["default"])(this, (0, _getPrototypeOf2["default"])(LineSearchError).apply(this, arguments));
  }

  return LineSearchError;
}( /*#__PURE__*/(0, _wrapNativeSuper2["default"])(Error));

exports.LineSearchError = LineSearchError;

var LineSearchNoProgressError = /*#__PURE__*/function (_LineSearchError) {
  (0, _inherits2["default"])(LineSearchNoProgressError, _LineSearchError);

  function LineSearchNoProgressError() {
    (0, _classCallCheck2["default"])(this, LineSearchNoProgressError);
    return (0, _possibleConstructorReturn2["default"])(this, (0, _getPrototypeOf2["default"])(LineSearchNoProgressError).apply(this, arguments));
  }

  return LineSearchNoProgressError;
}(LineSearchError);

exports.LineSearchNoProgressError = LineSearchNoProgressError;