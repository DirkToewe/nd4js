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
exports.SingularMatrixSolveError = void 0;

var _classCallCheck2 = _interopRequireDefault(require("@babel/runtime/helpers/classCallCheck"));

var _possibleConstructorReturn2 = _interopRequireDefault(require("@babel/runtime/helpers/possibleConstructorReturn"));

var _getPrototypeOf3 = _interopRequireDefault(require("@babel/runtime/helpers/getPrototypeOf"));

var _inherits2 = _interopRequireDefault(require("@babel/runtime/helpers/inherits"));

var _wrapNativeSuper2 = _interopRequireDefault(require("@babel/runtime/helpers/wrapNativeSuper"));

var _nd_array = require("../nd_array");

var SingularMatrixSolveError =
/*#__PURE__*/
function (_Error) {
  (0, _inherits2["default"])(SingularMatrixSolveError, _Error);

  function SingularMatrixSolveError(x) {
    var _getPrototypeOf2;

    var _this;

    (0, _classCallCheck2["default"])(this, SingularMatrixSolveError);

    for (var _len = arguments.length, args = new Array(_len > 1 ? _len - 1 : 0), _key = 1; _key < _len; _key++) {
      args[_key - 1] = arguments[_key];
    }

    _this = (0, _possibleConstructorReturn2["default"])(this, (_getPrototypeOf2 = (0, _getPrototypeOf3["default"])(SingularMatrixSolveError)).call.apply(_getPrototypeOf2, [this].concat(args)));
    if (!(x instanceof _nd_array.NDArray)) throw new Error('Assertion failed.');
    _this.x = x;
    return _this;
  }

  return SingularMatrixSolveError;
}((0, _wrapNativeSuper2["default"])(Error));

exports.SingularMatrixSolveError = SingularMatrixSolveError;