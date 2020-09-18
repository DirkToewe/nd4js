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
exports.checked_array = checked_array;
exports.IndexOutOfBoundsError = void 0;

var _typeof2 = _interopRequireDefault(require("@babel/runtime/helpers/typeof"));

var _classCallCheck2 = _interopRequireDefault(require("@babel/runtime/helpers/classCallCheck"));

var _inherits2 = _interopRequireDefault(require("@babel/runtime/helpers/inherits"));

var _possibleConstructorReturn2 = _interopRequireDefault(require("@babel/runtime/helpers/possibleConstructorReturn"));

var _getPrototypeOf2 = _interopRequireDefault(require("@babel/runtime/helpers/getPrototypeOf"));

var _wrapNativeSuper2 = _interopRequireDefault(require("@babel/runtime/helpers/wrapNativeSuper"));

var _is_array = require("./is_array");

function _createSuper(Derived) { var hasNativeReflectConstruct = _isNativeReflectConstruct(); return function _createSuperInternal() { var Super = (0, _getPrototypeOf2["default"])(Derived), result; if (hasNativeReflectConstruct) { var NewTarget = (0, _getPrototypeOf2["default"])(this).constructor; result = Reflect.construct(Super, arguments, NewTarget); } else { result = Super.apply(this, arguments); } return (0, _possibleConstructorReturn2["default"])(this, result); }; }

function _isNativeReflectConstruct() { if (typeof Reflect === "undefined" || !Reflect.construct) return false; if (Reflect.construct.sham) return false; if (typeof Proxy === "function") return true; try { Date.prototype.toString.call(Reflect.construct(Date, [], function () {})); return true; } catch (e) { return false; } }

var IndexOutOfBoundsError = /*#__PURE__*/function (_Error) {
  (0, _inherits2["default"])(IndexOutOfBoundsError, _Error);

  var _super = _createSuper(IndexOutOfBoundsError);

  function IndexOutOfBoundsError() {
    (0, _classCallCheck2["default"])(this, IndexOutOfBoundsError);
    return _super.apply(this, arguments);
  }

  return IndexOutOfBoundsError;
}( /*#__PURE__*/(0, _wrapNativeSuper2["default"])(Error));

exports.IndexOutOfBoundsError = IndexOutOfBoundsError;
var ARRAY_BOUNDS_CHECKER = {
  get: function get(array, key, proxy) {
    if ((0, _typeof2["default"])(key) !== 'symbol' && key % 1 === 0) {
      key |= 0;
      if (!(0 <= key && key < array.length)) throw new IndexOutOfBoundsError();
    }

    var val = array[key];
    if (val instanceof Function) return function () {
      for (var _len = arguments.length, args = new Array(_len), _key = 0; _key < _len; _key++) {
        args[_key] = arguments[_key];
      }

      return val.apply(this === proxy ? array : this, args);
    };
    return val;
  },
  set: function set(array, key, val) {
    if ((0, _typeof2["default"])(key) !== 'symbol' && key % 1 === 0) {
      key |= 0;
      if (!(0 <= key && key < array.length)) throw new IndexOutOfBoundsError();
    }

    array[key] = val;
    return true;
  }
};

function checked_array(arr) {
  if (!(0, _is_array.is_array)(arr)) throw new Error('Assertion failed.');
  return new Proxy(arr, ARRAY_BOUNDS_CHECKER);
}