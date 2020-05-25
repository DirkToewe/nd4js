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
exports.strong_wolfe = void 0;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _defineProperty2 = _interopRequireDefault(require("@babel/runtime/helpers/defineProperty"));

var _albaali_fletcher = require("./albaali_fletcher");

function ownKeys(object, enumerableOnly) { var keys = Object.keys(object); if (Object.getOwnPropertySymbols) { var symbols = Object.getOwnPropertySymbols(object); if (enumerableOnly) symbols = symbols.filter(function (sym) { return Object.getOwnPropertyDescriptor(object, sym).enumerable; }); keys.push.apply(keys, symbols); } return keys; }

function _objectSpread(target) { for (var i = 1; i < arguments.length; i++) { var source = arguments[i] != null ? arguments[i] : {}; if (i % 2) { ownKeys(Object(source), true).forEach(function (key) { (0, _defineProperty2["default"])(target, key, source[key]); }); } else if (Object.getOwnPropertyDescriptors) { Object.defineProperties(target, Object.getOwnPropertyDescriptors(source)); } else { ownKeys(Object(source)).forEach(function (key) { Object.defineProperty(target, key, Object.getOwnPropertyDescriptor(source, key)); }); } } return target; }

var strong_wolfe = function strong_wolfe() {
  var opt = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : {};
  opt = _objectSpread({}, opt);

  for (var _i = 0, _arr = ['c1', 'fRed', 'c2', 'gRed', 'c3', 'grow']; _i < _arr.length; _i++) {
    var _arr$_i = (0, _slicedToArray2["default"])(_arr[_i], 2),
        old_name = _arr$_i[0],
        new_name = _arr$_i[1];

    if (!(new_name in opt) && old_name in opt) {
      console.warn("strong_wolfe(opt): opt.".concat(old_name, " is deprecated, use opt.").concat(new_name, " instead."));
      opt[new_name] = opt[old_name];
      delete opt[old_name];
    }
  }

  return (0, _albaali_fletcher.albaali_fletcher)(opt);
};

exports.strong_wolfe = strong_wolfe;