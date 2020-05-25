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
exports.KahanSum = void 0;

var _classCallCheck2 = _interopRequireDefault(require("@babel/runtime/helpers/classCallCheck"));

var _createClass2 = _interopRequireDefault(require("@babel/runtime/helpers/createClass"));

var KahanSum = /*#__PURE__*/function () {
  function KahanSum() {
    (0, _classCallCheck2["default"])(this, KahanSum);
    this.sum = this.rst = 0; // <- rest
  }

  (0, _createClass2["default"])(KahanSum, [{
    key: "set",
    value: function set(sum) {
      var rst = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0;
      if (isNaN(sum *= 1)) throw new Error('Assertion failed.');
      if (isNaN(rst *= 1)) throw new Error('Assertion failed.');
      this.sum = sum;
      this.rst = rst;
    }
  }, {
    key: "add",
    value: function add(val) {
      var cor = val - this.rst,
          sum = this.sum + cor;
      this.rst = sum - this.sum - cor;
      this.sum = sum;
    }
  }]);
  return KahanSum;
}();

exports.KahanSum = KahanSum;