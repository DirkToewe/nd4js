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
exports.Complex = void 0;

var _classCallCheck2 = _interopRequireDefault(require("@babel/runtime/helpers/classCallCheck"));

var _createClass2 = _interopRequireDefault(require("@babel/runtime/helpers/createClass"));

var _math = _interopRequireDefault(require("../math"));

var Complex = /*#__PURE__*/function () {
  function Complex(re, im) {
    (0, _classCallCheck2["default"])(this, Complex);
    if (null == re) re = 0.0;
    if (null == im) im = 0.0;

    if (re.re != null || re.im != null) {
      im = re.im || 0;
      re = re.re || 0;
    }

    this.re = re * 1;
    this.im = im * 1; // TODO: remove these assertions

    if (isNaN(re)) throw new Error('Real value is NaN.');
    if (isNaN(im)) throw new Error('Imaginary value is NaN.');
    Object.freeze(this);
  }

  (0, _createClass2["default"])(Complex, [{
    key: "add",
    value: function add(re, im) {
      if (null == re) re = 0.0;
      if (null == im) im = 0.0;

      if (re.re != null || re.im != null) {
        im = re.im || 0;
        re = re.re || 0;
      }

      return new Complex(this.re + re, this.im + im);
    }
  }, {
    key: "sub",
    value: function sub(re, im) {
      if (null == re) re = 0.0;
      if (null == im) im = 0.0;

      if (re.re != null || re.im != null) {
        im = re.im || 0;
        re = re.re || 0;
      }

      return new Complex(this.re - re, this.im - im);
    }
  }, {
    key: "mul",
    value: function mul(re, im) {
      if (null == re) re = 0.0;
      if (null == im) im = 0.0;

      if (re.re != null || re.im != null) {
        im = re.im || 0;
        re = re.re || 0;
      }

      return new Complex(this.re * re - this.im * im, this.re * im + this.im * re);
    }
  }, {
    key: "div",
    value: function div(re, im) {
      if (null == re) re = 0.0;
      if (null == im) im = 0.0;

      if (re.re != null || re.im != null) {
        im = re.im || 0;
        re = re.re || 0;
      }

      if (im == 0) return new Complex(this.re / re, this.im / re); //    return new Complex(
      //      (this.re*re + this.im*im) / (re*re + im*im),
      //      (this.im*re - this.re*im) / (re*re + im*im)
      //    );

      if (Math.abs(re) >= Math.abs(im)) {
        var R = im / re;
        return new Complex((this.re + this.im * R) / (re + im * R), (this.im - this.re * R) / (re + im * R));
      } else {
        var _R = re / im;

        return new Complex((this.re * _R + this.im) / (re * _R + im), (this.im * _R - this.re) / (re * _R + im));
      }
    }
  }, {
    key: "equals",
    value: function equals(re, im) {
      if (null == re) re = 0.0;
      if (null == im) im = 0.0;

      if (re.re != null || re.im != null) {
        im = re.im || 0;
        re = re.re || 0;
      }

      return this.re == re && this.im == im;
    }
  }, {
    key: "abs",
    value: function abs() {
      return Math.hypot(this.re, this.im);
    }
  }, {
    key: "arg",
    value: function arg() {
      return Math.atan2(this.im, this.re);
    }
  }, {
    key: "conj",
    value: function conj() {
      return new Complex(this.re, -this.im);
    }
  }, {
    key: "toFixed",
    value: function toFixed(digits) {
      return new Complex(this.re.toFixed(digits), this.im.toFixed(digits));
    }
  }, {
    key: "sqrt",
    value: function sqrt() {
      // https://en.wikipedia.org/wiki/Square_root#Algebraic_formula
      var abs = this.abs();
      return new Complex(_math["default"].sqrt(_math["default"].mul(_math["default"].add(abs, this.re), 0.5)), _math["default"].sqrt(_math["default"].mul(_math["default"].sub(abs, this.re), 0.5)));
    } // https://nodejs.org/api/util.html#util_util_inspect_custom

  }, {
    key: Symbol["for"]('nodejs.util.inspect.custom'),
    value: function value() {
      return this.toString();
    }
  }, {
    key: "toString",
    value: function toString() {
      if (this.im == 0) return this.re.toString();
      if (this.re == 0) return "".concat(this.im, "j");
      if (this.im < 0) return "".concat(this.re, " - ").concat(-this.im, "j");
      return "".concat(this.re, " + ").concat(this.im, "j");
    }
  }, {
    key: "valueOf",
    value: function valueOf() {
      if (this.im == 0) return this.re;
      return this;
    }
  }]);
  return Complex;
}();

exports.Complex = Complex;