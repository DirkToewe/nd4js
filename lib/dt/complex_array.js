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
exports.Complex128Array = void 0;

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _regenerator = _interopRequireDefault(require("@babel/runtime/regenerator"));

var _classCallCheck2 = _interopRequireDefault(require("@babel/runtime/helpers/classCallCheck"));

var _createClass2 = _interopRequireDefault(require("@babel/runtime/helpers/createClass"));

var _typeof2 = _interopRequireDefault(require("@babel/runtime/helpers/typeof"));

var _complex = require("./complex");

require("util");

function createComplexArrayType(FloatArray) {
  var HANDLER = {
    get: function get(complexArray, property) {
      if ((0, _typeof2["default"])(property) !== 'symbol' && property % 1 === 0) return complexArray.get(property);
      return complexArray[property];
    },
    set: function set(complexArray, property, value) {
      if ((0, _typeof2["default"])(property) !== 'symbol' && property % 1 === 0) {
        complexArray.set(property, value);
        return true; // <- FIXME bounds check?
      }

      complexArray[property] = value;
      return true;
    }
  };

  var ComplexArray =
  /*#__PURE__*/
  function () {
    (0, _createClass2["default"])(ComplexArray, null, [{
      key: "BYTES_PER_ELEMENT",
      get: function get() {
        return FloatArray.BYTES_PER_ELEMENT * 2;
      }
    }, {
      key: "name",
      get: function get() {
        return "Complex".concat(ComplexArray.BYTES_PER_ELEMENT * 8, "Array");
      }
    }]);

    function ComplexArray(buffer, byteOffset, length) {
      (0, _classCallCheck2["default"])(this, ComplexArray);

      if (buffer % 1 === 0) {
        length = buffer;
        this._array = new FloatArray(length * 2);
      } else {
        if (null != length) length *= 2;
        this._array = new FloatArray(buffer, byteOffset, length);
      }

      Object.seal(this);
      return new Proxy(this, HANDLER);
    }

    (0, _createClass2["default"])(ComplexArray, [{
      key: Symbol.iterator,
      value:
      /*#__PURE__*/
      _regenerator["default"].mark(function value() {
        var i;
        return _regenerator["default"].wrap(function value$(_context) {
          while (1) {
            switch (_context.prev = _context.next) {
              case 0:
                i = 0;

              case 1:
                if (!(i < this.length)) {
                  _context.next = 7;
                  break;
                }

                _context.next = 4;
                return this[i];

              case 4:
                i++;
                _context.next = 1;
                break;

              case 7:
              case "end":
                return _context.stop();
            }
          }
        }, value, this);
      })
    }, {
      key: "keys",
      value:
      /*#__PURE__*/
      _regenerator["default"].mark(function keys() {
        var i;
        return _regenerator["default"].wrap(function keys$(_context2) {
          while (1) {
            switch (_context2.prev = _context2.next) {
              case 0:
                i = 0;

              case 1:
                if (!(i < this.length)) {
                  _context2.next = 7;
                  break;
                }

                _context2.next = 4;
                return i;

              case 4:
                i++;
                _context2.next = 1;
                break;

              case 7:
              case "end":
                return _context2.stop();
            }
          }
        }, keys, this);
      })
    }, {
      key: "values",
      value:
      /*#__PURE__*/
      _regenerator["default"].mark(function values() {
        return _regenerator["default"].wrap(function values$(_context3) {
          while (1) {
            switch (_context3.prev = _context3.next) {
              case 0:
                return _context3.delegateYield(this[Symbol.iterator](), "t0", 1);

              case 1:
              case "end":
                return _context3.stop();
            }
          }
        }, values, this);
      })
    }, {
      key: "entries",
      value:
      /*#__PURE__*/
      _regenerator["default"].mark(function entries() {
        var i;
        return _regenerator["default"].wrap(function entries$(_context4) {
          while (1) {
            switch (_context4.prev = _context4.next) {
              case 0:
                i = 0;

              case 1:
                if (!(i < this.length)) {
                  _context4.next = 7;
                  break;
                }

                _context4.next = 4;
                return [i, this[i]];

              case 4:
                i++;
                _context4.next = 1;
                break;

              case 7:
              case "end":
                return _context4.stop();
            }
          }
        }, entries, this);
      })
    }, {
      key: "forEach",
      value: function forEach(callback, thisArg) {
        for (var i = 0; i < source.length; i++) {
          callback.call(thisArg, source[i], i);
        }
      }
    }, {
      key: "map",
      value: function map(mapFn, thisArg) {
        return ComplexArray.from(this, mapFn, thisArg);
      }
    }, {
      key: "reduce",
      value: function reduce(reduceFn, initialValue) {
        var i = 0;

        if (null == initialValue) {
          if (this.length == 0) throw new TypeError('TypeError: Reduce of empty array with no initial value.');
          initialValue = this[i++];
        }

        for (; i < this.length; i++) {
          initialValue = reduceFn(initialValue, this[i], i, this);
        }

        return initialValue;
      }
    }, {
      key: "slice",
      value: function slice(begin, end) {
        if (null != begin) begin *= 2;
        if (null != end) end *= 2;
        return new ComplexArray(this._array.slice(begin, end).buffer);
      }
    }, {
      key: "subarray",
      value: function subarray(begin, end) {
        if (null != begin) begin *= 2;
        if (null != end) end *= 2;

        var sub = this._array.subarray(begin, end);

        return new ComplexArray(sub.buffer, sub.byteOffset, sub.length / 2);
      }
    }, {
      key: "join",
      value: function join(separator) {
        if (null == separator) separator = ', ';
        var result = '';

        for (var i = 0; i < this.length; i++) {
          if (i > 0) result += separator;
          result += this[i].toString();
        }

        return result;
      }
    }, {
      key: "fill",
      value: function fill(value, start, end) {
        if (!(value instanceof _complex.Complex)) value = new _complex.Complex(value);
        if (null == start) start = 0;
        if (null == end) end = this.length;
        if (0 > end) end += this.length;
        if (0 > start) start += this.length;

        for (var i = this._array.length; (i -= 2) >= 0;) {
          this._array[i + 1] = value.im;
          this._array[i + 0] = value.re;
        }

        return this;
      }
    }, {
      key: "get",
      value: function get(index) {
        return new _complex.Complex(this._array[2 * index + 0], this._array[2 * index + 1]);
      }
    }, {
      key: "set",
      value: function set(index, value) {
        index *= 2;

        if (typeof value === 'number' || !isNaN(value)) {
          this._array[index + 0] = value;
          this._array[index + 1] = isNaN(value) ? NaN : 0;
          return;
        }

        if (!('re' in value) || !('im' in value)) value = new _complex.Complex(value);
        this._array[index + 0] = value.re;
        this._array[index + 1] = value.im;
      }
    }, {
      key: Symbol["for"]('nodejs.util.inspect.custom'),
      value: function value(options) {
        return 'Complex128Array ' + util.inspect(Array.from(this), options);
      }
    }, {
      key: "toString",
      value: function toString(max_len) {
        return this.join(',');
      }
    }, {
      key: "buffer",
      get: function get() {
        return this._array.buffer;
      }
    }, {
      key: "byteOffset",
      get: function get() {
        return this._array.byteOffset;
      }
    }, {
      key: "byteLength",
      get: function get() {
        return this._array.byteLength;
      }
    }, {
      key: "length",
      get: function get() {
        return this._array.length >>> 1;
      }
    }, {
      key: Symbol.toStringTag,
      get: function get() {
        throw new Error('not yet implemented.');
      }
    }], [{
      key: "from",
      value: function from(source, mapFn, thisArg) {
        if (!source.hasOwnProperty('length')) source = (0, _toConsumableArray2["default"])(source); // fast-copy complex arrays

        if (null == mapFn) {
          if (source instanceof ComplexArray) return new ComplexArray(FloatArray.from(source._array).buffer, 0, source.length);

          if (source instanceof Int32Array || source instanceof Float32Array || source instanceof Float64Array) {
            var array = new FloatArray(source.length * 2);

            for (var i = source.length; i-- > 0;) {
              array[i << 1] = source[i];
            }

            return new ComplexArray(array.buffer, 0, source.length);
          }
        }

        var result = new ComplexArray(source.length);

        mapFn = mapFn || function (x) {
          return x;
        };

        for (var _i = 0; _i < source.length; _i++) {
          result[_i] = mapFn.call(thisArg, source[_i], _i);
        }

        return result;
      }
    }, {
      key: "of",
      value: function of() {
        var result = new ComplexArray(arguments.length);

        for (var i = 0; i < arguments.length; i++) {
          result[i] = i < 0 || arguments.length <= i ? undefined : arguments[i];
        }

        return result;
      }
    }]);
    return ComplexArray;
  }();

  ;
  return ComplexArray;
}

var Complex128Array = createComplexArrayType(Float64Array);
exports.Complex128Array = Complex128Array;