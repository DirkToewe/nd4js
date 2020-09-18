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
exports.istr_parse = istr_parse;
exports.istr_stringify = istr_stringify;
exports.istr_stringify_gen = istr_stringify_gen;

var _regenerator = _interopRequireDefault(require("@babel/runtime/regenerator"));

var _dt = require("../dt");

var _b = require("./b64");

var _ = require(".");

var _nd_array = require("../nd_array");

var _marked = /*#__PURE__*/_regenerator["default"].mark(istr_stringify_gen);

function _createForOfIteratorHelper(o, allowArrayLike) { var it; if (typeof Symbol === "undefined" || o[Symbol.iterator] == null) { if (Array.isArray(o) || (it = _unsupportedIterableToArray(o)) || allowArrayLike && o && typeof o.length === "number") { if (it) o = it; var i = 0; var F = function F() {}; return { s: F, n: function n() { if (i >= o.length) return { done: true }; return { done: false, value: o[i++] }; }, e: function e(_e) { throw _e; }, f: F }; } throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method."); } var normalCompletion = true, didErr = false, err; return { s: function s() { it = o[Symbol.iterator](); }, n: function n() { var step = it.next(); normalCompletion = step.done; return step; }, e: function e(_e2) { didErr = true; err = _e2; }, f: function f() { try { if (!normalCompletion && it["return"] != null) it["return"](); } finally { if (didErr) throw err; } } }; }

function _unsupportedIterableToArray(o, minLen) { if (!o) return; if (typeof o === "string") return _arrayLikeToArray(o, minLen); var n = Object.prototype.toString.call(o).slice(8, -1); if (n === "Object" && o.constructor) n = o.constructor.name; if (n === "Map" || n === "Set") return Array.from(o); if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n)) return _arrayLikeToArray(o, minLen); }

function _arrayLikeToArray(arr, len) { if (len == null || len > arr.length) len = arr.length; for (var i = 0, arr2 = new Array(len); i < len; i++) { arr2[i] = arr[i]; } return arr2; }

function istr_parse(b64_chars) {
  b64_chars = b64_chars[Symbol.iterator]();
  var dtype = [];

  for (;;) {
    var _b64_chars$next = b64_chars.next(),
        value = _b64_chars$next.value,
        done = _b64_chars$next.done;

    if (done) throw new Error('b64_parse(b64_chars): Invalid b64_chars.');
    if (value === '[') break;
    dtype.push(value);
  }

  dtype = dtype.join('').trim();
  if (dtype === '') throw new Error('b64_parse(b64_chars): dtype=object not (yet) supported.');
  (0, _dt._check_dtype)(dtype);
  var shape = [];

  for (var s = [];;) {
    var _b64_chars$next2 = b64_chars.next(),
        _value = _b64_chars$next2.value,
        _done = _b64_chars$next2.done;

    if (_done) throw new Error('b64_parse(b64_chars): Invalid b64_chars.');
    if (_value === ']' && s.length === 0) break;

    if (_value === ']' || _value === ',') {
      var d = s.join('');
      if (isNaN(d)) throw new Error("b64_parse(b64_chars): Invalid shape entry: \"".concat(d, "\"."));
      shape.push(1 * d);
      if (_value === ']') break;
      s = [];
      continue;
    }

    s.push(_value);
  }

  shape = Int32Array.from(shape);
  b64_chars = (0, _b.b64_decode_gen)(b64_chars);
  var DTypeArray = _dt.ARRAY_TYPES[dtype],
      data = new Uint8Array(shape.reduce(function (m, n) {
    return m * n;
  }, DTypeArray.BYTES_PER_ELEMENT)),
      word = new Uint8Array(DTypeArray.BYTES_PER_ELEMENT);

  for (var i = 0; i < data.length;) {
    for (var j = 0; j < word.length; j++) {
      var _b64_chars$next3 = b64_chars.next(),
          _value2 = _b64_chars$next3.value,
          _done2 = _b64_chars$next3.done;

      if (_done2) throw new Error('b64_parse(b64_chars): b64_chars too short.');
      word[j] = _value2;
    }

    if (!_.IS_LITTLE_ENDIAN) word.reverse();

    for (var _j = 0; _j < word.length; _j++, i++) {
      data[i] = word[_j];
    }
  }

  return new _nd_array.NDArray(shape, new DTypeArray(data.buffer));
}

function istr_stringify(ndarray) {
  var options = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : {};
  var result = '';

  var _iterator = _createForOfIteratorHelper(istr_stringify_gen(ndarray, options)),
      _step;

  try {
    for (_iterator.s(); !(_step = _iterator.n()).done;) {
      var c = _step.value;
      result += c;
    }
  } catch (err) {
    _iterator.e(err);
  } finally {
    _iterator.f();
  }

  return result;
}

function istr_stringify_gen(ndarray) {
  var options,
      data,
      _args = arguments;
  return _regenerator["default"].wrap(function istr_stringify_gen$(_context) {
    while (1) {
      switch (_context.prev = _context.next) {
        case 0:
          options = _args.length > 1 && _args[1] !== undefined ? _args[1] : {};
          ndarray = (0, _nd_array.asarray)(ndarray);
          data = ndarray.data;

          if (!(ndarray.dtype === 'object')) {
            _context.next = 5;
            break;
          }

          throw new Error("b64_format(A): A.dtype='".concat(ndarray.dtype, "' not supported."));

        case 5:
          if (_.IS_LITTLE_ENDIAN) {
            _context.next = 7;
            break;
          }

          throw new Error('Big endianness not (yet) supported.');

        case 7:
          return _context.delegateYield("".concat(ndarray.dtype, "[").concat(ndarray.shape.join(','), "]\n"), "t0", 8);

        case 8:
          return _context.delegateYield((0, _b.b64_encode_gen)(new Uint8Array(data.buffer, data.byteOffset, data.byteLength), options), "t1", 9);

        case 9:
        case "end":
          return _context.stop();
      }
    }
  }, _marked);
}