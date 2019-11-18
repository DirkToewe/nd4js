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
exports.npy_serialize = npy_serialize;
exports.npy_serialize_gen = npy_serialize_gen;
exports.npy_deserialize = npy_deserialize;

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _regenerator = _interopRequireDefault(require("@babel/runtime/regenerator"));

var _pyon = require("./pyon");

var _dt = require("../dt");

var _nd_array = require("../nd_array");

var _ = require(".");

var _marked =
/*#__PURE__*/
_regenerator["default"].mark(npy_serialize_gen);

var MAGIC_STRING = "\x93NUMPY";

function npy_serialize(A) {
  return Uint8Array.from(npy_serialize_gen(A));
}

function npy_serialize_gen(A) {
  var dt, header, i, headerLen, _i, _i2;

  return _regenerator["default"].wrap(function npy_serialize_gen$(_context) {
    while (1) {
      switch (_context.prev = _context.next) {
        case 0:
          A = (0, _nd_array.asarray)(A);

          dt = function () {
            switch (A.dtype) {
              default:
                throw new Error("nd_to_npy: A.dtype=".concat(A.dtype, " not yet supported."));

              case 'int32':
                return 'i4';

              case 'float32':
                return 'f4';

              case 'float64':
                return 'f8';

              case 'complex128':
                return 'c16';
            }
          }();

          if (_.IS_LITTLE_ENDIAN) dt = '<' + dt;else dt = '>' + dt;
          header = "{\"descr\": \"".concat(dt, "\", \"fortran_order\": False, \"shape\": (").concat(A.shape).concat(A.ndim > 0 ? ',' : '', ")}"); // MAGIC STRING

          i = 0;

        case 5:
          if (!(i < MAGIC_STRING.length)) {
            _context.next = 11;
            break;
          }

          _context.next = 8;
          return MAGIC_STRING.codePointAt(i);

        case 8:
          i++;
          _context.next = 5;
          break;

        case 11:
          _context.next = 13;
          return 1;

        case 13:
          _context.next = 15;
          return 0;

        case 15:
          headerLen = header.length + 11 + 63 >>> 6 << 6;

          if (!(headerLen > 0xFFFF)) {
            _context.next = 18;
            break;
          }

          throw new Error('nd_to_npy: Header too large.');

        case 18:
          _context.next = 20;
          return headerLen - 10 >>> 0 & 255;

        case 20:
          _context.next = 22;
          return headerLen - 10 >>> 8 & 255;

        case 22:
          if (!(headerLen % 64 !== 0)) {
            _context.next = 24;
            break;
          }

          throw new Error('Assertion failed.');

        case 24:
          _i = 0;

        case 25:
          if (!(_i < header.length)) {
            _context.next = 31;
            break;
          }

          _context.next = 28;
          return header.codePointAt(_i);

        case 28:
          _i++;
          _context.next = 25;
          break;

        case 31:
          _i2 = header.length + 11;

        case 32:
          if (!(_i2 < headerLen)) {
            _context.next = 38;
            break;
          }

          _context.next = 35;
          return ' '.codePointAt(0);

        case 35:
          _i2++;
          _context.next = 32;
          break;

        case 38:
          _context.next = 40;
          return '\n'.codePointAt(0);

        case 40:
          return _context.delegateYield(new Uint8Array(A.data.buffer, A.data.byteOffset, A.data.byteLength), "t0", 41);

        case 41:
        case "end":
          return _context.stop();
      }
    }
  }, _marked);
}

function npy_deserialize(npy_bytes) {
  var _result;

  var nRead = 0;

  var next = function () {
    var iter = npy_bytes[Symbol.iterator]();
    return function () {
      var _iter$next = iter.next(),
          value = _iter$next.value,
          done = _iter$next.done;

      ++nRead;
      if (done) throw new Error('npy_to_nd: byte sequence ended unexpectedly.');
      return value;
    };
  }();

  npy_bytes = undefined;
  var magic_str = '';

  while (magic_str.length < 6) {
    magic_str += String.fromCodePoint(next());
  }

  if (magic_str !== MAGIC_STRING) throw new Error("npy_to_nd: byte sequence does not start with '\x93NUMPY'.");
  var version = "".concat(next(), ".").concat(next());

  switch (version) {
    default:
      throw new Error("npy_bytes: npy-file version ".concat(version, " not supported."));

    case '1.0':
    case '2.0':
  }

  var header = function () {
    var headerLen = next() * (1 << 0) + next() * (1 << 8);

    if (version === '2.0') {
      headerLen += next() * (1 << 16) + next() * (1 << 24);
    }

    if (headerLen < 0) throw new Error('Assertion failed.');
    return (0, _pyon.pyon_parse)(
    /*#__PURE__*/
    _regenerator["default"].mark(function _callee() {
      return _regenerator["default"].wrap(function _callee$(_context2) {
        while (1) {
          switch (_context2.prev = _context2.next) {
            case 0:
              if (!(headerLen-- > 0)) {
                _context2.next = 5;
                break;
              }

              _context2.next = 3;
              return String.fromCodePoint(next());

            case 3:
              _context2.next = 0;
              break;

            case 5:
            case "end":
              return _context2.stop();
          }
        }
      }, _callee);
    })());
  }();

  var LE = function () {
    switch (header.descr[0]) {
      default:
        throw new Error("npy_to_nd: dtype '".concat(header.descr, "' not yet supported."));

      case '<':
        return true;

      case '>':
        return false;
    }
  }();

  var dtype = function () {
    switch (header.descr.slice(1)) {
      default:
        throw new Error("npy_to_nd: dtype '".concat(header.descr, "' not yet supported."));

      case 'i4':
        return 'int32';

      case 'f4':
        return 'float32';

      case 'f8':
        return 'float64';

      case 'c16':
        return 'complex128';
    }
  }();

  var DTypeArray = _dt.ARRAY_TYPES[dtype];
  var data = Uint8Array.from({
    length: header.shape.reduce(function (m, n) {
      return m * n;
    }, DTypeArray.BYTES_PER_ELEMENT)
  }, next);

  if (LE !== _.IS_LITTLE_ENDIAN) {
    var wordLen = DTypeArray.BYTES_PER_ELEMENT >>> dtype.startsWith('complex');

    for (var off = 0; off < data.length; off += wordLen) {
      for (var i = off, j = off + wordLen; i < --j; i++) {
        var data_i = data[i];
        data[i] = data[j];
        data[j] = data_i;
      }
    }
  }

  var shape = Int32Array.from(header.shape);
  if (header.fortran_order) shape.reverse();
  var result = new _nd_array.NDArray(shape, new DTypeArray(data.buffer));
  if (header.fortran_order && result.ndim > 1) result = (_result = result).transpose.apply(_result, (0, _toConsumableArray2["default"])(
  /*#__PURE__*/
  _regenerator["default"].mark(function _callee2() {
    var _i3;

    return _regenerator["default"].wrap(function _callee2$(_context3) {
      while (1) {
        switch (_context3.prev = _context3.next) {
          case 0:
            _i3 = result.ndim;

          case 1:
            if (!(_i3-- > 0)) {
              _context3.next = 6;
              break;
            }

            _context3.next = 4;
            return _i3;

          case 4:
            _context3.next = 1;
            break;

          case 6:
          case "end":
            return _context3.stop();
        }
      }
    }, _callee2);
  })()));
  return result;
}