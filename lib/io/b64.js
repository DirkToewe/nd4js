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
exports.b64_decode = b64_decode;
exports.b64_decode_gen = b64_decode_gen;
exports.b64_encode = b64_encode;
exports.b64_encode_gen = b64_encode_gen;

var _regenerator = _interopRequireDefault(require("@babel/runtime/regenerator"));

function _createForOfIteratorHelper(o, allowArrayLike) { var it; if (typeof Symbol === "undefined" || o[Symbol.iterator] == null) { if (Array.isArray(o) || (it = _unsupportedIterableToArray(o)) || allowArrayLike && o && typeof o.length === "number") { if (it) o = it; var i = 0; var F = function F() {}; return { s: F, n: function n() { if (i >= o.length) return { done: true }; return { done: false, value: o[i++] }; }, e: function e(_e) { throw _e; }, f: F }; } throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method."); } var normalCompletion = true, didErr = false, err; return { s: function s() { it = o[Symbol.iterator](); }, n: function n() { var step = it.next(); normalCompletion = step.done; return step; }, e: function e(_e2) { didErr = true; err = _e2; }, f: function f() { try { if (!normalCompletion && it["return"] != null) it["return"](); } finally { if (didErr) throw err; } } }; }

function _unsupportedIterableToArray(o, minLen) { if (!o) return; if (typeof o === "string") return _arrayLikeToArray(o, minLen); var n = Object.prototype.toString.call(o).slice(8, -1); if (n === "Object" && o.constructor) n = o.constructor.name; if (n === "Map" || n === "Set") return Array.from(o); if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n)) return _arrayLikeToArray(o, minLen); }

function _arrayLikeToArray(arr, len) { if (len == null || len > arr.length) len = arr.length; for (var i = 0, arr2 = new Array(len); i < len; i++) { arr2[i] = arr[i]; } return arr2; }

var _marked = /*#__PURE__*/_regenerator["default"].mark(b64_decode_gen),
    _marked2 = /*#__PURE__*/_regenerator["default"].mark(b64_encode_gen);

var BYTE_TO_CHAR = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/',
    WHITESPACES = '\f\n\r\t\v ',
    CHAR_TO_BYTE = Int32Array.from({
  length: 256
}, function () {
  return -2;
});
if (BYTE_TO_CHAR.length !== 64) throw new Error('Assertion failed.');

for (var i = BYTE_TO_CHAR.length; i-- > 0;) {
  CHAR_TO_BYTE[BYTE_TO_CHAR.charCodeAt(i)] = i;
}

for (var _i = WHITESPACES.length; _i-- > 0;) {
  CHAR_TO_BYTE[WHITESPACES.charCodeAt(_i)] = -1;
}

function b64_decode(b64_chars) {
  return Uint8Array.from(b64_decode_gen(b64_chars));
}

function b64_decode_gen(b64_chars) {
  var done, next6, warn, bit0to5, bit6to11, bit12to17, bit18to23;
  return _regenerator["default"].wrap(function b64_decode_gen$(_context) {
    while (1) {
      switch (_context.prev = _context.next) {
        case 0:
          warn = function _warn() {
            console.warn("_b64_decode(b64_chars): Base64 sequence ended unexpectedly.");
          };

          next6 = function _next() {
            while (true) {
              var value = void 0;

              var _b64_chars$next = b64_chars.next();

              value = _b64_chars$next.value;
              done = _b64_chars$next.done;
              done = done || '=' === value;
              if (done) return NaN;
              if ('string' !== typeof value || value.length !== 1) throw new Error("_b64_decode(b64_chars): Invalid token in character sequence: ".concat(value, "."));
              value = CHAR_TO_BYTE[value.charCodeAt(0)]; // skip whitespaces

              if (value === -1) continue;
              if (!(value >= 0)) // <- handles undefined
                throw new Error("_b64_decode(b64_chars): Invalid base64 character '".concat(value, "'."));
              return value;
            }
          };

          b64_chars = b64_chars[Symbol.iterator]();

        case 3:
          if (!true) {
            _context.next = 24;
            break;
          }

          bit0to5 = next6();

          if (!done) {
            _context.next = 7;
            break;
          }

          return _context.abrupt("return");

        case 7:
          bit6to11 = next6();

          if (!done) {
            _context.next = 10;
            break;
          }

          return _context.abrupt("return", warn());

        case 10:
          _context.next = 12;
          return bit0to5 << 2 | bit6to11 >> 4;

        case 12:
          bit12to17 = next6();

          if (!done) {
            _context.next = 15;
            break;
          }

          return _context.abrupt("return");

        case 15:
          _context.next = 17;
          return (bit6to11 & 15) << 4 | bit12to17 >> 2;

        case 17:
          bit18to23 = next6();

          if (!done) {
            _context.next = 20;
            break;
          }

          return _context.abrupt("return");

        case 20:
          _context.next = 22;
          return (bit12to17 & 3) << 6 | bit18to23 & 63;

        case 22:
          _context.next = 3;
          break;

        case 24:
        case "end":
          return _context.stop();
      }
    }
  }, _marked);
}

function b64_encode(bytes) {
  var options = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : {};
  var result = '';

  var _iterator = _createForOfIteratorHelper(b64_encode_gen(bytes, options)),
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

function b64_encode_gen(bytes) {
  var _ref,
      _ref$pad,
      pad,
      _ref$linewidth,
      linewidth,
      rawBytes,
      sixPacks,
      i,
      _iterator3,
      _step3,
      six,
      _args4 = arguments;

  return _regenerator["default"].wrap(function b64_encode_gen$(_context4) {
    while (1) {
      switch (_context4.prev = _context4.next) {
        case 0:
          _ref = _args4.length > 1 && _args4[1] !== undefined ? _args4[1] : {}, _ref$pad = _ref.pad, pad = _ref$pad === void 0 ? true : _ref$pad, _ref$linewidth = _ref.linewidth, linewidth = _ref$linewidth === void 0 ? 128 : _ref$linewidth;

          if (0 < linewidth) {
            _context4.next = 3;
            break;
          }

          throw new Error("_b64_encode(bytes,{linewidth}): Invalid linewidth: ".concat(linewidth, "."));

        case 3:
          rawBytes = bytes[Symbol.iterator]();
          bytes = /*#__PURE__*/_regenerator["default"].mark(function _callee() {
            var _iterator2, _step2, _byte;

            return _regenerator["default"].wrap(function _callee$(_context2) {
              while (1) {
                switch (_context2.prev = _context2.next) {
                  case 0:
                    _iterator2 = _createForOfIteratorHelper(rawBytes);
                    _context2.prev = 1;

                    _iterator2.s();

                  case 3:
                    if ((_step2 = _iterator2.n()).done) {
                      _context2.next = 12;
                      break;
                    }

                    _byte = _step2.value;
                    _byte *= 1;

                    if (0 <= _byte && _byte < 256) {
                      _context2.next = 8;
                      break;
                    }

                    throw new Error("_b64_encode(bytes): Invalid token in byte sequence: ".concat(_byte));

                  case 8:
                    _context2.next = 10;
                    return _byte;

                  case 10:
                    _context2.next = 3;
                    break;

                  case 12:
                    _context2.next = 17;
                    break;

                  case 14:
                    _context2.prev = 14;
                    _context2.t0 = _context2["catch"](1);

                    _iterator2.e(_context2.t0);

                  case 17:
                    _context2.prev = 17;

                    _iterator2.f();

                    return _context2.finish(17);

                  case 20:
                  case "end":
                    return _context2.stop();
                }
              }
            }, _callee, null, [[1, 14, 17, 20]]);
          })();
          sixPacks = /*#__PURE__*/_regenerator["default"].mark(function _callee2() {
            var value, done, six, _bytes$next, _bytes$next2, _bytes$next3;

            return _regenerator["default"].wrap(function _callee2$(_context3) {
              while (1) {
                switch (_context3.prev = _context3.next) {
                  case 0:
                    if (!true) {
                      _context3.next = 33;
                      break;
                    }

                    six = void 0;
                    _bytes$next = bytes.next();
                    value = _bytes$next.value;
                    done = _bytes$next.done;

                    if (!done) {
                      _context3.next = 7;
                      break;
                    }

                    return _context3.abrupt("return");

                  case 7:
                    _context3.next = 9;
                    return value >> 2;

                  case 9:
                    six = (value & 3) << 4;
                    _bytes$next2 = bytes.next();
                    value = _bytes$next2.value;
                    done = _bytes$next2.done;

                    if (!done) {
                      _context3.next = 17;
                      break;
                    }

                    _context3.next = 16;
                    return six;

                  case 16:
                    return _context3.abrupt("return");

                  case 17:
                    _context3.next = 19;
                    return value >> 4 | six;

                  case 19:
                    six = (value & 15) << 2;
                    _bytes$next3 = bytes.next();
                    value = _bytes$next3.value;
                    done = _bytes$next3.done;

                    if (!done) {
                      _context3.next = 27;
                      break;
                    }

                    _context3.next = 26;
                    return six;

                  case 26:
                    return _context3.abrupt("return");

                  case 27:
                    _context3.next = 29;
                    return value >> 6 | six;

                  case 29:
                    _context3.next = 31;
                    return value & 63;

                  case 31:
                    _context3.next = 0;
                    break;

                  case 33:
                  case "end":
                    return _context3.stop();
                }
              }
            }, _callee2);
          })();
          i = 0;
          _iterator3 = _createForOfIteratorHelper(sixPacks);
          _context4.prev = 8;

          _iterator3.s();

        case 10:
          if ((_step3 = _iterator3.n()).done) {
            _context4.next = 23;
            break;
          }

          six = _step3.value;

          if (!('number' !== typeof six)) {
            _context4.next = 14;
            break;
          }

          throw new Error('Assertion failed.');

        case 14:
          if (0 <= six && six < 64) {
            _context4.next = 16;
            break;
          }

          throw new Error('Assertion failed.');

        case 16:
          _context4.next = 18;
          return BYTE_TO_CHAR[six];

        case 18:
          if (!(0 === ++i % linewidth)) {
            _context4.next = 21;
            break;
          }

          _context4.next = 21;
          return '\n';

        case 21:
          _context4.next = 10;
          break;

        case 23:
          _context4.next = 28;
          break;

        case 25:
          _context4.prev = 25;
          _context4.t0 = _context4["catch"](8);

          _iterator3.e(_context4.t0);

        case 28:
          _context4.prev = 28;

          _iterator3.f();

          return _context4.finish(28);

        case 31:
          if (!pad) {
            _context4.next = 37;
            break;
          }

        case 32:
          if (!(0 !== i++ % 4)) {
            _context4.next = 37;
            break;
          }

          _context4.next = 35;
          return '=';

        case 35:
          _context4.next = 32;
          break;

        case 37:
        case "end":
          return _context4.stop();
      }
    }
  }, _marked2, null, [[8, 25, 28, 31]]);
}