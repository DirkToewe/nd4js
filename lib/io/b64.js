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
          warn = function _ref2() {
            console.warn("_b64_decode(b64_chars): Base64 sequence ended unexpectedly.");
          };

          next6 = function _ref() {
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
  var _iteratorNormalCompletion = true;
  var _didIteratorError = false;
  var _iteratorError = undefined;

  try {
    for (var _iterator = b64_encode_gen(bytes, options)[Symbol.iterator](), _step; !(_iteratorNormalCompletion = (_step = _iterator.next()).done); _iteratorNormalCompletion = true) {
      var c = _step.value;
      result += c;
    }
  } catch (err) {
    _didIteratorError = true;
    _iteratorError = err;
  } finally {
    try {
      if (!_iteratorNormalCompletion && _iterator["return"] != null) {
        _iterator["return"]();
      }
    } finally {
      if (_didIteratorError) {
        throw _iteratorError;
      }
    }
  }

  return result;
}

function b64_encode_gen(bytes) {
  var _ref3,
      _ref3$pad,
      pad,
      _ref3$linewidth,
      linewidth,
      rawBytes,
      sixPacks,
      i,
      _iteratorNormalCompletion3,
      _didIteratorError3,
      _iteratorError3,
      _iterator3,
      _step3,
      six,
      _args4 = arguments;

  return _regenerator["default"].wrap(function b64_encode_gen$(_context4) {
    while (1) {
      switch (_context4.prev = _context4.next) {
        case 0:
          _ref3 = _args4.length > 1 && _args4[1] !== undefined ? _args4[1] : {}, _ref3$pad = _ref3.pad, pad = _ref3$pad === void 0 ? true : _ref3$pad, _ref3$linewidth = _ref3.linewidth, linewidth = _ref3$linewidth === void 0 ? 128 : _ref3$linewidth;

          if (0 < linewidth) {
            _context4.next = 3;
            break;
          }

          throw new Error("_b64_encode(bytes,{linewidth}): Invalid linewidth: ".concat(linewidth, "."));

        case 3:
          rawBytes = bytes[Symbol.iterator]();
          bytes = /*#__PURE__*/_regenerator["default"].mark(function _callee() {
            var _iteratorNormalCompletion2, _didIteratorError2, _iteratorError2, _iterator2, _step2, _byte;

            return _regenerator["default"].wrap(function _callee$(_context2) {
              while (1) {
                switch (_context2.prev = _context2.next) {
                  case 0:
                    _iteratorNormalCompletion2 = true;
                    _didIteratorError2 = false;
                    _iteratorError2 = undefined;
                    _context2.prev = 3;
                    _iterator2 = rawBytes[Symbol.iterator]();

                  case 5:
                    if (_iteratorNormalCompletion2 = (_step2 = _iterator2.next()).done) {
                      _context2.next = 15;
                      break;
                    }

                    _byte = _step2.value;
                    _byte *= 1;

                    if (0 <= _byte && _byte < 256) {
                      _context2.next = 10;
                      break;
                    }

                    throw new Error("_b64_encode(bytes): Invalid token in byte sequence: ".concat(_byte));

                  case 10:
                    _context2.next = 12;
                    return _byte;

                  case 12:
                    _iteratorNormalCompletion2 = true;
                    _context2.next = 5;
                    break;

                  case 15:
                    _context2.next = 21;
                    break;

                  case 17:
                    _context2.prev = 17;
                    _context2.t0 = _context2["catch"](3);
                    _didIteratorError2 = true;
                    _iteratorError2 = _context2.t0;

                  case 21:
                    _context2.prev = 21;
                    _context2.prev = 22;

                    if (!_iteratorNormalCompletion2 && _iterator2["return"] != null) {
                      _iterator2["return"]();
                    }

                  case 24:
                    _context2.prev = 24;

                    if (!_didIteratorError2) {
                      _context2.next = 27;
                      break;
                    }

                    throw _iteratorError2;

                  case 27:
                    return _context2.finish(24);

                  case 28:
                    return _context2.finish(21);

                  case 29:
                  case "end":
                    return _context2.stop();
                }
              }
            }, _callee, null, [[3, 17, 21, 29], [22,, 24, 28]]);
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
          _iteratorNormalCompletion3 = true;
          _didIteratorError3 = false;
          _iteratorError3 = undefined;
          _context4.prev = 10;
          _iterator3 = sixPacks[Symbol.iterator]();

        case 12:
          if (_iteratorNormalCompletion3 = (_step3 = _iterator3.next()).done) {
            _context4.next = 26;
            break;
          }

          six = _step3.value;

          if (!('number' !== typeof six)) {
            _context4.next = 16;
            break;
          }

          throw new Error('Assertion failed.');

        case 16:
          if (0 <= six && six < 64) {
            _context4.next = 18;
            break;
          }

          throw new Error('Assertion failed.');

        case 18:
          _context4.next = 20;
          return BYTE_TO_CHAR[six];

        case 20:
          if (!(0 === ++i % linewidth)) {
            _context4.next = 23;
            break;
          }

          _context4.next = 23;
          return '\n';

        case 23:
          _iteratorNormalCompletion3 = true;
          _context4.next = 12;
          break;

        case 26:
          _context4.next = 32;
          break;

        case 28:
          _context4.prev = 28;
          _context4.t0 = _context4["catch"](10);
          _didIteratorError3 = true;
          _iteratorError3 = _context4.t0;

        case 32:
          _context4.prev = 32;
          _context4.prev = 33;

          if (!_iteratorNormalCompletion3 && _iterator3["return"] != null) {
            _iterator3["return"]();
          }

        case 35:
          _context4.prev = 35;

          if (!_didIteratorError3) {
            _context4.next = 38;
            break;
          }

          throw _iteratorError3;

        case 38:
          return _context4.finish(35);

        case 39:
          return _context4.finish(32);

        case 40:
          if (!pad) {
            _context4.next = 46;
            break;
          }

        case 41:
          if (!(0 !== i++ % 4)) {
            _context4.next = 46;
            break;
          }

          _context4.next = 44;
          return '=';

        case 44:
          _context4.next = 41;
          break;

        case 46:
        case "end":
          return _context4.stop();
      }
    }
  }, _marked2, null, [[10, 28, 32, 40], [33,, 35, 39]]);
}