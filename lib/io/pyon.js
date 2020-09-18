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
exports.pyon_parse = pyon_parse;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

function _createForOfIteratorHelper(o, allowArrayLike) { var it; if (typeof Symbol === "undefined" || o[Symbol.iterator] == null) { if (Array.isArray(o) || (it = _unsupportedIterableToArray(o)) || allowArrayLike && o && typeof o.length === "number") { if (it) o = it; var i = 0; var F = function F() {}; return { s: F, n: function n() { if (i >= o.length) return { done: true }; return { done: false, value: o[i++] }; }, e: function e(_e) { throw _e; }, f: F }; } throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method."); } var normalCompletion = true, didErr = false, err; return { s: function s() { it = o[Symbol.iterator](); }, n: function n() { var step = it.next(); normalCompletion = step.done; return step; }, e: function e(_e2) { didErr = true; err = _e2; }, f: function f() { try { if (!normalCompletion && it["return"] != null) it["return"](); } finally { if (didErr) throw err; } } }; }

function _unsupportedIterableToArray(o, minLen) { if (!o) return; if (typeof o === "string") return _arrayLikeToArray(o, minLen); var n = Object.prototype.toString.call(o).slice(8, -1); if (n === "Object" && o.constructor) n = o.constructor.name; if (n === "Map" || n === "Set") return Array.from(o); if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n)) return _arrayLikeToArray(o, minLen); }

function _arrayLikeToArray(arr, len) { if (len == null || len > arr.length) len = arr.length; for (var i = 0, arr2 = new Array(len); i < len; i++) { arr2[i] = arr[i]; } return arr2; }

function pyon_parse(char_seq) {
  // next: returns the next character
  // skip: returns the next non-whitespace character
  var _ref = function () {
    var end = false,
        counter = 0;

    var iter = char_seq[Symbol.iterator](),
        next = function next() {
      var _iter$next = iter.next(),
          value = _iter$next.value,
          done = _iter$next.done;

      ++counter;

      if (done) {
        if (end) throw new Error('pyon_parse: Character sequence ended unexpectedly.');
        end = true;
        return null;
      }

      return value;
    },
        skip = function skip() {
      for (var _char;;) {
        switch (_char = next()) {
          default:
            return _char;

          case '\f':
          case '\n':
          case '\r':
          case '\t':
          case '\v':
          case ' ':
        }
      }
    },
        err = function err(encountered, expected) {
      var prefix = expected == null ? 'Invalid character' : "Expected ".concat(expected, ", but ");
      throw new Error("pyon_parse: ".concat(prefix, " '").concat(encountered, "' encountered as ").concat(counter, "-th character."));
    };

    return [next, skip, err];
  }(),
      _ref2 = (0, _slicedToArray2["default"])(_ref, 3),
      next = _ref2[0],
      skip = _ref2[1],
      err = _ref2[2];

  function any() {
    switch (_char2) {
      default:
        err(_char2);

      case '{':
        return dict();

      case '"':
      case "'":
        return str();

      case '(':
      case '[':
        return list();

      case 'T':
        return True();

      case 'F':
        return False();

      case 'N':
        return None();

      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '.':
      case '8':
      case '+':
      case '9':
      case '-':
        return num();
    }
  }

  function num() {
    var num = '';

    loop: for (;; _char2 = next()) {
      switch (_char2) {
        case '\f':
        case '\n':
        case '\r':
        case '\t':
        case '\v':
        case ' ':
          _char2 = skip();

        default:
          break loop;

        case '0':
        case '1':
        case '2':
        case '3':
        case 'a':
        case 'A':
        case '4':
        case 'b':
        case 'B':
        case '5':
        case 'c':
        case 'C':
        case '6':
        case 'd':
        case 'D':
        case '7':
        case 'e':
        case 'E':
        case '.':
        case '8':
        case 'f':
        case 'F':
        case '+':
        case '9':
        case 'e':
        case 'E':
        case '-':
          num += _char2;
      }
    }

    num = num.trim() * 1;
    if (isNaN(num)) err(num);
    return num;
  }

  function dict() {
    var result = {};
    if ('{' !== _char2) throw new Error('Assertion failed.');

    for (;;) {
      if ('}' === (_char2 = skip())) {
        _char2 = skip();
        return result;
      }

      var key = any();
      if (':' !== _char2) err(_char2, end);
      _char2 = skip();
      result[key] = any();

      switch (_char2) {
        default:
          err(_char2, ",' or '}");

        case ',':
          continue;

        case '}':
          _char2 = skip();
          return result;
      }
    }
  }

  function list() {
    var end = function () {
      switch (_char2) {
        case '(':
          return ')';

        case '[':
          return ']';

        default:
          throw new Error('Assertion failed.');
      }
    }();

    var result = [];

    for (;;) {
      switch (_char2 = skip()) {
        case end:
          _char2 = skip();
          return result;

        default:
      }

      result.push(any());

      switch (_char2) {
        default:
          err(_char2, ",' or '".concat(end));

        case ',':
          continue;

        case end:
          _char2 = skip();
          return result;
      }
    }
  }

  function True() {
    if ('T' !== _char2) throw new Error('Assertion failed.');

    var _iterator = _createForOfIteratorHelper('rue'),
        _step;

    try {
      for (_iterator.s(); !(_step = _iterator.n()).done;) {
        var c = _step.value;
        _char2 = next();
        if (c !== _char2) err(_char2, c);
      }
    } catch (err) {
      _iterator.e(err);
    } finally {
      _iterator.f();
    }

    _char2 = skip();
    return true;
  }

  function False() {
    if ('F' !== _char2) throw new Error('Assertion failed.');

    var _iterator2 = _createForOfIteratorHelper('alse'),
        _step2;

    try {
      for (_iterator2.s(); !(_step2 = _iterator2.n()).done;) {
        var c = _step2.value;
        _char2 = next();
        if (c !== _char2) err(_char2, c);
      }
    } catch (err) {
      _iterator2.e(err);
    } finally {
      _iterator2.f();
    }

    _char2 = skip();
    return false;
  }

  function str() {
    var END = _char2;

    switch (END) {
      default:
        throw new Error('Assertion failed.');

      case "'":
      case '"':
    }

    var result = '';

    for (;;) {
      switch (_char2 = next()) {
        case END:
          _char2 = skip();
          return result;

        case '\\':
          switch (_char2 = next()) {
            case 'u':
              var codePoint = 0;

              for (var i = 0; i < 4; i++) {
                codePoint <<= 8;

                switch (_char2 = next()) {
                  default:
                    err(_char2, 'Hexadecimal digit');

                  case 'f':
                  case 'F':
                    ++codePoint;

                  case 'e':
                  case 'E':
                    ++codePoint;

                  case 'd':
                  case 'D':
                    ++codePoint;

                  case 'c':
                  case 'C':
                    ++codePoint;

                  case 'b':
                  case 'B':
                    ++codePoint;

                  case 'a':
                  case 'A':
                    ++codePoint;

                  case '9':
                    ++codePoint;

                  case '8':
                    ++codePoint;

                  case '7':
                    ++codePoint;

                  case '6':
                    ++codePoint;

                  case '5':
                    ++codePoint;

                  case '4':
                    ++codePoint;

                  case '3':
                    ++codePoint;

                  case '2':
                    ++codePoint;

                  case '1':
                    ++codePoint;

                  case '0':
                    ++codePoint;
                }
              }

              result += String.fromCodePoint(codePoint);
              continue;

            case 'b':
              result += '\b';
              continue;

            case 'f':
              result += '\f';
              continue;

            case 'n':
              result += '\n';
              continue;

            case 'r':
              result += '\r';
              continue;

            case 't':
              result += '\t';
              continue;

            case '"':
            case "'":
            case "\\": // FALL THROUGH

          }

        default:
          result += _char2;
      }
    }
  }

  function None() {
    if ('N' !== _char2) throw new Error('Assertion failed.');

    var _iterator3 = _createForOfIteratorHelper('one'),
        _step3;

    try {
      for (_iterator3.s(); !(_step3 = _iterator3.n()).done;) {
        var c = _step3.value;
        _char2 = next();
        if (c !== _char2) err(_char2, c);
      }
    } catch (err) {
      _iterator3.e(err);
    } finally {
      _iterator3.f();
    }

    _char2 = skip();
    return null;
  }

  var _char2 = skip();

  return any();
}