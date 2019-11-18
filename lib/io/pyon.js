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
    var _iteratorNormalCompletion = true;
    var _didIteratorError = false;
    var _iteratorError = undefined;

    try {
      for (var _iterator = 'rue'[Symbol.iterator](), _step; !(_iteratorNormalCompletion = (_step = _iterator.next()).done); _iteratorNormalCompletion = true) {
        var c = _step.value;
        _char2 = next();
        if (c !== _char2) err(_char2, c);
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

    _char2 = skip();
    return true;
  }

  function False() {
    if ('F' !== _char2) throw new Error('Assertion failed.');
    var _iteratorNormalCompletion2 = true;
    var _didIteratorError2 = false;
    var _iteratorError2 = undefined;

    try {
      for (var _iterator2 = 'alse'[Symbol.iterator](), _step2; !(_iteratorNormalCompletion2 = (_step2 = _iterator2.next()).done); _iteratorNormalCompletion2 = true) {
        var c = _step2.value;
        _char2 = next();
        if (c !== _char2) err(_char2, c);
      }
    } catch (err) {
      _didIteratorError2 = true;
      _iteratorError2 = err;
    } finally {
      try {
        if (!_iteratorNormalCompletion2 && _iterator2["return"] != null) {
          _iterator2["return"]();
        }
      } finally {
        if (_didIteratorError2) {
          throw _iteratorError2;
        }
      }
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
    var _iteratorNormalCompletion3 = true;
    var _didIteratorError3 = false;
    var _iteratorError3 = undefined;

    try {
      for (var _iterator3 = 'one'[Symbol.iterator](), _step3; !(_iteratorNormalCompletion3 = (_step3 = _iterator3.next()).done); _iteratorNormalCompletion3 = true) {
        var c = _step3.value;
        _char2 = next();
        if (c !== _char2) err(_char2, c);
      }
    } catch (err) {
      _didIteratorError3 = true;
      _iteratorError3 = err;
    } finally {
      try {
        if (!_iteratorNormalCompletion3 && _iterator3["return"] != null) {
          _iterator3["return"]();
        }
      } finally {
        if (_didIteratorError3) {
          throw _iteratorError3;
        }
      }
    }

    _char2 = skip();
    return null;
  }

  var _char2 = skip();

  return any();
}