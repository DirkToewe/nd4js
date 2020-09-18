'use strict';
/* This file is part of ND4JS.
 *
 * ND4JS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ND4JS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ND4JS. If not, see <http://www.gnu.org/licenses/>.
 */

var _interopRequireDefault = require("@babel/runtime/helpers/interopRequireDefault");

Object.defineProperty(exports, "__esModule", {
  value: true
});
var _exportNames = {
  linspace: true,
  range: true,
  cartesian_prod: true,
  enumerate: true,
  zip: true,
  repeat: true
};
exports.linspace = linspace;
exports.range = range;
exports.cartesian_prod = cartesian_prod;
exports.enumerate = enumerate;
exports.zip = zip;
exports.repeat = repeat;

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _regenerator = _interopRequireDefault(require("@babel/runtime/regenerator"));

var _min_max = require("./min_max");

Object.keys(_min_max).forEach(function (key) {
  if (key === "default" || key === "__esModule") return;
  if (Object.prototype.hasOwnProperty.call(_exportNames, key)) return;
  Object.defineProperty(exports, key, {
    enumerable: true,
    get: function get() {
      return _min_max[key];
    }
  });
});

function _createForOfIteratorHelper(o, allowArrayLike) { var it; if (typeof Symbol === "undefined" || o[Symbol.iterator] == null) { if (Array.isArray(o) || (it = _unsupportedIterableToArray(o)) || allowArrayLike && o && typeof o.length === "number") { if (it) o = it; var i = 0; var F = function F() {}; return { s: F, n: function n() { if (i >= o.length) return { done: true }; return { done: false, value: o[i++] }; }, e: function e(_e) { throw _e; }, f: F }; } throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method."); } var normalCompletion = true, didErr = false, err; return { s: function s() { it = o[Symbol.iterator](); }, n: function n() { var step = it.next(); normalCompletion = step.done; return step; }, e: function e(_e2) { didErr = true; err = _e2; }, f: function f() { try { if (!normalCompletion && it["return"] != null) it["return"](); } finally { if (didErr) throw err; } } }; }

function _unsupportedIterableToArray(o, minLen) { if (!o) return; if (typeof o === "string") return _arrayLikeToArray(o, minLen); var n = Object.prototype.toString.call(o).slice(8, -1); if (n === "Object" && o.constructor) n = o.constructor.name; if (n === "Map" || n === "Set") return Array.from(o); if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n)) return _arrayLikeToArray(o, minLen); }

function _arrayLikeToArray(arr, len) { if (len == null || len > arr.length) len = arr.length; for (var i = 0, arr2 = new Array(len); i < len; i++) { arr2[i] = arr[i]; } return arr2; }

var _marked = /*#__PURE__*/_regenerator["default"].mark(linspace),
    _marked2 = /*#__PURE__*/_regenerator["default"].mark(range),
    _marked4 = /*#__PURE__*/_regenerator["default"].mark(cartesian_prod),
    _marked5 = /*#__PURE__*/_regenerator["default"].mark(enumerate),
    _marked6 = /*#__PURE__*/_regenerator["default"].mark(zip),
    _marked7 = /*#__PURE__*/_regenerator["default"].mark(repeat);

function linspace(start, end, num) {
  var i, s;
  return _regenerator["default"].wrap(function linspace$(_context) {
    while (1) {
      switch (_context.prev = _context.next) {
        case 0:
          if (num % 1 === 0) {
            _context.next = 2;
            break;
          }

          throw new Error('linspace(start,end,num): num must be an integer greater than 1.');

        case 2:
          if (num > 1) {
            _context.next = 4;
            break;
          }

          throw new Error('linspace(start,end,num): num must be an integer.greater than 1.');

        case 4:
          if (isFinite(start)) {
            _context.next = 6;
            break;
          }

          throw new Error('linspace(start,end,num): start must be a finite number.');

        case 6:
          if (isFinite(end)) {
            _context.next = 8;
            break;
          }

          throw new Error('linspace(start,end,num): end must be a finite number.');

        case 8:
          i = 0;

        case 9:
          if (!(i < num)) {
            _context.next = 16;
            break;
          }

          s = i / (num - 1);
          _context.next = 13;
          return start * (1 - s) + s * end;

        case 13:
          i++;
          _context.next = 9;
          break;

        case 16:
        case "end":
          return _context.stop();
      }
    }
  }, _marked);
}

function range(start, stop) {
  var step,
      i,
      _i,
      _args2 = arguments;

  return _regenerator["default"].wrap(function range$(_context2) {
    while (1) {
      switch (_context2.prev = _context2.next) {
        case 0:
          step = _args2.length > 2 && _args2[2] !== undefined ? _args2[2] : 1;

          if (!(start % 1 !== 0)) {
            _context2.next = 3;
            break;
          }

          throw new Error('range(start, stop, step=1): start not a valid integer.');

        case 3:
          if (!(stop % 1 !== 0)) {
            _context2.next = 5;
            break;
          }

          throw new Error('range(start, stop, step=1): stop not a valid integer.');

        case 5:
          if (!(step % 1 !== 0)) {
            _context2.next = 7;
            break;
          }

          throw new Error('range(start, stop, step=1): step not a valid integer.');

        case 7:
          if (!(step === 0)) {
            _context2.next = 9;
            break;
          }

          throw new Error('range(start, stop, step=1): step must not be 0.');

        case 9:
          if (!(0 < step)) {
            _context2.next = 19;
            break;
          }

          i = start;

        case 11:
          if (!(i < stop)) {
            _context2.next = 17;
            break;
          }

          _context2.next = 14;
          return i;

        case 14:
          i += step;
          _context2.next = 11;
          break;

        case 17:
          _context2.next = 26;
          break;

        case 19:
          _i = start;

        case 20:
          if (!(_i > stop)) {
            _context2.next = 26;
            break;
          }

          _context2.next = 23;
          return _i;

        case 23:
          _i += step;
          _context2.next = 20;
          break;

        case 26:
        case "end":
          return _context2.stop();
      }
    }
  }, _marked2);
}

function cartesian_prod() {
  var _marked3,
      _len,
      seqs,
      _key,
      result,
      iter,
      _args4 = arguments;

  return _regenerator["default"].wrap(function cartesian_prod$(_context4) {
    while (1) {
      switch (_context4.prev = _context4.next) {
        case 0:
          iter = function _iter(i) {
            var _iterator, _step, x;

            return _regenerator["default"].wrap(function iter$(_context3) {
              while (1) {
                switch (_context3.prev = _context3.next) {
                  case 0:
                    if (!(i < seqs.length)) {
                      _context3.next = 20;
                      break;
                    }

                    _iterator = _createForOfIteratorHelper(seqs[i]);
                    _context3.prev = 2;

                    _iterator.s();

                  case 4:
                    if ((_step = _iterator.n()).done) {
                      _context3.next = 10;
                      break;
                    }

                    x = _step.value;
                    result[i] = x;
                    return _context3.delegateYield(iter(i + 1), "t0", 8);

                  case 8:
                    _context3.next = 4;
                    break;

                  case 10:
                    _context3.next = 15;
                    break;

                  case 12:
                    _context3.prev = 12;
                    _context3.t1 = _context3["catch"](2);

                    _iterator.e(_context3.t1);

                  case 15:
                    _context3.prev = 15;

                    _iterator.f();

                    return _context3.finish(15);

                  case 18:
                    _context3.next = 22;
                    break;

                  case 20:
                    _context3.next = 22;
                    return result.slice();

                  case 22:
                  case "end":
                    return _context3.stop();
                }
              }
            }, _marked3, null, [[2, 12, 15, 18]]);
          };

          _marked3 = /*#__PURE__*/_regenerator["default"].mark(iter);

          for (_len = _args4.length, seqs = new Array(_len), _key = 0; _key < _len; _key++) {
            seqs[_key] = _args4[_key];
          }

          seqs = seqs.map(function (seq) {
            return seq instanceof Array || seq instanceof Float32Array || seq instanceof Float64Array || seq instanceof Int8Array || seq instanceof Int16Array || seq instanceof Int32Array || seq instanceof Uint8Array || seq instanceof Uint16Array || seq instanceof Uint32Array ? seq.slice() : (0, _toConsumableArray2["default"])(seq);
          });
          result = new Array(seqs.length);
          return _context4.delegateYield(iter(0), "t0", 6);

        case 6:
        case "end":
          return _context4.stop();
      }
    }
  }, _marked4);
}

function enumerate(seq) {
  var i, _iterator2, _step2, x;

  return _regenerator["default"].wrap(function enumerate$(_context5) {
    while (1) {
      switch (_context5.prev = _context5.next) {
        case 0:
          i = 0;
          _iterator2 = _createForOfIteratorHelper(seq);
          _context5.prev = 2;

          _iterator2.s();

        case 4:
          if ((_step2 = _iterator2.n()).done) {
            _context5.next = 10;
            break;
          }

          x = _step2.value;
          _context5.next = 8;
          return [i++, x];

        case 8:
          _context5.next = 4;
          break;

        case 10:
          _context5.next = 15;
          break;

        case 12:
          _context5.prev = 12;
          _context5.t0 = _context5["catch"](2);

          _iterator2.e(_context5.t0);

        case 15:
          _context5.prev = 15;

          _iterator2.f();

          return _context5.finish(15);

        case 18:
        case "end":
          return _context5.stop();
      }
    }
  }, _marked5, null, [[2, 12, 15, 18]]);
}

function zip() {
  var _len2,
      seqs,
      _key2,
      N,
      i,
      next,
      _i2,
      item,
      _args6 = arguments;

  return _regenerator["default"].wrap(function zip$(_context6) {
    while (1) {
      switch (_context6.prev = _context6.next) {
        case 0:
          for (_len2 = _args6.length, seqs = new Array(_len2), _key2 = 0; _key2 < _len2; _key2++) {
            seqs[_key2] = _args6[_key2];
          }

          if (0 < seqs.length) {
            _context6.next = 3;
            break;
          }

          throw new Error('zip(...seqs): seqs.length must be at least 1.');

        case 3:
          N = seqs.length;

          for (i = N; i-- > 0;) {
            seqs[i] = seqs[i][Symbol.iterator]();
          }

        case 5:
          next = [];
          _i2 = 0;

        case 7:
          if (!(_i2 < N)) {
            _context6.next = 15;
            break;
          }

          item = seqs[_i2].next();

          if (!item.done) {
            _context6.next = 11;
            break;
          }

          return _context6.abrupt("return");

        case 11:
          next.push(item.value);

        case 12:
          _i2++;
          _context6.next = 7;
          break;

        case 15:
          _context6.next = 17;
          return next;

        case 17:
          _context6.next = 5;
          break;

        case 19:
        case "end":
          return _context6.stop();
      }
    }
  }, _marked6);
}

function repeat(n, seq) {
  var i;
  return _regenerator["default"].wrap(function repeat$(_context7) {
    while (1) {
      switch (_context7.prev = _context7.next) {
        case 0:
          if (undefined === seq) {
            seq = n;
            n = Infinity;
          }

          if (!(n !== Infinity)) {
            _context7.next = 6;
            break;
          }

          if (!(n % 1 !== 0)) {
            _context7.next = 4;
            break;
          }

          throw new Error('repeat(n,seq): n must be a non-negative int or Infinity.');

        case 4:
          if (n >= 0) {
            _context7.next = 6;
            break;
          }

          throw new Error('repeat(n,seq): n must be a non-negative int or Infinity.');

        case 6:
          seq = seq instanceof Array || seq instanceof Float32Array || seq instanceof Float64Array || seq instanceof Int8Array || seq instanceof Int16Array || seq instanceof Int32Array || seq instanceof Uint8Array || seq instanceof Uint16Array || seq instanceof Uint32Array ? seq.slice() : (0, _toConsumableArray2["default"])(seq);
          i = n;

        case 8:
          if (!(i-- > 0)) {
            _context7.next = 12;
            break;
          }

          return _context7.delegateYield(seq, "t0", 10);

        case 10:
          _context7.next = 8;
          break;

        case 12:
        case "end":
          return _context7.stop();
      }
    }
  }, _marked7);
}