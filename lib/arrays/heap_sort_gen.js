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
exports.heap_sort_gen = heap_sort_gen;

var _regenerator = _interopRequireDefault(require("@babel/runtime/regenerator"));

var _marked = /*#__PURE__*/_regenerator["default"].mark(heap_sort_gen);

function heap_sort_gen(items) {
  var compare_fn,
      i,
      len,
      less,
      swap,
      parent,
      child,
      siftDown,
      j,
      _args = arguments;
  return _regenerator["default"].wrap(function heap_sort_gen$(_context) {
    while (1) {
      switch (_context.prev = _context.next) {
        case 0:
          compare_fn = _args.length > 1 && _args[1] !== undefined ? _args[1] : function (x, y) {
            return (x > y) - (x < y);
          };

          if (!(items.length % 1 !== 0)) {
            _context.next = 3;
            break;
          }

          throw new Error('Assertion failed.');

        case 3:
          if (items.length >= 0) {
            _context.next = 5;
            break;
          }

          throw new Error('Assertion failed.');

        case 5:
          if (compare_fn instanceof Function) {
            _context.next = 7;
            break;
          }

          throw new Error('Assertion failed.');

        case 7:
          i = 0;
          len = items.length;
          less = function less(i, j) {
            return compare_fn(items[i], items[j]) < 0;
          }, swap = function swap(i, j) {
            var x = items[i];
            items[i] = items[j];
            items[j] = x;
          };
          parent = function parent(child) {
            return len - 1 - (len - 2 - child >> 1);
          }, child = function child(parent) {
            return len - 2 - (len - 1 - parent << 1);
          };

          siftDown = function siftDown(p) {
            var c = child(p);

            if (c > i) {
              if (!less(c, c - 1)) --c;

              if (less(c, p)) {
                swap(c, p);
                siftDown(c);
              }
            } else if (c === i && less(c, p)) swap(c, p);
          }; // HEAPIFY


          for (j = parent(0); j < len; j++) {
            siftDown(j);
          } // EXTRACT MINIMA FROM HEAP


        case 13:
          if (!(i < len)) {
            _context.next = 20;
            break;
          }

          swap(i, len - 1);
          _context.next = 17;
          return items[i++];

        case 17:
          siftDown(len - 1); // <- reinstate heap property

          _context.next = 13;
          break;

        case 20:
        case "end":
          return _context.stop();
      }
    }
  }, _marked);
}