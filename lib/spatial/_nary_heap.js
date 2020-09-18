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
exports.NAryHeap = void 0;

var _classCallCheck2 = _interopRequireDefault(require("@babel/runtime/helpers/classCallCheck"));

var _createClass2 = _interopRequireDefault(require("@babel/runtime/helpers/createClass"));

var ARITY = 8,
    child = function child(parent) {
  return 1 + parent * ARITY;
},
    parent = function parent(child) {
  return (child - 1) / ARITY | 0;
};

var NAryHeap = /*#__PURE__*/function () {
  function NAryHeap() {
    (0, _classCallCheck2["default"])(this, NAryHeap);
    this._heap = [];
    Object.seal(this);
  }

  (0, _createClass2["default"])(NAryHeap, [{
    key: "add",
    value: function add(item) {
      if (!('key' in item)) throw new Error('Assertion failed.');
      var heap = this._heap;
      heap.push(item); // SIFT UP

      var i = heap.length - 1;

      for (;;) {
        var p = parent(i);
        if (i === 0 || item.key >= heap[p].key) break;
        heap[i] = heap[p];
        i = p;
      }

      heap[i] = item;
    }
  }, {
    key: "popMin",
    value: function popMin() {
      var heap = this._heap,
          result = heap[0],
          filler = heap.pop(); // <- fill the root with the rightmost entry (then we sift it down again)
      // SIFT DOWN

      if (0 < heap.length) for (var i = 0;;) {
        var p = i;
        heap[i] = filler;
        var c = child(i); // find largest child that is larger than the filler

        var C = Math.min(c + ARITY, heap.length);

        for (; c < C; c++) {
          if (heap[c].key < heap[i].key) i = c;
        }

        if (i === p) break;
        heap[p] = heap[i]; // <- move smallest value to root/parent
      }
      return result;
    }
  }, {
    key: "size",
    get: function get() {
      return this._heap.length;
    }
  }, {
    key: "min",
    get: function get() {
      return this._heap[0];
    }
  }]);
  return NAryHeap;
}();

exports.NAryHeap = NAryHeap;