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
exports.KDTree = void 0;

var _regenerator = _interopRequireDefault(require("@babel/runtime/regenerator"));

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _createClass2 = _interopRequireDefault(require("@babel/runtime/helpers/createClass"));

var _classCallCheck2 = _interopRequireDefault(require("@babel/runtime/helpers/classCallCheck"));

var _test_data_generators = require("../_test_data_generators");

var _is_array = require("../arrays/is_array");

var _norm = require("../la/norm");

var _nary_heap = require("./_nary_heap");

// TODO:
//   - Add NDArray support
var EUCLIDEAN_DISTANCE = function EUCLIDEAN_DISTANCE(x, y) {
  var len = x.length;
  if (len !== y.length) throw new Error('Assertion failed.');
  var norm = new _norm.FrobeniusNorm();

  for (var i = len; i-- > 0;) {
    norm.include(x[i] - y[i]);
  }

  return norm.result;
};

var Branch = function Branch(axis, threshold, c0, c1) {
  (0, _classCallCheck2["default"])(this, Branch);
  Object.assign(this, {
    axis: axis,
    threshold: threshold,
    c0: c0,
    c1: c1
  });
  Object.seal(this);
};

var Leaf = function Leaf(point) {
  (0, _classCallCheck2["default"])(this, Leaf);
  this.point = point;
  Object.seal(this);
};

var HeapedBranch = function HeapedBranch(distance, branch, nearest) {
  (0, _classCallCheck2["default"])(this, HeapedBranch);
  this.key = distance;
  this.nearest = nearest;
  Object.assign(this, branch);
  Object.seal(this);
};

var HeapedLeaf = function HeapedLeaf(distance, leaf) {
  (0, _classCallCheck2["default"])(this, HeapedLeaf);
  this.key = distance;
  this.point = leaf.point;
  Object.seal(this);
};

var KDTree = /*#__PURE__*/function () {
  function KDTree(points) {
    (0, _classCallCheck2["default"])(this, KDTree);
    points = (0, _toConsumableArray2["default"])(points);
    var ndim = points[0].length;
    if (points.some(function (pt) {
      return pt.length !== ndim;
    })) throw new Error('new KDTree(points): points must be iterable of (typed) arrays of same length (NDArrays are not yet supported).');
    if (points.some(function (pt) {
      return !(0, _is_array.is_array)(pt);
    })) throw new Error('new KDTree(points): points must be iterable of (typed) arrays of same length (NDArrays are not yet supported).');

    this.root = function build_tree(axis, start, stop) {
      //*DEBUG*/      if( !(start >= 0            ) ) throw new Error('Assertion failed.');
      //*DEBUG*/      if( !(start <  points.length) ) throw new Error('Assertion failed.');
      //*DEBUG*/      if( !(stop  >  start        ) ) throw new Error('Assertion failed.');
      //*DEBUG*/      if( !(stop  <= points.length) ) throw new Error('Assertion failed.');
      //*DEBUG*/      if( !(axis >= 0  ) ) throw new Error('Assertion failed.');
      //*DEBUG*/      if( !(axis < ndim) ) throw new Error('Assertion failed.');
      var swap = function swap(i, j) {
        if (!(i >= start)) throw new Error('Assertion failed.');
        if (!(j >= start)) throw new Error('Assertion failed.');
        if (!(i < stop)) throw new Error('Assertion failed.');
        if (!(j < stop)) throw new Error('Assertion failed.');
        var magnum_pi = points[i];
        points[i] = points[j];
        points[j] = magnum_pi;
      };

      if (1 === stop - start) return new Leaf(points[start]);

      for (;;) {
        // random split very similar to the quick select algorithm
        var threshold = points[(0, _test_data_generators._rand_int)(start, stop)][axis]; // partition using the threshold

        var l = start,
            r = start;

        for (var i = start; i < stop; i++) {
          var pi = points[i][axis];

          if (pi <= threshold) {
            swap(i, r);
            if (pi < threshold) swap(l++, r);
            r++;
          }
        }

        if (!(l >= start)) throw new Error('Assertion failed.');
        if (!(l <= r)) throw new Error('Assertion failed.');
        if (!(r <= stop)) throw new Error('Assertion failed.'); //*DEBUG*/        for( let i=start; i <  l   ; i++ ) if( ! (points[i][axis] < threshold) ) throw new Error('Assertion failed.');
        //*DEBUG*/        for( let i=r    ; i <  stop; i++ ) if( ! (points[i][axis] > threshold) ) throw new Error('Assertion failed.');
        //*DEBUG*/        for( let i=l    ; i <  r   ; i++ ) if( ! (points[i][axis]===threshold) ) throw new Error('Assertion failed.');

        var mid = start + stop >>> 1;
        mid = Math.max(mid, l);
        mid = Math.min(mid, r);
        var ax = (axis + 1) % ndim;
        return new Branch(axis, threshold, build_tree(ax, start, mid), build_tree(ax, mid, stop));
      }
    }(0, 0, points.length);

    Object.seal(this);
  }

  (0, _createClass2["default"])(KDTree, [{
    key: "nearest_gen",
    value: /*#__PURE__*/_regenerator["default"].mark(function nearest_gen(queryPoint) {
      var heap, enqueue, item, nearest, axis, threshold, c0, c1, n0, n1;
      return _regenerator["default"].wrap(function nearest_gen$(_context) {
        while (1) {
          switch (_context.prev = _context.next) {
            case 0:
              heap = new _nary_heap.NAryHeap();

              enqueue = function enqueue(nearest, node) {
                return heap.add(node instanceof Branch ? new HeapedBranch(EUCLIDEAN_DISTANCE(queryPoint, nearest), node, nearest) : new HeapedLeaf(EUCLIDEAN_DISTANCE(queryPoint, node.point), node));
              };

              enqueue(queryPoint, this.root);

            case 3:
              if (!(heap.size > 0)) {
                _context.next = 18;
                break;
              }

              item = heap.popMin();

              if (!(item instanceof HeapedBranch)) {
                _context.next = 14;
                break;
              }

              nearest = item.nearest, axis = item.axis, threshold = item.threshold, c0 = item.c0, c1 = item.c1;
              n0 = nearest, n1 = nearest;

              if (nearest[axis] > threshold) {
                n0 = n0.slice();
                n0[axis] = threshold;
              }

              if (nearest[axis] < threshold) {
                n1 = n1.slice();
                n1[axis] = threshold;
              }

              enqueue(n0, c0);
              enqueue(n1, c1);
              _context.next = 16;
              break;

            case 14:
              _context.next = 16;
              return item.point;

            case 16:
              _context.next = 3;
              break;

            case 18:
            case "end":
              return _context.stop();
          }
        }
      }, nearest_gen, this);
    })
  }]);
  return KDTree;
}();

exports.KDTree = KDTree;