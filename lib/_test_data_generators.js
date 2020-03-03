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

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports._rand_int = _rand_int;
exports._rand_rows0 = _rand_rows0;
exports._rand_cols0 = _rand_cols0;
exports._rand_rankdef = _rand_rankdef;

var _tabulate = require("./tabulate");

var _matmul = require("./la/matmul");

var _rand_ortho = require("./la/rand_ortho");

function _rand_int(from, until) {
  if (0 !== from % 1) throw new Error('Assertion failed.');
  if (0 !== until % 1) throw new Error('Assertion failed.');
  if (!(from < until)) throw new Error('Assertion failed.');
  return Math.floor(Math.random() * (until - from)) + from;
} ///** Generates a series of shapes that are broadcast compatible except for their tails.
// */
//export function _rand_shapes( ...tails ) // <- TODO
//{
//  throw new Error('Not yet implemented.');
//}


function _rand_rows0() {
  for (var _len = arguments.length, shape = new Array(_len), _key = 0; _key < _len; _key++) {
    shape[_key] = arguments[_key];
  }

  if (!(2 <= shape.length)) throw new Error('Assertion failed.');
  var AA = (0, _tabulate.tabulate)(shape, 'float64', function () {
    return Math.random() * 8 - 4;
  });
  var N = shape.pop(),
      M = shape.pop();
  var A = AA.data;
  var rows = new Int32Array(M);

  for (var A_off = A.length; (A_off -= M * N) >= 0;) {
    for (var i = M; i-- > 0;) {
      rows[i] = i;
    }

    var rank = _rand_int(0, M);

    for (var m = M; m > rank;) {
      var k = _rand_int(0, m--),
          _i = rows[k];

      rows[k] = rows[m];

      for (var j = N; j-- > 0;) {
        A[A_off + N * _i + j] = 0;
      }
    }
  }

  Object.freeze(A.buffer);
  Object.freeze(AA);
  return AA;
}

function _rand_cols0() {
  for (var _len2 = arguments.length, shape = new Array(_len2), _key2 = 0; _key2 < _len2; _key2++) {
    shape[_key2] = arguments[_key2];
  }

  if (!(2 <= shape.length)) throw new Error('Assertion failed.');
  var AA = (0, _tabulate.tabulate)(shape, 'float64', function () {
    return Math.random() * 8 - 4;
  });
  var N = shape.pop(),
      M = shape.pop();
  var A = AA.data;
  var cols = new Int32Array(N);

  for (var A_off = A.length; (A_off -= M * N) >= 0;) {
    for (var i = N; i-- > 0;) {
      cols[i] = i;
    }

    var rank = _rand_int(0, N);

    for (var n = N; n > rank;) {
      var k = _rand_int(0, n--),
          j = cols[k];

      cols[k] = cols[n];

      for (var _i2 = M; _i2-- > 0;) {
        A[A_off + N * _i2 + j] = 0;
      }
    }
  }

  Object.freeze(A.buffer);
  Object.freeze(AA);
  return AA;
}

function _rand_rankdef() {
  for (var _len3 = arguments.length, shape = new Array(_len3), _key3 = 0; _key3 < _len3; _key3++) {
    shape[_key3] = arguments[_key3];
  }

  if (!(2 <= shape.length)) throw new Error('Assertion failed.');
  var N = shape.pop(),
      M = shape.pop(),
      L = Math.min(M, N); // use random, mostly rank-deficient SVD to generate test matrix A

  var ranks = (0, _tabulate.tabulate)(shape, 'int32', function () {
    return _rand_int(0, L + 1);
  }),
      // <- ranks
  U = _rand_ortho.rand_ortho.apply(void 0, ['float64'].concat(shape, [M, L])),
      V = _rand_ortho.rand_ortho.apply(void 0, ['float64'].concat(shape, [L, N])),
      S = (0, _tabulate.tabulate)([].concat(shape, [L, L]), 'float64', function () {
    for (var _len4 = arguments.length, idx = new Array(_len4), _key4 = 0; _key4 < _len4; _key4++) {
      idx[_key4] = arguments[_key4];
    }

    var j = idx.pop(),
        i = idx.pop(),
        rank = ranks.apply(void 0, idx);
    if (i !== j || rank <= i) return 0;
    return Math.random() * 8 + 0.1;
  }),
      A = (0, _matmul.matmul)(U, S, V); // add random scaling


  for (var S_off = S.data.length; (S_off -= L * L) >= 0;) {
    var scale = Math.random() * 1e9 + 1;

    for (var i = L; i-- > 0;) {
      S.data[S_off + L * i + i] *= scale;
    }
  }

  Object.freeze(ranks.data.buffer);
  Object.freeze(A.data.buffer);
  return [A, ranks];
}