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

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.regular_simplex = regular_simplex;

var _nd_array = require("../nd_array");

var _dt = require("../dt");

function regular_simplex(dtype, N) {
  if (null == N) {
    N = dtype;
    dtype = 'float64';
  } else {
    if (!(dtype in _dt.ARRAY_TYPES)) throw new Error("regular_simplex(dtype,N): dtype '".concat(dtype, "' not supported."));
    if (dtype.startsWith('int')) throw new Error("regular_simplex(dtype,N): dtype must not be integral.");
  }

  if (N % 1 !== 0) throw new Error('regular_simplex(N): N must be a valid integer.');
  if (!(0 < N)) throw new Error('regular_simplex(N): N must be greater than 0.');
  N |= 0;
  var V = new _dt.ARRAY_TYPES[dtype]((N + 1) * N),
      shape = Int32Array.of(N + 1, N); // The idea of the algorithm is to "recursively" build an i-simplex from an (i-1)
  // simplex. To do that we find the center of the (i-1)-simplex and move from there
  // along a new coordinate axis (x[i-1]-axis).

  for (var i = 1; i < N + 1; i++) {
    for (var j = 0; j < i - 1; j++) {
      V[N * i + j] = V[N * (i - 1) + j];
    }

    V[N * i + (i - 2)] /= i; // <- there is (i-1) points at x[i-2]=0 and 1 point at x[i-2]=V[N*i+(i-2)], i.e. the midpoint is at 1/i along the x[i-2]-axis

    V[N * i + (i - 1)] = Math.sqrt((i + 1) / (2 * i)); // <- height of the i-simplex, see: https://math.stackexchange.com/questions/1697870/height-of-n-simplex
  }

  return new _nd_array.NDArray(shape, V);
}