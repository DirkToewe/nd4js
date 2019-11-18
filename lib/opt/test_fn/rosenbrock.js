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
exports.rosenbrock = rosenbrock;
exports.rosenbrock_grad = rosenbrock_grad;

var _nd_array = require("../../nd_array");

var _tabulate = require("../../tabulate");

function rosenbrock(x) {
  x = (0, _nd_array.asarray)(x);
  if (x.ndim < 1) throw new Error('rosenbrock(x): x.ndim must be at least 1.');
  if (x.shape[x.ndim - 1] < 2) throw new Error('rosenbrock(x): x.shape[-1] must be at least 2.');
  var N = x.shape[x.ndim - 1],
      dtype = x.dtype === 'float32' ? 'float32' : 'float64';
  return (0, _tabulate.tabulate)(x.shape.slice(0, -1), dtype, function () {
    var sum = 0;

    for (var _len = arguments.length, i = new Array(_len), _key = 0; _key < _len; _key++) {
      i[_key] = arguments[_key];
    }

    for (var j = N - 1; j-- > 0;) {
      sum += 100 * Math.pow(x.apply(void 0, i.concat([j + 1])) - Math.pow(x.apply(void 0, i.concat([j])), 2), 2) + Math.pow(1 - x.apply(void 0, i.concat([j])), 2);
    }

    return sum;
  });
}

function rosenbrock_grad(x) {
  x = (0, _nd_array.asarray)(x);
  if (x.ndim < 1) throw new Error('rosenbrock(x): x.ndim must be at least 1.');
  if (x.shape[x.ndim - 1] < 2) throw new Error('rosenbrock(x): x.shape[-1] must be at least 2.');
  var N = x.shape[x.ndim - 1],
      dtype = x.dtype === 'float32' ? 'float32' : 'float64';
  return x.mapElems(dtype, function (xj) {
    for (var _len2 = arguments.length, i = new Array(_len2 > 1 ? _len2 - 1 : 0), _key2 = 1; _key2 < _len2; _key2++) {
      i[_key2 - 1] = arguments[_key2];
    }

    var j = i.pop();
    var result = 0;
    if (j < N - 1) result += 400 * (Math.pow(xj, 2) - x.apply(void 0, i.concat([j + 1]))) * xj - 2 * (1 - xj);
    if (j > 0) result += 200 * (xj - Math.pow(x.apply(void 0, i.concat([j - 1])), 2));
    return result;
  });
}