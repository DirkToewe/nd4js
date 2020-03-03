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
 * along with ND.JS. If not, see <http://www.gnu.org/licenses/>.
 */

var _interopRequireDefault = require("@babel/runtime/helpers/interopRequireDefault");

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.rosenbrock = rosenbrock;
exports.rosenbrock_grad = rosenbrock_grad;
exports.rosenbrock_hess = rosenbrock_hess;
exports.rosenbrock_lsq = rosenbrock_lsq;
exports.rosenbrock_lsq_jac = rosenbrock_lsq_jac;

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _nd_array = require("../../nd_array");

var _tabulate = require("../../tabulate");

function rosenbrock(x) {
  // https://en.wikipedia.org/wiki/Rosenbrock_function
  x = (0, _nd_array.asarray)(x);
  if (x.ndim < 1) throw new Error('rosenbrock(x): x.ndim must be at least 1.');
  if (x.shape[x.ndim - 1] < 2) throw new Error('rosenbrock(x): x.shape[-1] must be at least 2.');
  var N = x.shape[x.ndim - 1];
  var F_shape = x.shape.slice(0, -1),
      F = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length / N);
  x = x.data;

  for (var x_off = x.length, F_off = F.length; F_off-- > 0;) {
    x_off -= N;

    for (var i = N - 1; i-- > 0;) {
      var xj = x[x_off + i + 1],
          xi = x[x_off + i],
          u = xj - xi * xi,
          v = 1 - xi;
      F[F_off] += 100 * u * u + v * v;
    }
  }

  return new _nd_array.NDArray(F_shape, F);
}

function rosenbrock_grad(x) {
  x = (0, _nd_array.asarray)(x);
  if (x.ndim < 1) throw new Error('rosenbrock(x): x.ndim must be at least 1.');
  if (x.shape[x.ndim - 1] < 2) throw new Error('rosenbrock(x): x.shape[-1] must be at least 2.');
  var N = x.shape[x.ndim - 1];
  var G_shape = x.shape,
      G = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length);
  x = x.data;

  for (var off = x.length; off >= 0;) {
    off -= N;

    for (var i = N; i-- > 0;) {
      var xh = x[off + i - 1],
          xi = x[off + i],
          xj = x[off + i + 1];
      if (i < N - 1) G[off + i] = 400 * (xi * xi - xj) * xi - 2 * (1 - xi);
      if (0 < i) G[off + i] += 200 * (xi - xh * xh);
    }
  }

  return new _nd_array.NDArray(G_shape, G);
}

function rosenbrock_hess(x) {
  x = (0, _nd_array.asarray)(x);
  if (x.ndim < 1) throw new Error('rosenbrock(x): x.ndim must be at least 1.');
  if (x.shape[x.ndim - 1] < 2) throw new Error('rosenbrock(x): x.shape[-1] must be at least 2.');
  var N = x.shape[x.ndim - 1];
  var H_shape = Int32Array.of.apply(Int32Array, (0, _toConsumableArray2["default"])(x.shape).concat([N])),
      H = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length * N);
  x = x.data;

  for (var H_off = H.length, x_off = x.length; x_off >= 0;) {
    x_off -= N;
    H_off -= N * N;

    for (var i = N; i-- > 0;) {
      if (0 < i) {
        H[H_off + N * i + i] = 200;
        H[H_off + N * (i - 1) + i] = H[H_off + N * i + (i - 1)] = -400 * x[x_off + i - 1];
      }

      if (i < N - 1) {
        var xi = x[x_off + i];
        H[H_off + N * i + i] -= 400 * x[x_off + i + 1] - 1200 * xi * xi - 2;
      }
    }
  }

  return new _nd_array.NDArray(H_shape, H);
}

function rosenbrock_lsq(x) {
  x = (0, _nd_array.asarray)(x);
  if (x.ndim < 1) throw new Error('rosenbrock_lsq(x): x.ndim must be at least 1.');
  if (x.shape[x.ndim - 1] < 2) throw new Error('rosenbrock_lsq(x): x.shape[-1] must be at least 2.');
  var N = x.shape[x.ndim - 1],
      M = (N - 1) * 2;
  var F_shape = x.shape.slice(),
      F = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length / N * M);
  F_shape[F_shape.length - 1] = M;
  x = x.data;

  for (var x_off = x.length, F_off = F.length; (F_off -= M) >= 0;) {
    x_off -= N;

    for (var i = N - 1; i-- > 0;) {
      var xj = x[x_off + i + 1],
          xi = x[x_off + i],
          u = xj - xi * xi,
          v = 1 - xi;
      F[F_off + 2 * i + 0] = u * 10;
      F[F_off + 2 * i + 1] = v;
    }
  }

  return new _nd_array.NDArray(F_shape, F);
}

function rosenbrock_lsq_jac(x) {
  x = (0, _nd_array.asarray)(x);
  if (x.ndim < 1) throw new Error('rosenbrock(x): x.ndim must be at least 1.');
  if (x.shape[x.ndim - 1] < 2) throw new Error('rosenbrock(x): x.shape[-1] must be at least 2.');
  var N = x.shape[x.ndim - 1],
      M = (N - 1) * 2;
  var J_shape = Int32Array.of.apply(Int32Array, (0, _toConsumableArray2["default"])(x.shape.slice(0, -1)).concat([M, N])),
      J = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length * M);
  x = x.data;

  for (var x_off = x.length, J_off = J.length; (J_off -= M * N) >= 0;) {
    x_off -= N;

    for (var i = N - 1; i-- > 0;) {
      var j = i + 1,
          xj = x[x_off + j],
          xi = x[x_off + i];
      J[J_off + N * (2 * i + 0) + i] += -20 * xi;
      J[J_off + N * (2 * i + 0) + j] += 10;
      J[J_off + N * (2 * i + 1) + i] -= 1;
    }
  }

  return new _nd_array.NDArray(J_shape, J);
}