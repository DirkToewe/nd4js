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
exports.freudenstein_roth = void 0;

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _nd_array = require("../../nd_array");

// REFERENCES
// ----------
// .. [1] "Testing Unconstrained Optimization Software"
//         Jorge J. Moré, Burton S. Garbow, Kenneth E. Hillstrom
//         ACM Transactions on Mathematical Software, Vol 7, No. 1, March 1982, pp. 17-41
// .. [2] "Numerical solutions of systems of nonlinear equations"
//         F. Freudenstein, B. Roth,
//         ACM 10, 4 (Oct. 1963), pp. 550-556
var freudenstein_roth = function freudenstein_roth(x) {
  x = (0, _nd_array.asarray)(x);
  if (x.ndim < 1) throw new Error('freudenstein_roth(x): x.ndim must be at least 1.');
  if (x.shape[x.ndim - 1] !== 2) throw new Error('freudenstein_roth(x): x.shape[-1] must be 2.');
  var F_shape = x.shape.slice(0, -1),
      F = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length / 2);
  x = x.data;

  for (var x_off = x.length, F_off = F.length; F_off-- > 0;) {
    x_off -= 2;
    var x1 = x[x_off + 0],
        x2 = x[x_off + 1];
    var f1 = -13 + x1 + ((5 - x2) * x2 - 2) * x2,
        f2 = -29 + x1 + ((1 + x2) * x2 - 14) * x2;
    F[F_off] = f1 * f1 + f2 * f2;
  }

  return new _nd_array.NDArray(F_shape, F);
};

exports.freudenstein_roth = freudenstein_roth;
freudenstein_roth.nIn = freudenstein_roth.nOut = 2;
freudenstein_roth.roots = [[5, 4]];
freudenstein_roth.minima_global = [[5, 4]];
freudenstein_roth.minima = [[5, 4], [(53 - Math.sqrt(352)) / 3, (2 - Math.sqrt(22)) / 3]];

freudenstein_roth.grad = function (x) {
  x = (0, _nd_array.asarray)(x);
  if (x.ndim < 1) throw new Error('freudenstein_roth.grad(x): x.ndim must be at least 1.');
  if (x.shape[x.ndim - 1] !== 2) throw new Error('freudenstein_roth.grad(x): x.shape[-1] must be 2.');
  var G_shape = x.shape,
      G = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length);
  x = x.data;

  for (var off = x.length; (off -= 2) >= 0;) {
    var x1 = x[off + 0],
        x2 = x[off + 1]; //    const f1 = -13 + x1 + ((5-x2)*x2 - 2)*x2,
    //          f2 = -29 + x1 + ((1+x2)*x2 -14)*x2,
    //          d1 = (10-3*x2)*x2 - 2,
    //          d2 = ( 2+3*x2)*x2 -14
    //    G[off+0] = 2*(f1+f2);
    //    G[off+1] = 2*(f1*d1 + f2*d2);

    G[off + 0] = -84 + 4 * x1 + (-32 + 12 * x2) * x2;
    G[off + 1] = 864 - (32 - 24 * x2) * x1 + (24 + (-240 + (8 + (-40 + 12 * x2) * x2) * x2) * x2) * x2;
  }

  return new _nd_array.NDArray(G_shape, G);
};

freudenstein_roth.hess = function (x) {
  x = (0, _nd_array.asarray)(x);
  if (x.ndim < 1) throw new Error('freudenstein_roth.hess(x): x.ndim must be at least 1.');
  if (x.shape[x.ndim - 1] !== 2) throw new Error('freudenstein_roth.hess(x): x.shape[-1] must be 2.');
  var H_shape = Int32Array.of.apply(Int32Array, (0, _toConsumableArray2["default"])(x.shape).concat([2])),
      H = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length * 2);
  x = x.data;

  for (var H_off = H.length, x_off = x.length; (x_off -= 2) >= 0;) {
    H_off -= 4;
    var x1 = x[x_off + 0],
        x2 = x[x_off + 1];
    var f1 = -13 + x1 + ((5 - x2) * x2 - 2) * x2,
        f2 = -29 + x1 + ((1 + x2) * x2 - 14) * x2,
        d1 = (10 - 3 * x2) * x2 - 2,
        d2 = (2 + 3 * x2) * x2 - 14;
    H[H_off + 0] = 4;
    H[H_off + 1] = H[H_off + 2] = -32 + 24 * x2;
    H[H_off + 3] = f1 * (20 - 12 * x2) + 2 * d1 * d1 + f2 * (4 + 12 * x2) + 2 * d2 * d2;
  }

  return new _nd_array.NDArray(H_shape, H);
};

freudenstein_roth.lsq = function (x) {
  x = (0, _nd_array.asarray)(x);
  if (x.ndim < 1) throw new Error('freudenstein_roth.lsq(x): x.ndim must be at least 1.');
  if (x.shape[x.ndim - 1] !== 2) throw new Error('freudenstein_roth.lsq(x): x.shape[-1] must be 2.');
  var F_shape = x.shape.slice(),
      F = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length);
  x = x.data;

  for (var off = x.length; (off -= 2) >= 0;) {
    var x1 = x[off + 0],
        x2 = x[off + 1];
    F[off + 0] = -13 + x1 + ((5 - x2) * x2 - 2) * x2;
    F[off + 1] = -29 + x1 + ((1 + x2) * x2 - 14) * x2;
  }

  return new _nd_array.NDArray(F_shape, F);
};

freudenstein_roth.lsq_jac = function (x) {
  x = (0, _nd_array.asarray)(x);
  if (x.ndim < 1) throw new Error('freudenstein_roth.lsq_jac(x): x.ndim must be at least 1.');
  if (x.shape[x.ndim - 1] < 2) throw new Error('freudenstein_roth.lsq_jac(x): x.shape[-1] must be 2.');
  var N = x.shape[x.ndim - 1],
      M = (N - 1) * 2;
  var J_shape = Int32Array.of.apply(Int32Array, (0, _toConsumableArray2["default"])(x.shape).concat([2])),
      J = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length * M);
  x = x.data;

  for (var x_off = x.length, J_off = J.length; (J_off -= 4) >= 0;) {
    x_off -= 2;
    var x1 = x[x_off + 0],
        x2 = x[x_off + 1];
    J[J_off + 0] = 1;
    J[J_off + 1] = (10 - 3 * x2) * x2 - 2;
    J[J_off + 2] = 1;
    J[J_off + 3] = (2 + 3 * x2) * x2 - 14;
  }

  return new _nd_array.NDArray(J_shape, J);
};

Object.freeze(freudenstein_roth);