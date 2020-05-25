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
exports.powell_badscale = void 0;

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _nd_array = require("../../nd_array");

var _root1d_bisect = require("../root1d_bisect");

// REFERENCES
// ----------
// .. [1] "Testing Unconstrained Optimization Software"
//         Jorge J. Moré, Burton S. Garbow, Kenneth E. Hillstrom
//         ACM Transactions on Mathematical Software, Vol 7, No. 1, March 1982, pp. 17-41
// .. [2] "A Hybrid Method for Nonlinear Equations"
//         M.J.D. Powell
//         In "Numerical Methods for Nonlinear Algebratic Equations", Gordon & Breach, New York, 1970, pp. 87-114
var powell_badscale = function powell_badscale(x) {
  x = (0, _nd_array.asarray)(x);
  if (x.ndim < 1) throw new Error('powell_badscale(x): x.ndim must be at least 1.');
  if (x.shape[x.ndim - 1] !== 2) throw new Error('powell_badscale(x): x.shape[-1] must be 2.');
  var F_shape = x.shape.slice(0, -1),
      F = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length / 2);
  x = x.data;

  for (var x_off = x.length, F_off = F.length; F_off-- > 0;) {
    x_off -= 2;
    var x1 = x[x_off + 0],
        x2 = x[x_off + 1];
    var f1 = 1e4 * x1 * x2 - 1,
        f2 = Math.exp(-x1) + Math.exp(-x2) - 1.0001;
    F[F_off] = f1 * f1 + f2 * f2;
  }

  return new _nd_array.NDArray(F_shape, F);
};

exports.powell_badscale = powell_badscale;
powell_badscale.nIn = powell_badscale.nOut = 2;

powell_badscale.minima_global = powell_badscale.roots = function () {
  var x1 = 0.0000109815932969981745569,
      x2 = 9.1061467398665240109;
  return [[x1, x2], [x2, x1]];
}();

powell_badscale.minima = function () {
  var x = (0, _root1d_bisect.root1d_bisect)( //    x => 1e4*(1e4*x*x - 1)*x - 2*Math.exp(-2*x) + 1.0001*Math.exp(-x),
  //    x => 1e4*(1e4*x*x - 1)*x - Math.exp(-x) * (2*Math.exp(-x) - 1.0001),
  function (x) {
    return 1e4 * (1e4 * x * x - 1) * x + Math.exp(-x) * (0.0001 - Math.expm1(-x) - Math.exp(-x));
  }, -0.02, -0.0005);
  return [].concat((0, _toConsumableArray2["default"])(powell_badscale.minima_global), [[x, x]]);
}();

powell_badscale.grad = function (x) {
  x = (0, _nd_array.asarray)(x);
  if (x.ndim < 1) throw new Error('powell_badscale.grad(x): x.ndim must be at least 1.');
  if (x.shape[x.ndim - 1] !== 2) throw new Error('powell_badscale.grad(x): x.shape[-1] must be 2.');
  var G_shape = x.shape,
      G = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length);
  x = x.data;

  for (var off = x.length; (off -= 2) >= 0;) {
    var x1 = x[off + 0],
        x2 = x[off + 1];
    var f1 = 1e4 * x1 * x2 - 1,
        f2 = Math.exp(-x1) + Math.exp(-x2) - 1.0001;
    G[off + 0] = 2 * (f1 * 1e4 * x2 - f2 * Math.exp(-x1));
    G[off + 1] = 2 * (f1 * 1e4 * x1 - f2 * Math.exp(-x2));
  }

  return new _nd_array.NDArray(G_shape, G);
};

powell_badscale.hess = function (x) {
  x = (0, _nd_array.asarray)(x);
  if (x.ndim < 1) throw new Error('powell_badscale.hess(x): x.ndim must be at least 1.');
  if (x.shape[x.ndim - 1] !== 2) throw new Error('powell_badscale.hess(x): x.shape[-1] must be 2.');
  var H_shape = Int32Array.of.apply(Int32Array, (0, _toConsumableArray2["default"])(x.shape).concat([2])),
      H = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length * 2);
  x = x.data;

  for (var H_off = H.length, x_off = x.length; (x_off -= 2) >= 0;) {
    H_off -= 4;
    var x1 = x[x_off + 0],
        x2 = x[x_off + 1];
    var f1 = 1e4 * x1 * x2 - 1,
        f2 = Math.exp(-x1) + Math.exp(-x2) - 1.0001;
    H[H_off + 0] = 2 * (1e8 * x2 * x2 + f2 * Math.exp(-x1) + Math.exp(-2 * x1));
    H[H_off + 1] = H[H_off + 2] = 2 * (f1 * 1e4 + 1e8 * x1 * x2 + Math.exp(-x1 - x2));
    H[H_off + 3] = 2 * (1e8 * x1 * x1 + f2 * Math.exp(-x2) + Math.exp(-2 * x2));
  }

  return new _nd_array.NDArray(H_shape, H);
};

powell_badscale.lsq = function (x) {
  x = (0, _nd_array.asarray)(x);
  if (x.ndim < 1) throw new Error('powell_badscale.lsq(x): x.ndim must be at least 1.');
  if (x.shape[x.ndim - 1] !== 2) throw new Error('powell_badscale.lsq(x): x.shape[-1] must be 2.');
  var F_shape = x.shape.slice(),
      F = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length);
  x = x.data;

  for (var off = x.length; (off -= 2) >= 0;) {
    var x1 = x[off + 0],
        x2 = x[off + 1];
    F[off + 0] = 1e4 * x1 * x2 - 1;
    F[off + 1] = Math.exp(-x1) + Math.exp(-x2) - 1.0001;
  }

  return new _nd_array.NDArray(F_shape, F);
};

powell_badscale.lsq_jac = function (x) {
  x = (0, _nd_array.asarray)(x);
  if (x.ndim < 1) throw new Error('powell_badscale.lsq_jac(x): x.ndim must be at least 1.');
  if (x.shape[x.ndim - 1] !== 2) throw new Error('powell_badscale.lsq_jac(x): x.shape[-1] must be 2.');
  var J_shape = Int32Array.of.apply(Int32Array, (0, _toConsumableArray2["default"])(x.shape).concat([2])),
      J = new (x.dtype === 'float32' ? Float32Array : Float64Array)(x.data.length * 2);
  x = x.data;

  for (var x_off = x.length, J_off = J.length; (J_off -= 4) >= 0;) {
    x_off -= 2;
    var x1 = x[x_off + 0],
        x2 = x[x_off + 1];
    J[J_off + 0] = 1e4 * x2;
    J[J_off + 1] = 1e4 * x1;
    J[J_off + 2] = -Math.exp(-x1);
    J[J_off + 3] = -Math.exp(-x2);
  }

  return new _nd_array.NDArray(J_shape, J);
};

Object.freeze(powell_badscale);