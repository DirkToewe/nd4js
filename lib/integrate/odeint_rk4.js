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
exports.odeint_rk4 = odeint_rk4;

var _nd_array = require("../nd_array");

var _zip_elems = require("../zip_elems");

function odeint_rk4(dy, y0, t0, dt) {
  y0 = (0, _nd_array.asarray)(y0); // https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods

  var k1 = dy(y0, t0).mapElems(function (dy) {
    return dy * dt;
  }),
      k2 = dy((0, _zip_elems.zip_elems)([y0, k1], function (y0, k1) {
    return y0 + k1 / 2;
  }), t0 + dt / 2).mapElems(function (dy) {
    return dy * dt;
  }),
      k3 = dy((0, _zip_elems.zip_elems)([y0, k2], function (y0, k2) {
    return y0 + k2 / 2;
  }), t0 + dt / 2).mapElems(function (dy) {
    return dy * dt;
  }),
      k4 = dy((0, _zip_elems.zip_elems)([y0, k3], function (y0, k3) {
    return y0 + k3;
  }), t0 + dt).mapElems(function (dy) {
    return dy * dt;
  });
  return (0, _zip_elems.zip_elems)([y0, k1, k2, k3, k4], function (y0, k1, k2, k3, k4) {
    return y0 + (k1 + 2 * (k2 + k3) + k4) / 6;
  });
}