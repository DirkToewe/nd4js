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

var _interopRequireDefault = require("@babel/runtime/helpers/interopRequireDefault");

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.solve = solve;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _rrqr = require("./rrqr");

function solve(A, y) {
  var _rrqr_decomp = (0, _rrqr.rrqr_decomp)(A),
      _rrqr_decomp2 = (0, _slicedToArray2["default"])(_rrqr_decomp, 3),
      Q = _rrqr_decomp2[0],
      R = _rrqr_decomp2[1],
      P = _rrqr_decomp2[2];

  return (0, _rrqr.rrqr_solve)(Q, R, P, y);
}