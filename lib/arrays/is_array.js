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
exports.is_array = void 0;

var _complex_array = require("../dt/complex_array");

var TYPESET = new Set([Array, Float32Array, Float64Array, Int8Array, Int16Array, Int32Array, Uint8Array, Uint16Array, Uint32Array, Uint8ClampedArray, _complex_array.Complex128Array].map(function (clazz) {
  return clazz.prototype;
}));

var is_array = function is_array(obj) {
  return TYPESET.has(Object.getPrototypeOf(obj));
};

exports.is_array = is_array;