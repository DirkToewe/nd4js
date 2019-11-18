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
var _exportNames = {
  IS_LITTLE_ENDIAN: true,
  WHITESPACES: true
};
exports.WHITESPACES = exports.IS_LITTLE_ENDIAN = void 0;

var _b = require("./b64");

Object.keys(_b).forEach(function (key) {
  if (key === "default" || key === "__esModule") return;
  if (Object.prototype.hasOwnProperty.call(_exportNames, key)) return;
  Object.defineProperty(exports, key, {
    enumerable: true,
    get: function get() {
      return _b[key];
    }
  });
});

var _istr = require("./istr");

Object.keys(_istr).forEach(function (key) {
  if (key === "default" || key === "__esModule") return;
  if (Object.prototype.hasOwnProperty.call(_exportNames, key)) return;
  Object.defineProperty(exports, key, {
    enumerable: true,
    get: function get() {
      return _istr[key];
    }
  });
});

var _npy = require("./npy");

Object.keys(_npy).forEach(function (key) {
  if (key === "default" || key === "__esModule") return;
  if (Object.prototype.hasOwnProperty.call(_exportNames, key)) return;
  Object.defineProperty(exports, key, {
    enumerable: true,
    get: function get() {
      return _npy[key];
    }
  });
});
var IS_LITTLE_ENDIAN = new Int16Array(Uint8Array.of(1, 0).buffer)[0] === 1;
exports.IS_LITTLE_ENDIAN = IS_LITTLE_ENDIAN;
var WHITESPACES = '\f\n\r\t\v ';
exports.WHITESPACES = WHITESPACES;