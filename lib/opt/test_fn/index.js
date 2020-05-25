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

var _beale = require("./beale");

Object.keys(_beale).forEach(function (key) {
  if (key === "default" || key === "__esModule") return;
  Object.defineProperty(exports, key, {
    enumerable: true,
    get: function get() {
      return _beale[key];
    }
  });
});

var _brown_badscale = require("./brown_badscale");

Object.keys(_brown_badscale).forEach(function (key) {
  if (key === "default" || key === "__esModule") return;
  Object.defineProperty(exports, key, {
    enumerable: true,
    get: function get() {
      return _brown_badscale[key];
    }
  });
});

var _freudenstein_roth = require("./freudenstein_roth");

Object.keys(_freudenstein_roth).forEach(function (key) {
  if (key === "default" || key === "__esModule") return;
  Object.defineProperty(exports, key, {
    enumerable: true,
    get: function get() {
      return _freudenstein_roth[key];
    }
  });
});

var _helical_valley = require("./helical_valley");

Object.keys(_helical_valley).forEach(function (key) {
  if (key === "default" || key === "__esModule") return;
  Object.defineProperty(exports, key, {
    enumerable: true,
    get: function get() {
      return _helical_valley[key];
    }
  });
});

var _jennrich_sampson = require("./jennrich_sampson");

Object.keys(_jennrich_sampson).forEach(function (key) {
  if (key === "default" || key === "__esModule") return;
  Object.defineProperty(exports, key, {
    enumerable: true,
    get: function get() {
      return _jennrich_sampson[key];
    }
  });
});

var _powell_badscale = require("./powell_badscale");

Object.keys(_powell_badscale).forEach(function (key) {
  if (key === "default" || key === "__esModule") return;
  Object.defineProperty(exports, key, {
    enumerable: true,
    get: function get() {
      return _powell_badscale[key];
    }
  });
});

var _rastrigin = require("./rastrigin");

Object.keys(_rastrigin).forEach(function (key) {
  if (key === "default" || key === "__esModule") return;
  Object.defineProperty(exports, key, {
    enumerable: true,
    get: function get() {
      return _rastrigin[key];
    }
  });
});

var _rosenbrock = require("./rosenbrock");

Object.keys(_rosenbrock).forEach(function (key) {
  if (key === "default" || key === "__esModule") return;
  Object.defineProperty(exports, key, {
    enumerable: true,
    get: function get() {
      return _rosenbrock[key];
    }
  });
});