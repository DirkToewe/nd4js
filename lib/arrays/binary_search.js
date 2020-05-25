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
exports.binary_rangesearch = binary_rangesearch;
exports.binary_search = binary_search;

function binary_rangesearch(from, until, compass_fn) {
  if (!(from <= until)) throw new Error('binary_rangesearch(from, until, compass_fn): from must not be greater than until.');
  until -= 1;

  while (from <= until) {
    var mid = from + until >>> 1,
        c = compass_fn(mid);
    if (c < 0) from = mid + 1;else if (c > 0) until = mid - 1;else return mid;
  }

  return ~from;
}

function binary_search(array, key) {
  var compare_fn = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : function (x, y) {
    return (x > y) - (x < y);
  };
  var from = 0,
      to = array.length - 1;

  while (from <= to) {
    var mid = from + to >>> 1,
        c = compare_fn(array[mid], key);
    if (c < 0) from = mid + 1;else if (c > 0) to = mid - 1;else return mid;
  }

  return ~from;
}