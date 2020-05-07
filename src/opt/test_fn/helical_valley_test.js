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

import {generic_test_test_fn} from './_generic_test_test_fn'
import {helical_valley} from './helical_valley'


const E = Number.EPSILON**0.25; // <- margin to avoid bad numerical gradients


const test_ranges = {
  'x1 > 0': [ [+E,+4], [-4,+4], [-4,+4] ],
  'x2 > 0': [ [-4,+4], [+E,+4], [-4,+4] ],
  'x1 < 0': [ [-4,-E], [-4,+4], [-4,+4] ]
};


for( const [test,range] of Object.entries(test_ranges) )
{
  const fn = x => helical_valley(x);
  Object.assign(fn, helical_valley);
  Object.defineProperty(fn, 'name', {value: `${helical_valley.name} (${test})`, writable: false})
  generic_test_test_fn(fn, range);
}
