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

import {asarray} from '../nd_array';
import {tabulate} from '../tabulate';

import {generic_test_min_gen        } from './_generic_test_min';
import {generic_test_min_gen_bounded} from './_generic_test_min_bounded';
import {min_lbfgsb_gen} from './lbfgsb';


const min_lbfgsb_gen_unbounded = (fg,x0) => {
  x0 = asarray('float64', x0);

  const bounds = tabulate(
    [x0.shape[0],2],
    'float64',
    (i,j) => [-Infinity, +Infinity][j]
  );

  return min_lbfgsb_gen(fg, x0, bounds);
}
Object.defineProperty(min_lbfgsb_gen_unbounded, 'name', {value: `${min_lbfgsb_gen.name} [unbounded]`, writable: false});


generic_test_min_gen        (min_lbfgsb_gen_unbounded);
generic_test_min_gen_bounded(min_lbfgsb_gen);
