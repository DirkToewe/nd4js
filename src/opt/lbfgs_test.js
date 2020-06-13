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

import {lsq_lbfgs_gen,
        min_lbfgs_gen} from './lbfgs';
import {generic_test_lsq_gen} from "./_generic_test_lsq";
import {generic_test_min_gen} from './_generic_test_min';


generic_test_lsq_gen(lsq_lbfgs_gen);
generic_test_min_gen(min_lbfgs_gen);
