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

import {generic_test_fit_gen    } from "./_generic_test_fit";
import {generic_test_lsq_gen    } from "./_generic_test_lsq";
import {generic_test_min_gen    } from "./_generic_test_min";
import {generic_test_fit_odr_gen} from "./_generic_test_fit_odr";
import {fit_dogleg_gen,
    fit_odr_dogleg_gen,
        lsq_dogleg_gen,
        min_dogleg_gen} from "./dogleg";

generic_test_fit_odr_gen(fit_odr_dogleg_gen);
generic_test_min_gen    (    min_dogleg_gen);
generic_test_lsq_gen    (    lsq_dogleg_gen);
generic_test_fit_gen    (    fit_dogleg_gen);
