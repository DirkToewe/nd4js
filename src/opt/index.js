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

// AGENDA:
//   - Quadratic Programming
//   - Orthogonal Distance Regression
//   - Derivative-Free Multivariate Minimizers

import * as line_search from './line_search'
import * as test_fn from './test_fn'

export {
  line_search,
  test_fn
}

export * from './dogleg'
export * from './fit_lin'
export * from './gss'
export * from './lbfgs'
export * from './lbfgsb'
export * from './line_search' // <- FIXME: this should not be here remove in next major version
export * from './lm'
export * from './num_grad'
export * from './optimization_error'
export * from './root1d_bisect'
export * from './root1d_brent'
export * from './root1d_illinois'
