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

import {generic_test_line_search        } from './_generic_test_line_search'
import {generic_test_line_search_bounded} from './_generic_test_line_search_bounded'
import {strong_wolfe as __strong_wolfe} from './strong_wolfe'


const strong_wolfe = opt =>
{
  const                                     line_search = __strong_wolfe(opt),
                        strong_wolfe = f => line_search(f);
  Object.defineProperty(strong_wolfe, 'name', {value: `strong_wolfe(${JSON.stringify({...line_search})})`, writable: false});
  Object.assign(        strong_wolfe, line_search);
  Object.freeze(        strong_wolfe);
  return                strong_wolfe;
}


generic_test_line_search_bounded( strong_wolfe() );
generic_test_line_search_bounded( strong_wolfe({fRed: 0.2           }) );
generic_test_line_search_bounded( strong_wolfe({           gRed: 0.7}) );
generic_test_line_search_bounded( strong_wolfe({fRed: 0.4, gRed: 0.6}) );


generic_test_line_search( strong_wolfe() );
generic_test_line_search( strong_wolfe({fRed: 0.2           }) );
generic_test_line_search( strong_wolfe({           gRed: 0.7}) );
generic_test_line_search( strong_wolfe({fRed: 0.4, gRed: 0.6}) );
