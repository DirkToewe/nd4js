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

import {min_nelder_mead_gen} from "./nelder_mead";

import {generic_test_min_gradless_gen} from "./_generic_test_min_gradless";

const RESTARTS = 1;

function *nelder_mead_with_restarts( F, x )
{
  let f;
  for( let i=0; i++ <= RESTARTS; )
  {
    const iter = min_nelder_mead_gen(F,x);

    const {value: [nextX,nextF]} = iter.next();

    if( i===1 ||   nextF < f )
      yield [nextX,nextF];
    f = nextF;

    for( const [x,f] of iter )
      yield[x,f];
  }
}
Object.defineProperty(nelder_mead_with_restarts, 'name', {value: `${min_nelder_mead_gen.name} with ${RESTARTS} restart${RESTARTS===1 ? '' : 's'}`, writable: false});
Object.freeze(nelder_mead_with_restarts);

generic_test_min_gradless_gen(nelder_mead_with_restarts)
