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


export function shuffle(
  array,
  rand_int = (from,until) => {
    if( 0  !==  from %1 ) throw new Error('Assertion failed.');
    if( 0  !==  until%1 ) throw new Error('Assertion failed.');
    if( !(from < until) ) throw new Error('Assertion failed.');

    return Math.floor( Math.random() * (until-from) ) + from;
  },
  from, until
)
{
  if( array.length%1 !== 0 )
    throw new Error('Assertion failed.');

  if( null == from ) from = 0;
  if( null == until) until= array.length;

  if(0!==from %1        ) throw new Error('Assertion failed.');
  if(0!==        until%1) throw new Error('Assertion failed.');
  if( ! (from <= until) ) throw new Error('Assertion failed.');

  // https://en.wikipedia.org/wiki/Fisher-Yates_shuffle
  for( let i=from; i < until-1; i++ )
  { const j = rand_int(i,until),
         aj = array[j];
              array[j] = array[i];
                         array[i] = aj;
  }
}
