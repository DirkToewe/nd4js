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

import {forEachItemIn} from '../jasmine_utils'
import {eye} from './eye'


describe('eye', () => {
  forEachItemIn(
    function*(){
      for( let i=1; i <= 8; i++ ) { yield [i]
      for( let j=1; j <= 8; j++ ) { yield [i,j]
      for( let k=1; k <= 8; k++ ) { yield [i,j,k]
      for( let l=1; l <= 8; l++ ) { yield [i,j,k,l] }}}}
    }()
  ).it('works on random shapes', shape => {
    const I = eye(...shape)

    if( shape.length === 1 )
      shape = [...shape, ...shape]
    shape = Int32Array.from(shape)

    expect(I.shape).toEqual(shape)

    for( const [idx,I_idx] of I.elems() ) {
      const [i,j] = idx.slice(-2)
      expect(I_idx).toBe( 1*(i===j) )
    }
  })
})
