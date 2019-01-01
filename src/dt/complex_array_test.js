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

import {Complex} from './complex'
import {Complex128Array} from './complex_array'

const randInt = (from,until) => Math.floor( Math.random()*(until-from) ) + from;

describe('Complex128Array', () => {

  beforeEach( () => {
    jasmine.addCustomEqualityTester( (a,b) => {
      if( a instanceof Complex ) return a.equals(b);
      if( b instanceof Complex ) return b.equals(a);
    })
  })

  it('of() should work correctly', () => {
    for( let i=4; i-- > 0; )
    {
      const objArr = Array.from(
              { length: randInt(0,32) },
              () => new Complex(
                Math.random()*128 - 64,
                Math.random()*128 - 64
              )
            ),
            cmpArr = Complex128Array.of(...objArr);
  
      for( let i=objArr.length; i-- > 0; )
        expect(objArr[i]).toEqual(cmpArr[i]);
    }
  })
})
