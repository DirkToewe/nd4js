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

import {forEachItemIn, CUSTOM_MATCHERS} from '../jasmine_utils'
import {tabulate} from '../tabulate'
import {diag, diag_mat} from './diag'


describe('diag', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  forEachItemIn(
    function*(){
      function* shapes() {
        for( let i=0; i++ < 7; )
        for( let j=0; j++ < 7; ) { yield [i,j]
        for( let k=0; k++ < 7; ) { yield [i,j,k]
        for( let l=0; l++ < 7; ) { yield [i,j,k,l] }}}
      }

      for( const shape of shapes() )
      {
        const A = tabulate( shape, 'int32', (...idx) => idx.reduce((flat,s) => 10*flat + s+1, 0) ),
           [M,N]= shape.slice(-2)
        Object.freeze(A.data.buffer)
        for( let off=-M; ++off < N; )
          yield [A,off]
      }
    }()
  ).it('diag works given generated shapes', ([A,off]) => {
    const D = diag(A,off)

    for( const [idx,D_idx] of D.elems() )
    {
      const [k] = idx.slice(-1)
      expect(D_idx).toEqual(A(
        ...idx.slice(0,-1),
        Math.max(k,k-off),
        Math.max(k,k+off)
      ))
    }
  })

  forEachItemIn(
    function*(){
      function* shapes() {
        for( let i=1; i <= 8; i++ ) { yield [i]
        for( let j=1; j <= 8; j++ ) { yield [i,j]
        for( let k=1; k <= 8; k++ ) { yield [i,j,k] }}}
      }

      for( const shape of shapes() )
        yield tabulate( shape, 'int32', (...idx) => idx.reduce((flat,s) => 10*flat+s+1, 0) )
    }()
  ).it('diag_mat works given generated examples', D => {
    const A = diag_mat(D)
    expect(A).toBeDiagonal()
    for( const [idx,D_idx] of D.elems() )
    {
      const [k] = idx.slice(-1)
      expect( A(...idx,k) ).toEqual(D_idx)
    }
  })
})
