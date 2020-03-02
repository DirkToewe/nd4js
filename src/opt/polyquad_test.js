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

import {roots1d_polyquad} from './polyquad'


describe('polyquad', () => {
  forEachItemIn(
    function*() {
      const N = 24;

      function* xRange() {
        for( let i=-N; i <= +N; i++ )
          yield Math.PI * i/42;

        for( let i=-N; i <= +N; i++ )
          yield Number.EPSILON*i;

        for( let run=N; run-- > 0; )
          yield Math.random()*8-4;
      }

      function* cRange() {
        for( let run=N; run-- > 0; )
          yield Math.random()*8-4;
      }

      for( const x1 of xRange() )
      for( const x2 of xRange() )
      for( const c  of cRange() )
        yield [c,x1,x2];
    }()
  ).it('roots1d_polyquad(a,b,c) satisfies element-wise mixed stability', ([c, x1, x2]) => {
    let a = +(x1*x2) * c,
        b = -(x1+x2) * c;

    const [y1,y2] = roots1d_polyquad(a,b,c),
              tol = 1e-7;

    if( x1 > x2 ) [x1,x2] = [x2,x1];

    expect(y1).not.toBeGreaterThan(y2);

    expect(y1).not.toBeNaN();
    expect(y2).not.toBeNaN();
    expect( Math.abs(y1-x1) ).not.toBeGreaterThan( Math.abs(x1)*tol );
    expect( Math.abs(y2-x2) ).not.toBeGreaterThan( Math.abs(x2)*tol );
  })
})
