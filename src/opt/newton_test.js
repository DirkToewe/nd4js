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

import {CUSTOM_MATCHERS, forEachItemIn} from '../jasmine_utils'
import {array, NDArray} from '../nd_array'
import {zip_elems} from '../zip_elems'

import {norm} from '../la/norm'

import {root_newton_gen} from './newton'

import {rosenbrock, rosenbrock_grad, rosenbrock_hess} from './test_fn/rosenbrock'


describe('newton', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  forEachItemIn(
    function*(){

      for( let run=1024; run-- > 0; )
      {
        const X  = Math.random()*2 - 1,
              c  = Math.random()*2 + 1,
              x0 = Math.random()*4 - 2;

        const fH = ({data: [x,y]}) => [
           [ (x-X)*(x*x+c) ],
          [[ (x-X)*2*x + (x*x+c) ]]
        ].map( arr => {
          arr = array(arr);
          Object.freeze(arr.data.buffer);
          return arr;
        });

        yield [fH, X, x0];
      }

    }()
  ).it('root_newton_gen solves generated univariate example.', ([fH, X, x0]) => {
    let x,f,H, nIter = -1,
          F = Infinity;

    for( [x,f,H] of root_newton_gen(fH,[x0]) )
    {
      expect(x).toEqual( jasmine.any(NDArray) )
      expect(f).toEqual( jasmine.any(NDArray) )
      expect(H).toEqual( jasmine.any(NDArray) )

      expect(x.ndim).toBe(1)
      expect(f.ndim).toBe(1)
      expect(H.ndim).toBe(2)

      expect(x.shape).toEqual( Int32Array.of(1) )
      expect(f.shape).toEqual( Int32Array.of(1) )
      expect(H.shape).toEqual( Int32Array.of(1,1) )

      const df = norm( zip_elems([F,f], (x,y) => x-y) );
      F = f;

      if( df <= 1e-8 )
        break;
      expect(++nIter).toBeLessThan(32);
    }

    expect( x(0) ).toBeAllCloseTo(X);
    expect(f).toBeAllCloseTo(0);
  })


  forEachItemIn(
    function*(){

      for( let run=1024; run-- > 0; )
      {
        const X  = Math.random()*2 - 1,
              Y  = Math.random()*2 - 1,
              x0 = Math.random()*4 - 2,
              y0 = Math.random()*4 - 2;

        const fH = ({data: [x,y]}) => [
           [ x-X,
             x*y - X*Y ],
          [[1, 0],
           [y, x]]
        ].map( arr => {
          arr = array(arr);
          Object.freeze(arr.data.buffer);
          return arr;
        });

        yield [fH, X,Y, x0,y0];
      }

    }()
  ).it('root_newton_gen solves generated bivariate example.', ([fH, X,Y, x0,y0]) => {
    let xy,f,H, nIter = -1,
           F = Infinity;

    for( [xy,f,H] of root_newton_gen(fH,[x0,y0]) )
    {
      expect(xy).toEqual( jasmine.any(NDArray) )
      expect(f ).toEqual( jasmine.any(NDArray) )
      expect(H ).toEqual( jasmine.any(NDArray) )

      expect(xy.ndim).toBe(1)
      expect(f .ndim).toBe(1)
      expect(H .ndim).toBe(2)

      expect(xy.shape).toEqual( Int32Array.of(2) )
      expect(f .shape).toEqual( Int32Array.of(2) )
      expect(H .shape).toEqual( Int32Array.of(2,2) )

      const df = norm( zip_elems([F,f], (x,y) => x-y) );
      F = f;

      if( df <= 1e-8 )
        break;
      expect(++nIter).toBeLessThan(32);
    }

    expect( xy(0) ).toBeAllCloseTo(X);
    expect( xy(1) ).toBeAllCloseTo(Y);
    expect(f).toBeAllCloseTo(0);
  })
})
