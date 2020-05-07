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

import {NDArray} from '../../nd_array'
import {cartesian_prod,
        linspace,
        range} from '../../iter'
import {forEachItemIn, CUSTOM_MATCHERS} from '../../jasmine_utils'
import {_rand_int} from '../../_test_data_generators'
import {zip_elems} from '../../zip_elems'

import {norm} from '../../la/norm'

import {beale}             from '../test_fn/beale';
import {brown_badscale}    from '../test_fn/brown_badscale';
import {freudenstein_roth} from '../test_fn/freudenstein_roth';
import {JennrichSampson}   from '../test_fn/jennrich_sampson';
import {powell_badscale}   from '../test_fn/powell_badscale';
import {Rastrigin}         from '../test_fn/rastrigin'
import {Rosenbrock}        from '../test_fn/rosenbrock'

import {LineSearchNoProgressError} from './line_search_error'


export function generic_test_line_search_with_test_fn( line_search, test_fn, x_range )
{
  if( x_range.length !== test_fn.nIn )
    throw new Error('Assertion failed.');
  for( const lo_hi of x_range )
  {
    if( lo_hi.length !== 2 )
      throw new Error('Assertion failed.');

    const [lo,hi] = lo_hi;
    if( ! (lo < hi))
      throw new Error('Assertion failed.');
  }


  forEachItemIn(
    function(){
      const N = Math.round( 2**(16/test_fn.nIn) );

      return cartesian_prod(
        ...x_range.map( r => linspace(...r,N) )
      )
    }()
  ).it(`works given generated ${test_fn.name} examples`, x0 => {
    x0 = Object.freeze(x0);

    let nCalls = 0
    const fg = x => {
      ++nCalls
      return [
        test_fn(x),
        test_fn.grad(x)
      ]
    };

    const [f0,g0] = fg(x0), negDir = g0.mapElems('float64', x => x*0.1),
       [c1,c2,c3] = [0.4,0.8,1.6]

    let x,f,g;
    try {
      ( [x,f,g] = line_search({c1,c2,c3})(fg)(x0,f0,g0, negDir) );
    }
    catch(err) {
      if( err instanceof LineSearchNoProgressError )
        return;
      throw err
    }

    expect(x).toEqual(jasmine.any(NDArray))
    expect(f).toEqual(jasmine.any(Number))
    expect(g).toEqual(jasmine.any(NDArray))

    expect( test_fn.grad(x) ).toBeAllCloseTo(g, {atol:0, gtol:0})

    expect(x.shape).toEqual( Int32Array.of(x0.length) )
    expect(g.shape).toEqual( Int32Array.of(x0.length) )

    expect(nCalls).toBeLessThan(64)

    const [p,p0] = [g,g0].map( g => {
      let pg = 0
      for( let [i]=negDir.shape; i-- > 0; )
        pg -= negDir.data[i] * g.data[i]
      return pg
    })

    expect( Math.abs(p) ).not.toBeGreaterThan( -c2*p0 )

    let α = norm( zip_elems([x,x0], (x,y) => x-y) ) / norm(negDir)
    expect( f - f0 ).not.toBeGreaterThan( c1*α*p0 )
  });
}


export function generic_test_line_search( line_search )
{
  describe(`${line_search.name} [generic tests]`, () => {
    beforeEach( () => {
      jasmine.addMatchers(CUSTOM_MATCHERS);
    });

    generic_test_line_search_with_test_fn( line_search, beale, [[-4.5,+4.5],
                                                                [-4.5,+4.5]])

    generic_test_line_search_with_test_fn( line_search, brown_badscale, [[1e+6 - 2e+5, 1e+6 + 2e+5],
                                                                         [2e-6 - 4e-7, 2e-6 + 4e-7]]);

    generic_test_line_search_with_test_fn( line_search, freudenstein_roth, [[-Math.PI*5, +Math.PI*5],
                                                                            [-Math.PI*5, +Math.PI*5]] )

    for( const m of range(2,11) )
    {
      const R = Math.log(1+m) / m;

      generic_test_line_search_with_test_fn( line_search, new JennrichSampson(m), [[-R, +R],
                                                                                   [-R, +R]] );
    }

    generic_test_line_search_with_test_fn( line_search, powell_badscale, [[1e-6, 9.5],
                                                                          [1e-6, 9.5]])

    for( const length of range(1,9) )
      generic_test_line_search_with_test_fn( line_search, new Rastrigin(length), Array.from({length}, () => [-Math.PI*11,+33]) )

    for( const length of range(2,9) )
      generic_test_line_search_with_test_fn( line_search, new Rosenbrock(length), Array.from({length}, () => [-Math.PI*3,+Math.PI*3]) )
  });
}
