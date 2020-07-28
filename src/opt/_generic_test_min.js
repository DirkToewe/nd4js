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
 * along with ND.JS. If not, see <http://www.gnu.org/licenses/>.
 */

import {forEachItemIn, CUSTOM_MATCHERS} from '../jasmine_utils';
import {cartesian_prod,
        linspace,
        range} from "../iter";
import {NDArray} from '../nd_array';

import {norm} from '../la/norm'

import {KDTree} from "../spatial/kd_tree";

import {LineSearchError,
        LineSearchNoProgressError} from "./line_search/line_search_error";

import {beale}             from './test_fn/beale';
import {brown_badscale}    from './test_fn/brown_badscale';
import {freudenstein_roth} from './test_fn/freudenstein_roth';
import {helical_valley}    from "./test_fn/helical_valley";
import {JennrichSampson}   from './test_fn/jennrich_sampson';
import {powell_badscale}   from './test_fn/powell_badscale';
import {Rastrigin}         from './test_fn/rastrigin'
import {Rosenbrock}        from './test_fn/rosenbrock'
import { OptimizationNoProgressError } from './optimization_error';


export function generic_test_min_gen_with_test_fn( minimize_gen, test_fn, x_range )
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


  const kdTree = new KDTree(test_fn.minima);


  const closest_min = x => {
    const [nearest] = kdTree.nearest_gen(x.data);
    return nearest;
  };


  const test_body = x0 =>
  {
    Object.freeze(x0);

    let nCalls = 0
    const fg = x => {
      if( ++nCalls > 512*1024 )
        throw new Error('Too many iterations.');
      return [
        test_fn(x),
        test_fn.grad(x)
      ]
    }

    let x,f,g, nIter = -1
    try
    {
      for( [x,f,g] of minimize_gen(fg, x0) )
      {
        expect(x).toEqual( jasmine.any(NDArray) )
        expect(f).toEqual( jasmine.any(Number ) )
        expect(g).toEqual( jasmine.any(NDArray) )

        expect(x.ndim).toBe(1)
        expect(g.ndim).toBe(1)

        expect(x.shape).toEqual( Int32Array.of(x0.length) )
        expect(g.shape).toEqual( Int32Array.of(x0.length) )

        expect(f).toBeAllCloseTo(test_fn     (x), {rtol:0, atol:0})
        expect(g).toBeAllCloseTo(test_fn.grad(x), {rtol:0, atol:0})

        const gNorm = norm(g) / Math.sqrt(g.data.length);
        if(   gNorm <= 1e-8 )
          break
        expect(++nIter).toBeLessThan(128*1024)
      }
    }
    catch( err ) {
      if(    ! (err instanceof   LineSearchNoProgressError)
          && ! (err instanceof OptimizationNoProgressError) )
        throw err;
    }

    const TOL = {rtol: 1e-4, atol: 2e-4};

    expect(g).toBeAllCloseTo(0, TOL)
    expect(x).toBeAllCloseTo(closest_min(x), TOL); // <- TODO: generalize
  }


  forEachItemIn(
    function(){
      const N = Math.ceil( 2**(10/test_fn.nIn) );

      return cartesian_prod(
        ...x_range.map( r => linspace(...r,N) )
      )
    }()  
  ).it(`works with ${test_fn.name} given generated starting points`, test_body)


  forEachItemIn(
    function*(){
      for( let run=0; run++ < 1337; )
        yield x_range.map( ([lo,hi]) => {
          const s = Math.random();
          return lo*(1-s) + s*hi;
        })
    }()
  ).it(`works with ${test_fn.name} given random starting points`, test_body)
}


export function generic_test_min_gen( minimize_gen )
{
  describe(`${minimize_gen.name} [generic tests]`, () => {

    beforeEach( () => {
      jasmine.addMatchers(CUSTOM_MATCHERS)
    })

    generic_test_min_gen_with_test_fn( minimize_gen, beale, [[+0.1, +5.0],
                                                             [-0.5, +1.0]] );

    generic_test_min_gen_with_test_fn( minimize_gen, brown_badscale, [[1e+6 - 2e+5, 1e+6 + 2e+5],
                                                                      [2e-6 - 4e-7, 2e-6 + 4e-7]] ); // <- TODO: increase range of starting points once line search improved

    generic_test_min_gen_with_test_fn( minimize_gen, freudenstein_roth, [[ 4, Math.PI*4],
                                                                         [-Math.PI/2, 5]] );

    generic_test_min_gen_with_test_fn( minimize_gen, helical_valley, [[-Math.PI/2, +2],
                                                                      [-Math.PI/2, +2],
                                                                      [-Math.PI/2, +2]] );

    generic_test_min_gen_with_test_fn( minimize_gen, new JennrichSampson(10), [[0.15, 0.32],
                                                                               [0.15, 0.32]] );
/*
    generic_test_min_gen_with_test_fn( minimize_gen, powell_badscale, [[-9.6, +9.5], // <- avoids starting at x1=x2 which leads to a saddle point
                                                                       [-9.5, +9.6]] );
*/

    for( const length of range(2,4) )
      generic_test_min_gen_with_test_fn(
        minimize_gen, new Rosenbrock(length), Array.from({length}, () => [-Math.PI*3,+Math.PI*3])
      );
/*
    for( const length of range(1,4) )
      generic_test_min_gen_with_test_fn(
        minimize_gen, new Rastrigin(length), Array.from({length}, () => [-Math.PI*11,+Math.PI*11])
      );
*/
  });
}
