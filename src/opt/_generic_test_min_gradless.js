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

import {nextUp} from "../dt/float64_utils";

import {KDTree} from "../spatial/kd_tree";

import {OptimizationNoProgressError} from "./optimization_error";

import {beale}             from './test_fn/beale';
import {brown_badscale}    from './test_fn/brown_badscale';
import {freudenstein_roth} from './test_fn/freudenstein_roth';
import {helical_valley}    from "./test_fn/helical_valley";
import {JennrichSampson}   from './test_fn/jennrich_sampson';
import {powell_badscale}   from './test_fn/powell_badscale';
import {Rastrigin}         from './test_fn/rastrigin'
import {Rosenbrock}        from './test_fn/rosenbrock'


export function generic_test_min_gradless_gen_with_test_fn( minimize_gen, test_fn, x_range )
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

    const MAX_CALLS = 1e8;

    let nCalls = 0

    const F = x => {
      if( ++nCalls > MAX_CALLS )
        throw new Error('Too many function calls.');
      return test_fn(x);
    };

    let f, x, nIter = -1,
              fPrev = nextUp( test_fn(x0) );

    try
    {
      for( [x,f] of minimize_gen(F, x0) )
      {
        // Jasmine create huge memory overhead so we have to use these hacky assertions instead... yikes
        if( !(x instanceof NDArray) ) throw new Error('Assertion failed.');
        if(  'number' !== typeof f  ) throw new Error('Assertion failed.');

        if( ! isFinite(f) ) throw new Error('Assertion failed.');

        if( f != test_fn(x) )
          throw new Error('Assertion failed.');

        if( ! (fPrev > f) ) throw new Error('Assertion failed.');
        else   fPrev = f;

        if( x.ndim     !== 1         ) throw new Error('Assertion failed.');
        if( x.shape[0] !== x0.length ) throw new Error('Assertion failed.');

        if( !(++nIter <= MAX_CALLS) ) throw new Error('Too many iterations.');
      }
    }
    catch(    err ) {
      if( ! ( err instanceof OptimizationNoProgressError ) )
        throw err;
    }

    const TOL = {rtol: 5e-3, atol: 3e-3};

    expect(x).toBeAllCloseTo(closest_min(x), TOL); // <- TODO: generalize
  }


  forEachItemIn(
    function(){
      const N = Math.round( 2**(5.5/test_fn.nIn) );
      if( !(N > 0) )
        throw new Error('Assertion failed.');

      return cartesian_prod(
        ...x_range.map( r => linspace(...r,N) )
      )
    }()  
  ).it(`works with ${test_fn.name} given generated x0`, test_body);


  forEachItemIn(
    function*(rng){
      for( let run=0; run++ < 42; )
        yield x_range.map( ([lo,hi]) => rng.uniform(lo,hi) );
    }
  ).it(`works with ${test_fn.name} given random x0`, test_body);
}


export function generic_test_min_gradless_gen( minimize_gen )
{
  describe(`${minimize_gen.name} [generic gradient-free tests]`, () => {

    beforeEach( () => {
      jasmine.addMatchers(CUSTOM_MATCHERS)
    })

    generic_test_min_gradless_gen_with_test_fn( minimize_gen, beale, [[+0.1, +5.0],
                                                                      [-0.5, +1.0]] );

    // generic_test_min_gradless_gen_with_test_fn( minimize_gen, brown_badscale, [[1e+6 - 2e+5, 1e+6 + 2e+5],
    //                                                                            [2e-6 - 4e-7, 2e-6 + 4e-7]] ); // <- TODO: increase range of starting points once line search improved

    generic_test_min_gradless_gen_with_test_fn( minimize_gen, freudenstein_roth, [[ 4, Math.PI*4],
                                                                                  [-Math.PI/2, 5]]);

    // generic_test_min_gradless_gen_with_test_fn( minimize_gen, helical_valley, [[-Math.PI/2, +2],
    //                                                                            [-Math.PI/2, +2],
    //                                                                            [-Math.PI/2, +2]] );

    generic_test_min_gradless_gen_with_test_fn( minimize_gen, new JennrichSampson(10), [[0.15, 0.32],
                                                                                        [0.15, 0.32]] );

    // generic_test_min_gradless_gen_with_test_fn( minimize_gen, powell_badscale, [[-10.1, +10.0], // <- avoids starting at x1=x2 which leads to a saddle point
    //                                                                             [-10.0, +10.1]] );

    for( const length of range(1,4) )
      generic_test_min_gradless_gen_with_test_fn(
        minimize_gen, new Rastrigin(length), Array.from({length}, () => [-Math.PI*11,+Math.PI*11])
      );

    for( const length of range(2,4) )
      generic_test_min_gradless_gen_with_test_fn(
        minimize_gen, new Rosenbrock(length), Array.from({length}, () => [-Math.PI*1,+Math.PI*2])
      );

  });
}
