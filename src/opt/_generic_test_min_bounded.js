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
import {NDArray, asarray} from '../nd_array';

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


export function generic_test_min_gen_bounded_with_test_fn( minimize_bounded_gen, test_fn, x_range )
{
  if( x_range.length !== test_fn.nIn )
    throw new Error('Assertion failed.');
  for( const lo_hi of x_range )
  {
    if( lo_hi.length !== 2 )
      throw new Error('Assertion failed.');

    const  [lo,hi] = lo_hi;
    if( ! (lo < hi))
      throw new Error('Assertion failed.');
  }


  function rand_float( min, max )
  {
    const s = Math.random();
    return min*(1-s) + s*max;
  }


  function rand_bounds( X )
  {
    const N = X.length;

    const bounds = x_range.map( ([x_min,x_max],i) => {
      const x = X[i];

      if( ! (x_min <= x_max) )
        throw new Error('Assertion failed.');

      const dx = x_max - x_min;

      x_min -= dx/2;
      x_max += dx/2;

           if( Math.random() < 0.05 ) x_min = -Infinity;
      else if( Math.random() < 0.05 ) x_min = x;
      else                            x_min = rand_float(x,x_min);

           if( Math.random() < 0.05 ) x_max = +Infinity;
      else if( Math.random() < 0.05 ) x_max = x;
      else                            x_max = rand_float(x,x_max);

      if( ! (x_min <= x) ) throw new Error('Assertion failed.');
      if( ! (x_max >= x) ) throw new Error('Assertion failed.');

      // return [-Infinity, +Infinity];
      return [x_min,x_max];
    });

    return Object.freeze(bounds);
  }


  const kdTree = new KDTree(test_fn.minima);


  const closest_min = x => {
    const [nearest] = kdTree.nearest_gen(x.data);
    return nearest;
  };


  const test_body = ([x0,bounds]) =>
  {
    const N = x0.length;

    let nCalls = 0
    const fg = x => {
      x = asarray('float64', x);
      expect(x.shape).toEqual( Int32Array.of(N) );

      x.data.forEach( (xi,i) => {
        const [x_min,x_max] = bounds[i];
        if( ! (xi <= x_max) ) throw new Error('Assertion failed.');
        if( ! (xi >= x_min) ) throw new Error('Assertion failed.');
      })

      if( ++nCalls > 64*1024 )
        throw new Error('Too many function calls.');
      return [
        test_fn(x),
        test_fn.grad(x)
      ]
    }

    let x,f,g,G, nIter = -1
    try
    {
      for( [x,f,g,G] of minimize_bounded_gen(fg, x0, bounds) )
      {
        expect(x).toEqual( jasmine.any(NDArray) )
        expect(f).toEqual( jasmine.any(Number ) )
        expect(g).toEqual( jasmine.any(NDArray) )
        expect(G).toEqual( jasmine.any(NDArray) )
        expect(f).not.toBeNaN();
        expect(f).not.toBePositiveInfinity();
        expect(f).not.toBeNegativeInfinity();

        expect(x.ndim).toBe(1)
        expect(g.ndim).toBe(1)
        expect(G.ndim).toBe(1)

        expect(x.shape).toEqual( Int32Array.of(x0.length) )
        expect(g.shape).toEqual( Int32Array.of(x0.length) )
        expect(G.shape).toEqual( Int32Array.of(x0.length) )

        expect(f).toBeAllCloseTo(test_fn     (x), {rtol:0, atol:0})
        expect(G).toBeAllCloseTo(test_fn.grad(x), {rtol:0, atol:0})

        for( let i=N; i-- > 0; )
        {
          const xi = x.data[i],
                gi = g.data[i],
                Gi = G.data[i];

          const [x_min,x_max] = bounds[i];

               if( 0 < Gi && xi===x_min ) expect(gi).toBe(0);
          else if( 0 > Gi && xi===x_max ) expect(gi).toBe(0);
          else                            expect(gi).toBe(Gi);
        }

        const gNorm = norm(g) / Math.sqrt(g.data.length);
        if(   gNorm <= 1e-8 )
          break
        expect(++nIter).toBeLessThan(16*1024)
      }
    }
    catch( err ) {
      if( ! (err instanceof LineSearchNoProgressError) )
        throw err;
    }
    // console.log({nIter, nCalls})

    const internalSolution = x.data.every( (xi,i) => {
      const [x_min,x_max] = bounds[i];
      return x_min < xi && xi < x_max;
    });

    const TOL = {rtol: 1e-4, atol: 1e-3};

    expect(g).toBeAllCloseTo(0, TOL);
    if( internalSolution )
      expect(x).toBeAllCloseTo(closest_min(x), TOL); // <- TODO: generalize
  }


  forEachItemIn(
    function*(){
      const N = Math.ceil( 2**(9.5/test_fn.nIn) );

      const seq = cartesian_prod(
        ...x_range.map( r => linspace(...r,N) )
      );

      for( const x of seq ) {
        Object.freeze(x);
        yield Object.freeze([x, rand_bounds(x)]);
      }
    }()  
  ).it(`works with ${test_fn.name} given generated starting points`, test_body)


  forEachItemIn(
    function*(){
      for( let run=0; run++ < 733; )
      {
        const x = x_range.map( ([lo,hi]) => {
          const s = Math.random();
          return lo*(1-s) + s*hi;
        })

        yield Object.freeze([x, rand_bounds(x)]);
      }
    }()
  ).it(`works with ${test_fn.name} given random starting points`, test_body)
}


export function generic_test_min_gen_bounded( minimize_bounded_gen )
{
  describe(`${minimize_bounded_gen.name} [bounded] [generic tests]`, () => {

    beforeEach( () => {
      jasmine.addMatchers(CUSTOM_MATCHERS)
    })

    generic_test_min_gen_bounded_with_test_fn( minimize_bounded_gen, beale, [[+0.1, +5.0],
                                                                             [-0.5, +1.0]] );
/*
    generic_test_min_gen_bounded_with_test_fn( minimize_bounded_gen, brown_badscale, [[1e+6 - 2e+5, 1e+6 + 2e+5],
                                                                                      [2e-6 - 4e-7, 2e-6 + 4e-7]] ); // <- TODO: increase range of starting points once line search improved
*/
    generic_test_min_gen_bounded_with_test_fn( minimize_bounded_gen, freudenstein_roth, [[ Math.PI  , +14],
                                                                                         [-Math.PI/2,  +6]] );

    generic_test_min_gen_bounded_with_test_fn( minimize_bounded_gen, helical_valley, [[-Math.PI/2, +2],
                                                                                      [-Math.PI/2, +2],
                                                                                      [-Math.PI/2, +1]] );

    generic_test_min_gen_bounded_with_test_fn( minimize_bounded_gen, new JennrichSampson(10), [[-1, 0.35],
                                                                                               [-1, 0.35]] );
/*
    generic_test_min_gen_bounded_with_test_fn( minimize_bounded_gen, powell_badscale, [[-10.1, +10.0], // <- avoids starting at x1=x2 which leads to a saddle point
                                                                                       [-10.0, +10.1]] );
*/
    for( const length of range(1,4) )
      generic_test_min_gen_bounded_with_test_fn(
        minimize_bounded_gen, new Rastrigin(length), Array.from({length}, () => [-Math.PI*11,+Math.PI*11])
      );

    for( const length of range(2,4) )
      generic_test_min_gen_bounded_with_test_fn(
        minimize_bounded_gen, new Rosenbrock(length), Array.from({length}, () => [-Math.PI*1,+Math.PI*2])
      );

  });
}
