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

import {cartesian_prod,
        linspace,
        range} from "../iter";
import {forEachItemIn, CUSTOM_MATCHERS} from '../jasmine_utils';
import {NDArray} from '../nd_array';

import {norm} from "../la/norm";

import {KDTree} from "../spatial/kd_tree";

import {beale}             from './test_fn/beale';
import {brown_badscale}    from './test_fn/brown_badscale';
import {freudenstein_roth} from './test_fn/freudenstein_roth';
import {helical_valley}    from "./test_fn/helical_valley";
import {JennrichSampson}   from './test_fn/jennrich_sampson';
import {powell_badscale}   from './test_fn/powell_badscale';
import {Rastrigin}         from './test_fn/rastrigin'
import {Rosenbrock}        from './test_fn/rosenbrock'

import {OptimizationNoProgressError} from "./optimization_error";


export function generic_test_lsq_gen_with_test_fn( lsq_gen, test_fn, x_range )
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
    const fJ = x => {
      expect(++nCalls).toBeLessThan(8192);
      const f = test_fn.lsq(x),
            J = test_fn.lsq_jac(x);
      return [f,J];
    };

    const N = test_fn.nIn,
          M = test_fn.nOut;

    let x, mse, mse_grad;
    try
    {
      for( [x, mse, mse_grad] of lsq_gen(fJ, x0) )
      {
        expect(x       ).toEqual( jasmine.any(NDArray) )
        expect(mse     ).toEqual( jasmine.any(NDArray) )
        expect(mse_grad).toEqual( jasmine.any(NDArray) )

        expect(x       .ndim).toBe(1)
        expect(mse     .ndim).toBe(0)
        expect(mse_grad.ndim).toBe(1)

        expect(x       .shape).toEqual( Int32Array.of(N) )
        expect(mse_grad.shape).toEqual( Int32Array.of(N) )

        expect(mse     ).toBeAllCloseTo( test_fn(x) / M )
        expect(mse_grad).toBeAllCloseTo( test_fn.grad(x).mapElems(x => x/M) )

        const gNorm = norm(mse_grad);
        if(   gNorm <= 1e-8 )
          break;
      }
    }
    catch(    onpe ) {
      if( ! ( onpe instanceof OptimizationNoProgressError ) )
        throw onpe;
    }

    const TOL = {rtol: 1e-3, atol: 1e-3};

    expect(mse_grad).toBeAllCloseTo(0, TOL)
    expect(x       ).toBeAllCloseTo(closest_min(x), TOL)
  }


  forEachItemIn(
    function(){
      const N = Math.round( 2**(13/test_fn.nIn) );

      return cartesian_prod(
        ...x_range.map( r => linspace(...r,N) )
      );
    }()
  ).it(`works with ${test_fn.name} given generated starting points`, test_body)


  forEachItemIn(
    function*(){
      for( let run=0; run++ < 3*1337; )
        yield x_range.map( ([lo,hi]) => {
          const s = Math.random();
          return lo*(1-s) + s*hi;
        })
    }()
  ).it(`works with ${test_fn.name} given random starting points`, test_body)
}


export function generic_test_lsq_gen( lsq_gen )
{
  describe(`${lsq_gen.name} [generic tests]`, () => {
    beforeEach( () => {
      jasmine.addMatchers(CUSTOM_MATCHERS)
    })

    generic_test_lsq_gen_with_test_fn( lsq_gen, beale, [[+0.1, +5.0],
                                                        [-0.5, +0.8]] );

    generic_test_lsq_gen_with_test_fn( lsq_gen, brown_badscale, [[-1e+6, +2e+6],
                                                                 [-1e+6, +1e+6]] );

    generic_test_lsq_gen_with_test_fn( lsq_gen, freudenstein_roth, [[-16, +16],
                                                                    [-16, +16]] );

    generic_test_lsq_gen_with_test_fn( lsq_gen, helical_valley, [[-Math.PI/2, +2],
                                                                 [-Math.PI/2, +2],
                                                                 [-Math.PI/2, +2]] );

    generic_test_lsq_gen_with_test_fn( lsq_gen, new JennrichSampson(10), [[0, 0.6],
                                                                          [0, 0.8]] );

    generic_test_lsq_gen_with_test_fn( lsq_gen, powell_badscale, [[-12.1, +12.0], // <- avoids starting at x1=x2 which leads to a saddle point
                                                                  [-12.0, +12.1]] );
/* TODO: get the following test to work
    for( const length of range(1,4) )
      generic_test_lsq_gen_with_test_fn(
        lsq_gen, new Rastrigin(length), Array.from({length}, () => [-Math.PI*11,+Math.PI*11])
      );
*/
    for( const length of range(2,4) )
      generic_test_lsq_gen_with_test_fn(
        lsq_gen, new Rosenbrock(length), Array.from({length}, () => [-Math.PI*3,+Math.PI*3])
      );
  })
}
