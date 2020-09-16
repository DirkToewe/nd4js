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

import {norm} from "../la/norm";

import {OptimizationNoProgressError} from "./optimization_error";


export function generic_test_tls_gen( tls_gen )
{
  describe(`${tls_gen.name} [generic tests]`, () => {
    beforeEach( () => {
      jasmine.addMatchers(CUSTOM_MATCHERS)
    })


    forEachItemIn(
      function*(rng){ 
        for( let run=0; run++ < 173; )
        { 
          // Let's test with the following minimization problem
          //
          //     min. f(a,b,x,y) = ( u / (1 + (a*x)²) )² +
          //                       ( v / (1 + (b*y)²) )² + x² + y²
          //
          // Which can be expressed as Total Least Squares (TLS) problem.
          const u = rng.uniform(-4,+4),
                v = rng.uniform(-4,+4),
                a0= rng.uniform(-2,+2),
                b0= rng.uniform(-2,+2),
                x0= rng.uniform(-1,+1) / 64,
                y0= rng.uniform(-1,+1) / 64;
          yield [u,v, a0,b0, x0,y0];
        }
      }
    ).it(`works given a hand-crafted example`, ([u,v, a0,b0, x0,y0]) => {

      const fgg = ([a,b],[x,y]) => [
        [ u / (1 + (a*x)**2),
          v / (1 + (b*y)**2) ],

        [[   u / (1 + (a*x)**2)**2 * -2 * a*x*x, 0],
         [0, v / (1 + (b*y)**2)**2 * -2 * b*y*y   ]],

        [ u / (1 + (a*x)**2)**2 * -2 * a*a*x,
          v / (1 + (b*y)**2)**2 * -2 * b*b*y]
      ];

      const loss = ([a,b],[x,y]) => {
        const z = u / (1 + (a*x)**2),
              w = v / (1 + (b*y)**2);
        return (x*x + y*y + z*z + w*w) / 4;
      };

      const dloss_dab = ([a,b],[x,y]) => [
        -u*u / (1 + (a*x)**2)**3 * a*x*x,
        -v*v / (1 + (b*y)**2)**3 * b*y*y
      ];

      const dloss_dxy = ([a,b],[x,y]) => [
        x/2 - u*u / (1 + (a*x)**2)**3 * a*a*x,
        y/2 - v*v / (1 + (b*y)**2)**3 * b*b*y
      ];

      let     ab,xy,  mse,dmse_dp,dmse_dx, dy, nIter=0;
      try {
        for( [ab,xy, mse,dmse_dp,dmse_dx, dy] of tls_gen(fgg, [a0,b0],[x0,y0]) )
        {
          if( !(++nIter <= 64*1024) )
            throw new Error('Assertion failed.');

          expect( mse   ).toBeAllCloseTo(  loss    (ab,xy) );
          expect(dmse_dp).toBeAllCloseTo( dloss_dab(ab,xy) );
          expect(dmse_dx).toBeAllCloseTo( dloss_dxy(ab,xy) );

          if(  norm(dmse_dp) <= 1e-8
            && norm(dmse_dx) <= 1e-8 ) break;
        }
      }
      catch(err) {
        if( !(err instanceof OptimizationNoProgressError) )
          throw err;
      }

      expect(dmse_dp).toBeAllCloseTo(0);
      expect(dmse_dx).toBeAllCloseTo(0);

      expect(nIter).toBeGreaterThan(0);
    });
  });
}
