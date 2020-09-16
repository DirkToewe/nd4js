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

import {forEachItemIn, CUSTOM_MATCHERS} from '../jasmine_utils'
import {tabulate} from "../tabulate";

import {matmul2} from "./matmul";


export function generic_test_solve( solve )
{
  describe(`${solve.name} [generic SOLVE tests]`, () => {

    beforeEach( () => {
      jasmine.addMatchers(CUSTOM_MATCHERS)
    });


    forEachItemIn(
      function*(rng){
        for( let run=512; run-- > 0; )
        {
          let ndim = rng.int(0,3),
            shapes = [ Array.from({length: ndim}, () => rng.int(1,4)) ]
            shapes.splice( rng.int(0,2), 0, shapes[0].slice( rng.int(0,ndim+1) ) )

          // add random broadcasts
          for( let d=ndim; d > 0; d-- )
          for( let i=rng.int(0,2); i-- > 0; ) {
            const     shape = shapes[rng.int(0,2)],
                  j = shape.length - d
            if(0<=j)  shape[j] = 1
          }

          const M = rng.int(1,24); shapes[0].push(M,M)
          const N = rng.int(1,24); shapes[1].push(M,N)

          yield shapes.map(
            s => tabulate(s, 'float64', () => rng.uniform(-4,+4))
          )
        }
      }
    ).it('solves random square examples', ([A,y]) => {
      const x = solve(A,y);
      expect( matmul2(A,x) ).toBeAllCloseTo(y);
    });


    forEachItemIn(
      function*(rng){
        function* sizes() {
          const steps_per_binade = 3;

          for( let N=0; N++ < 16; )
            yield  N;

          for( let run=4*steps_per_binade; run++ < 7*steps_per_binade; )
            yield Math.round( 2 ** (run/steps_per_binade) );
        }

        for( const M of sizes() )
        {    const N = rng.int(1,8);
          yield [
            tabulate([M,M], 'float64', () => rng.uniform(-4,+4)),
            tabulate([M,N], 'float64', () => rng.uniform(-4,+4))
          ];
        }
      }
    ).it('solves generated large matrices', ([A,y]) => {
      const x = solve(A,y);
      expect( matmul2(A,x) ).toBeAllCloseTo(y);
    });

  });
}
