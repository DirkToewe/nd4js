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
import {tabulate} from '../tabulate'
import {_rand_int,
        _shuffle} from '../_test_data_generators'

import {matmul,
        matmul2} from '../la/matmul'
import {norm} from '../la/norm'
import {solve} from '../la/solve'

import {LBFGS_Solver       } from './_lbfgs_solver'
import {LBFGS_SolverRefImpl,
        _rand_updates      } from './_lbfgs_solver_test_utils'


describe('LBFGS_Solver', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  for( const [suffix, rand_scale] of [
    [''           , () => 1                  ],
    [' and scales', () => Math.random() + 0.5]
  ])
  {
    forEachItemIn(
      function*(){
        for( let run=0; run++ < 2; )
        for( let M=0; M++ <  8; )
        for( let N=0; N++ < 16; )
          yield [M,N];
      }()
    ).it('compute_Hv works given random updates' + suffix, ([M,N]) => {
      const ref = new LBFGS_SolverRefImpl(M,N),
            tst = new LBFGS_Solver       (M,N);

      expect(M).toBeGreaterThan(0);
      expect(N).toBeGreaterThan(0);

      const shape = Int32Array.of(N,1);

      let i=0;
      for( let [dx,dg] of _rand_updates(N) )
      {
        ref.update(dx.data, dg.data);
        tst.update(dx.data, dg.data);

        const       scale = rand_scale();
        ref.scale = scale;

        const compute_Hv = v => {
          v = Float64Array.from(v.data);
          tst.compute_Hv_phase1(v); for( let i=v.length; i-- > 0; ) v[i] /= scale;
          tst.compute_Hv_phase2(v);
          return v;
        };

        const indices = Int32Array.from({length: N}, (_,i) => i);
        _shuffle(indices);
        Object.freeze(indices.buffer);

        const {B} = ref,
               Z  = tabulate([N,N], (i,j) => 1*(i === indices[j]) ),
               y  = tabulate([N,1], () => Math.random()*8 - 4);

        const b = matmul(Z.T, B, Z),
              x = solve(b,y);

        const X = compute_Hv( matmul2(Z,y) );

        const Zx = matmul2(Z,x).reshape(-1),
            tol = {
              rtol: 0,
              atol: 1e-4 * Math.max(norm(X), norm(x))
            };
        expect(X).toBeAllCloseTo(Zx, tol);

        if( ++i >= 64 ) break;
      }
    })
  }
});