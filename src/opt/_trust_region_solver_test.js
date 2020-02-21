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

import {concat} from '../concat'
import {forEachItemIn, CUSTOM_MATCHERS} from '../jasmine_utils'
import {NDArray} from '../nd_array'
import {tabulate} from '../tabulate'
import {_rand_rankdef} from '../_test_data_generators'
import {zip_elems} from '../zip_elems'

import {diag_mat} from '../la/diag'
import {lstsq} from '../la/lstsq'
import {matmul2} from '../la/matmul'
import {FrobeniusNorm} from '../la/norm'
import {solve} from '../la/solve'

import {TrustRegionSolverLSQ} from './_trust_region_solver'


const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from;


describe('TrustRegionSolverLSQ', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  forEachItemIn(
    function*(){
      for( let run=1024; run-- > 0; )
      { const N = randInt(1,64),
              M = randInt(N,64);

        function* fJ(){
          for( let run=randInt(1,64); run-- > 0; )
          {
            const f = tabulate([M,1],'float64', () => Math.random()*4-2),
                  J = tabulate([M,N],'float64', () => Math.random()*4-2);

            Object.freeze(f);
            Object.freeze(J);
            Object.freeze(f.data.buffer);
            Object.freeze(J.data.buffer);

            yield [f,J];
          }
        }

        yield [M,N, fJ()];
      }
    }()
  ).it(`computeMinGlobal() works for random over-determined examples`, ([M,N, fJ]) => {
    const solver = new TrustRegionSolverLSQ(M,N),
         X_shape = Int32Array.of(N,1);

    for( const [f,J] of fJ )
    {
      solver.update(f,J);
      solver.computeMinGlobal();

      const x = new NDArray(X_shape, solver.X.slice()),
           Jx = matmul2(J,x);

      // every least square solution satisfies the normal equaltion Aᵀ(Ax - y) = 0
      expect( matmul2(J.T, zip_elems([Jx,f], (x,y) => x-y) ) ).toBeAllCloseTo(0);
    }
  })

  forEachItemIn(
    function*(){
      for( let run=1024; run-- > 0; )
      { const N = randInt(1,64),
              M = randInt(1,64);

        function* fJ(){
          for( let run=randInt(1,64); run-- > 0; )
          {
            const f = tabulate([M,1],'float64', () => Math.random()*4-2),
                  J = tabulate([M,N],'float64', () => Math.random()*4-2);

            Object.freeze(f);
            Object.freeze(J);
            Object.freeze(f.data.buffer);
            Object.freeze(J.data.buffer);

            yield [f,J];
          }
        }

        yield [M,N, fJ()];
      }
    }()
  ).it(`computeMinGlobal() works for random examples`, ([M,N, fJ]) => {
    const solver = new TrustRegionSolverLSQ(M,N),
         X_shape = Int32Array.of(N,1);

    const D = new Float64Array(N); // <- USED TO CHECK DIAGONAL SCALING
    D.fill(Number.MIN_VALUE);

    const NORM = new FrobeniusNorm();

    for( const [f,J] of fJ )
    {
      for( let j=N; j-- > 0; )
      {
        NORM.reset();
        for( let i=M; i-- > 0; )
          NORM.include( J(i,j) );
        D[j] = Math.max(D[j], NORM.result);
      }

      solver.update(f,J);        expect(solver.D).toBeAllCloseTo(D, {atol:0, rtol:Number.EPSILON*5});
      solver.computeMinGlobal(); expect(solver.D).toBeAllCloseTo(D, {atol:0, rtol:Number.EPSILON*5});

      const x = new NDArray(X_shape, solver.X.slice()),
           Jx = matmul2(J,x);

      // every least square solution satisfies the normal equaltion Aᵀ(Ax - y) = 0
      expect( matmul2(J.T, zip_elems([Jx,f], 'float64', (x,y) => x-y) ) ).toBeAllCloseTo(0);
    }
  })


  forEachItemIn(
    function*(){
      for( let run=1024; run-- > 0; )
      { const M = randInt(1,64),
              N = randInt(M,64);

        function* fJ(){
          for( let run=randInt(1,32); run-- > 0; )
          {
            const f = tabulate([M,1],'float64', () => Math.random()*4-2),
                  J = tabulate([M,N],'float64', () => Math.random()*4-2);

            Object.freeze(f);
            Object.freeze(J);
            Object.freeze(f.data.buffer);
            Object.freeze(J.data.buffer);

            yield [f,J];
          }
        }

        yield [M,N, fJ()];
      }
    }()
  ).it(`computeMinGlobal() works for random under-determined examples`, ([M,N, fJ]) => {
    const solver = new TrustRegionSolverLSQ(M,N),
         X_shape = Int32Array.of(N,1);

    for( const [f,J] of fJ )
    {
      solver.update(f,J);
      solver.computeMinGlobal();

      const x = new NDArray(X_shape, solver.X.slice());

      // check that x is a solution
      expect( matmul2(J,x) ).toBeAllCloseTo(f);

      const D = new NDArray(X_shape, solver.D);

      const JD = zip_elems([J, D.T], 'float64', (j,d) => j/d),
            Dx = zip_elems([x, D  ], 'float64', (x,d) => x*d);

      expect( matmul2(JD,Dx) ).toBeAllCloseTo(f);

      // check that it's the minimum norm solution
      // (taking into account diagonal scaling of the Trust-Region method)
      // https://www.math.usm.edu/lambers/mat419/lecture15.pdf
      const JJT = matmul2(JD,JD.T);
      expect(Dx).toBeAllCloseTo(matmul2(JD.T, solve(JJT,f)), {rtol: 1e-4});
    }
  })


  forEachItemIn(
    function*(){
      for( let run=1024; run-- > 0; )
      { const N = randInt(1,64),
              M = randInt(1,64);

        function* fJ(){
          for( let run=randInt(1,32); run-- > 0; )
          {
            const f = tabulate([M,1],'float64', () => Math.random()*4-2),
                 [J] = _rand_rankdef(M,N);

            Object.freeze(f);
            Object.freeze(J);
            Object.freeze(f.data.buffer);
            Object.freeze(J.data.buffer);

            yield [f,J];
          }
        }

        yield [M,N, fJ()];
      }
    }()
  ).it(`computeMinGlobal() works for random rank-deficient examples`, ([M,N, fJ]) => {
    const solver = new TrustRegionSolverLSQ(M,N),
         X_shape = Int32Array.of(N,1);

    for( const [f,J] of fJ )
    {
      solver.update(f,J);
      solver.computeMinGlobal();

      const x = new NDArray(X_shape, solver.X.slice()),
           Jx = matmul2(J,x);

      // check that x is a solution
      expect( matmul2(J.T, zip_elems([Jx,f], 'float64', (x,y) => x-y) ) ).toBeAllCloseTo(0);

      const D = new NDArray(X_shape, solver.D),
           JD = zip_elems([J, D.T], 'float64', (j,d) => j/d),
            Dx= zip_elems([x, D  ], 'float64', (x,d) => x*d);

      expect(Dx).toBeAllCloseTo( lstsq(JD,f) );
    }
  })


  forEachItemIn(
    function*(){
      for( let run=512; run-- > 0; )
      { const N = randInt(1,32),
              M = randInt(1,32); // <- TODO: undo this after debugging

        function* fJ(){
          for( let run=randInt(1,16); run-- > 0; )
          {
            const f = tabulate([M,1],'float64', () => Math.random()*8-4),
                  J = tabulate([M,N],'float64', () => Math.random()*8-4);
//                 [J] = _rand_rankdef(M,N);

            Object.freeze(f);
            Object.freeze(J);
            Object.freeze(f.data.buffer);
            Object.freeze(J.data.buffer);

            function* lambdas() {
              for( let trial=randInt(1,16); trial-- > 0; )
                yield Math.random()*8 + 1;
            }

            yield [f,J, lambdas()];
          }
        }

        yield [M,N, fJ()];
      }
    }()
  ).it(`computeMinRegularized(λSqrt) works for random rank-deficient examples`, ([M,N, fJ]) => {
    const solver = new TrustRegionSolverLSQ(M,N),
         X_shape = Int32Array.of(N,1);

    for( const [f,J, lambdas] of fJ )
    {
      solver.update(f,J);
      solver.computeMinGlobal();

      for( const λ of lambdas )
      {
        solver.computeMinRegularized(λ);

        const x = new NDArray(X_shape, solver.X.slice()),
             λD = diag_mat(solver.D).mapElems(d => d*λ);

        const B = concat([J,λD]),
              z = concat([f,tabulate([N,1], 'float64', () => 0)]);

//        const A = zip_elems([matmul2(J.T,J), λD], (j,d) => j + d*d),
//              y = matmul2(J.T, f);
//        expect( solve(A,y) ).toBeAllCloseTo( lstsq(B,z) );
//        expect(x).toBeAllCloseTo( solve(A,y) );

        expect(x).toBeAllCloseTo( lstsq(B,z) );
      }
    }
  })
})
