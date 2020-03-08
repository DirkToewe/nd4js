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
import {_rand_cols0,
        _rand_rankdef} from '../_test_data_generators'
import {zip_elems} from '../zip_elems'

import {diag_mat} from '../la/diag'
import {lstsq} from '../la/lstsq'
import {matmul,
        matmul2} from '../la/matmul'
import {FrobeniusNorm} from '../la/norm'
import {rand_ortho} from '../la/rand_ortho'
import {solve} from '../la/solve'
import {svd_lstsq} from '../la/svd'
import {svd_jac_2sided} from '../la/svd_jac_2sided'

import {num_grad} from './num_grad'
import {TrustRegionSolverLSQ} from './_trust_region_solver'

import {computeMinGlobal_overdet_gen,
        computeMinGlobal_underdet_gen} from './_trust_region_solver_test_data'


const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from;


describe('TrustRegionSolverLSQ', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  forEachItemIn(
    function*(){
      for( let run=512; run-- > 0; )
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
         X_shape = Int32Array.of(N,1),
         f_shape = Int32Array.of(M);

    for( const [f,J] of fJ )
    {
      const neg_f = new NDArray(f_shape, f.data.map(x => -x));

      solver.update(neg_f,J);

      for( let i=2; i-- > 0; )
      {
        solver.computeMinGlobal();

        const x = new NDArray(X_shape, solver.X.slice()),
             Jx = matmul2(J,x);

        // every least square solution satisfies the normal equaltion Aᵀ(Ax - y) = 0
        expect( matmul2(J.T, zip_elems([Jx,f], (x,y) => x-y) ) ).toBeAllCloseTo(0);
      }
    }
  })


  forEachItemIn(
    function*(){
      for( let run=768; run-- > 0; )
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
         X_shape = Int32Array.of(N,1),
         f_shape = Int32Array.of(M);

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

      const neg_f = new NDArray(f_shape, f.data.map(x => -x));

      solver.update(neg_f,J);
      expect(solver.D).toBeAllCloseTo(D, {atol:0, rtol:Number.EPSILON*5});

      for( let i=2; i-- > 0; )
      {
        solver.computeMinGlobal();
        expect(solver.D).toBeAllCloseTo(D, {atol:0, rtol:Number.EPSILON*5});

        const x = new NDArray(X_shape, solver.X.slice()),
             Jx = matmul2(J,x);

        // every least square solution satisfies the normal equaltion Aᵀ(Ax - y) = 0
        expect( matmul2(J.T, zip_elems([Jx,f], 'float64', (x,y) => x-y) ) ).toBeAllCloseTo(0);
      }
    }
  })


  forEachItemIn(
    function*(){
      for( let run=384; run-- > 0; )
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
         X_shape = Int32Array.of(N,1),
         f_shape = Int32Array.of(M);

    for( const [f,J] of fJ )
    {
      const neg_f = new NDArray(f_shape, f.data.map(x => -x));

      solver.update(neg_f,J);
      for( let i=2; i-- > 0; )
      {
        solver.computeMinGlobal();

        const x = new NDArray(X_shape, solver.X.slice());

        // check that x is a solution
        expect( matmul2(J,x) ).toBeAllCloseTo(f);

        const D = new NDArray(X_shape, solver.D);

        const JD = zip_elems([J, D.T], 'float64', (j,d) => d===0 ? j : j/d),
              Dx = zip_elems([x, D  ], 'float64', (x,d) => x*d);

        expect( matmul2(JD,Dx) ).toBeAllCloseTo(f);

        // check that it's the minimum norm solution
        // (taking into account diagonal scaling of the Trust-Region method)
        // https://www.math.usm.edu/lambers/mat419/lecture15.pdf
        const JJT = matmul2(JD,JD.T);
        expect(Dx).toBeAllCloseTo(matmul2(JD.T, solve(JJT,f)), {rtol: 5e-4});
      }
    }
  })


  forEachItemIn(
    function*(){
      for( let run=256; run-- > 0; )
      { const N = randInt(1,48),
              M = randInt(1,48);

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
         X_shape = Int32Array.of(N,1),
         f_shape = Int32Array.of(M);

    for( const [f,J] of fJ )
    {
      const neg_f = new NDArray(f_shape, f.data.map(x => -x));

      solver.update(neg_f,J);
      for( let i=2; i-- > 0; )
      {
        solver.computeMinGlobal();

        const x = new NDArray(X_shape, solver.X.slice()),
             Jx = matmul2(J,x);

        // check that x is a solution
        expect( matmul2(J.T, zip_elems([Jx,f], 'float64', (x,y) => x-y) ) ).toBeAllCloseTo(0);

        const D = new NDArray(X_shape, solver.D),
             JD = zip_elems([J, D.T], 'float64', (j,d) => d===0 ? j : j/d),
              Dx= zip_elems([x, D  ], 'float64', (x,d) => x*d);

        expect(Dx).toBeAllCloseTo( svd_lstsq(...svd_jac_2sided(JD),f) );
      }
    }
  })


  forEachItemIn(
    function*(){
      function* shapes() {
        for( let M=8; M > 1; M-- )
        for( let N=8; N > 1; N-- )
          yield [M,N];

        for( let run=192; run-- > 0; )
        {
          const M = randInt(1,64),
                N = randInt(1,64); // <- TODO remove after testing
          yield [M,N];
        }
      }

      for( const [M,N] of shapes() )
      {
        function* fJ()
        {
          for( let run=randInt(1,32); run-- > 0; )
          {
            const f =  tabulate([M,1],'float64', () => Math.random()*4-2),
                  J =_rand_cols0(M,N);

            Object.freeze(f);
            Object.freeze(f.data.buffer);

            yield [f,J];
          }
        }

        yield [M,N, fJ()];
      }
    }()
  ).it('computeMinGlobal() works for random examples with zero columns in J', ([M,N, fJ]) => {
    const solver = new TrustRegionSolverLSQ(M,N),
         X_shape = Int32Array.of(N,1),
         f_shape = Int32Array.of(M);

    for( const [f,J] of fJ )
    {
      const neg_f = new NDArray(f_shape, f.data.map(x => -x));

      solver.update(neg_f,J);
      for( let i=2; i-- > 0; )
      {
        solver.computeMinGlobal();

        const x = new NDArray(X_shape, solver.X.slice()),
             Jx = matmul2(J,x);

        // check that x is a solution
        expect( matmul2(J.T, zip_elems([Jx,f], 'float64', (x,y) => x-y) ) ).toBeAllCloseTo(0);

        const D = new NDArray(X_shape, solver.D),
             JD = zip_elems([J, D.T], 'float64', (j,d) => d===0 ? j : j/d),
              Dx= zip_elems([x, D  ], 'float64', (x,d) => x*d);

        expect(Dx).toBeAllCloseTo( svd_lstsq(...svd_jac_2sided(JD),f) );
      }
    }
  })


  forEachItemIn(
    computeMinGlobal_overdet_gen()
  ).it(`computeMinRegularized(0) works for pre-generated over-determined examples`, ([f,J, R,DR]) => {
    const [M,N] = J.shape,
         solver = new TrustRegionSolverLSQ(M,N)

    for( let i=4; i-- > 0; )
    {
      solver.update(f,J);
      for( let l=4; l-- > 0; )
      {
        solver.computeMinGlobal();

        expect( solver.scaledNorm(solver.X) ).toBeAllCloseTo(R);

        for( let j=4; j-- > 0; )
        {
          for( let k=4; k-- > 0; )
          {
            const [r,dr] = solver.computeMinRegularized(0);

            expect( r).toBeAllCloseTo( R);
            expect(dr).toBeAllCloseTo(DR);

            expect( solver.scaledNorm(solver.X) ).toBeAllCloseTo(R);
          }
          solver.computeMinRegularized( Math.random() );
        }
      }
    }
  })


  forEachItemIn(
    computeMinGlobal_underdet_gen()
  ).it(`computeMinRegularized(0) works for pre-generated under-determined examples`, ([f,J, R,DR]) => {
    const [M,N] = J.shape,
         solver = new TrustRegionSolverLSQ(M,N)

    for( let i=4; i-- > 0; )
    {
      solver.update(f,J);
      for( let l=4; l-- > 0; )
      {
        solver.computeMinGlobal();

        expect( solver.scaledNorm(solver.X) ).toBeAllCloseTo(R);

        for( let j=4; j-- > 0; )
        {
          for( let k=4; k-- > 0; )
          {
            const [r,dr] = solver.computeMinRegularized(0);

            expect( r).toBeAllCloseTo( R);
            expect(dr).toBeAllCloseTo(DR);

            expect( solver.scaledNorm(solver.X) ).toBeAllCloseTo(R);
          }
          solver.computeMinRegularized( Math.random() );
        }
      }
    }
  })


  forEachItemIn(
    function*(){
      function* shapes() {
        for( let N=0  ; N++ < 4; )
        for( let M=N-1; M++ < 4; )
          yield [M,N];

        for( let run=512; run-- > 0; ) {
          const N = randInt(1,16),
                M = randInt(N,16);
          yield[M,N];
        }
      }

      for( const [M,N] of shapes() )
      {
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
              for( let trial=randInt(1,16); trial-- > 0; ) {
                yield Math.random()*8 + 0.1;
              }
            }

            yield [f,J, lambdas()];
          }
        }

        yield [M,N, fJ()];
      }
    }()
  ).it(`computeMinRegularized(λ) works for random over-determined examples`, ([M,N, fJ]) => {
    const solver = new TrustRegionSolverLSQ(M,N),
         X_shape = Int32Array.of(N,1),
         f_shape = Int32Array.of(M);

    const g = num_grad( λ => {
      const [r,dr] = solver.computeMinRegularized(λ);
      return r;
    });

    for( const [f,J, lambdas] of fJ )
    {
      const neg_f = new NDArray(f_shape, f.data.map(x => -x));

      const X = svd_lstsq(...svd_jac_2sided(J), f);
      Object.freeze(X.data.buffer);
      Object.freeze(X);

      solver.update(neg_f,J);
      for( let l=2; l-- > 0; )
      {
        solver.computeMinGlobal();

        const x = new NDArray(X_shape, solver.X.slice());

        expect(x).toBeAllCloseTo(X);

        for( const λ of lambdas )
        {
          const [r,dr] = solver.computeMinRegularized(λ),
                 λSqrt = Math.sqrt(λ);

          const x = new NDArray(X_shape, solver.X.slice()),
               λD = diag_mat(solver.D).mapElems(d => d*λSqrt),
                B = concat([J,λD]),
                z = concat([f,tabulate([N,1], 'float64', () => 0)]);

//          const A = zip_elems([matmul2(J.T,J), λD], (j,d) => j + d*d),
//                y = matmul2(J.T, f);
//          expect( solve(A,y) ).toBeAllCloseTo( lstsq(B,z) );
//          expect(x).toBeAllCloseTo( solve(A,y) );

          expect(x).toBeAllCloseTo( svd_lstsq(...svd_jac_2sided(B), z) );

          expect(dr).toBeAllCloseTo( g(λ), {rtol: 1e-3} );
        }
      }
    }
  })


  forEachItemIn(
    function*(){
      function* shapes() {
        for( let run=384; run-- > 0; ) {
          const M = randInt(1,16),
                N = randInt(M,16);
          yield[M,N];
        }
      }

      for( const [M,N] of shapes() )
      {
        const L = Math.min(M,N);

        function* fJ(){
          for( let run=randInt(1,16); run-- > 0; )
          {
            const rank = randInt(1,1+L);

            const f = tabulate([M,1],'float64', () => Math.random()*8-4),
                  J = tabulate([M,N],'float64', () => Math.random()*8-4);

            Object.freeze(f);
            Object.freeze(J);
            Object.freeze(f.data.buffer);
            Object.freeze(J.data.buffer);

            function* lambdas() {
              for( let trial=randInt(1,16); trial-- > 0; ) {
                yield Math.random()*8 + 0.1;
              }
            }

            yield [f,J, lambdas()];
          }
        }

        yield [M,N, fJ()];
      }
    }()
  ).it(`computeMinRegularized(λ) works for random under-determined examples`, ([M,N, fJ]) => {
    const solver = new TrustRegionSolverLSQ(M,N),
         X_shape = Int32Array.of(N,1),
         f_shape = Int32Array.of(M);

    const g = num_grad( λ => {
      const [r,dr] = solver.computeMinRegularized(λ);
      return r;
    });

    for( const [f,J, lambdas] of fJ )
    {
      const neg_f = new NDArray(f_shape, f.data.map(x => -x));

      solver.update(neg_f,J);
      for( let l=2; l-- > 0; )
      {
        solver.computeMinGlobal();

        const x = new NDArray(X_shape, solver.X.slice()),
              D = new NDArray(X_shape, solver.D.slice()),
             JD = zip_elems([J,D.T], (j,d) => d===0 ? j : j/d),
             dx = zip_elems([x,D  ], (x,d) => x*d),
             DX = svd_lstsq(...svd_jac_2sided(JD), f);

        expect(dx).toBeAllCloseTo(DX);

        for( const λ of lambdas )
        {
          const [r,dr] = solver.computeMinRegularized(λ),
                 λSqrt = Math.sqrt(λ);

          if( isNaN( r) ) throw new Error('Assertion failed.');
          if( isNaN(dr) ) throw new Error('Assertion failed.');

          const x = new NDArray(X_shape, solver.X.slice()),
                D = new NDArray(X_shape, solver.D.slice()),
               JD = zip_elems([J,D.T], (j,d) => d===0 ? j : j/d),
               dx = zip_elems([x,D  ], (x,d) => x*d),
                I = tabulate([N,N], (i,j) => (i===j)*λSqrt);

          const B = concat([JD,I]),
                z = concat([f,tabulate([N,1], 'float64', () => 0)]);

//          const A = zip_elems([matmul2(J.T,J), λD], (j,d) => j + d*d),
//                y = matmul2(J.T, f);
//          expect( solve(A,y) ).toBeAllCloseTo( lstsq(B,z) );
//          expect(x).toBeAllCloseTo( solve(A,y) );

          const DX = svd_lstsq(...svd_jac_2sided(B), z);

          expect(dx).toBeAllCloseTo(DX);

          expect(dr).toBeAllCloseTo( g(λ), {rtol: 1e-3, atol: 1e-6} );
        }
      }
    }
  })


  forEachItemIn(
    function*(){
      function* shapes() {
        for( let N=0  ; N++ < 4; )
        for( let M=N-1; M++ < 4; )
          yield [M,N];

        for( let run=768; run-- > 0; ) {
          const N = randInt(1,16),
                M = randInt(1,16);
          yield[M,N];
        }
      }

      for( const [M,N] of shapes() )
      {
        const L = Math.min(M,N);

        function* fJ(){
          for( let run=randInt(1,16); run-- > 0; )
          {
            const rank = randInt(1,1+L);

            const U = rand_ortho('float64', M,L),
                  V = rand_ortho('float64', L,N),
                  S = tabulate([L,L], 'float64', (i,j) => {
                    if( i !== j || rank <= i ) return 0;

                    return Math.random() + 0.1
                  });

            const f = tabulate([M,1],'float64', () => Math.random()*8-4),
                  J = matmul(U,S,V);

            Object.freeze(f);
            Object.freeze(J);
            Object.freeze(f.data.buffer);
            Object.freeze(J.data.buffer);

            function* lambdas() {
              for( let trial=randInt(1,16); trial-- > 0; ) {
                yield Math.random()*8 + 0.1;
              }
            }

            yield [f,J, lambdas()];
          }
        }

        yield [M,N, fJ()];
      }
    }()
  ).it(`computeMinRegularized(λ) works for random rank-deficient examples`, ([M,N, fJ]) => {
    const solver = new TrustRegionSolverLSQ(M,N),
         X_shape = Int32Array.of(N,1),
         f_shape = Int32Array.of(M);

    const g = num_grad( λ => {
      const [r,dr] = solver.computeMinRegularized(λ);
      return r;
    });

    for( const [f,J, lambdas] of fJ )
    {
      const neg_f = new NDArray(f_shape, f.data.map(x => -x));

      solver.update(neg_f,J);
      for( let l=2; l-- > 0; )
      {
        solver.computeMinGlobal();

        const x = new NDArray(X_shape, solver.X.slice()),
              D = new NDArray(X_shape, solver.D.slice()),
             JD = zip_elems([J,D.T], (j,d) => d===0 ? j : j/d),
             dx = zip_elems([x,D  ], (x,d) => x*d),
             DX = svd_lstsq(...svd_jac_2sided(JD), f);

        expect(dx).toBeAllCloseTo(DX);

        for( const λ of lambdas )
        {
          const [r,dr] = solver.computeMinRegularized(λ),
                 λSqrt = Math.sqrt(λ);

          if( isNaN( r) ) throw new Error('Assertion failed.');
          if( isNaN(dr) ) throw new Error('Assertion failed.');

          const x = new NDArray(X_shape, solver.X.slice()),
                D = new NDArray(X_shape, solver.D.slice()),
                I = tabulate([N,N], (i,j) => (i===j)*λSqrt),
               JD = zip_elems([J,D.T], (j,d) => d===0 ? j : j/d),
               dx = zip_elems([x,D  ], (x,d) => x*d);

          const B = concat([JD,I]),
                z = concat([f,tabulate([N,1], 'float64', () => 0)]);

//          const A = zip_elems([matmul2(J.T,J), λD], (j,d) => j + d*d),
//                y = matmul2(J.T, f);
//          expect( solve(A,y) ).toBeAllCloseTo( lstsq(B,z) );
//          expect(x).toBeAllCloseTo( solve(A,y) );

          const DX = svd_lstsq(...svd_jac_2sided(B), z);

          expect(dx).toBeAllCloseTo(DX);

          expect(dr).toBeAllCloseTo( g(λ), {rtol: 1e-3, atol: 1e-6} );
        }
      }
    }
  })


  forEachItemIn(
    function*(){
      function* shapes() {
        for( let N=0  ; N++ < 4; )
        for( let M=N-1; M++ < 4; )
          yield [M,N];

        for( let run=768; run-- > 0; ) {
          const N = randInt(1,16),
                M = randInt(1,16);
          yield[M,N];
        }
      }

      for( const [M,N] of shapes() )
      {
        const L = Math.min(M,N);

        function* fJ(){
          for( let run=randInt(1,16); run-- > 0; )
          {
            const f = tabulate([M,1],'float64', () => Math.random()*8-4),
                  J = _rand_cols0(M,N);

            Object.freeze(f);
            Object.freeze(f.data.buffer);

            function* lambdas() {
              for( let trial=randInt(1,16); trial-- > 0; ) {
                yield Math.random()*8 + 0.1;
              }
            }

            yield [f,J, lambdas()];
          }
        }

        yield [M,N, fJ()];
      }
    }()
  ).it('computeMinRegularized(λ) works for random examples with zero columns in J', ([M,N, fJ]) => {
    const solver = new TrustRegionSolverLSQ(M,N),
         X_shape = Int32Array.of(N,1),
         f_shape = Int32Array.of(M);

    const g = num_grad( λ => {
      const [r,dr] = solver.computeMinRegularized(λ);
      return r;
    });

    for( const [f,J, lambdas] of fJ )
    {
      const neg_f = new NDArray(f_shape, f.data.map(x => -x));

      solver.update(neg_f,J);
      for( let l=2; l-- > 0; )
      {
        solver.computeMinGlobal();

        const x = new NDArray(X_shape, solver.X.slice()),
              D = new NDArray(X_shape, solver.D.slice()),
             JD = zip_elems([J,D.T], (j,d) => d===0 ? j : j/d),
             dx = zip_elems([x,D  ], (x,d) => x*d),
             DX = svd_lstsq(...svd_jac_2sided(JD), f);

        expect(dx).toBeAllCloseTo(DX);

        for( const λ of lambdas )
        {
          const [r,dr] = solver.computeMinRegularized(λ),
                 λSqrt = Math.sqrt(λ);

          if( ! isFinite( r) ) throw new Error('Assertion failed.');
          if( ! isFinite(dr) ) throw new Error('Assertion failed.');

          const x = new NDArray(X_shape, solver.X.slice()),
                D = new NDArray(X_shape, solver.D.slice()),
                I = tabulate([N,N], (i,j) => (i===j)*λSqrt),
               JD = zip_elems([J,D.T], (j,d) => d===0 ? j : j/d),
               dx = zip_elems([x,D  ], (x,d) => x*d);

          const B = concat([JD,I]),
                z = concat([f,tabulate([N,1], 'float64', () => 0)]);

//          const A = zip_elems([matmul2(J.T,J), λD], (j,d) => j + d*d),
//                y = matmul2(J.T, f);
//          expect( solve(A,y) ).toBeAllCloseTo( lstsq(B,z) );
//          expect(x).toBeAllCloseTo( solve(A,y) );

          const DX = svd_lstsq(...svd_jac_2sided(B), z);

          expect(dx).toBeAllCloseTo(DX);

          expect(dr).toBeAllCloseTo( g(λ), {rtol: 1e-3, atol: 1e-6} );
        }
      }
    }
  })


  forEachItemIn(
    function*(){
      for( let run=768; run-- > 0; )
      { const N = randInt(1,32),
              M = randInt(1,32);

        function* fJ(){
          for( let run=randInt(1,32); run-- > 0; )
          {
            const f = tabulate([M  ],'float64', () => Math.random()*4-2),
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
  ).it(`G is the correct gradient`, ([M,N, fJ]) => {
    const solver = new TrustRegionSolverLSQ(M,N);

    const NORM = new FrobeniusNorm();

    for( const [f,J] of fJ )
    {
      const [M,N] = J.shape;

      const g = num_grad( x => {
        if( x.ndim     !== 1 ) throw new Error('Assertion failed.');
        if( x.shape[0] !== N ) throw new Error('Assertion failed.');
//        x = x.data;

        let result = 0.0;

        for( let i=M; i-- > 0; )
        {
          let s = 0;

          for( let j=N; j-- > 0; )
            s += J(i,j) * x(j);

          s += f(i);
          result += s*s;
        }

        return result/2;
      });

      solver.update(f,J);

      expect( solver.G ).toBeAllCloseTo( g( new Float64Array(N) ) );

//      if( Math.random() < 0.5 )
//        solver.computeMinGlobal();
    }
  })


  forEachItemIn(
    function*(){
      for( let run=512; run-- > 0; )
      { const N = randInt(1,64),
              M = randInt(1,64);

        function* fJ(){
          for( let run=randInt(1,64); run-- > 0; )
          {
            const f = tabulate([M  ],'float64', () => Math.random()*4-2),
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
  ).it(`cauchyPointTravel() works for random examples`, ([M,N, fJ]) => {
    const solver = new TrustRegionSolverLSQ(M,N),
         X_shape = Int32Array.of(N,1),
         f_shape = Int32Array.of(M);

    for( const [f,J] of fJ )
    {
      const g = num_grad(c => { // <- linearized approximation
        const {M,N, G} = solver;

        let result = 0;

        for( let i=M; i-- > 0; )
        {
          let s = f.data[i];

          for( let j=N; j-- > 0; )
            s += J(i,j) * G[j] * c;

          result += s*s;
        }
        
        return result/2;
      });

      solver.update(f,J);

      const  c = solver.cauchyPointTravel();
      expect(c).not.toBeGreaterThan(0);

      if( isFinite(c) ) {
        expect(  g(c) ).toBeAllCloseTo(0);
      }
      else {
        expect( g(+1) ).toBeAllCloseTo( g(0) );
        expect( g(-1) ).toBeAllCloseTo( g(0) );
      }

      if( Math.random() < 0.5 )
        solver.computeMinGlobal();
    }
  })
})
