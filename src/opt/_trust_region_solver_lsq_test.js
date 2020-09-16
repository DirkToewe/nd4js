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
import {NDArray, array} from '../nd_array'
import {tabulate} from '../tabulate'
import {_rand_cols0,
        _rand_rankdef} from '../_test_data_generators'
import {zip_elems} from '../zip_elems'

import {diag_mat} from '../la/diag'
import {matmul2} from '../la/matmul'
import {FrobeniusNorm, norm} from '../la/norm'
import {solve} from '../la/solve'
import {svd_lstsq} from '../la/svd'
import {svd_jac_2sided} from '../la/svd_jac_2sided'

import {cartesian_prod,
        linspace,
        range} from "../iter";

import {beale}             from './test_fn/beale';
import {brown_badscale}    from './test_fn/brown_badscale';
import {freudenstein_roth} from './test_fn/freudenstein_roth';
import {helical_valley}    from "./test_fn/helical_valley";
import {JennrichSampson}   from './test_fn/jennrich_sampson';
import {powell_badscale}   from './test_fn/powell_badscale';
import {Rastrigin}         from './test_fn/rastrigin'
import {Rosenbrock}        from './test_fn/rosenbrock'


import {num_grad} from './num_grad'
import {TrustRegionSolverLSQ} from './_trust_region_solver_lsq'

import {computeNewtonRegularized0_overdet_gen,
        computeNewtonRegularized0_underdet_gen} from './_trust_region_solver_lsq_test_data'

const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from;


describe('TrustRegionSolverLSQ', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  forEachItemIn(
    function*(){
      for( let run=128; run-- > 0; )
      { const N = randInt(1,32),
              M = randInt(N,32);

        function* samples(){
          for( let run=randInt(1,32); run-- > 0; )
          {
            const dx = tabulate([  N],'float64', () => Math.random()*4-2),
                   f = tabulate([M  ],'float64', () => Math.random()*4-2),
                   J = tabulate([M,N],'float64', () => Math.random()*4-2);

            Object.freeze(dx);
            Object.freeze(f);
            Object.freeze(J);
            Object.freeze(dx.data.buffer);
            Object.freeze(f.data.buffer);
            Object.freeze(J.data.buffer);

            yield [dx,f,J];
          }
        }

        yield [M,N, samples()];
      }
    }()
  ).it(`computeNewton() works for random  over-determined examples`, ([M,N, samples]) => {

    let solver, dx,f,J, x0 = 0;

    let nCalls = 0;
    function fJ( x )
    {
      expect(nCalls++).toBe(0);
      expect(x).toBeAllCloseTo(x0);
      return [f,J];
    }

    const   D = new Float64Array(N), // <- USED TO CHECK DIAGONAL SCALING
      x_shape = Int32Array.of(N,1);

    let iter=0;
    for( [dx,f,J] of samples )
    {
      nCalls = 0;

      x0 = zip_elems([x0,dx], (x,dx) => x+dx);

      if( undefined === solver ) {
        solver = new TrustRegionSolverLSQ(fJ, x0);
        expect(iter).toBe(0);
      }
      else {
        solver.considerMove(dx.data);
        solver.makeConsideredMove();
        expect(iter).toBeGreaterThan(0);
      }

      expect(nCalls).toBe(1);

      ++iter;

      expect(solver.M).toBe(M);
      expect(solver.N).toBe(N);
      expect( solver.X0 ).toBeAllCloseTo( x0.data );
      expect( solver.F0 ).toBeAllCloseTo(  f.data );
      expect( solver.J0 ).toBeAllCloseTo(  J.data );

      for( let j=N; j-- > 0; )
      {
        const norm = new FrobeniusNorm();
        for( let i=M; i-- > 0; )
          norm.include( J(i,j) );
        D[j] = Math.max(D[j], norm.result);
      }

      expect(solver.D).toBeAllCloseTo(D, {rtol: Number.EPSILON*5, atol: 0});

      solver.computeNewton();
      expect(solver.rank).toBe(N);

      expect(solver.D).toBeAllCloseTo(D, {rtol: Number.EPSILON*5, atol: 0});

      f = f.reshape(-1,1);
      const x = new NDArray(x_shape, solver.newton_dX),
           Jx = matmul2(J,x);

      // Newton point has to satisfy normal equaltion Jᵀ(Jx + f) = 0
      expect( matmul2(J.T, zip_elems([Jx,f], (x,y) => x+y) ) ).toBeAllCloseTo(0, {atol: 1e-7});
    }
  })


  forEachItemIn(
    function*(rng){
      for( let run=96; run-- > 0; )
      { const M = rng.int(1,42),
              N = rng.int(M,42);

        function* samples(){
          for( let run=rng.int(1,42); run-- > 0; )
          {
            const dx = tabulate([  N],'float64', () => rng.uniform(-2,+2) ),
                   f = tabulate([M  ],'float64', () => rng.uniform(-2,+2) ),
                   J = tabulate([M,N],'float64', () => rng.uniform(-2,+2) );

            Object.freeze(dx);
            Object.freeze(f);
            Object.freeze(J);
            Object.freeze(dx.data.buffer);
            Object.freeze(f.data.buffer);
            Object.freeze(J.data.buffer);

            yield [dx,f,J];
          }
        }

        yield [M,N, samples()];
      }
    }
  ).it(`computeNewton() works for random under-determined examples`, ([M,N, samples]) => {

    let solver, dx,f,J, x0 = 0;

    let nCalls = 0;
    function fJ( x )
    {
      expect(nCalls++).toBe(0);
      expect(x).toBeAllCloseTo(x0);
      return [f,J];
    }

    const D_test = new Float64Array(N), // <- USED TO CHECK DIAGONAL SCALING
         x_shape = Int32Array.of(N,1);

    let iter=0;
    for( [dx,f,J] of samples )
    {
      nCalls = 0;

      x0 = zip_elems([x0,dx], (x,dx) => x+dx);

      if( undefined === solver ) {
        solver = new TrustRegionSolverLSQ(fJ, x0);
        expect(iter).toBe(0);
      }
      else {
        solver.considerMove(dx.data);
        solver.makeConsideredMove();
        expect(iter).toBeGreaterThan(0);
      }

      expect(nCalls).toBe(1);

      ++iter;

      expect(solver.M).toBe(M);
      expect(solver.N).toBe(N);
      expect( solver.X0 ).toBeAllCloseTo( x0.data );
      expect( solver.F0 ).toBeAllCloseTo(  f.data );
      expect( solver.J0 ).toBeAllCloseTo(  J.data );

      for( let j=N; j-- > 0; )
      {
        const norm = new FrobeniusNorm();
        for( let i=M; i-- > 0; )
          norm.include( J(i,j) );
        D_test[j] = Math.max(D_test[j], norm.result);
      }

      expect(solver.D).toBeAllCloseTo(D_test, {rtol: Number.EPSILON*5, atol: 0});

      solver.computeNewton();
      expect(solver.rank).toBe(M);

      expect(solver.D).toBeAllCloseTo(D_test, {rtol: Number.EPSILON*5, atol: 0});

      f = f.reshape(-1,1);
      const x = new NDArray(x_shape, solver.newton_dX),
           Jx = matmul2(J,x);

      const neg_f = f.mapElems(x => -x);

      // check that x is a solution
      expect(Jx).toBeAllCloseTo(neg_f);

      const D = new NDArray(x_shape, solver.D);

      const JD = zip_elems([J, D.T], 'float64', (j,d) => d===0 ? j : j/d),
            Dx = zip_elems([x, D  ], 'float64', (x,d) => x*d);

      expect( matmul2(JD,Dx) ).toBeAllCloseTo(neg_f);

      // check that it's the minimum norm solution
      // (taking into account diagonal scaling of the Trust-Region method)
      // https://www.math.usm.edu/lambers/mat419/lecture15.pdf
      const JJT = matmul2(JD,JD.T);
      expect(Dx).toBeAllCloseTo(matmul2(JD.T, solve(JJT,neg_f)), {rtol: 5e-4});
    }
  })


  forEachItemIn(
    function*(){
      for( let run=128; run-- > 0; )
      { const N = randInt(1,32),
              M = randInt(1,32);

        function* samples(){
          for( let run=randInt(1,32); run-- > 0; )
          {
            const dx = tabulate([  N],'float64', () => Math.random()*4-2),
                   f = tabulate([M  ],'float64', () => Math.random()*4-2),
                  [J,rank] = _rand_rankdef(M,N);

            expect(rank%1).toBe(0);

            Object.freeze(dx);
            Object.freeze(f);
            Object.freeze(J);
            Object.freeze(dx.data.buffer);
            Object.freeze(f.data.buffer);
            Object.freeze(J.data.buffer);

            yield [dx,f,J, rank|0];
          }
        }

        yield [M,N, samples()];
      }
    }()
  ).it(`computeNewton() works for random  rank-deficient  examples`, ([M,N, samples]) => {

    let solver, dx,f,J,rank, x0 = 0;

    let nCalls = 0;
    function fJ( x )
    {
      expect(nCalls++).toBe(0);
      expect(x).toBeAllCloseTo(x0);
      return [f,J];
    }

    const D_test = new Float64Array(N), // <- USED TO CHECK DIAGONAL SCALING
         x_shape = Int32Array.of(N,1);

    let iter=0;
    for( [dx,f,J,rank] of samples )
    {
      nCalls = 0;

      x0 = zip_elems([x0,dx], (x,dx) => x+dx);

      if( undefined === solver ) {
        solver = new TrustRegionSolverLSQ(fJ, x0);
        expect(iter).toBe(0);
      }
      else {
        solver.considerMove(dx.data);
        solver.makeConsideredMove();
        expect(iter).toBeGreaterThan(0);
      }

      expect(nCalls).toBe(1);

      ++iter;
      expect(solver.M).toBe(M);
      expect(solver.N).toBe(N);
      expect( solver.X0 ).toBeAllCloseTo( x0.data );
      expect( solver.F0 ).toBeAllCloseTo(  f.data );
      expect( solver.J0 ).toBeAllCloseTo(  J.data );

      for( let j=N; j-- > 0; )
      {
        const norm = new FrobeniusNorm();
        for( let i=M; i-- > 0; )
          norm.include( J(i,j) );
        D_test[j] = Math.max(D_test[j], norm.result);
      }

      expect(solver.D).toBeAllCloseTo(D_test, {rtol: Number.EPSILON*5, atol: 0});

      solver.computeNewton();
      expect(solver.rank).toBe(rank);

      expect(solver.D).toBeAllCloseTo(D_test, {rtol: Number.EPSILON*5, atol: 0});

      f = f.reshape(-1,1);
      const x = new NDArray(x_shape, solver.newton_dX),
           Jx = matmul2(J,x);

      // Newton point has to satisfy normal equaltion Jᵀ(Jx + f) = 0
      expect( matmul2(J.T, zip_elems([Jx,f], (x,y) => x+y) ) ).toBeAllCloseTo(0, {atol: 1e-7});

      const D = new NDArray(x_shape, solver.D),
           JD = zip_elems([J, D.T], 'float64', (j,d) => d===0 ? j : j/d),
            Dx= zip_elems([x, D  ], 'float64', (x,d) => x*d);

      const neg_f = f.mapElems(x => -x);

      expect(Dx).toBeAllCloseTo( svd_lstsq(...svd_jac_2sided(JD),neg_f) );
    }
  })


  forEachItemIn(
    function*(){
      function* shapes() {
        for( let M=8; M > 1; M-- )
        for( let N=8; N > 1; N-- )
          yield [M,N];

        for( let run=128; run-- > 0; )
        {
          const M = randInt(1,32),
                N = randInt(1,32); // <- TODO remove after testing
          yield [M,N];
        }
      }

      for( const [M,N] of shapes() )
      {
        function* samples()
        {
          for( let run=randInt(1,32); run-- > 0; )
          {
            const dx =  tabulate([  N],'float64', () => Math.random()*4-2),
                   f =  tabulate([M  ],'float64', () => Math.random()*4-2),
                   J =_rand_cols0(M,N);

            Object.freeze(dx);
            Object.freeze(f);
            Object.freeze(J);
            Object.freeze(dx.data.buffer);
            Object.freeze(f.data.buffer);
            Object.freeze(J.data.buffer);

            yield [dx,f,J];
          }
        }

        yield [M,N, samples()];
      }
    }()
  ).it('computeNewton() works for random examples with zero columns in J', ([M,N, samples]) => {

    let solver, dx,f,J, x0 = 0;

    let nCalls = 0;
    function fJ( x )
    {
      expect(nCalls++).toBe(0);
      expect(x).toBeAllCloseTo(x0);
      return [f,J];
    }

    const D_test = new Float64Array(N), // <- USED TO CHECK DIAGONAL SCALING
         x_shape = Int32Array.of(N,1);

    let iter=0;
    for( [dx,f,J] of samples )
    {
      nCalls = 0;

      x0 = zip_elems([x0,dx], (x,dx) => x+dx);

      if( undefined === solver ) {
        solver = new TrustRegionSolverLSQ(fJ, x0);
        expect(iter).toBe(0);
      }
      else {
        solver.considerMove(dx.data);
        solver.makeConsideredMove();
        expect(iter).toBeGreaterThan(0);
      }

      expect(nCalls).toBe(1);

      ++iter;
      expect(solver.M).toBe(M);
      expect(solver.N).toBe(N);
      expect( solver.X0 ).toBeAllCloseTo( x0.data );
      expect( solver.F0 ).toBeAllCloseTo(  f.data );
      expect( solver.J0 ).toBeAllCloseTo(  J.data );

      for( let j=N; j-- > 0; )
      {
        const norm = new FrobeniusNorm();
        for( let i=M; i-- > 0; )
          norm.include( J(i,j) );
        D_test[j] = Math.max(D_test[j], norm.result);
      }

      expect(solver.D).toBeAllCloseTo(D_test, {rtol: Number.EPSILON*5, atol: 0});

      solver.computeNewton();

      expect(solver.D).toBeAllCloseTo(D_test, {rtol: Number.EPSILON*5, atol: 0});

      f = f.reshape(-1,1);
      const x = new NDArray(x_shape, solver.newton_dX),
           Jx = matmul2(J,x);

      // Newton point has to satisfy normal equaltion Jᵀ(Jx + f) = 0
      expect( matmul2(J.T, zip_elems([Jx,f], (x,y) => x+y) ) ).toBeAllCloseTo(0, {atol: 1e-7});

      const D = new NDArray(x_shape, solver.D),
           JD = zip_elems([J, D.T], 'float64', (j,d) => d===0 ? j : j/d),
            Dx= zip_elems([x, D  ], 'float64', (x,d) => x*d);

      const neg_f = f.mapElems(x => -x);

      expect(Dx).toBeAllCloseTo( svd_lstsq(...svd_jac_2sided(JD),neg_f) );
    }
  })


  forEachItemIn(
    computeNewtonRegularized0_overdet_gen()
  ).it(`computeNewtonRegularized(0) works for pre-generated  over-determined examples`, ([f,J, R,DR]) => {

    const fJ = function(){
      const _f = f.mapElems(),
            _J = J.mapElems();

      return x => {
        return [_f,_J];
      };
    }();

    const [M,N] = J.shape,                     x0 = tabulate([N], 'float64', () => Math.random()*4-2),
         solver = new TrustRegionSolverLSQ(fJ, x0);

    for( let i=4; i-- > 0; )
    {
      for( let l=4; l-- > 0; )
      {
        solver.computeNewton();

        expect( solver.scaledNorm(solver.newton_dX) ).toBeAllCloseTo(R);

        for( let j=4; j-- > 0; )
        {
          for( let k=4; k-- > 0; )
          {
            const [r,dr] = solver.computeNewtonRegularized(0);

            expect( r).toBeAllCloseTo( R);
            expect(dr).toBeAllCloseTo(DR);

            expect( solver.scaledNorm(solver.regularized_dX) ).toBeAllCloseTo(R);
          }
          // solver.computeNewtonRegularized( Math.random() );
        }
      }

      const dx = Float64Array.from({length: N}, () => Math.random()*4-2);
      solver.considerMove(dx);
      solver.makeConsideredMove();
    }
  })


  forEachItemIn(
    computeNewtonRegularized0_underdet_gen()
  ).it(`computeNewtonRegularized(0) works for pre-generated under-determined examples`, ([f,J, R,DR]) => {

    const fJ = function(){
      const _f = f.mapElems(),
            _J = J.mapElems();

      return x => {
        return [_f,_J];
      };
    }();

    const [M,N] = J.shape,                     x0 = tabulate([N], 'float64', () => Math.random()*4-2),
         solver = new TrustRegionSolverLSQ(fJ, x0);

    for( let i=4; i-- > 0; )
    {
      for( let l=4; l-- > 0; )
      {
        solver.computeNewton();

        expect( solver.scaledNorm(solver.newton_dX) ).toBeAllCloseTo(R);

        for( let j=4; j-- > 0; )
        {
          for( let k=4; k-- > 0; )
          {
            const [r,dr] = solver.computeNewtonRegularized(0);

            expect( r).toBeAllCloseTo( R);
            expect(dr).toBeAllCloseTo(DR);

            expect( solver.scaledNorm(solver.regularized_dX) ).toBeAllCloseTo(R);
          }
          // solver.computeMinRegularized( Math.random() );
        }
      }

      const dx = Float64Array.from({length: N}, () => Math.random()*4-2);
      solver.considerMove(dx);
      solver.makeConsideredMove();
    }
  })


  forEachItemIn(
    function*(){
      function* shapes() {
        for( let N=0  ; N++ < 4; )
        for( let M=N-1; M++ < 4; )
          yield [M,N];

        for( let run=96; run-- > 0; ) {
          const N = randInt(1,16),
                M = randInt(N,16);
          yield[M,N];
        }
      }

      for( const [M,N] of shapes() )
      {
        function* samples(){
          for( let run=randInt(1,16); run-- > 0; )
          {
            const dx = tabulate([  N],'float64', () => Math.random()*8-4),
                   f = tabulate([M  ],'float64', () => Math.random()*8-4),
                   J = tabulate([M,N],'float64', () => Math.random()*8-4);

            Object.freeze(dx);
            Object.freeze( f);
            Object.freeze( J);
            Object.freeze(dx.data.buffer);
            Object.freeze( f.data.buffer);
            Object.freeze( J.data.buffer);

            function* lambdas() {
              for( let trial=randInt(1,16); trial-- > 0; ) {
                yield Math.random()*8 + 0.1;
              }
            }

            yield [dx, f,J, lambdas()];
          }
        }

        yield [M,N, samples()];
      }
    }()
  ).it(`computeNewtonRegularized(λ) works for random  over-determined examples`, ([M,N, samples]) => {
    expect(M).toBeGreaterThanOrEqual(N);

    let solver, lambdas, dx,f,J, x0 = 0;

    const g = num_grad( λ => {
      const [r,dr] = solver.computeNewtonRegularized(λ);
      return r;
    });

    let nCalls = 0;
    function fJ( x )
    {
      expect(nCalls++).toBe(0);
      expect(x).toBeAllCloseTo(x0);
      return [f,J];
    }

    const x_shape = Int32Array.of(N,1),
          f_shape = Int32Array.of(M,1);

    let iter=0
    for( [dx,f,J, lambdas] of samples )
    {
      nCalls = 0;

      x0 = zip_elems([x0,dx], (x,dx) => x+dx);

      if( undefined === solver ) {
        solver = new TrustRegionSolverLSQ(fJ, x0);
        expect(iter).toBe(0);
      }
      else {
        solver.considerMove(dx.data);
        solver.makeConsideredMove();
        expect(iter).toBeGreaterThan(0);
      }

      expect(nCalls).toBe(1);

      ++iter;
      expect(solver.M).toBe(M);
      expect(solver.N).toBe(N);
      expect( solver.X0 ).toBeAllCloseTo( x0.data );
      expect( solver.F0 ).toBeAllCloseTo(  f.data );
      expect( solver.J0 ).toBeAllCloseTo(  J.data );

      const neg_f = new NDArray(f_shape, f.data.map(x => -x));

      const X = svd_lstsq(...svd_jac_2sided(J), neg_f);
      Object.freeze(X.data.buffer);
      Object.freeze(X);

      for( let l=2; l-- > 0; )
      {
        solver.computeNewton();

        const x = new NDArray(x_shape, solver.newton_dX.slice());

        expect(x).toBeAllCloseTo(X);

        for( const λ of lambdas )
        {
          const [r,dr] = solver.computeNewtonRegularized(λ),
                 λSqrt = Math.sqrt(λ);

          const x = new NDArray(x_shape, solver.regularized_dX.slice()),
               λD = diag_mat(solver.D).mapElems(d => d*λSqrt),
                B = concat([J,λD]),
                z = concat([neg_f,tabulate([N,1], 'float64', () => 0)]);

          expect(x).toBeAllCloseTo( svd_lstsq(...svd_jac_2sided(B), z) );

          expect(dr).toBeAllCloseTo( g(λ) );
        }
      }
    }
  })


  forEachItemIn(
    function*(){
      function* shapes() {
        for( let M=0  ; M++ < 4; )
        for( let N=M-1; N++ < 4; )
          yield [M,N];

        for( let run=73; run-- > 0; ) {
          const M = randInt(1,16),
                N = randInt(M,16);
          yield[M,N];
        }
      }

      for( const [M,N] of shapes() )
      {
        function* samples(){
          for( let run=randInt(1,16); run-- > 0; )
          {
            const dx = tabulate([  N],'float64', () => Math.random()*8-4),
                   f = tabulate([M  ],'float64', () => Math.random()*8-4),
                   J = tabulate([M,N],'float64', () => Math.random()*8-4);

            Object.freeze(dx);
            Object.freeze( f);
            Object.freeze( J);
            Object.freeze(dx.data.buffer);
            Object.freeze( f.data.buffer);
            Object.freeze( J.data.buffer);

            function* lambdas() {
              for( let trial=randInt(1,16); trial-- > 0; ) {
                yield Math.random()*8 + 0.1;
              }
            }

            yield [dx, f,J, lambdas()];
          }
        }

        yield [M,N, samples()];
      }
    }()
  ).it(`computeNewtonRegularized(λ) works for random under-determined examples`, ([M,N, samples]) => {
    expect(M).toBeLessThanOrEqual(N);

    let solver, lambdas, dx,f,J, x0 = 0;

    const g = num_grad( λ => {
      const [r,dr] = solver.computeNewtonRegularized(λ);
      return r;
    });

    let nCalls = 0;
    function fJ( x )
    {
      expect(nCalls++).toBe(0);
      expect(x).toBeAllCloseTo(x0);
      return [f,J];
    }

    const x_shape = Int32Array.of(N,1),
          f_shape = Int32Array.of(M,1);

    let iter=0
    for( [dx,f,J, lambdas] of samples )
    {
      nCalls = 0;

      x0 = zip_elems([x0,dx], (x,dx) => x+dx);

      if( undefined === solver ) {
        solver = new TrustRegionSolverLSQ(fJ, x0);
        expect(iter).toBe(0);
      }
      else {
        solver.considerMove(dx.data);
        solver.makeConsideredMove();
        expect(iter).toBeGreaterThan(0);
      }

      expect(nCalls).toBe(1);

      ++iter;
      expect(solver.M).toBe(M);
      expect(solver.N).toBe(N);
      expect( solver.X0 ).toBeAllCloseTo( x0.data );
      expect( solver.F0 ).toBeAllCloseTo(  f.data );
      expect( solver.J0 ).toBeAllCloseTo(  J.data );

      const neg_f = new NDArray(f_shape, f.data.map(x => -x));

      for( let l=2; l-- > 0; )
      {
        solver.computeNewton();

        const x = new NDArray(x_shape, solver.newton_dX.slice()),
              D = new NDArray(x_shape, solver.D.slice()),
             JD = zip_elems([J,D.T], (j,d) => d===0 ? j : j/d),
             dx = zip_elems([x,D  ], (x,d) => x*d),
             DX = svd_lstsq(...svd_jac_2sided(JD), neg_f);

        expect(dx).toBeAllCloseTo(DX);

        for( const λ of lambdas )
        {
          const [r,dr] = solver.computeNewtonRegularized(λ),
                 λSqrt = Math.sqrt(λ);

          if( isNaN( r) ) throw new Error('Assertion failed.');
          if( isNaN(dr) ) throw new Error('Assertion failed.');

          const x = new NDArray(x_shape, solver.regularized_dX.slice()),
                D = new NDArray(x_shape, solver.D.slice()),
               JD = zip_elems([J,D.T], (j,d) => d===0 ? j : j/d),
               dx = zip_elems([x,D  ], (x,d) => x*d),
                I = tabulate([N,N], (i,j) => (i===j)*λSqrt);

          const B = concat([JD,I]),
                z = concat([neg_f,tabulate([N,1], 'float64', () => 0)]);

          const DX = svd_lstsq(...svd_jac_2sided(B), z);

          expect(dx).toBeAllCloseTo(DX);

          expect(dr).toBeAllCloseTo( g(λ) );//, {rtol: 1e-3, atol: 1e-6} );
        }
      }
    }
  })


  forEachItemIn(
    function*(){
      function* shapes() {
        for( let M=0; M++ < 4; )
        for( let N=0; N++ < 4; )
          yield [M,N];

        for( let run=73; run-- > 0; ) {
          const M = randInt(1,16),
                N = randInt(M,16);
          yield[M,N];
        }
      }

      for( const [M,N] of shapes() )
      {
        function* samples(){
          for( let run=randInt(1,16); run-- > 0; )
          {
            const dx =     tabulate([  N],'float64', () => Math.random()*8-4),
                   f =     tabulate([M  ],'float64', () => Math.random()*8-4),
                  [J]= _rand_rankdef(M,N);

            Object.freeze(dx);
            Object.freeze( f);
            Object.freeze( J);
            Object.freeze(dx.data.buffer);
            Object.freeze( f.data.buffer);
            Object.freeze( J.data.buffer);

            function* lambdas() {
              for( let trial=randInt(1,16); trial-- > 0; ) {
                yield Math.random()*8 + 0.1;
              }
            }

            yield [dx, f,J, lambdas()];
          }
        }

        yield [M,N, samples()];
      }
    }()
  ).it(`computeNewtonRegularized(λ) works for random  rank-deficient  examples`, ([M,N, samples]) => {

    let solver, lambdas, dx,f,J, x0 = 0;

    const g = num_grad( λ => {
      const [r,dr] = solver.computeNewtonRegularized(λ);
      return r;
    });

    let nCalls = 0;
    function fJ( x )
    {
      expect(nCalls++).toBe(0);
      expect(x).toBeAllCloseTo(x0);
      return [f,J];
    }

    const x_shape = Int32Array.of(N,1),
          f_shape = Int32Array.of(M,1);

    let iter=0
    for( [dx,f,J, lambdas] of samples )
    {
      nCalls = 0;

      x0 = zip_elems([x0,dx], (x,dx) => x+dx);

      if( undefined === solver ) {
        solver = new TrustRegionSolverLSQ(fJ, x0);
        expect(iter).toBe(0);
      }
      else {
        solver.considerMove(dx.data);
        solver.makeConsideredMove();
        expect(iter).toBeGreaterThan(0);
      }

      expect(nCalls).toBe(1);

      ++iter;
      expect(solver.M).toBe(M);
      expect(solver.N).toBe(N);
      expect( solver.X0 ).toBeAllCloseTo( x0.data );
      expect( solver.F0 ).toBeAllCloseTo(  f.data );
      expect( solver.J0 ).toBeAllCloseTo(  J.data );

      const neg_f = new NDArray(f_shape, f.data.map(x => -x));

      for( let l=2; l-- > 0; )
      {
        solver.computeNewton();

        const x = new NDArray(x_shape, solver.newton_dX.slice()),
              D = new NDArray(x_shape, solver.D.slice()),
             JD = zip_elems([J,D.T], (j,d) => d===0 ? j : j/d),
             dx = zip_elems([x,D  ], (x,d) => x*d),
             DX = svd_lstsq(...svd_jac_2sided(JD), neg_f);

        expect(dx).toBeAllCloseTo(DX);

        for( const λ of lambdas )
        {
          const [r,dr] = solver.computeNewtonRegularized(λ),
                 λSqrt = Math.sqrt(λ);

          if( isNaN( r) ) throw new Error('Assertion failed.');
          if( isNaN(dr) ) throw new Error('Assertion failed.');

          const x = new NDArray(x_shape, solver.regularized_dX.slice()),
                D = new NDArray(x_shape, solver.D.slice()),
                I = tabulate([N,N], (i,j) => (i===j)*λSqrt),
               JD = zip_elems([J,D.T], (j,d) => d===0 ? j : j/d),
               dx = zip_elems([x,D  ], (x,d) => x*d);

          const B = concat([JD,I]),
                z = concat([neg_f,tabulate([N,1], 'float64', () => 0)]);

          const DX = svd_lstsq(...svd_jac_2sided(B), z);

          expect(dx).toBeAllCloseTo(DX);

          expect(dr).toBeAllCloseTo( g(λ) );//, {rtol: 1e-3, atol: 1e-6} );
        }
      }
    }
  })


  forEachItemIn(
    function*(){
      function* shapes() {
        for( let M=0; M++ < 4; )
        for( let N=0; N++ < 4; )
          yield [M,N];

        for( let run=128; run-- > 0; ) {
          const M = randInt(1,16),
                N = randInt(M,16);
          yield[M,N];
        }
      }

      for( const [M,N] of shapes() )
      {
        function* samples(){
          for( let run=randInt(1,14); run-- > 0; )
          {
            const dx =   tabulate([  N],'float64', () => Math.random()*8-4),
                   f =   tabulate([M  ],'float64', () => Math.random()*8-4),
                   J = _rand_cols0(M,N);

            Object.freeze(dx);
            Object.freeze( f);
            Object.freeze( J);
            Object.freeze(dx.data.buffer);
            Object.freeze( f.data.buffer);
            Object.freeze( J.data.buffer);

            function* lambdas() {
              for( let trial=randInt(1,16); trial-- > 0; ) {
                yield Math.random()*8 + 0.1;
              }
            }

            yield [dx, f,J, lambdas()];
          }
        }

        yield [M,N, samples()];
      }
    }()
  ).it('computeNewtonRegularized(λ) works for random examples with zero columns in J', ([M,N, samples]) => {

    let solver, lambdas, dx,f,J, x0 = 0;

    const g = num_grad( λ => {
      const [r,dr] = solver.computeNewtonRegularized(λ);
      return r;
    });

    let nCalls = 0;
    function fJ( x )
    {
      expect(nCalls++).toBe(0);
      expect(x).toBeAllCloseTo(x0);
      return [f,J];
    }

    const x_shape = Int32Array.of(N,1),
          f_shape = Int32Array.of(M,1);

    let iter=0
    for( [dx,f,J, lambdas] of samples )
    {
      nCalls = 0;

      x0 = zip_elems([x0,dx], (x,dx) => x+dx);

      if( undefined === solver ) {
        solver = new TrustRegionSolverLSQ(fJ, x0);
        expect(iter).toBe(0);
      }
      else {
        solver.considerMove(dx.data);
        solver.makeConsideredMove();
        expect(iter).toBeGreaterThan(0);
      }

      expect(nCalls).toBe(1);

      ++iter;
      expect(solver.M).toBe(M);
      expect(solver.N).toBe(N);
      expect( solver.X0 ).toBeAllCloseTo( x0.data );
      expect( solver.F0 ).toBeAllCloseTo(  f.data );
      expect( solver.J0 ).toBeAllCloseTo(  J.data );

      const neg_f = new NDArray(f_shape, f.data.map(x => -x));

      for( let l=2; l-- > 0; )
      {
        solver.computeNewton();

        const x = new NDArray(x_shape, solver.newton_dX.slice()),
              D = new NDArray(x_shape, solver.D.slice()),
             JD = zip_elems([J,D.T], (j,d) => d===0 ? j : j/d),
             dx = zip_elems([x,D  ], (x,d) => x*d),
             DX = svd_lstsq(...svd_jac_2sided(JD), neg_f);

        expect(dx).toBeAllCloseTo(DX);

        for( const λ of lambdas )
        {
          const [r,dr] = solver.computeNewtonRegularized(λ),
                 λSqrt = Math.sqrt(λ);

          if( isNaN( r) ) throw new Error('Assertion failed.');
          if( isNaN(dr) ) throw new Error('Assertion failed.');

          const x = new NDArray(x_shape, solver.regularized_dX.slice()),
                D = new NDArray(x_shape, solver.D.slice()),
                I = tabulate([N,N], (i,j) => (i===j)*λSqrt),
               JD = zip_elems([J,D.T], (j,d) => d===0 ? j : j/d),
               dx = zip_elems([x,D  ], (x,d) => x*d);

          const B = concat([JD,I]),
                z = concat([neg_f,tabulate([N,1], 'float64', () => 0)]);

          const DX = svd_lstsq(...svd_jac_2sided(B), z);

          expect(dx).toBeAllCloseTo(DX);

          expect(dr).toBeAllCloseTo( g(λ) );//, {rtol: 1e-3, atol: 1e-6} );
        }
      }
    }
  })


  function test_x0f0J0g0( test_fn, x_range )
  {
    const {nIn: N, nOut: M} = test_fn;

    const fJ = x => [
      test_fn.lsq(x),
      test_fn.lsq_jac(x)
    ];

    const data_gens = {
      [`[X0,F0,J0,G0] and the report are consistent with ${test_fn.name} at generated sampling points`]: function(){
        const N = Math.ceil( 2**(11.1/test_fn.nIn) );

        return cartesian_prod(
          ...x_range.map( r => linspace(...r,N) )
        )
      },
      [`[X0,F0,J0,G0] and the report are consistent with ${test_fn.name} at random sampling points`]: function*(){
        for( let run=0; run++ < 3*1337; )
          yield x_range.map( ([lo,hi]) => {
            const s = Math.random();
            return lo*(1-s) + s*hi;
          })
      }
    };

    for( const [title, data_gen_fn] of Object.entries(data_gens) )
      forEachItemIn(
        function*(){
          for( let repeat=0; ++repeat < 4; )
          {
            let x0 = new Float64Array(N);

            const solver = new TrustRegionSolverLSQ(fJ, x0);

            for( const x of data_gen_fn() )
            {
              const {data: dx} = zip_elems([x,x0], (x,x0) => x-x0);

              solver.considerMove(dx);

              expect( solver.report_x         ).toBeAllCloseTo(                 x  );
              expect( solver.report_f         ).toBeAllCloseTo( test_fn.lsq    (x) );
              expect( solver.report_J         ).toBeAllCloseTo( test_fn.lsq_jac(x) );
              expect( solver.report_loss      ).toBeAllCloseTo( test_fn        (x) / M );
              expect( solver.report_loss_grad ).toBeAllCloseTo( test_fn.grad   (x).mapElems(x => x/M) );

              solver.makeConsideredMove();

              x0.set(solver.X0);
              yield [solver,x0];
            }
          }
        }()
      ).it(title, ([solver,x]) => {
        const J0 = new NDArray(Int32Array.of(M,N), solver.J0);
        expect(solver.X0).toBeAllCloseTo(x);
        expect(solver.F0).toBeAllCloseTo( test_fn.lsq    (x) );
        expect(       J0).toBeAllCloseTo( test_fn.lsq_jac(x) );
        expect(solver.G0).toBeAllCloseTo( test_fn.grad   (x).mapElems(x => x/2) );
      });
  }


  test_x0f0J0g0( beale, [[+0.1, +5.0],
                         [-0.5, +1.0]] );

  test_x0f0J0g0( brown_badscale, [[1e+6 - 2e+5, 1e+6 + 2e+5],
                                  [2e-6 - 4e-7, 2e-6 + 4e-7]] ); // <- TODO: increase range of starting points once line search improved

  test_x0f0J0g0( freudenstein_roth, [[-Math.PI/2, +16],
                                     [-Math.PI/2,  +8]] );

  test_x0f0J0g0( helical_valley, [[-Math.PI/2, +2],
                                  [-Math.PI/2, +2],
                                  [-Math.PI/2, +2]] );

  test_x0f0J0g0( new JennrichSampson(10), [[-1, 0.35],
                                           [-1, 0.35]] );

  test_x0f0J0g0( powell_badscale, [[-10.1, +10.0],
                                   [-10.0, +10.1]] );

  for( const length of range(1,4) )
    test_x0f0J0g0(
      new Rastrigin(length), Array.from({length}, () => [-Math.PI*11,+Math.PI*11])
    );

  for( const length of range(2,4) )
    test_x0f0J0g0(
      new Rosenbrock(length), Array.from({length}, () => [-Math.PI*3,+Math.PI*3])
    );


  forEachItemIn(
    function*(){
      function* shapes() {
        for( let M=8; M > 1; M-- )
        for( let N=8; N > 1; N-- )
          yield [M,N];

        for( let run=128; run-- > 0; )
        {
          const M = randInt(1,32),
                N = randInt(1,32); // <- TODO remove after testing
          yield [M,N];
        }
      }

      for( const [M,N] of shapes() )
      {
        function* samples()
        {
          for( let run=randInt(1,32); run-- > 0; )
          {
            const dx = tabulate([  N],'float64', () => Math.random()*4-2),
                   f = tabulate([M  ],'float64', () => Math.random()*4-2),
                   J = tabulate([M,N],'float64', () => Math.random()*4-2);

            Object.freeze(dx);
            Object.freeze(f);
            Object.freeze(J);
            Object.freeze(dx.data.buffer);
            Object.freeze(f.data.buffer);
            Object.freeze(J.data.buffer);

            yield [dx,f,J];
          }
        }

        yield [M,N, samples()];
      }
    }()
  ).it('[X0,F0,J0,G0] are computed correctly', ([M,N, samples]) => { // <- TODO this could also be well tested using the different test_fn

    let solver, dx,f,J, x0 = 0, zero = new Float64Array(N);

    const fJ = x => [f,J];

    const g = num_grad( x => {
      if( x.ndim     !== 1 ) throw new Error('Assertion failed.');
      if( x.shape[0] !== N ) throw new Error('Assertion failed.');

      const [predict_loss] = solver.considerMove(x.data);
      return predict_loss*M/2;
    });

    let iter=0;
    for( [dx,f,J] of samples )
    {
      x0 = zip_elems([x0,dx], (x,dx) => x+dx);

      if( undefined === solver ) {
        solver = new TrustRegionSolverLSQ(fJ, x0);
        expect(iter).toBe(0);
      }
      else {
        solver.considerMove(dx.data);
        solver.makeConsideredMove();
        expect(iter).toBeGreaterThan(0);
      }

      ++iter;
      expect(solver.M).toBe(M);
      expect(solver.N).toBe(N);

      expect( solver.X0 ).toBeAllCloseTo( x0.data );
      expect( solver.F0 ).toBeAllCloseTo(  f.data );
      expect( solver.J0 ).toBeAllCloseTo(  J.data );
      expect( solver.G0 ).toBeAllCloseTo( g(zero) );
    }
  })


  forEachItemIn(
    function*(){
      function* shapes() {
        for( let M=8; M > 1; M-- )
        for( let N=8; N > 1; N-- )
          yield [M,N];

        for( let run=337; run-- > 0; )
        {
          const M = randInt(1,48),
                N = randInt(1,48); // <- TODO remove after testing
          yield [M,N];
        }
      }

      for( const [M,N] of shapes() )
      {
        function* samples()
        {
          for( let run=randInt(1,32); run-- > 0; )
          {
            const dx = tabulate([  N],'float64', () => Math.random()*4-2),
                   f = tabulate([M  ],'float64', () => Math.random()*4-2),
                   J = tabulate([M,N],'float64', () => Math.random()*4-2);

            Object.freeze(dx);
            Object.freeze(f);
            Object.freeze(J);
            Object.freeze(dx.data.buffer);
            Object.freeze(f.data.buffer);
            Object.freeze(J.data.buffer);

            yield [dx,f,J];
          }
        }

        yield [M,N, samples()];
      }
    }()
  ).it(`cauchyTravel() works for random examples`, ([M,N, samples]) => { // <- TODO this could also be well tested using the different test_fn

    let solver, dx,f,J, x0 = 0, zero = new Float64Array(N);

    const fJ = x => [f,J];

    const g = num_grad(c => { // <- linearized approximation
      const                                      dx = solver.G0.map(g => g*c),
            [predict_loss] = solver.considerMove(dx);
      return predict_loss;
    });

    let iter=0;
    for( [dx,f,J] of samples )
    {
      x0 = zip_elems([x0,dx], (x,dx) => x+dx);

      if( undefined === solver ) {
        solver = new TrustRegionSolverLSQ(fJ, x0);
        expect(iter).toBe(0);
      }
      else {
        solver.considerMove(dx.data);
        solver.makeConsideredMove();
        expect(iter).toBeGreaterThan(0);
      }

      ++iter;
      expect(solver.M).toBe(M);
      expect(solver.N).toBe(N);

      expect( solver.X0 ).toBeAllCloseTo( x0.data );
      expect( solver.F0 ).toBeAllCloseTo(  f.data );
      expect( solver.J0 ).toBeAllCloseTo(  J.data );

      const  c = solver.cauchyTravel();
      expect(c).not.toBeGreaterThan(0);

      if( isFinite(c) ) {
        expect(  g(c) ).toBeAllCloseTo(0);
      }
      else {
        expect( g(+1) ).toBeAllCloseTo( g(0) );
        expect( g(-1) ).toBeAllCloseTo( g(0) );
      }
    }
  })
})
