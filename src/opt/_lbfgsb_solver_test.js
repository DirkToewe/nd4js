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
import {NDArray} from '../nd_array'
import {tabulate} from '../tabulate'
import {_shuffle} from '../_test_data_generators'

import {heap_sort_gen} from '../arrays/heap_sort_gen'

import {diag_mat} from '../la/diag'
import {matmul,
        matmul2} from '../la/matmul'
import {norm} from '../la/norm'
import {solve} from '../la/solve'

import {LBFGSB_Solver       } from './_lbfgsb_solver'
import { LBFGS_SolverRefImpl,
        _rand_updates_v2    } from './_lbfgs_solver_test_utils'


// References
// ----------
// .. [1] "REPRESENTATIONS OF QUASI-NEWTON MATRICES AND THEIR USE IN LIMITED MEMORY METHODS"
//         by Richard H. Byrd, Jorge Nocedal and Robert B. Schnabel
//         Technical Report NAM-03, June 1992; revised January 21, 1996
//         http://users.iems.northwestern.edu/~nocedal/PDFfiles/representations.pdf
//
// .. [2] "A LIMITED MEMORY ALGORITHM FOR BOUND CONSTRAINED OPTIMIZATION"
//         by Richard H. Byrd, Peihuang Lu, Jorge Nocedal and Ciyou Zhu
//         Technical Report NAM-08; Revised: May 1994
//         http://users.iems.northwestern.edu/~nocedal/PDFfiles/limited.pdf


/** Computes the dot product of two arrays.
 */
const dot = (u,v) => {
  if( u.length%1 !== 0       ) throw new Error('Assertion failed.');
  if( v.length%1 !== 0       ) throw new Error('Assertion failed.');
  if( u.length   !== v.length) throw new Error('Assertion failed.');

  let uv = 0;
  for( let i=u.length; i-- > 0; )
    uv += u[i]*v[i];

  return uv;
};


/** A slow reference implementation of the generalized Cauchy
 *  point computation (see [2]).
 */
function* generalized_cauchy(G,B, X, bounds, indices, complete=true)
{
  const N = X.length;

  if(      B.ndim     !==2  ) throw new Error('Assertion failed.');
  if(      B.shape[0] !==  N) throw new Error('Assertion failed.');
  if(      B.shape[1] !==  N) throw new Error('Assertion failed.');
  if(      G.length   !==  N) throw new Error('Assertion failed.');
  if( bounds.length   !==2*N) throw new Error('Assertion failed.');
  if(indices.length   !==  N) throw new Error('Assertion failed.');
  if( ! X.every((x,i) => bounds[2*i+0] <= x) ) throw new Error('Assertion failed.');
  if( ! X.every((x,i) => bounds[2*i+1] >= x) ) throw new Error('Assertion failed.');

  if( ! indices.slice().sort((i,j) => i-j).every((i,I) => i===I) )
    throw new Error('Assertion failed.');

  for( let  i=N; i-- > 0; ) {
    if( !(X[i] >= bounds[2*i+0]) ) throw new Error('Assertion failed');
    if( !(X[i] <= bounds[2*i+1]) ) throw new Error('Assertion failed');
  }

  const dX = Float64Array.from(X, () => 0),
        G0 = G.slice(), // <- memoize initial gradient
         D = G.slice(); // <- current cauchy direction

  const shape = Int32Array.of(N,1),
           dx = new NDArray(shape, dX),
           d  = new NDArray(shape, D);

  // travel along the gradient at which the bound is hit for each dimension
  const travels = Float64Array.from(X, (x,i) => {
    if( 0 === G0[i] )
      return Infinity;
    const          end = bounds[2*i + (G0[i]<0)];
    return (X[i] - end) / G0[i];
  });

  let n_fix = 0,
         df = 0;

  const order = heap_sort_gen(indices, (i,j) => {
    const ti = travels[i],
          tj = travels[j];
    return (ti > tj) - (ti < tj);
  });

  // work off the dimension in the order in which the bounds are hit
  loop:for( const i of order )
    if( 0 !== G0[i] )
    { if( !(X[i] >= bounds[2*i+0]) ) throw new Error('Assertion failed');
      if( !(X[i] <= bounds[2*i+1]) ) throw new Error('Assertion failed');

      // travel along remaining gradient at which we hit the (i+1)-th bounds
      const         end = bounds[2*i + (G0[i]<0)],
        t = (X[i] - end) / G0[i];

      if( !(0 <= t) )
        throw new Error('Assertion failed.');

      if( 0===t ) {
        D[i] = 0;
        ++n_fix;
        continue loop;
      }

      yield [n_fix, X.slice(), df];

      const fp = dot(G,D), // <- univariate polynomial approximation along remaining gradient
            fpp= matmul(d.T, B, d).data[0];

      // B is supposed to be positive definite
      if( !(0 <= fpp) )
        throw new Error('Assertion failed.');

      // Cauchy point along the remaining gradient
      const cp = Math.min( t, fp / (2*fpp) ); // <- only consider Cauchy point in direction of remaining gradient

      if( !complete && !(cp >= t) )
        break loop;

      if( 0 < cp )
      {
        for( let i=N; i-- > 0; )
        {  X[i] -= cp * D[i];
          dX[i] -= cp * D[i];
          const [lo,hi] = bounds.subarray(2*i,2*i+2);
          if( X[i] < lo ) { if(X[i] < lo-1e-8) throw new Error('Assertion failed'); X[i]=lo; }
          if( X[i] > hi ) { if(X[i] > hi+1e-8) throw new Error('Assertion failed'); X[i]=hi; }
        }

        // recompute gradient at the current point
        matmul2(B, dx).data.forEach(
          (dg,i) => G[i] = G0[i] + dg*2
        );

        df = matmul(dx.T, B, dx).data[0] + dot(G0,dX);

        if( cp < t && !( Math.abs(dot(D,G)) <= 1e-6) )
          throw new Error('Assertion failed: ' + dot(G0,G));
      }

      if( complete && !(cp >= t) )
        break loop;

      D[i] = 0;
      X[i] = end; // <- we hit the bounds so let's set the X value exactly to it
      ++n_fix;
    }

  yield [n_fix, X.slice(), df];
}


for( let M=0; M++ <  7; )
for( let N=0; N++ < 13; )
  describe(`LBFGSB_Solver(${M},${N.toString().padStart(2)})`, () => {
    beforeEach( () => {
      jasmine.addMatchers(CUSTOM_MATCHERS)
    })


    if( M === N )
      for( const complete of [[],[true],[false]] )
        forEachItemIn(
          function*(rng){
            for( let run=0; run++ < 1337; )
            {
              const G =           Float64Array.from({length: N}, () => rng.uniform( -2,  +2)),
                    D = diag_mat( Float64Array.from({length: N}, () => rng.uniform(0.2, 1.8)) ),
                    Q = rng.ortho(N),
                    B = matmul(Q, D, Q.T),
               bounds = Float64Array.from({length: 2*N}, () => rng.uniform(-2,+2) );

              for( let i=bounds.length; (i-=2) >= 0; )
                if( bounds[i] > bounds[i+1] )
                  [ bounds[i],  bounds[i+1] ] = [ bounds[i+1], bounds[i] ];

              for( let i=bounds.length; (i-=2) >= 0; )
                if( bounds[i] > bounds[i+1] )
                  throw new Error('Assertion failed.');

              const X = Float64Array.from({length: N}, (_,i) => {
                i *= 2;     const [lo,hi] = [bounds[i],bounds[i+1]];
                return rng.uniform(lo,hi)
              });

              Object.freeze(     G.buffer);
              Object.freeze(     X.buffer);
              Object.freeze(bounds.buffer);
              Object.freeze(B.data.buffer);
              Object.freeze(B);

              yield [G,B,X,bounds];
            }
          }
        ).it(`generalized_cauchy(complete=${complete}) reference implementation works given random examples`, ([G,B,X,bounds]) => {
          // THE TESTED POLYNOMIAL IS: F(X) = G.T @ X + X.T @ B @ X
          const                 N = G.length;
          expect(X.length).toBe(N);

          const g = G.slice(),
                x = X.slice();

          const compute_G = x => {
            const dx = new NDArray(
              Int32Array.of(N,1),
              x.map((x,i) => x - X[i])
            );
            return matmul2(B, dx).data.map((dg,i) => G[i] + 2*dg);
          };

          const compute_F = x => {
            const dx = new NDArray(
              Int32Array.of(N,1),
              x.map((x,i) => x - X[i])
            );
            return dot(G,dx.data) + matmul(dx.T, B, dx).data[0];
          };

          const F_init = compute_F(x);

          const     indices = Int32Array.from({length: N}, (_,i) => i), // <- TODO shuffle
                      rest = generalized_cauchy(g,B,x,bounds,indices, ...complete);
          let [_,X0] = rest.next().value;
          expect(X0.length).toBe(N);
          expect(X0).toEqual(x);

          let proj_grad1 = NaN;

          let   n_fix, X1, df;
          for( [n_fix, X1, df] of rest )
          {
            expect(X1.length).toBe(N);
            expect(X1).toEqual(x);

            // FROZEN/FIXED DIMENSION MUST NOT CHANGE ANYMORE
            if( 0 < n_fix )
                expect( Float64Array.from( indices.subarray(0,n_fix-1), i => X1[i] ) )
              .toEqual( Float64Array.from( indices.subarray(0,n_fix-1), i => X0[i] ) );

            // CHANGE IN X MUST ALWAYS BE ALONG REMAINING (INITIAL) GRADIENT
            if( n_fix < N )
            {
              const DX = Float64Array.from( indices.subarray(n_fix), i => X1[i] - X0[i] ),
                    DG = Float64Array.from( indices.subarray(n_fix), i => G[i] );

              const norm_dx = norm(DX),
                    norm_dg = norm(DG);

              if( norm_dx !== 0 )
              {
                for( let i=DX.length; i-- > 0; ) {
                  DX[i] /= -norm_dx;
                  DG[i] /=  norm_dg;
                }

                expect(DX).toBeAllCloseTo(DG);
              }
            }

            // x must always be within bounds
            expect(
              X1.every( (x,i) => bounds[2*i+0] <= x && x <= bounds[2*i+1] )
            ).toBe(true)

            // indices must always be the same full set of integer values
            expect(
              indices.slice().sort((i,j) => i-j).every((i,I) => i===I)
            ).toBe(true);

            // for fixed indices i, X[i] must be touching one of the bounds
            for( const i of indices.subarray(0,n_fix) )
              expect( bounds.subarray(2*i, 2*i+2) ).toContain(X1[i]);

            // g should alway reflect the gradient at current x
            expect(g).toBeAllCloseTo( compute_G(X1) );

            const                   dx = X0.map( (x0,i) => X1[i] - x0 );
                  proj_grad1 = dot( dx, compute_G(X1) );
            const proj_grad0 = dot( dx, compute_G(X0) );

            // traveling along the generalized cauchy path, gradient has to be negative
            // continously decrease until (generalized) Cauchy point is reached.
            expect(proj_grad0).toBeLessThan(1e-8);
            expect(proj_grad1).toBeLessThan(1e-8);
            // as more bounds are hit and we travel down the quadratic hypersurface,
            // the gradient has to become less negative (closer to 0).
            expect(proj_grad1).toBeGreaterThan(proj_grad0 - 1e-8);

            const F0 = compute_F(X0),
                  F1 = compute_F(X1);

            // traveling along the generalized cauchy path, function values have to
            //  continously decrease until (generalized) Cauchy point is reached.
            expect(F1).toBeLessThan(F0 + 1e-8);

            expect(df).toBeAllCloseTo(F1-F_init);

            // set previous x for next iteration
            X0 = X1;
          }

          if( complete[0] !== false && n_fix < N )
            expect(proj_grad1).toBeAllCloseTo(0);
        });


    for( const [suffix, rand_scale] of [
      [''           , rng => () => 1                    ],
      [' and scales', rng => () => rng.uniform(0.5, 1.5)]
    ])
    {
      forEachItemIn(
        function*(rng){
          const rand_scale_rng = rand_scale(rng);
          for( let run=0; run++ < 17; )
            yield function*(){
              let len = rng.int(32,64);

              for( const [dx,dg] of _rand_updates_v2(rng,N) )
                if( --len < 0 )
                  break;
                else
                  yield [
                    dx,dg,
                    rand_scale_rng(),
                    tabulate([N,1], 'float64', () => rng.uniform(-4,+4))
                  ];
            }();
        }
      ).it('compute_Hv works given all free variables (n_fix=0) given random updates' + suffix, updates => {
        const ref = new LBFGS_SolverRefImpl(M,N),
              tst = new LBFGSB_Solver      (M,N);

        expect(M).toBeGreaterThan(0);
        expect(N).toBeGreaterThan(0);

        for( const [dx,dg, scale, y] of updates )
        {
          ref.update(dx.data, dg.data);
          tst.update(dx.data, dg.data);

          ref.scale = scale;
          tst.scale = scale;

          const {B} = ref;

          const x = solve(B,y),
                X = new Float64Array(N);

          tst.compute_Hv(y.data, X);

          const tol = {
            rtol: 0,
            atol: 1e-4 * Math.max(norm(X), norm(x))
          };
          expect(X).toBeAllCloseTo(x.data, tol);
        }
      });


      forEachItemIn(
        function*(rng){
          const rand_scale_rng = rand_scale(rng);
          for( let run=0; run++ < 8; )
            yield function*(){
              let len = rng.int(32,64);

              for( const [dx,dg] of _rand_updates_v2(rng,N) )
              {
                if( --len < 0 )
                  break;

                const         indices = Int32Array.from({length: N}, (_,i) => i);
                  rng.shuffle(indices);
                Object.freeze(indices.buffer);

                yield [
                  dx,dg, indices,
                  rand_scale_rng(),
                  tabulate([N,1], 'float64', () => rng.uniform(-4,+4))
                ];
              }
            }();
        }
      ).it('compute_subspace_Hv works with all free variables (n_fix=0) given random updates' + suffix, updates => {
        const ref = new LBFGS_SolverRefImpl(M,N),
              tst = new LBFGSB_Solver      (M,N);

        expect(M).toBeGreaterThan(0);
        expect(N).toBeGreaterThan(0);

        for( const [dx,dg, indices, scale, y] of updates )
        {
          ref.update(dx.data, dg.data);
          tst.update(dx.data, dg.data);

          ref.scale = scale;
          tst.scale = scale;

          const {B} = ref,
                 Z  = tabulate([N,N], (i,j) => 1*(i === indices[j]) );

          const b = matmul(Z.T, B, Z),
                x = solve(b,y);

          const X = new Float64Array(N);
          tst.compute_subspace_Hv(matmul2(Z,y).data, X, indices,0);

          const Zx = matmul2(Z,x).reshape(-1),
               tol = {
                 rtol: 0,
                 atol: 1e-4 * Math.max(norm(X), norm(x))
               };
          expect(X).toBeAllCloseTo(Zx, tol);
        }
      });

      forEachItemIn(
        function*(rng){
          const rand_scale_rng = rand_scale(rng);
          for( let run=0; run++ < 17; )
            yield function*(){
              let len = rng.int(32,64);

              for( const [dx,dg] of _rand_updates_v2(rng,N) )
              {
                if( --len < 0 )
                  break;

                const         indices = Int32Array.from({length: N}, (_,i) => i);
                  rng.shuffle(indices);
                Object.freeze(indices.buffer);

                const n_fix = rng.int(0,N);

                yield [
                  dx,dg, indices, n_fix,
                  rand_scale_rng(),
                  tabulate([N-n_fix,1], 'float64', () => rng.uniform(-4,+4))
                ];
              }
            }()
        }
      ).it('compute_subspace_Hv works given random updates' + suffix, updates => {
        const ref = new LBFGS_SolverRefImpl(M,N),
              tst = new LBFGSB_Solver      (M,N);

        expect(M).toBeGreaterThan(0);
        expect(N).toBeGreaterThan(0);

        for( const [dx,dg, indices, n_fix, scale, y] of updates )
        {
          ref.update(dx.data, dg.data);
          tst.update(dx.data, dg.data);

          ref.scale = scale;
          tst.scale = scale;

          const {B} = ref,
                 Z  = tabulate([N,N-n_fix], (i,j) => 1*(i === indices[n_fix+j]) );

          const b = matmul(Z.T, B, Z),
                x = solve(b,y);

          const X = new Float64Array(N);
          tst.compute_subspace_Hv(matmul2(Z,y).data, X, indices,n_fix);

          const Zx = matmul2(Z,x).reshape(-1),
              tol = {
                rtol: 0,
                atol: 1e-4 * Math.max(norm(X), norm(x))
              };
          expect(X).toBeAllCloseTo(Zx, tol);
        }
      });

      // FIXME: rework these tests
      if( M===N )
      for( const args of [[],[true],[false]] )
        forEachItemIn(
          function*(rng){
            const rand_scale_rng = rand_scale(rng);

            for( let run=0; run++ < 2048; )
            {
              const   G = Float64Array.from({length: N}, () => rng.uniform( -2, +2) ),
                      d = Float64Array.from({length: N}, () => rng.uniform(0.2,1.8) ),
                      D = diag_mat(d),
                      Q = rng.ortho(N),
                      B = matmul(Q.T, D, Q),
                 bounds = Float64Array.from({length: 2*N}, () => rng.uniform(-2,+2) ),
                updates = [];

              let i=0;
              for( const {data: row} of Q )
              {
                const dxScale = rng.uniform(0.2,1.8),
                      dgScale = d[i];
                i++;
                const dx = row.map(dx => dx*dxScale),
                      dg =  dx.map(dg => dg*dgScale);
                Object.freeze(dx.buffer);
                Object.freeze(dg.buffer);
                updates.push([dx,dg]);
              }

              Object.freeze(updates);
              Object.freeze(B.data.buffer);
              Object.freeze(B);
      
              for( let i=bounds.length; (i-=2) >= 0; )
                if( bounds[i] > bounds[i+1] )
                  [ bounds[i],  bounds[i+1] ] = [ bounds[i+1], bounds[i] ];
      
              for( let i=bounds.length; (i-=2) >= 0; )
                if( bounds[i] > bounds[i+1] )
                  throw new Error('Assertion failed.');
      
              const X = Float64Array.from({length: N}, (_,i) => {
                i *= 2;     const [lo,hi] = [bounds[i],bounds[i+1]];
                return rng.uniform(lo,hi);
              });

              const scale = rand_scale_rng();
      
              Object.freeze(     G.buffer);
              Object.freeze(     X.buffer);
              Object.freeze(bounds.buffer);
              Object.freeze(B.data.buffer);
              Object.freeze(B);

              const indices_ref = Int32Array.from({length: N}, (_,i) => i),
                    indices_tst = Int32Array.from({length: N}, (_,i) => i);
              rng.shuffle(indices_ref);
              rng.shuffle(indices_tst);
      
              yield {G,B,X,bounds,updates,scale,indices_ref,indices_tst};
            }
          }
        ).it(`compute_cauchyGeneralized(complete=${args}) works given random examples${suffix}`, ({G,B,X,bounds,updates,scale,indices_ref,indices_tst}) => {
          const                 N = G.length;
          expect(X.length).toBe(N);
      
          const   G_ref = G.slice(),
                  X_ref = X.slice();
      
          const [n_ref,df_ref] = function(){
            let   n, df;
            for( [n,,df] of generalized_cauchy(
              G_ref,
              B.mapElems(x => x/2), // <- LBFGS B is actually Hessian times two
              X_ref,
              bounds,
              indices_ref,
              ...args
            )); // PASS
            return [n,df]
          }();
          

          const lbfgs = new LBFGSB_Solver(N,N);

          for( const [dx,dg] of updates )
            lbfgs.update(dx,dg);
          lbfgs.scale = scale;

          ;{
            const be = Array.from({length: N}, (_,i) => {
              const bei = new Float64Array(N*2);
              lbfgs.compute_be(i, bei);
              return bei;
            });
            
            const b = tabulate([N,N], (i,j) =>
              lbfgs.compute_ubbv(be[i], 1*(i===j), be[j])
            );
            
            expect(b).toBeAllCloseTo(B);
          };

          const   G_tst = G.slice(),
                  X_tst = X.slice();

          const [n_tst,df_tst] = lbfgs.compute_cauchyGeneralized(G_tst, X_tst, bounds, indices_tst, ...args);

          expect(
            indices_tst.slice().sort((i,j) => i-j).every((i,I) => i===I)
          ).toBe(true);

          expect(n_tst).toBe(n_ref);
          expect(X_tst).toBeAllCloseTo(X_ref);
          expect(G_tst).toBeAllCloseTo(G_ref, {atol: 1e-6});

            expect( indices_tst.subarray(0,n_tst) )
          .toEqual( indices_ref.subarray(0,n_ref) )

          expect(df_tst).toBeAllCloseTo(df_ref);
        });


      forEachItemIn(
        function*(rng){
          const rand_scale_rng = rand_scale(rng);

          for( let run=0; run++ < 8; )
            yield function*(){
              let len = rng.int(32,64);

              for( const [dx,dg] of _rand_updates_v2(rng,N) )
              {
                if( --len < 0 )
                  break;
                const uv = [[ 
                    Float64Array.from({length: N}, () => rng.uniform(-4,+4) ),
                    Float64Array.from({length: N}, () => rng.uniform(-4,+4) ), rand_scale_rng()
                ],[ Float64Array.from({length: N}, () => rng.uniform(-4,+4) ),
                    Float64Array.from({length: N}, () => rng.uniform(-4,+4) ), rand_scale_rng()
                ]];
                Object.freeze(uv[0][0].buffer);
                Object.freeze(uv[0][1].buffer);
                Object.freeze(uv[1][0].buffer);
                Object.freeze(uv[1][1].buffer);
                Object.freeze(uv);
                yield [dx,dg, uv];
              }
            }();
        }
      ).it('compute_bv+compute_ubbv work given random updates' + suffix, updates => {
        const ref = new LBFGS_SolverRefImpl(M,N),
              tst = new LBFGSB_Solver      (M,N);

        expect(M).toBeGreaterThan(0);
        expect(N).toBeGreaterThan(0);

        const shape = Int32Array.of(N,1);

        for( const [dx,dg, uv] of updates )
        {
          ref.update(dx.data, dg.data);
          tst.update(dx.data, dg.data);
          for( const [u,v, scale] of uv )
          {
            ref.scale = scale;
            tst.scale = scale;

            const {B} = ref,
                   bu = new Float64Array(2*M),
                   bv = new Float64Array(2*M);

            const U = new NDArray(shape, u),
                  V = new NDArray(shape, v);

            const UBU = matmul(U.T, B, U),
                  UBV = matmul(U.T, B, V),
                  VBV = matmul(V.T, B, V);
            Object.freeze(UBU); Object.freeze(UBU.data.buffer);
            Object.freeze(UBV); Object.freeze(UBV.data.buffer);
            Object.freeze(VBV); Object.freeze(VBV.data.buffer);
            
            for( let outer=0; outer++ < 2; )
            {
              tst.compute_bv(u, bu);

              const    uBu = tst.compute_ubbv(bu, dot(u,u), bu);
              expect([[uBu]]).toBeAllCloseTo(UBU);

              for( let inner=0; inner++ < 2; )
              {
                tst.compute_bv(v, bv);

                const    uBv = tst.compute_ubbv(bu, dot(u,v), bv);
                expect([[uBv]]).toBeAllCloseTo(UBV);

                const    vBu = tst.compute_ubbv(bv, dot(v,u), bu);
                expect([[vBu]]).toBeAllCloseTo(UBV);

                const    vBv = tst.compute_ubbv(bv, dot(v,v), bv);
                expect([[vBv]]).toBeAllCloseTo(VBV);
              }
            }
          }
        }
      });


      forEachItemIn(
        function*(rng){
          for( let run=0; run++ < 24; )
            yield [
              function*(){
                let len = rng.int(32,64);

                for( const [dx,dg] of _rand_updates_v2(rng,N) )
                {
                  if( --len < 0 )
                    break;

                  const         indices = Int32Array.from({length: N}, (_,i) => i);
                    rng.shuffle(indices);
                  Object.freeze(indices.buffer);

                  const n_fix = rng.int(0,N);

                  yield [
                    dx,dg, indices, n_fix,
                    tabulate([N-n_fix,1], 'float64', () => rng.uniform(-4,+4))
                  ];
                }
              }(),
              rand_scale(rng)
            ];
        }
      ).it('compute_be works given random updates' + suffix, ([updates,rand_scale]) => {
        const tst = new LBFGSB_Solver(M,N);

        expect(M).toBeGreaterThan(0);
        expect(N).toBeGreaterThan(0);

        const shape = Int32Array.of(N,1);

        let m=0;
        for( let [dx,dg] of updates )
        {
          tst.update(dx.data, dg.data);

          m = Math.min(M,m+1);

          const BEI = new Float64Array(2*M),
                bei = new Float64Array(2*M);

          for( let i=0; i < N; i++ )
          {
            const EI = Float64Array.from({length: N}, (_,j) => (i===j)*1);
            Object.freeze(EI.buffer);

            for( let run=2; run-- > 0; )
            {
              for( let k=2*M; k-- > 0; )
              {
                BEI[k] = Math.random()*8-4;
                bei[k] = Math.random()*8-4;
              }
              tst.scale = rand_scale();
              tst.compute_bv(EI, BEI);
              tst.compute_be( i, bei);

                expect( bei.subarray(0,m*2) )
              .toEqual( BEI.subarray(0,m*2) );
            }
          }
        }
      });


      if( M===N )
        forEachItemIn(
          function*(rng){
            const rand_scale_rng = rand_scale(rng);
            for( let run=0; run++ < 1337; )
            {
              const   d = Float64Array.from({length: N}, () => rng.uniform(0.2,1.8) ),
                      D = diag_mat(d),
                      Q = rng.ortho(N),
                      B = matmul(Q.T, D, Q),
                updates = [];

              let i=0;
              for( const {data: row} of Q )
              {
                const dxScale = rng.uniform(0.2,1.8),
                      dgScale = d[i];
                i++;
                const dx = row.map(dx => dx*dxScale),
                      dg =  dx.map(dg => dg*dgScale);
                Object.freeze(dx.buffer);
                Object.freeze(dg.buffer);
                updates.push([dx,dg]);
              }

              Object.freeze(updates);
              Object.freeze(B.data.buffer);
              Object.freeze(B);

              yield [updates, B, rand_scale_rng()];
            }
          }
        ).it('approximates Hessian given random multivariate polynomials' + suffix, ([updates, B, scale]) => {
          const N = B.shape[0];

          const lbfgs = new LBFGSB_Solver(N,N);

          for( const [dx,dg] of updates )
            lbfgs.update(dx,dg);
          lbfgs.scale = scale;

          const be = Array.from({length: N}, (_,i) => {
            const bei = new Float64Array(N*2);
            lbfgs.compute_be(i, bei);
            return bei;
          });

          const b = tabulate([N,N], (i,j) =>
            lbfgs.compute_ubbv(be[i], 1*(i===j), be[j])
          );

          expect(b).toBeAllCloseTo(B);
        });


      forEachItemIn(
        function*(rng){
          const rand_scale_rng = rand_scale(rng);

          for( let run=0; run++ < 7; )
            yield function*(){
              let len = rng.int(32,64);

              for( const [dx,dg] of _rand_updates_v2(rng,N) )
              {
                if( --len < 0 )
                  break;
                const vs = [[ Float64Array.from({length: N}, () => rng.uniform(-4,+4) ), rand_scale_rng() ],
                            [ Float64Array.from({length: N}, () => rng.uniform(-4,+4) ), rand_scale_rng() ]];
                Object.freeze(vs[0][0].buffer);
                Object.freeze(vs[1][0].buffer);
                Object.freeze(vs);
                yield [dx,dg, vs];
              }
            }();
        }
      ).it('compute_Bv work given random updates' + suffix, updates => {
        const ref = new LBFGS_SolverRefImpl(M,N),
              tst = new LBFGSB_Solver      (M,N);

        expect(M).toBeGreaterThan(0);
        expect(N).toBeGreaterThan(0);

        const shape = Int32Array.of(N,1);

        for( let [dx,dg, vs] of updates )
        {
          ref.update(dx.data, dg.data);
          tst.update(dx.data, dg.data);
          for( const [v,scale] of vs )
          {
            ref.scale = scale;
            tst.scale = scale;

            const {B} = ref;

            const V = new NDArray(shape, v .slice()),
                 BV = matmul2(B,V);

            Object.freeze( V ); Object.freeze( V .data.buffer);
            Object.freeze(BV ); Object.freeze(BV .data.buffer);

            for( let outer=0; outer++ < 2; )
            {
              const bv = Float64Array.from({length: N}, () => Math.random()*8-4);
              Object.freeze(bv.buffer);
              tst.compute_Bv(v, bv);

              expect( new NDArray(shape,bv) ).toBeAllCloseTo(BV);
            }
          }
        }
      })
    }
  });
