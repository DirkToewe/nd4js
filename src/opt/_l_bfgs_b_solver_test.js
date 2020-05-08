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
import {asarray, NDArray} from '../nd_array'
import {tabulate} from '../tabulate'
import {_rand_int,
        _shuffle} from '../_test_data_generators'
import {zip_elems} from '../zip_elems'

import {diag_mat} from '../la/diag'
import {eye} from '../la/eye'
import {matmul,
        matmul2} from '../la/matmul'
import {norm} from '../la/norm'
import {rand_ortho} from '../la/rand_ortho'
import {solve} from '../la/solve'

import {_heap_sort} from './_opt_utils'
import {L_BFGS_B_Solver} from './_l_bfgs_b_solver'

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


/** A slow reference implementation of the LBFGS Hessian approximation
 *  used to test the much more efficient `LBFGSSolver`.
 */
class L_BFGS_SolverRef
{
  constructor( M, N )
  {
    if( M%1 !== 0 ) throw new Error('Assertion failed.');
    if( N%1 !== 0 ) throw new Error('Assertion failed.');
    if( !(0 < M) ) throw new Error('Assertion failed.');
    if( !(0 < N) ) throw new Error('Assertion failed.');
    this.M = M|0;
    this.N = N|0;
    this.dx = [];
    this.dg = [];
    this.scale = 1;
  }

  update( dx, dg )
  {
    if( dx.length%1 !== 0 ) throw new Error('Assertion failed.');
    if( dg.length%1 !== 0 ) throw new Error('Assertion failed.');

    const {dx: dX, dg: dG, N} = this;

    dx = Float64Array.from(dx);
    dg = Float64Array.from(dg);

    if( dx.length !== N ) throw new Error('Assertion failed: ' + dx.length);
    if( dg.length !== N ) throw new Error('Assertion failed: ' + dg.length);

    const shape = Int32Array.of(N,1);

    dx = new NDArray(shape, dx);
    dg = new NDArray(shape, dg);
    Object.freeze(dx.buffer);
    Object.freeze(dg.buffer);

    if( dX.length != dG.length ) throw new Error('Assertion failed.');
    dX.push(dx);
    dG.push(dg);
    if( dX.length > this.M ) {
        dX.shift();
        dG.shift();
    }
    if( dX.length > this.M ) throw new Error('Assertion failed.');
  }

  *_history()
  {
    const {dx,dg} = this;
    if( dx.length != dg.length ) throw new Error('Assertion failed.');
    for( let i=0; i < dx.length; i++ )
      yield [dx[i], dg[i]];
  }

  get B()
  {
    const {N} = this;

    let B = eye(N).mapElems(x => x * this.scale);

    for( const [dx,dg] of this._history() )
      B = zip_elems(
        [ B,
          matmul2(dg, dg.T),
          matmul2(dg.T, dx),
          matmul(B,dx, dx.T, B.T),
          matmul(dx.T, B, dx) ],
        (B, gg, gx, BxxB, xBx) => B + gg/gx - BxxB/xBx
      );

    Object.freeze(B.data.buffer);
    Object.freeze(B);
    return B;
  }

  get H()
  {
    const {N} = this;

    const I = eye(N);

    let H = I.mapElems(x => x / this.scale);

    for( const [dx,dg] of this._history() )
    {
      const gx = matmul2(dg.T, dx),
             W = zip_elems(
               [I, matmul2(dx, dg.T), gx],
               (I,xg,gx) => I - xg/gx
             );

      H = zip_elems(
        [ matmul(W, H, W.T),
          matmul2(dx, dx.T),
          gx ],
        (WHW,xx,gx) => WHW + xx/gx
      );
    }

    Object.freeze(H.data.buffer);
    Object.freeze(H);
    return        H;
  }
}


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
function* generalized_cauchy(G,B, X, bounds, indices)
{
  if(      B.ndim     !==2         ) throw new Error('Assertion failed.');
  if(      B.shape[0] !==  X.length) throw new Error('Assertion failed.');
  if(      B.shape[1] !==  X.length) throw new Error('Assertion failed.');
  if(      G.length   !==  X.length) throw new Error('Assertion failed.');
  if( bounds.length   !==2*X.length) throw new Error('Assertion failed.');
  if(indices.length   !==  X.length) throw new Error('Assertion failed.');
  if( ! X.every((x,i) => bounds[2*i+0] <= x) ) throw new Error('Assertion failed.');
  if( ! X.every((x,i) => bounds[2*i+1] >= x) ) throw new Error('Assertion failed.');

  if( ! indices.slice().sort((i,j) => i-j).every((i,I) => i===I) )
    throw new Error('Assertion failed.');

  for( let i=X.length; i-- > 0; )
  {
    if( !(X[i] >= bounds[2*i+0]) ) throw new Error('Assertion failed');
    if( !(X[i] <= bounds[2*i+1]) ) throw new Error('Assertion failed');
  }

  B = asarray(B);
  const dX = Float64Array.from(X, () => 0);

  const shape = Int32Array.of(X.length,1),
            g = new NDArray(shape,  G),
           dx = new NDArray(shape, dX);

  // travel along the gradient at which the bound is hit for each dimension
  const travels = Float64Array.from(X, (x,i) => {
    if( 0 === G[i] )
      return Infinity;
    const          end = bounds[2*i + (G[i]<0)];
    return (X[i] - end) / G[i];
  });

  let n_fix = 0;

  // work off the dimension in the order in which the bounds are hit
  loop:for( const i of _heap_sort(indices, (i,j) => travels[i] < travels[j]) )
    if( 0!==G[i] )
    {
      if( !(X[i] >= bounds[2*i+0]) ) throw new Error('Assertion failed');
      if( !(X[i] <= bounds[2*i+1]) ) throw new Error('Assertion failed');

      // travel along remaining gradient at which the hit the (i+1)-th bounds
      const         end = bounds[2*i + (G[i]<0)],
        t = (X[i] - end) / G[i];

      if( !(0 <= t) )
        throw new Error('Assertion failed.');

      if( 0===t ) {
        G[i] = 0;
        ++n_fix;
        continue loop;
      }

      yield [n_fix, X.slice()];

      const      γ = matmul2(B, dx).data.map((dg,i) => G[i] + dg*2), // <- gradient at the current point
        fp = dot(γ,G), // <- univariate polynomial approximation along remaining gradient
        fpp= matmul(g.T, B, g).data[0];

      // Cauchy point along the remaining gradient
      let cp = fp / (2*fpp); // <- only consider Cauchy point in direction of remaining gradient

      // B is supposed to be positive definite
      if( !(0 <= fpp) )
        throw new Error('Assertion failed.');

      if( cp < t )
      {
        if( 0 < cp )
        {
          for( let i=X.length; i-- > 0; )
          {  X[i] -= cp * G[i];
            dX[i] -= cp * G[i];
            const [lo,hi] = bounds.subarray(2*i,2*i+2);
            if( X[i] < lo ) { if(X[i] < lo-1e-8) throw new Error('Assertion failed'); X[i]=lo; }
            if( X[i] > hi ) { if(X[i] > hi+1e-8) throw new Error('Assertion failed'); X[i]=hi; }
          }
          // at the Cauchy point, the remaining gradient is supposed to be orthogonal to the gradient
          const g = matmul2(B, dx).data.map((dg,i) => G[i] + dg*2)
          if( !( Math.abs(dot(G,g)) <= 1e-6) )
            throw new Error('Assertion failed: ' + dot(G,g));
        }
        break loop;
      }

      for( let i=X.length; i-- > 0; )
      {  X[i] -= t * G[i];
        dX[i] -= t * G[i];
        const [lo,hi] = bounds.subarray(2*i, 2*i+2);
        if( X[i] < lo ) { if(X[i] < lo-1e-8) throw new Error('Assertion failed'); X[i]=lo; }
        if( X[i] > hi ) { if(X[i] > hi+1e-8) throw new Error('Assertion failed'); X[i]=hi; }
      }
      X[i] = end; // <- we hit the bounds so let's set the X value exactly to it
      G[i] = 0; // <- we hit the bounds along this axis so we can move any further along this dimension
      ++n_fix;
    }

  yield [n_fix, X.slice()];
}


/** Produces an infinite series of LBFGS update input pairs.
 *  The are filtered to result in a sufficiently well-conditioned
 *  (positive definite) Hessian approximation.
 */
function* _rand_updates( N )
{
  if( N%1 !== 0 ) throw new Error('Assertion failed.');
  if( ! (0 < N) ) throw new Error('Assertion failed.');

  const shape = Int32Array.of(N,1);

  loop:for(;;)
  {
    const dx = Float64Array.from({length: N}, () => Math.random()*2-1),
          dg = Float64Array.from({length: N}, () => Math.random()*2-1);

    const nx = Math.hypot(...dx),
          ng = Math.hypot(...dg);

    for( let i=N; i-- > 0; ) {
      dx[i] /= nx;
      dg[i] /= ng;
    }

    let dot_xg = dot(dx,dg);
    if( dot_xg < 0 ) {
      for( let i=N; i-- > 0; )
        dx[i] *= -1;
      dot_xg *= -1;
    }
    if( dot_xg < 0.05 )
      continue loop;

    const sx = Math.random()*1 + 0.5,
          sg = Math.random()*1 + 0.5;

    for( let i=N; i-- > 0; ) {
      dx[i] *= sx;
      dg[i] *= sg;
    }

    Object.freeze(dx.buffer);
    Object.freeze(dg.buffer);
    yield Object.freeze([
      new NDArray(shape,dx),
      new NDArray(shape,dg)
    ]);
  }
}


describe('L_BFGS_B_Solver', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  forEachItemIn(
    function*(){
      for( let run=0; run++ < 512; )
      for( let N=0; N++ < 24; )
      {
        const G =           Float64Array.from({length: N}, () => Math.random()*4.0 - 2.0),
              D = diag_mat( Float64Array.from({length: N}, () => Math.random()*1.8 + 0.2) ),
              Q = rand_ortho(N),
              B = matmul(Q, D, Q.T),
         bounds = Float64Array.from({length: 2*N}, () => Math.random()*4.0 - 2.0);

        for( let i=bounds.length; (i-=2) >= 0; )
          if( bounds[i] > bounds[i+1] )
            [ bounds[i],  bounds[i+1] ] = [ bounds[i+1], bounds[i] ];

        for( let i=bounds.length; (i-=2) >= 0; )
          if( bounds[i] > bounds[i+1] )
            throw new Error('Assertion failed.');

        const X = Float64Array.from({length: N}, (_,i) => {
          i *= 2;
          const [lo,hi] = [bounds[i],bounds[i+1]],
                      s = Math.random();
          return lo*(1-s) + s*hi;
        });

        Object.freeze(     G.buffer);
        Object.freeze(     X.buffer);
        Object.freeze(bounds.buffer);
        Object.freeze(B.data.buffer);
        Object.freeze(B);

        yield [G,B,X,bounds];
      }
    }()
  ).it('generalized_cauchy() reference implementation works given random examples', ([G,B,X,bounds]) => {
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
      return dot(G,dx.data) + matmul(dx.T, B, dx).data[0] / 2;
    };

    const     indices = Int32Array.from({length: N}, (_,i) => i), // <- TODO shuffle
                 rest = generalized_cauchy(g,B,x,bounds,indices);
    let [_,X0] = rest.next().value;
    expect(X0.length).toBe(N);
    expect(X0).toEqual(x);

    let proj_grad1 = NaN;

    for( const [n_fix,X1] of rest )
    {
      expect(X1.length).toBe(N);
      expect(X1).toEqual(x);

      expect(
        X1.every( (x,i) => bounds[2*i+0] <= x && x <= bounds[2*i+1] )
      ).toBe(true)

      expect(
        indices.slice().sort((i,j) => i-j).every((i,I) => i===I)
      ).toBe(true);

      for( let i=0    ; i < n_fix; i++ ) expect(g[indices[i]]).toBe(0);
      for( let i=n_fix; i < N    ; i++ ) expect(g[indices[i]]).toBe(G[indices[i]]);

      const                   dir = X0.map( (x0,i) => X1[i] - x0 );
            proj_grad1 = dot( dir, compute_G(X1) );
      const proj_grad0 = dot( dir, compute_G(X0) );

      expect(proj_grad0).toBeLessThan(1e-8);
      expect(proj_grad1).toBeLessThan(1e-8);
      expect(proj_grad1).toBeGreaterThan(proj_grad0 - 1e-8);

      const F0 = compute_F(X0),
            F1 = compute_F(X1);

      expect(F1).toBeLessThan(F0 + 1e-8);

      X0 = X1;
    }

    if( g.some(g => g!==0) )
      expect(proj_grad1).toBeAllCloseTo(0);
  });


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
    ).it('compute_subspace_Hv works with all free variables (n_fix=0) given random updates' + suffix, ([M,N]) => {
      const ref = new L_BFGS_SolverRef(M,N),
            tst = new L_BFGS_B_Solver (M,N);

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
        tst.scale = scale;

        const indices = Int32Array.from({length: N}, (_,i) => i);
        _shuffle(indices);
        Object.freeze(indices.buffer);

        const {B} = ref,
               Z  = tabulate([N,N], (i,j) => 1*(i === indices[j]) ),
               y  = tabulate([N,1], () => Math.random()*8 - 4);

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

        if( ++i >= 64 ) break;
      }
    })

    forEachItemIn(
      function*(){
        for( let run=0; run++ < 2; )
        for( let M=0; M++ <  8; )
        for( let N=0; N++ < 16; )
          yield [M,N];
      }()
    ).it('compute_subspace_Hv works given random updates' + suffix, ([M,N]) => {
      const ref = new L_BFGS_SolverRef(M,N),
            tst = new L_BFGS_B_Solver (M,N);

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
        tst.scale = scale;

        const indices = Int32Array.from({length: N}, (_,i) => i);
        _shuffle(indices);
        Object.freeze(indices.buffer);

        const n_fix = _rand_int(0,N);

        const {B} = ref,
               Z  = tabulate([N,N-n_fix], (i,j) => 1*(i === indices[n_fix+j]) ),
               y  = tabulate([N-n_fix,1], () => Math.random()*8 - 4);

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

        if( ++i >= 64 ) break;
      }
    })

    forEachItemIn(
      function*(){
        for( let run=0; run++ < 512; )
        for( let N=0; N++ < 24; )
        {
          const   G = Float64Array.from({length: N}, () => Math.random()*4.0 - 2.0),
                  d = Float64Array.from({length: N}, () => Math.random()*1.8 + 0.2),
                  D = diag_mat(d),
                  Q = rand_ortho(N),
                  B = matmul(Q.T, D, Q),
             bounds = Float64Array.from({length: 2*N}, () => Math.random()*4.0 - 2.0),
            updates = [];

          let i=0;
          for( const {data: row} of Q )
          {
            const dxScale = Math.random()*1.8 + 0.2,
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
            i *= 2;
            const [lo,hi] = [bounds[i],bounds[i+1]],
                        s = Math.random();
            return lo*(1-s) + s*hi;
          });
  
          Object.freeze(     G.buffer);
          Object.freeze(     X.buffer);
          Object.freeze(bounds.buffer);
          Object.freeze(B.data.buffer);
          Object.freeze(B);
  
          yield {G,B,X,bounds,updates};
        }
      }()
    ).it('compute_cauchyGeneralized() works given random examples' + suffix, ({G,B,X,bounds,updates}) => {
      const                 N = G.length;
      expect(X.length).toBe(N);
  
      const   G_ref = G.slice(),
              X_ref = X.slice(),
        indices_ref = Int32Array.from({length: N}, (_,i) => i); // TODO shuffle
  
      let   n_ref;
      for( [n_ref] of generalized_cauchy(
        G_ref,
        B.mapElems(x => x/2), // <- LBFGS B is actually Hessian times two
        X_ref,
        bounds,
        indices_ref
      ))
        /*pass*/;

      const lbfgs = new L_BFGS_B_Solver(N,N);

      for( const [dx,dg] of updates )
        lbfgs.update(dx,dg);
      lbfgs.scale = rand_scale();

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
              X_tst = X.slice(),
        indices_tst = Int32Array.from({length: N}, (_,i) => i); // TODO shuffle

      const n_tst = lbfgs.compute_cauchyGeneralized(G_tst, X_tst, bounds, indices_tst);

      expect(
        indices_tst.slice().sort((i,j) => i-j).every((i,I) => i===I)
      ).toBe(true);

      expect(n_tst).toBe(n_ref);
      expect(X_tst).toBeAllCloseTo(X_ref);
      expect(G_tst).toEqual(G_ref);
        expect( indices_tst.subarray(0,n_tst) )
      .toEqual( indices_ref.subarray(0,n_ref) )
    });


    forEachItemIn(
      function*(){
        for( let M=0; M++ <  8; )
        for( let N=0; N++ < 16; )
          yield [M,N];
      }()
    ).it('compute_bv and compute_ubbv work given random updates' + suffix, ([M,N]) => {
      const ref = new L_BFGS_SolverRef(M,N),
            tst = new L_BFGS_B_Solver (M,N);

      expect(M).toBeGreaterThan(0);
      expect(N).toBeGreaterThan(0);

      const shape = Int32Array.of(N,1);

      let i=0;
      for( let [dx,dg] of _rand_updates(N) )
      {
        ref.update(dx.data, dg.data);
        tst.update(dx.data, dg.data);
        for( let outmost=0; outmost++ < 2; )
        {
          const       scale = rand_scale();
          ref.scale = scale;
          tst.scale = scale;

          const {B} = ref;

          const u = Float64Array.from({length: N}, () => Math.random()*8-4),
                v = Float64Array.from({length: N}, () => Math.random()*8-4),
               bu = new Float64Array(2*M),
               bv = new Float64Array(2*M);
          Object.freeze(u.buffer);
          Object.freeze(v.buffer);

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

        if( ++i >= 64 ) break;
      }
    })


    forEachItemIn(
      function*(){
        for( let M=0; M++ < 12; )
        for( let N=0; N++ < 24; )
          yield [M,N];
      }()
    ).it('compute_be works given random updates' + suffix, ([M,N]) => {
      const tst = new L_BFGS_B_Solver(M,N);

      expect(M).toBeGreaterThan(0);
      expect(N).toBeGreaterThan(0);

      const shape = Int32Array.of(N,1);

      let i=0;
      for( let [dx,dg] of _rand_updates(N) )
      {
        tst.update(dx.data, dg.data);

        const m = Math.min(M,i+1);

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

        if( ++i >= 128 ) break;
      }
    })


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 512; )
        for( let N=0; N++ < 32; )
        {
          const   d = Float64Array.from({length: N}, () => Math.random()*1.8 + 0.2),
                  D = diag_mat(d),
                  Q = rand_ortho(N),
                  B = matmul(Q.T, D, Q),
            updates = [];

          let i=0;
          for( const {data: row} of Q )
          {
            const dxScale = Math.random()*1.8 + 0.2,
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

          yield [updates,B];
        }
      }()
    ).it('approximates Hessian given random multivariate polynomials' + suffix, ([updates,B]) => {
      const N = B.shape[0];

      const lbfgs = new L_BFGS_B_Solver(N,N);

      for( const [dx,dg] of updates )
        lbfgs.update(dx,dg);
      lbfgs.scale = rand_scale();

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
  }
})


describe('L_BFGS_SolverRef', () => {
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
        for( let M=0; M++ <  8; )
        for( let N=0; N++ < 16; )
          yield [M,N];
      }()
    ).it('computes B correctly given random updates' + suffix, ([M,N]) => {
      const lbfgs = new L_BFGS_SolverRef(M,N);
    
      expect(M).toBeGreaterThan(0);
      expect(N).toBeGreaterThan(0);
    
      let i=0;
      for( let [dx,dg] of _rand_updates(N) )
      {
        lbfgs.update(dx.data, dg.data);
        lbfgs.scale = rand_scale();
    
        const {B} = lbfgs;
        expect( matmul2(B,dx) ).toBeAllCloseTo(dg);
    
        if( ++i >= 128 ) break;
      }
    })
    
    
    forEachItemIn(
      function*(){
        for( let M=0; M++ <  8; )
        for( let N=0; N++ < 16; )
          yield [M,N];
      }()
    ).it('computes H correctly given random updates' + suffix, ([M,N]) => {
      const lbfgs = new L_BFGS_SolverRef(M,N);
    
      expect(M).toBeGreaterThan(0);
      expect(N).toBeGreaterThan(0);
    
      const shape = Int32Array.of(N,1);
      Object.freeze(shape.buffer);
    
      let i=0;
      for( let [dx,dg] of _rand_updates(N) )
      {
        lbfgs.update(dx.data, dg.data);
        lbfgs.scale = rand_scale();
    
        const {H} = lbfgs;
    
        const atol = Number.EPSILON * norm(H) * norm(dx);
        expect( matmul2(H,dg) ).toBeAllCloseTo(dx, {atol});
    
        if( ++i >= 128 ) break;
      }
    })


    forEachItemIn(
      function*(){
        for( let M=0; M++ <  8; )
        for( let N=0; N++ < 16; )
          yield [M,N];
      }()
    ).it('H=inv(B) given random updates' + suffix, ([M,N]) => {
      const lbfgs = new L_BFGS_SolverRef(M,N);
    
      expect(M).toBeGreaterThan(0);
      expect(N).toBeGreaterThan(0);
    
      const shape = Int32Array.of(N,1);
      Object.freeze(shape.buffer);
    
      const I = eye(N);
    
      let i=0;
      for( let [dx,dg] of _rand_updates(N) )
      {
        lbfgs.update(dx.data, dg.data);
        lbfgs.scale = rand_scale();
    
        const {H,B} = lbfgs;
    
        const tol = {
          rtol: 0,
          atol: Number.EPSILON * (1<<17) * norm(H) * norm(B)
        };
        expect( matmul2(B,H) ).toBeAllCloseTo(I, tol);
        expect( matmul2(H,B) ).toBeAllCloseTo(I, tol);
    
        if( ++i >= 64 ) break;
      }
    })


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 512; )
        for( let N=0; N++ < 16; )
        {
          const   d = Float64Array.from({length: N}, () => Math.random()*1.8 + 0.2),
                  D = diag_mat(d),
                  Q = rand_ortho(N),
                  B = matmul(Q.T, D, Q),
            updates = [];

          let i=0;
          for( const {data: row} of Q )
          {
            const dxScale = Math.random()*1.8 + 0.2,
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

          yield [updates,B];
        }
      }()
    ).it('approximates Hessian given random multivariate polynomials' + suffix, ([updates,B]) => {
      const N = B.shape[0];

      const lbfgs = new L_BFGS_SolverRef(N,N);

      for( const [dx,dg] of updates )
        lbfgs.update(dx,dg);
      lbfgs.scale = rand_scale();

      const  b = lbfgs.B;
      expect(b).toBeAllCloseTo(B);
    });
  }
})
