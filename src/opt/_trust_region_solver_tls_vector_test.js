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

import {concat} from "../concat";
import {forEachItemIn, CUSTOM_MATCHERS} from '../jasmine_utils'
import {asarray, NDArray} from "../nd_array";
import {stack} from "../stack";
import {tabulate} from "../tabulate";
import {zip_elems} from '../zip_elems';

import {checked_array} from "../arrays/_checked_array";
import {     is_array} from "../arrays/is_array";

import {AleaRNG} from "../rand/alea_rng";

import {eye} from "../la/eye";
import {lstsq} from "../la/lstsq";
import {matmul2} from "../la/matmul";
import {qr_decomp_full} from "../la/qr";
import {svd_lstsq} from "../la/svd";
import {svd_jac_2sided} from "../la/svd_jac_2sided";

import {num_grad} from "./num_grad";
import {TrustRegionSolverLSQ} from "./_trust_region_solver_lsq";
import {TrustRegionSolverTLS} from "./_trust_region_solver_tls";


//
// CREATE VECTOR TEST FUNCTIONS
// ----------------------------
function rand_func( rng, NX, NY, NP )
{
  const atoms = [
    function(){
      const f = x => 1;
      f.grad  = x => 0;
      f.str   = sym => `1`;
      return f;
    }(),
    function(){
      const f = x => x;
      f.grad  = x => 1;
      f.str   = sym => `${sym}`;
      return f;
    }(),
    function(){
      const f = x => x*x;
      f.grad  = x => 2*x;
      f.str   = sym => `${sym}²`;
      return f;
    }(),
    function(){
      const f = x => { const exp = Math.exp(x); return -0.5 + +1 / (1+exp); };
      f.grad  = x => { const exp = Math.exp(x); return        -1 / (1/exp + 2 + exp); };
      f.str   = sym => `( 1 / (1 + exp(${sym}) )`;
      return f;
    }()
  ];

  const NC = 4,
    coeffs = Float64Array.from({length: NY*NC}, () => rng.uniform(-1.5,+1.5) );

  const composition = [];

  for( let i=0; i < NY*NC; i++ )
  {
    let ip = rng.int(0,NP),
        ix = rng.int(0,NX),
        fp = atoms[rng.int(0,atoms.length)],
        fx = atoms[rng.int(0,atoms.length)];

    composition.push([ip,ix,fp,fx]);
  }

  const f = (p,x) => 
  {
    p = asarray('float64', p);
    if( p.ndim     !== 1  ) throw new Error('Assertion failed.');
    if( p.shape[0] !== NP ) throw new Error('Assertion failed.');
    p = p.data;

    x = asarray('float64', x);
    if(    2 !== x.ndim     ) throw new Error('Assertion failed.');
    if(   NX !== x.shape[1] ) throw new Error('Assertion failed.');
    const MX  =  x.shape[0];
    x = x.data;

    const f = new Float64Array(MX*NY);

    for( let i=0; i < NY; i++ )
    for( let j=0; j < NC; j++ )
    {
      const [ip,ix,fp,fx] = composition[NC*i+j];
      for( let k=0; k < MX; k++ )
        f[NY*k+i] += coeffs[NC*i+j] * fp( p[ip] ) * fx( x[NX*k+ix] );
    }

    return new NDArray( Int32Array.of(MX,NY), f );
  }

  const fgg = (p,x) =>
  {
    p = asarray('float64', p);
    if( p.ndim     !== 1  ) throw new Error('Assertion failed.');
    if( p.shape[0] !== NP ) throw new Error('Assertion failed.');
    p = p.data;

    x = asarray('float64', x);
    if(    2 !== x.ndim     ) throw new Error('Assertion failed.');
    if(   NX !== x.shape[1] ) throw new Error('Assertion failed.');
    const MX  =  x.shape[0];
    x = x.data;

    const f = new Float64Array(MX*NY),
         gp = new Float64Array(MX*NY*NP),
         gx = new Float64Array(MX*NY*NX);

    for( let i=0; i < NY; i++ )
    for( let j=0; j < NC; j++ )
    {
      const [ip,ix,fp,fx] = composition[NC*i+j];

      for( let k=0; k < MX; k++ )
      {
             f[NY*k+i]     += coeffs[NC*i+j] * fp     (p[ip]) * fx     (x[NX*k+ix]);
        gp[NP*(NY*k+i)+ip] += coeffs[NC*i+j] * fp.grad(p[ip]) * fx     (x[NX*k+ix]);
        gx[NX*(NY*k+i)+ix] += coeffs[NC*i+j] * fp     (p[ip]) * fx.grad(x[NX*k+ix]);
      }
    }

    return [
      new NDArray( Int32Array.of(MX,NY),    f  ),
      new NDArray( Int32Array.of(MX,NY,NP), gp ),
      new NDArray( Int32Array.of(MX,NY,NX), gx )
    ];
  }

  return {
    NX, NY, NP,
    f, fgg
  };
}


const vector_funcs = {};

;{
  const rng = new AleaRNG('TrustRegionSolverTLS vector functions');

  for( let NX=0; NX++ < 4; )
  for( let NY=0; NY++ < 4; )
    vector_funcs[`rand_func${NX}${NY}`] = rand_func(rng, NX,NY, (NX+NY)*2);
}


vector_funcs.poly22 = {
  NP: 4,
  NX: 2,
  NY: 2,
  f: ([a,b,c,d], xy) => {
    xy = asarray('float64', xy);
    const MX = xy.shape[0];
    if(  2 !== xy.shape[1] ) throw new Error('Assertion failed.');
    xy = xy.data;

    const f = new Float64Array(MX*2);

    for( let i=0; i < MX; i++ )
    { const x = xy[2*i+0],
            y = xy[2*i+1];
      f[2*i+0] = a     + b*x + c*x*y + d*y;
      f[2*i+1] = a*x*y + b   + c*y   + d*x;
    }

    return new NDArray(Int32Array.of(MX,2), f);
  },
  fgg: ([a,b,c,d], xy) => {
    xy = asarray('float64', xy);
    const MX = xy.shape[0];
    if(  2 !== xy.shape[1] ) throw new Error('Assertion failed.');
    xy = xy.data;

    const f = new Float64Array(MX*2),
         gp = new Float64Array(MX*2*4),
         gx = new Float64Array(MX*2*2);

    for( let i=0; i < MX; i++ )
    { const x = xy[2*i+0],
            y = xy[2*i+1];
      f[2*i+0] = a     + b*x + c*x*y + d*y;
      f[2*i+1] = a*x*y + b   + c*y   + d*x;
    }

    for( let i=0; i < MX; i++ )
    { const x = xy[2*i+0],
            y = xy[2*i+1];
      gp[8*i+0] = 1;
      gp[8*i+1] = x;
      gp[8*i+2] = x*y;
      gp[8*i+3] =   y;
      gp[8*i+4] = x*y;
      gp[8*i+5] = 1;
      gp[8*i+6] =   y;
      gp[8*i+7] = x;
    }

    for( let i=0; i < MX; i++ )
    { const x = xy[2*i+0],
            y = xy[2*i+1];
      gx[4*i+0] = b + c*y;
      gx[4*i+1] = d + c*x;
      gx[4*i+2] = d + a*y;
      gx[4*i+3] = c + a*x;
    }

    return [
      new NDArray(Int32Array.of(MX,2  ),  f),
      new NDArray(Int32Array.of(MX,2,4), gp),
      new NDArray(Int32Array.of(MX,2,2), gx)
    ];
  }
};


//
// TEST VECTOR TEST FUNCTIONS
// --------------------------
describe('TrustRegionSolverTLS test vector functions', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  });


  for( const [name,{NP,NX,NY,f,fgg}] of Object.entries(vector_funcs) )
    forEachItemIn(
      function*(rng){
        for( let run=0; run++ < 256; )
        {
          const MX = rng.int(1,64);
          yield [
            tabulate([NP],    'float64', () => rng.uniform(-2,+2) ),
            tabulate([MX,NX], 'float64', () => rng.uniform(-2,+2) )
          ];
        }
      }
    ).it(`${name} derivatives are correct`, ([p, x]) => {
      const GP = num_grad( p => f(p,x) )(p),
            GX = stack( [...x].map( num_grad( x => f(p,x.reshape(1,NX)).reshape(NY) ) ) );

      const [F,gp,gx] = fgg(p,x);

      expect(F).toBeAllCloseTo( f(p,x) );
      expect(gp).toBeAllCloseTo(GP);
      expect(gx).toBeAllCloseTo(GX);
    });
})


//
// TEST TrustRegionSolverTLS USING VECTOR TEST FUNCTIONS
// -----------------------------------------------------
for( const [name,{NP,NX,NY,f,fgg}] of Object.entries(vector_funcs) )
  describe(`TrustRegionSolverTLS [vec. func.: ${name}]`, () => {
    beforeEach( () => {
      jasmine.addMatchers(CUSTOM_MATCHERS)
    });


    forEachItemIn(
      function*(rng){
        for( let run=0; run++ < 337; )
        {
          const             MX = rng.int( 1, Math.ceil(42/Math.max(NX,NY)) ),
             p0 = tabulate(   [NP], 'float64', () => rng.uniform(-2,+2) ),
            dx0 = tabulate([MX,NX], 'float64', () => rng.uniform(-2,+2) );
          yield [p0,dx0];
        }
      }
    ).it(`initial X0,F0,G0,D,J is correct given random examples`, ([p0, dx0]) => {
      const     solver = new TrustRegionSolverTLS(fgg, p0,dx0),
        {M,N} = solver;

      const MX = dx0.shape[0];

      const                            X0 = concat([dx0.reshape(-1), p0]);
      expect(solver.X0).toBeAllCloseTo(X0);

      const F = X => {
        const p = new NDArray( p0.shape, X.data.slice(  -NP) ),
             dx = new NDArray(dx0.shape, X.data.slice(0,-NP) ),
             dy = f(p,dx);
        return concat([dx.reshape(-1), dy.reshape(-1)]);
      };

      const                            F0 = F(X0);
      expect(solver.F0).toBeAllCloseTo(F0);

      const         loss = x => F(x).data.reduce((loss,r) => loss + r*r/2, 0.0);
      expect(solver.loss).toBeAllCloseTo( loss(X0) / F0.shape[0] * 2 );

      expect(solver.G0).toBeAllCloseTo(num_grad(loss)(X0), {atol:1e-5});

      const  J = tabulate([M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j));
      expect(J).toBeAllCloseTo( num_grad(F)(X0) );

      const                            G0 = matmul2(F0.reshape(1,-1), J);
      expect(solver.G0).toBeAllCloseTo(G0);

      const D = J.reduceElems(0, 'float64', (x,y) => Math.hypot(x,y));
      expect(solver.D).toBeAllCloseTo(D);
    });


    forEachItemIn(
      function*(rng){
        for( let run=0; run++ < 137; )
        {
          const             MX = rng.int( 1, Math.ceil(42/Math.max(NX,NY)) ),
             p0 = tabulate(   [NP], 'float64', () => rng.uniform(-2,+2) ),
            dx0 = tabulate([MX,NX], 'float64', () => rng.uniform(-2,+2) );

          function* seq() {
            for( let step=0; step++ < 6; )
              yield tabulate([MX*NX+NP,1], 'float64', () => rng.uniform(-2,+2) );
          }

          yield [p0,dx0, seq()]
        }
      }
    ).it(`considerMove(dX) returns correct loss and loss prediction`, ([p0, dx0, seq]) => {
      const     solver = new TrustRegionSolverTLS(fgg, p0,dx0),
        {M,N} = solver;

      const                            X0 = concat([dx0.reshape(-1), p0]);
      expect(solver.X0).toBeAllCloseTo(X0);

      const F = X => {
        const p = new NDArray( p0.shape, X.data.slice(  -NP) ),
             dx = new NDArray(dx0.shape, X.data.slice(0,-NP) ),
             dy = f(p,dx);
        return concat([dx.reshape(-1), dy.reshape(-1)]);
      };

      const F0 = F(X0),
            G0 = num_grad( dX => solver.considerMove(dX.data)[0] * M / 2 )( new Float64Array(X0.shape[0]) ),
          loss = F => F.data.reduce((loss,r) => loss + r*r/M, 0.0);

      expect(solver.F0).toBeAllCloseTo(F0);
      expect(solver.loss).toBeAllCloseTo( loss(F0) );
      expect(solver.G0).toBeAllCloseTo(G0, {atol: 1e-5});

      const J = tabulate([M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) );

      for( const dX of seq )
      {
        const X1 = X0.mapElems((x0,i) => x0+dX.data[i]);

        const [loss_predict,
               loss_new] = solver.considerMove(dX.data);
        expect(loss_new).toBeAllCloseTo( loss( F(X1) ) );

        const                                     F_predict = zip_elems( [F0.reshape(-1,1), matmul2(J,dX)], (f,Jdx) => f + Jdx );
        expect(loss_predict).toBeAllCloseTo( loss(F_predict) );
      }

      expect(solver.F0).toBeAllCloseTo(F0);
      expect(solver.loss).toBeAllCloseTo( loss(F0) );
      expect(solver.G0).toBeAllCloseTo(G0, {atol: 1e-5});
    });


    forEachItemIn(
      function*(rng){
        for( let run=0; run++ < 1733; )
        {
          const             MX = rng.int(1,96),
             p0 = tabulate(   [NP], 'float64', () => rng.uniform(-2,+2) ),
            dx0 = tabulate([MX,NX], 'float64', () => rng.uniform(-2,+2) );
          yield [p0,dx0]
        }
      }
    ).it(`cauchyTravel() finds min. along gradient direction`, ([p0, dx0]) => {
      const solver = new TrustRegionSolverTLS(fgg, p0,dx0);

      const g = num_grad( cp => {
        const                                      dx = solver.G0.map( g => g*cp ),
              [loss_predict] = solver.considerMove(dx);
        return loss_predict;
      });

      const  cp = solver.cauchyTravel();
      expect(cp).toBeLessThanOrEqual(0);

      expect( g(cp) ).toBeAllCloseTo(0, {atol: 1e-5});
    });


    forEachItemIn(
      function*(rng){
        for( let run=0; run++ < 96; )
        {
          const             MX = rng.int( 1, Math.ceil(42/Math.max(NX,NY)) ),
             p0 = tabulate(   [NP], 'float64', () => rng.uniform(-2,+2) ),
            dx0 = tabulate([MX,NX], 'float64', () => rng.uniform(-2,+2) );
          yield [rng, p0,dx0]
        }
      }
    ).it(`initial computeNewton() produces no out-of-bounds errors given random examples`, ([rng, p0,dx0]) => {
      const        solver = new TrustRegionSolverTLS(fgg, p0,dx0),
        {M,N,MX} = solver;

      for( const [k,v] of Object.entries(solver) )
        if( is_array(v) )
          solver[k] = checked_array(v);

      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = rng.uniform(0.1,1.9);
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = rng.uniform(-2,+2);
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = rng.uniform(-2,+2);
      for( let i=M       ; i-- > 0; ) solver. F0[i] = rng.uniform(-2,+2);
      for( let i=N       ; i-- > 0; ) solver.  D[i] = rng.uniform(0.01, 2 ) * (rng.uniform(0,1) < 0.99 || i < MX*NX);

      const J = tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) );

      const              F_shape = Int32Array.of(M,1),
        F = new NDArray( F_shape, solver.F0.slice() );

      solver.computeNewton();

      // check that J,F is unmodified by computeNewton
      expect(J).toBeAllCloseTo( tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ), {rtol:0, atol:0} );
      expect(F).toBeAllCloseTo(  new NDArray( F_shape, solver.F0.slice() ) );

      expect(solver.rank).toBe( Math.min(M,N) );

      if( J.data.some(isNaN) ) throw new Error('Assertion failed.');

      const [r,V,D] = solver.__DEBUG_RVD,
              [Q,R] = qr_decomp_full( matmul2( zip_elems([J,D], (j,d) => j/d), V.T ) );
      // QR decomp. may differ by sign
      ;{
        const L = Math.min(M,N),
              RR = R.data,
              rr = r.data,
              QQ = Q.data;
        for( let i=L; i-- > 0; )
        {
          let dot = 0.0;
          for( let j=N; j-- > 0; )
            dot += RR[N*i+j] * rr[N*i+j];

          if( dot < 0 )
          {
            for( let j=N; j-- > 0; ) RR[N*i+j] *= -1;
            for( let j=M; j-- > 0; ) QQ[M*j+i] *= -1;
          }
        }
      }

      const                                   I = eye(N);
      expect( matmul2(V,V.T) ).toBeAllCloseTo(I);
      expect( matmul2(V.T,V) ).toBeAllCloseTo(I);

      if( r.data.some(isNaN) ) throw new Error('Assertion failed.');
      if( V.data.some(isNaN) ) throw new Error('Assertion failed.');

      expect(r).toBeAllCloseTo(R);

      const                    {rank} = solver;
      expect( solver.QF.slice(0,rank) ).toBeAllCloseTo( matmul2(Q.T,F).data.slice(0,rank) );
    });


    forEachItemIn(
      function*(rng){
        for( let run=0; run++ < 173; )
        {
          const             MX = rng.int( Math.ceil(NP/NY), Math.ceil(42/Math.max(NX,NY)) ),
             p0 = tabulate(   [NP], 'float64', () => rng.uniform(-2,+2) ),
            dx0 = tabulate([MX,NX], 'float64', () => rng.uniform(-2,+2) );
          yield [rng, p0,dx0]
        }
      }
    ).it(`initial computeNewton() QR decomposes J correctly given random over-determined examples`, ([rng, p0,dx0]) => {
      const        solver = new TrustRegionSolverTLS(fgg, p0,dx0),
        {M,N,MX} = solver;

      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = rng.uniform(0.1,1.9);
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = rng.uniform(-2,+2);
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = rng.uniform(-2,+2);
      for( let i=M       ; i-- > 0; ) solver. F0[i] = rng.uniform(-2,+2);

      const J = tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) );

      const              F_shape = Int32Array.of(M,1),
        F = new NDArray( F_shape, solver.F0.slice() );

      solver.computeNewton();

      // check that J,F is unmodified by computeNewton
      expect(J).toBeAllCloseTo( tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ), {rtol:0, atol:0} );
      expect(F).toBeAllCloseTo(  new NDArray( F_shape, solver.F0.slice() ) );

      expect(solver.rank).toBe(N);

      if( J.data.some(isNaN) ) throw new Error('Assertion failed.');

      const [r,V,D] = solver.__DEBUG_RVD,
              [Q,R] = qr_decomp_full( matmul2( zip_elems([J,D], (j,d) => j/d), V.T ) );
      // QR decomp. may differ by sign
      ;{
        const L = Math.min(M,N),
              RR = R.data,
              rr = r.data,
              QQ = Q.data;
        for( let i=L; i-- > 0; )
        {
          let dot = 0.0;
          for( let j=N; j-- > 0; )
            dot += RR[N*i+j] * rr[N*i+j];

          if( dot < 0 )
          {
            for( let j=N; j-- > 0; ) RR[N*i+j] *= -1;
            for( let j=M; j-- > 0; ) QQ[M*j+i] *= -1;
          }
        }
      }

      const                                   I = eye(N);
      expect( matmul2(V,V.T) ).toBeAllCloseTo(I);
      expect( matmul2(V.T,V) ).toBeAllCloseTo(I);

      if( r.data.some(isNaN) ) throw new Error('Assertion failed.');
      if( V.data.some(isNaN) ) throw new Error('Assertion failed.');

      expect(r).toBeAllCloseTo(R);

      const                    {rank} = solver;
      expect( solver.QF.slice(0,rank) ).toBeAllCloseTo( matmul2(Q.T,F).data.slice(0,rank) );
    });


    forEachItemIn(
      function*(rng){
        for( let run=0; run++ < 173; )
        {
          const             MX = rng.int( 1, Math.ceil(42/Math.max(NX,NY)) ),
             p0 = tabulate(   [NP], 'float64', () => rng.uniform(-2,+2) ),
            dx0 = tabulate([MX,NX], 'float64', () => rng.uniform(-2,+2) );
          yield [rng, p0,dx0]
        }
      }
    ).it(`initial computeNewton() QR decomposes J correctly given random examples`, ([rng, p0,dx0]) => {
      const        solver = new TrustRegionSolverTLS(fgg, p0,dx0),
        {M,N,MX} = solver;

      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = rng.uniform(0.1,1.9);
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = rng.uniform(-2,+2);
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = rng.uniform(-2,+2);
      for( let i=M       ; i-- > 0; ) solver. F0[i] = rng.uniform(-2,+2);
      for( let i=N       ; i-- > 0; ) solver.  D[i] = rng.uniform(0.01, 2 ) * (rng.uniform(0,1) < 0.99 || i < MX*NX);

      const J = tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) );

      const              F_shape = Int32Array.of(M,1),
        F = new NDArray( F_shape, solver.F0.slice() );

      solver.computeNewton();

      // check that J,F is unmodified by computeNewton
      expect(J).toBeAllCloseTo( tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ), {rtol:0, atol:0} );
      expect(F).toBeAllCloseTo(  new NDArray( F_shape, solver.F0.slice() ) );

      expect(solver.rank).toBe( Math.min(M,N) );

      if( J.data.some(isNaN) ) throw new Error('Assertion failed.');

      const [r,V,D] = solver.__DEBUG_RVD,
              [Q,R] = qr_decomp_full( matmul2( zip_elems([J,D], (j,d) => j/d), V.T ) );
      // QR decomp. may differ by sign
      ;{
        const L = Math.min(M,N),
              RR = R.data,
              rr = r.data,
              QQ = Q.data;
        for( let i=L; i-- > 0; )
        {
          let dot = 0.0;
          for( let j=N; j-- > 0; )
            dot += RR[N*i+j] * rr[N*i+j];

          if( dot < 0 )
          {
            for( let j=N; j-- > 0; ) RR[N*i+j] *= -1;
            for( let j=M; j-- > 0; ) QQ[M*j+i] *= -1;
          }
        }
      }

      const                                   I = eye(N);
      expect( matmul2(V,V.T) ).toBeAllCloseTo(I);
      expect( matmul2(V.T,V) ).toBeAllCloseTo(I);

      if( r.data.some(isNaN) ) throw new Error('Assertion failed.');
      if( V.data.some(isNaN) ) throw new Error('Assertion failed.');

      expect(r).toBeAllCloseTo(R);

      const                    {rank} = solver;
      expect( solver.QF.slice(0,rank) ).toBeAllCloseTo( matmul2(Q.T,F).data.slice(0,rank) );
    });


    forEachItemIn(
      function*(rng){
        for( let run=0; run++ < 173; )
        {
          const             MX = rng.int( 1, Math.ceil(42/Math.max(NX,NY)) ),
             p0 = tabulate(   [NP], 'float64', () => rng.uniform(-2,+2) ),
            dx0 = tabulate([MX,NX], 'float64', () => rng.uniform(-2,+2) );
          yield [rng, p0,dx0]
        }
      }
    ).it(`initial computeNewton() QR decomposes J correctly given random rank-deficient examples`, ([rng, p0,dx0]) => {
      const        solver = new TrustRegionSolverTLS(fgg, p0,dx0),
        {M,N,MX} = solver;

      const [{data: J22},rnk] = rng.rankDef(MX*NY,NP);
      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = rng.uniform(0.1, 1.9);
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = rng.uniform( -2, +2 );
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = J22[i];
      for( let i=M       ; i-- > 0; ) solver. F0[i] = rng.uniform( -2, +2 );
      for( let i=N       ; i-- > 0; ) solver.  D[i] = rng.uniform(0.01, 2 ) * (rng.uniform(0,1) < 0.99 || i < MX*NX);

      const J = tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) );

      const              F_shape = Int32Array.of(M,1),
        F = new NDArray( F_shape, solver.F0.slice() );

      solver.computeNewton();

      expect(solver.rank).toBe( MX*NX + rnk );

      // check that J,F is unmodified by computeNewton
      expect(J).toBeAllCloseTo( tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ), {rtol:0, atol:0} );
      expect(F).toBeAllCloseTo(  new NDArray( F_shape, solver.F0.slice() ) );

      expect(solver.rank).toBe( MX*NX + rnk );

      if( J.data.some(isNaN) ) throw new Error('Assertion failed.');

      const [r,V,D] = solver.__DEBUG_RVD,
              [Q,R] = qr_decomp_full( matmul2( zip_elems([J,D], (j,d) => j/d), V.T ) );
      // QR decomp. may differ by sign
      ;{
        const L = Math.min(M,N),
              RR = R.data,
              rr = r.data,
              QQ = Q.data;
        for( let i=L; i-- > 0; )
        {
          let dot = 0.0;
          for( let j=N; j-- > 0; )
            dot += RR[N*i+j] * rr[N*i+j];

          if( dot < 0 )
          {
            for( let j=N; j-- > 0; ) RR[N*i+j] *= -1;
            for( let j=M; j-- > 0; ) QQ[M*j+i] *= -1;
          }
        }
      }

      const                                   I = eye(N);
      expect( matmul2(V,V.T) ).toBeAllCloseTo(I);
      expect( matmul2(V.T,V) ).toBeAllCloseTo(I);

      if( r.data.some(isNaN) ) throw new Error('Assertion failed.');
      if( V.data.some(isNaN) ) throw new Error('Assertion failed.');

      expect(r).toBeAllCloseTo(R);

      const                    {rank} = solver;
      expect( solver.QF.slice(0,rank) ).toBeAllCloseTo( matmul2(Q.T,F).data.slice(0,rank) );
    });


    forEachItemIn(
      function*(rng){
        for( let run=0; run++ < 173; )
        {
          const             MX = rng.int( 1, Math.ceil(42/Math.max(NX,NY)) ),
             p0 = tabulate(   [NP], 'float64', () => rng.uniform(-2,+2) ),
            dx0 = tabulate([MX,NX], 'float64', () => rng.uniform(-2,+2) );
          yield [rng, p0,dx0]
        }
      }
    ).it(`initial computeNewton() QR decomposes J correctly given random sparse examples`, ([rng, p0,dx0]) => {
      const        solver = new TrustRegionSolverTLS(fgg, p0,dx0),
        {M,N,MX} = solver;

      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = rng.uniform(0.1, 1.9);
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = rng.uniform( -2, +2 ) * (rng.uniform(0,1) < 0.95);
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = rng.uniform( -2, +2 ) * (rng.uniform(0,1) < 0.95);
      for( let i=M       ; i-- > 0; ) solver. F0[i] = rng.uniform( -2, +2 );
      for( let i=N       ; i-- > 0; ) solver.  D[i] = rng.uniform(0.01, 2 ) * (rng.uniform(0,1) < 0.99 || i < MX*NX);

      const J = tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) );

      const              F_shape = Int32Array.of(M,1),
        F = new NDArray( F_shape, solver.F0.slice() );

      solver.computeNewton();

      // check that J,F is unmodified by computeNewton
      expect(J).toBeAllCloseTo( tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ), {rtol:0, atol:0} );
      expect(F).toBeAllCloseTo(  new NDArray( F_shape, solver.F0.slice() ) );

      if( J.data.some(isNaN) ) throw new Error('Assertion failed.');

      const [r,V,D] = solver.__DEBUG_RVD,
              [Q,R] = qr_decomp_full( matmul2( zip_elems([J,D], (j,d) => j/d), V.T ) );
      // QR decomp. may differ by sign
      ;{
        const L = Math.min(M,N),
              RR = R.data,
              rr = r.data,
              QQ = Q.data;
        for( let i=L; i-- > 0; )
        {
          let dot = 0.0;
          for( let j=N; j-- > 0; )
            dot += RR[N*i+j] * rr[N*i+j];

          if( dot < 0 )
          {
            for( let j=N; j-- > 0; ) RR[N*i+j] *= -1;
            for( let j=M; j-- > 0; ) QQ[M*j+i] *= -1;
          }
        }
      }

      const                                   I = eye(N);
      expect( matmul2(V,V.T) ).toBeAllCloseTo(I);
      expect( matmul2(V.T,V) ).toBeAllCloseTo(I);

      if( r.data.some(isNaN) ) throw new Error('Assertion failed.');
      if( V.data.some(isNaN) ) throw new Error('Assertion failed.');

      expect(r).toBeAllCloseTo(R);

      const                    {rank} = solver;
      expect( solver.QF.slice(0,rank) ).toBeAllCloseTo( matmul2(Q.T,F).data.slice(0,rank) );
    });


    forEachItemIn(
      function*(rng){
        for( let run=0; run++ < 256; )
        {
          const             MX = rng.int( Math.ceil(NP/NY), Math.ceil(42/Math.max(NX,NY)) ),
             p0 = tabulate(   [NP], 'float64', () => rng.uniform(-2,+2) ),
            dx0 = tabulate([MX,NX], 'float64', () => rng.uniform(-2,+2) );
          yield [rng, p0,dx0]
        }
      }
    ).it(`initial computeNewton() solves random over-determined examples`, ([rng, p0,dx0]) => {
      const        solver = new TrustRegionSolverTLS(fgg, p0,dx0),
        {M,N,MX} = solver;

      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = rng.uniform(0.1, 1.9);
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = rng.uniform(-2,+2);
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = rng.uniform(-2,+2);
      for( let i=M       ; i-- > 0; ) solver. F0[i] = rng.uniform(-2,+2);

      const J = tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) );

      const              F_shape = Int32Array.of(M,1),
        F = new NDArray( F_shape, solver.F0.slice() );

      solver.computeNewton();

      expect(solver.rank).toBe(N);

      // check that J,F is unmodified by computeNewton
      expect(J).toBeAllCloseTo( tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ), {rtol:0, atol:0} );
      expect(F).toBeAllCloseTo(  new NDArray( F_shape, solver.F0.slice() ) );

      const  X = new NDArray( Int32Array.of(N,1), solver.newton_dX.slice() );
      expect(X).toBeAllCloseTo( lstsq(J,F).mapElems('float64', x => -x) );
    });


    forEachItemIn(
      function*(rng){
        for( let run=0; run++ < 512; )
        {
          const             MX = rng.int( 1, Math.ceil(42/Math.max(NX,NY)) ),
             p0 = tabulate(   [NP], 'float64', () => rng.uniform(-2,+2) ),
            dx0 = tabulate([MX,NX], 'float64', () => rng.uniform(-2,+2) );
          yield [rng, p0,dx0]
        }
      }
    ).it(`initial computeNewton() solves random examples with J21=0`, ([rng, p0,dx0]) => {
      const        solver = new TrustRegionSolverTLS(fgg, p0,dx0),
        {M,N,MX} = solver;

      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = rng.uniform(0.1, 1.9);
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = 0;
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = rng.uniform( -2, +2 );
      for( let i=M       ; i-- > 0; ) solver. F0[i] = rng.uniform( -2, +2 );
      for( let i=N       ; i-- > 0; ) solver.  D[i] = rng.uniform(0.01, 2 ) * (rng.uniform(0,1) < 0.99 || i < MX*NX);

      const J  = tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ),
             D = new NDArray( Int32Array.of(N,1), solver.D.slice() ),
            JD = zip_elems([J,D.T], (j,d) => d===0 ? j : j/d);

      const              F_shape = Int32Array.of(M,1),
        F = new NDArray( F_shape, solver.F0.slice() );

      solver.computeNewton();

      expect(solver.rank).toBe( Math.min(M,N) );

      // check that J,F is unmodified by computeNewton
      expect(J).toBeAllCloseTo( tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ), {rtol:0, atol:0} );
      expect(F).toBeAllCloseTo(  new NDArray( F_shape, solver.F0.slice() ) );

      const   X = new NDArray( Int32Array.of(N,1), solver.newton_dX.slice() ),
             DX = zip_elems([D,X], (d,x) => d===0 ? x : d*x );
      expect(DX).toBeAllCloseTo( lstsq(JD,F).mapElems('float64', x => -x) );
    });


    forEachItemIn(
      function*(rng){
        for( let run=0; run++ < 373; )
        {
          const             MX = rng.int( 1, Math.ceil(42/Math.max(NX,NY)) ),
             p0 = tabulate(   [NP], 'float64', () => rng.uniform(-2,+2) ),
            dx0 = tabulate([MX,NX], 'float64', () => rng.uniform(-2,+2) );
          yield [rng, p0,dx0];
        }
      }
    ).it(`initial computeNewton() solves random rank-deficient examples with J21=0`, ([rng, p0,dx0]) => {
      const        solver = new TrustRegionSolverTLS(fgg, p0,dx0),
        {M,N,MX} = solver;

      const [{data: J22},rnk] = rng.rankDef(MX*NY,NP);

      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = rng.uniform(0.1, 1.9);
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = 0;
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = J22[i];
      for( let i=M       ; i-- > 0; ) solver. F0[i] = rng.uniform( -2, +2 );
      for( let i=N       ; i-- > 0; ) solver.  D[i] = rng.uniform(0.01, 2 ) * (rng.uniform(0,1) < 0.99 || i < MX*NX);

      const J  = tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ),
             D = new NDArray( Int32Array.of(N,1), solver.D.slice() ),
            JD = zip_elems([J,D.T], (j,d) => d===0 ? j : j/d);

      const              F_shape = Int32Array.of(M,1),
        F = new NDArray( F_shape, solver.F0.slice() );

      solver.computeNewton();

      expect(solver.rank).toBe(MX*NX + rnk);

      // check that J,F is unmodified by computeNewton
      expect(J).toBeAllCloseTo( tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ), {rtol:0, atol:0} );
      expect(F).toBeAllCloseTo(  new NDArray( F_shape, solver.F0.slice() ) );

      const   X = new NDArray( Int32Array.of(N,1), solver.newton_dX.slice() ),
             DX = zip_elems([D,X], (d,x) => d===0 ? x : d*x );
      expect(DX).toBeAllCloseTo( svd_lstsq(...svd_jac_2sided(JD),F).mapElems('float64', x => -x) );
    });


    forEachItemIn(
      function*(rng){
        for( let run=0; run++ < 337; )
        {
          const             MX = rng.int( 1, Math.ceil(42/Math.max(NX,NY)) ),
             p0 = tabulate(   [NP], 'float64', () => rng.uniform(-2,+2) ),
            dx0 = tabulate([MX,NX], 'float64', () => rng.uniform(-2,+2) );
          yield [rng, p0,dx0]
        }
      }
    ).it(`initial computeNewton() solves random examples`, ([rng, p0,dx0]) => {
      const        solver = new TrustRegionSolverTLS(fgg, p0,dx0),
        {M,N,MX} = solver;

      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = rng.uniform(0.1, 1.9);
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = rng.uniform( -2, +2 );
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = rng.uniform( -2, +2 );
      for( let i=M       ; i-- > 0; ) solver. F0[i] = rng.uniform( -2, +2 );
      for( let i=N       ; i-- > 0; ) solver.  D[i] = rng.uniform(0.01, 2 ) * (rng.uniform(0,1) < 0.99 || i < MX*NX);

      const J  = tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ),
             D = new NDArray( Int32Array.of(N,1), solver.D.slice() ),
            JD = zip_elems([J,D.T], (j,d) => d===0 ? j : j/d);

      const              F_shape = Int32Array.of(M,1),
        F = new NDArray( F_shape, solver.F0.slice() );

      solver.computeNewton();

      expect(solver.rank).toBe( Math.min(M,N) );

      // check that J,F is unmodified by computeNewton
      expect(J).toBeAllCloseTo( tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ), {rtol:0, atol:0} );
      expect(F).toBeAllCloseTo(  new NDArray( F_shape, solver.F0.slice() ) );

      const   X = new NDArray( Int32Array.of(N,1), solver.newton_dX.slice() ),
             DX = zip_elems([D,X], (d,x) => d===0 ? x : d*x );
      expect(DX).toBeAllCloseTo( lstsq(JD,F).mapElems('float64', x => -x) );
    });


    forEachItemIn(
      function*(rng){
        for( let run=0; run++ < 73; )
        {
          const             MX = rng.int( 1, Math.ceil(42/Math.max(NX,NY)) ),
             p0 = tabulate(   [NP], 'float64', () => rng.uniform(-2,+2) ),
            dx0 = tabulate([MX,NX], 'float64', () => rng.uniform(-2,+2) );
          yield [rng, p0,dx0]
        }
      }
    ).it(`initial computeNewton() solves random rank-deficient examples`, ([rng, p0,dx0]) => {
      const        solver = new TrustRegionSolverTLS(fgg, p0,dx0),
        {M,N,MX} = solver;

      const [{data: J22},rnk] = rng.rankDef(MX*NY,NP);

      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = rng.uniform(0.1, 1.9);
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = rng.uniform( -2, +2 );
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = J22[i];
      for( let i=M       ; i-- > 0; ) solver. F0[i] = rng.uniform( -2, +2 );
      for( let i=N       ; i-- > 0; ) solver.  D[i] = rng.uniform(0.01, 2 ) * (rng.uniform(0,1) < 0.99 || i < MX*NX);

      const J  = tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ),
             D = new NDArray( Int32Array.of(N,1), solver.D.slice() ),
            JD = zip_elems([J,D.T], (j,d) => d===0 ? j : j/d);

      const              F_shape = Int32Array.of(M,1),
        F = new NDArray( F_shape, solver.F0.slice() );

      solver.computeNewton();

      expect(solver.rank).toBe(MX*NX + rnk);

      // check that J,F is unmodified by computeNewton
      expect(J).toBeAllCloseTo( tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ), {rtol:0, atol:0} );
      expect(F).toBeAllCloseTo(  new NDArray( F_shape, solver.F0.slice() ) );

      const   X = new NDArray( Int32Array.of(N,1), solver.newton_dX.slice() ),
             DX = zip_elems([D,X], (d,x) => d===0 ? x : d*x );
      expect(DX).toBeAllCloseTo( svd_lstsq(...svd_jac_2sided(JD),F).mapElems('float64', x => -x) );
    });


    forEachItemIn(
      function*(rng){
        for( let run=0; run++ < 73; )
        {
          const             MX = rng.int( 1, Math.ceil(42/Math.max(NX,NY)) ),
             p0 = tabulate(   [NP], 'float64', () => rng.uniform(-2,+2) ),
            dx0 = tabulate([MX,NX], 'float64', () => rng.uniform(-2,+2) );
          yield [rng, p0,dx0]
        }
      }
    ).it(`initial computeNewton() solves random sparse examples`, ([rng, p0,dx0]) => {
      const        solver = new TrustRegionSolverTLS(fgg, p0,dx0),
        {M,N,MX} = solver;

      for( let i=MX*N    ; i-- > 0; ) solver.J11[i] = rng.uniform(0.1, 1.9);
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = rng.uniform( -2, +2 ) * (rng.uniform(0,1) < 0.95);
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = rng.uniform( -2, +2 ) * (rng.uniform(0,1) < 0.95);
      for( let i=M       ; i-- > 0; ) solver. F0[i] = rng.uniform( -2, +2 );
      for( let i=N       ; i-- > 0; ) solver.  D[i] = rng.uniform(0.01, 2 ) * (rng.uniform(0,1) < 0.99 || i < MX*NX);

      const J  = tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ),
             D = new NDArray( Int32Array.of(N,1), solver.D.slice() ),
            JD = zip_elems([J,D.T], (j,d) => d===0 ? j : j/d);

      const              F_shape = Int32Array.of(M,1),
        F = new NDArray( F_shape, solver.F0.slice() );

      solver.computeNewton();

      // check that J,F is unmodified by computeNewton
      expect(J).toBeAllCloseTo( tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ), {rtol:0, atol:0} );
      expect(F).toBeAllCloseTo(  new NDArray( F_shape, solver.F0.slice() ) );

      const   X = new NDArray( Int32Array.of(N,1), solver.newton_dX.slice() ),
             DX = zip_elems([D,X], (d,x) => d===0 ? x : d*x );
      expect(DX).toBeAllCloseTo( svd_lstsq(...svd_jac_2sided(JD),F).mapElems('float64', x => -x) );
    });


    const data_computeNewtonRegularized = Object.freeze([[
      'random over-determined examples', function*(rng) {
        for( let run=0; run++ < 32; )
        {
          const             MX = rng.int( Math.ceil(NP/NY), Math.ceil(42/Math.max(NX,NY)) ),
             p0 = tabulate(   [NP], 'float64', () => rng.uniform(-2,+2) ),
            dx0 = tabulate([MX,NX], 'float64', () => rng.uniform(-2,+2) );

          const     solver = new TrustRegionSolverTLS(fgg, p0,dx0),
            {M,N} = solver;

          for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = rng.uniform(0.1, 1.9);
          for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = rng.uniform( -2, +2 );
          for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = rng.uniform( -2, +2 );
          for( let i=M       ; i-- > 0; ) solver. F0[i] = rng.uniform( -2, +2 );
          for( let i=N       ; i-- > 0; ) solver.  D[i] = rng.uniform(0.01, 2 ) * (rng.uniform(0,1) < 0.99 || i < MX*NX);

          const         lambdas = Array.from({length: 4}, () => rng.uniform(0,2) );
                        lambdas[rng.int(0,lambdas.length)] = 0;
          if(rng.bool())lambdas[rng.int(0,lambdas.length)] = 0;
          Object.freeze(lambdas);

          yield [solver, lambdas];
        }
      }
    ],[
      'random under-determined examples', function*(rng) {
        for( let run=0; run++ < 137; )
        {
          const             MX = rng.int(1,Math.ceil(NP/NY)),
             p0 = tabulate(   [NP], 'float64', () => rng.uniform(-2,+2) ),
            dx0 = tabulate([MX,NX], 'float64', () => rng.uniform(-2,+2) );

          const     solver = new TrustRegionSolverTLS(fgg, p0,dx0),
            {M,N} = solver;

          for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = rng.uniform(0.1, 1.9);
          for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = rng.uniform( -2, +2);
          for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = rng.uniform( -2, +2);
          for( let i=M       ; i-- > 0; ) solver. F0[i] = rng.uniform( -2, +2);
          for( let i=N       ; i-- > 0; ) solver.  D[i] = rng.uniform(0.01, 2) * (rng.uniform(0,1) < 0.99 || i < MX*NX);

          const         lambdas = Array.from({length: 4}, () => rng.uniform(0,2) );
                        lambdas[rng.int(0,lambdas.length)] = 0;
          if(rng.bool())lambdas[rng.int(0,lambdas.length)] = 0;
          Object.freeze(lambdas);

          yield [solver, lambdas];
        }
      }
    ],[
      'random examples', function*(rng) {
        for( let run=0; run++ < 37; )
        {
          const             MX = rng.int( 1, Math.ceil(42/Math.max(NX,NY)) ),
             p0 = tabulate(   [NP], 'float64', () => rng.uniform(-2,+2) ),
            dx0 = tabulate([MX,NX], 'float64', () => rng.uniform(-2,+2) );

          const     solver = new TrustRegionSolverTLS(fgg, p0,dx0),
            {M,N} = solver;

          for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = rng.uniform(0.1, 1.9);
          for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = rng.uniform( -2, +2);
          for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = rng.uniform( -2, +2);
          for( let i=M       ; i-- > 0; ) solver. F0[i] = rng.uniform( -2, +2);
          for( let i=N       ; i-- > 0; ) solver.  D[i] = rng.uniform(0.01, 2) * (rng.uniform(0,1) < 0.99 || i < MX*NX);

          const         lambdas = Array.from({length: 4}, () => rng.uniform(0,2) );
                        lambdas[rng.int(0,lambdas.length)] = 0;
          if(rng.bool())lambdas[rng.int(0,lambdas.length)] = 0;
          Object.freeze(lambdas);

          yield [solver, lambdas];
        }
      }
    ]]);


    for( const [example_type,items] of data_computeNewtonRegularized )
      forEachItemIn(
        function*(rng) {
          let n=0;
          for( const item of items(rng) ) {
            if(++n > 31) break;
            yield item;
          }
        }
      ).it(`computeNewtonRegularized(λ) does not produce out-of-bounds errors solving ${example_type}`, ([solver, lambdas]) => {
        const {M,N} = solver;
      
        for( const [k,v] of Object.entries(solver) )
          if( is_array(v) )
            solver[k] = checked_array(v);
      
        const J1 = tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ),
               D = new NDArray( Int32Array.of(N,1), solver.D.slice() );
      
        const               F_shape = Int32Array.of(M,1),
          F1 = new NDArray( F_shape, solver.F0.slice() );
      
        for( const λ of lambdas )
        {
          const λSqrt = Math.sqrt(λ);
      
          solver.computeNewtonRegularized(λ);
      
          const X = new NDArray( Int32Array.of(N,1), solver.regularized_dX.slice() );
      
          // check that J,F is unmodified by computeNewton
          expect(J1).toBeAllCloseTo( tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ), {rtol:0, atol:0} );
          expect(F1).toBeAllCloseTo( new NDArray( F_shape, solver.F0.slice() ) );
          expect( D).toBeAllCloseTo( new NDArray( Int32Array.of(N,1), solver.D.slice() ) );
      
          if( 0 === λ )
          {
            const  JD = zip_elems([J1,D.T], (j,d) => d===0 ? j : j/d),
                   DX = zip_elems([D ,X  ], (d,x) => d===0 ? x : d*x);
            expect(DX).toBeAllCloseTo( lstsq(JD,F1).mapElems('float64', x => -x) );
          }
          else
          {
            const J2 = tabulate( [N,N], 'float64', (i,j) => i!==j ? 0 : (D(i,0)*λSqrt || 1) ),
                  F2 = tabulate( [N,1], 'float64',    () => 0 ),
                   J = concat([J1,J2]),
                   F = concat([F1,F2]);
      
            const  X = new NDArray( Int32Array.of(N,1), solver.regularized_dX.slice() );
            expect(X).toBeAllCloseTo( lstsq(J,F).mapElems('float64', x => -x) );
          }
        }
      });


    for( const [example_type,items] of data_computeNewtonRegularized )
      forEachItemIn(
        items
      ).it(`computeNewtonRegularized(λ) solves ${example_type}`, ([solver, lambdas]) => {
        const {M,N} = solver;
  
        const J1 = tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ),
               D = new NDArray( Int32Array.of(N,1), solver.D.slice() );
  
        const               F_shape = Int32Array.of(M,1),
          F1 = new NDArray( F_shape, solver.F0.slice() );
  
        for( const λ of lambdas )
        {
          const λSqrt = Math.sqrt(λ);
  
          solver.computeNewtonRegularized(λ);
  
          const X = new NDArray( Int32Array.of(N,1), solver.regularized_dX.slice() );
  
          // check that J,F is unmodified by computeNewton
          expect(J1).toBeAllCloseTo( tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ), {rtol:0, atol:0} );
          expect(F1).toBeAllCloseTo( new NDArray( F_shape, solver.F0.slice() ) );
          expect( D).toBeAllCloseTo( new NDArray( Int32Array.of(N,1), solver.D.slice() ) );
  
          if( 0 === λ )
          {
            const  JD = zip_elems([J1,D.T], (j,d) => d===0 ? j : j/d),
                    DX = zip_elems([D ,X  ], (d,x) => d===0 ? x : d*x);
            expect(DX).toBeAllCloseTo( lstsq(JD,F1).mapElems('float64', x => -x) );
          }
          else
          {
            const J2 = tabulate( [N,N], 'float64', (i,j) => i!==j ? 0 : (D(i,0)*λSqrt || 1) ),
                  F2 = tabulate( [N,1], 'float64',    () => 0 ),
                    J = concat([J1,J2]),
                    F = concat([F1,F2]);
  
            const  X = new NDArray( Int32Array.of(N,1), solver.regularized_dX.slice() );
            expect(X).toBeAllCloseTo( lstsq(J,F).mapElems('float64', x => -x) );
          }
        }
      });


    for( const [example_type,items] of data_computeNewtonRegularized )
      forEachItemIn(
        items
      ).it(`computeNewtonRegularized(0) returns correct [r,dr] given ${example_type}`, ([solver]) => {
        const {M,N} = solver;
  
        const fJ = x => [
          new NDArray(Int32Array.of(M), solver.F0.slice()),
          tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) )
        ];
  
        const                                   x0 = new NDArray(Int32Array.of(N), solver.X0.slice()),
        reference = new TrustRegionSolverLSQ(fJ,x0);
        reference.D.set(solver.D);
  
        const [r,dr] =    solver.computeNewtonRegularized(0),
              [R,DR] = reference.computeNewtonRegularized(0);
  
        expect(solver.newton_dX).toBeAllCloseTo(reference.newton_dX);
        expect( r).toBeAllCloseTo( R);
        expect(dr).toBeAllCloseTo(DR);
      });


    for( const [example_type,items] of data_computeNewtonRegularized )
      forEachItemIn(
        items
      ).it(`computeNewtonRegularized(λ) returns correct [r,dr] given ${example_type}`, ([solver, lambdas]) => {
        const {M,N} = solver;
  
        const fJ = x => [
          new NDArray(Int32Array.of(M), solver.F0.slice()),
          tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) )
        ];
  
        const                                   x0 = new NDArray(Int32Array.of(N), solver.X0.slice()),
        reference = new TrustRegionSolverLSQ(fJ,x0);
        reference.D.set(solver.D);
  
        for( const λ of lambdas )
        {
          const [r,dr] =    solver.computeNewtonRegularized(λ);
          const [R,DR] = reference.computeNewtonRegularized(λ);
  
          if( 0 === λ )
          expect(solver.     newton_dX).toBeAllCloseTo(reference.     newton_dX);
          expect(solver.regularized_dX).toBeAllCloseTo(reference.regularized_dX);
          expect( r).toBeAllCloseTo( R);
          expect(dr).toBeAllCloseTo(DR);
        }
      });
  });
