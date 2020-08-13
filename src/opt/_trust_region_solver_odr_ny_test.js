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
import {array, NDArray} from "../nd_array";
import {rand_normal} from "../rand_normal";
import {stack} from "../stack";
import {tabulate} from "../tabulate";
import {_rand_int,
        _rand_rankdef} from "../_test_data_generators";
import {zip_elems} from '../zip_elems';

import {matmul2} from "../la/matmul";
import {qr_decomp,
        qr_decomp_full} from "../la/qr";
import {svd_lstsq} from "../la/svd";
import {svd_jac_2sided} from "../la/svd_jac_2sided";

import {num_grad} from "./num_grad";
import {TrustRegionSolverLSQ   } from "./_trust_region_solver_lsq";
import {TrustRegionSolverODR_NY} from "./_trust_region_solver_odr_ny";


function rand_func( NX, NY, NP )
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
    coeffs = Float64Array.from({length: NY*NC}, () => Math.random()*3 - 1.5);

  const composition = [];

  for( let i=0; i < NY*NC; i++ )
  {
    let ip = _rand_int(0,NP),
        ix = _rand_int(0,NX),
        fp = atoms[_rand_int(0,atoms.length)],
        fx = atoms[_rand_int(0,atoms.length)];

    composition.push([ip,ix,fp,fx]);
  }

  const f = p => x => {
    const f = new Float64Array(NY);

    for( let i=0; i < NY; i++ )
    for( let j=0; j < NC; j++ )
    {
      const [ip,ix,fp,fx] = composition[NC*i+j];
      f[i] += coeffs[NC*i+j] * fp( p(ip) ) * fx( x(ix) );
    }

    return new NDArray( Int32Array.of(NY), f );
  };

  const fgg = p => x => {
    const f = new Float64Array(NY),
         gp = new Float64Array(NY*NP),
         gx = new Float64Array(NY*NX);

    for( let i=0; i < NY; i++ )
    for( let j=0; j < NC; j++ )
    {
      const [ip,ix,fp,fx] = composition[NC*i+j];

      f[i]        += coeffs[NC*i+j] * fp     ( p(ip) ) * fx     ( x(ix) );
      gp[NP*i+ip] += coeffs[NC*i+j] * fp.grad( p(ip) ) * fx     ( x(ix) );
      gx[NX*i+ix] += coeffs[NC*i+j] * fp     ( p(ip) ) * fx.grad( x(ix) );
    }

    return [
      new NDArray( Int32Array.of(NY),    f  ),
      new NDArray( Int32Array.of(NY,NP), gp ),
      new NDArray( Int32Array.of(NY,NX), gx )
    ];
  };

  return {
    NX, NY, NP,
    f, fgg
  };
}


const funcs = {};


for( let NX=0; NX++ < 3; )
for( let NY=0; NY++ < 4; )
  funcs[`rand_func${NX}${NY}`] = rand_func(NX,NY, (NX+NY)*2);


funcs.poly22 = {
  NP: 4,
  NX: 2,
  NY: 2,
  f: ([a,b,c,d]) => ([x,y]) => array('float64', [
    a     + b*x + c*x*y + d*y,
    a*x*y + b   + c*y   + d*x
  ]),
  fgg: ([a,b,c,d]) => ([x,y]) => {
    const f = array('float64', [
            a     + b*x + c*x*y + d*y,
            a*x*y + b   + c*y   + d*x
          ]),
          gp = array('float64', [
            [1,   x,  x*y, y],
            [x*y, 1,  y,   x]
          ]),
          gx = array('float64', [
            [b + c*y, c*x + d],
            [d + a*y, a*x + c]
          ]);
    return [f,gp,gx];
  }
};


for( const [name,{NP,NX,NY,f,fgg}] of Object.entries(funcs) )
  describe(`TrustRegionSolverODR_NY [test_fn: ${name}]`, () => {
    beforeEach( () => {
      jasmine.addMatchers(CUSTOM_MATCHERS)
    });


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 32*1024; )
          yield [
            tabulate([NP], 'float64', () => Math.random()*4 - 2),
            tabulate([NX], 'float64', () => Math.random()*4 - 2)
          ];
      }()
    ).it(`test_fn derivatives are correct`, ([p, x]) => {
      const GP = num_grad( p => f(p)(x) )(p),
            GX = num_grad( x => f(p)(x) )(x);

      const [F,gp,gx] = fgg(p)(x);

      expect(F).toBeAllCloseTo( f(p)(x) );
      expect(gp).toBeAllCloseTo(GP);
      expect(gx).toBeAllCloseTo(GX);
    });


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 2048; )
        {
          const MX = _rand_int(1,16),
            shape = [MX,NX],
                p = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
               p0 = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
              dx0 = tabulate(shape, 'float64', () => Math.random()*4 - 2);

          let x = tabulate(shape,   'float64', () => Math.random()*8 - 4),
              y = stack( [...x].map(f(p)) );

          x = x.mapElems( 'float64', x => x + rand_normal() / 8 );
          y = y.mapElems( 'float64', y => y + rand_normal() / 8 );

          yield [x,y, p0,dx0]
        }
      }()
    ).it(`initial X0,F0,G0,D,J is correct given random examples`, ([x, y, p0, dx0]) => {
      const solver = new TrustRegionSolverODR_NY(x,y, fgg, p0,dx0);

      const MX = x.shape[0];

      const                            X0 = concat([dx0.reshape(-1), p0]);
      expect(solver.X0).toBeAllCloseTo(X0);

      const F = X => {
        const p = array('float64', X.data.slice(-NP)),
             dx = new NDArray(x.shape, X.data.slice(0,-NP)),
            xdx = zip_elems([x,dx], 'float64', (x,dx) => x+dx), fx = stack( [...xdx].map(f(p)) ),
             dy = zip_elems([fx,y], (fx,y) => fx-y);

        return concat([dx.reshape(-1), dy.reshape(-1)]);
      };

      const                            F0 = F(X0);
      expect(solver.F0).toBeAllCloseTo(F0);

      const         loss = x => F(x).data.reduce((loss,r) => loss + r*r/2, 0.0);
      expect(solver.loss).toBeAllCloseTo( loss(X0) / F0.shape[0] * 2 );

      expect(solver.G0).toBeAllCloseTo(num_grad(loss)(X0), {atol:1e-5});

      const {M,N} = solver,                        J = tabulate([M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j)),
                    G0 = matmul2(F0.reshape(1,-1), J);
      expect(solver.G0).toBeAllCloseTo(G0);

      const D = J.reduceElems(0, 'float64', (x,y) => Math.hypot(x,y));
      expect(solver.D).toBeAllCloseTo(D);
    });


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 512; )
        {
          const MX = _rand_int(1,16),
            shape = [MX,NX],
                p = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
               p0 = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
              dx0 = tabulate(shape, 'float64', () => Math.random()*4 - 2);

          let x = tabulate(shape, 'float64', () => Math.random()*8 - 4),
              y = stack( [...x].map(f(p)) );

          x = x.mapElems('float64', x => x + rand_normal() / 8 );
          y = y.mapElems('float64', y => y + rand_normal() / 8 );

          function* seq() {
            for( let step=0; step++ < 8; )
              yield tabulate([MX*NX+NP,1], 'float64', () => Math.random()*4 - 2);
          }

          yield [x,y, p0,dx0, seq()]
        }
      }()
    ).it(`considerMove(dX) returns correct loss and loss prediction`, ([x, y, p0, dx0, seq]) => {
      const solver = new TrustRegionSolverODR_NY(x,y, fgg, p0,dx0),
        {M,N} = solver;

      const                            X0 = concat([dx0.reshape(-1), p0]);
      expect(solver.X0).toBeAllCloseTo(X0);

      const F = X => {
        const p = array('float64', X.data.slice(-NP)),
             dx = new NDArray(x.shape, X.data.slice(0,-NP)),
            xdx = zip_elems([x,dx], 'float64', (x,dx) => x+dx), fx = stack( [...xdx].map(f(p)) ),
             dy = zip_elems([fx,y], (fx,y) => fx-y);

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
      function*(){
        for( let run=0; run++ < 1024; )
        {
          const MX = _rand_int(1,64),
            shape = [MX,NX],
                p = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
               p0 = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
              dx0 = tabulate(shape, 'float64', () => Math.random()*4 - 2);

          let x = tabulate(shape, 'float64', () => Math.random()*8 - 4),
              y = stack( [...x].map(f(p)) );

          x = x.mapElems('float64', x => x + rand_normal() / 8 );
          y = y.mapElems('float64', y => y + rand_normal() / 8 );

          yield [x,y, p0,dx0]
        }
      }()
    ).it(`cauchyTravel() finds min. along gradient direction`, ([x, y, p0, dx0]) => {
      const solver = new TrustRegionSolverODR_NY(x,y, fgg, p0,dx0);

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
      function*(){
        for( let run=0; run++ < 1024; )
        {
          const MX = _rand_int(NP,64), // <- TODO: replace with `_rand_int(1,64)` after rank-deficient implementation is finished
             shape = [MX,NX],
                 p = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
                p0 = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
               dx0 = tabulate(shape, 'float64', () => Math.random()*4 - 2);

          let x = tabulate(shape, 'float64', () => Math.random()*8 - 4),
              y = stack( [...x].map(f(p)) );

          x = x.mapElems('float64', x => x + rand_normal() / 8 );
          y = y.mapElems('float64', y => y + rand_normal() / 8 );

          yield [x,y, p0,dx0]
        }
      }()
    ).it(`[OLD!!!] initial computeNewton() QR decomposes sparse part of J correctly given random examples`, ([x, y, p0, dx0]) => {
      const solver = new TrustRegionSolverODR_NY(x,y, fgg, p0,dx0),
          {M,N,MX} = solver;

      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = Math.random()*1.8 + 0.1;
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = Math.random()*4 - 2;
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = Math.random()*4 - 2;
      for( let i=M       ; i-- > 0; ) solver. F0[i] = Math.random()*4 - 2;

      const J = tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) );

      const              F_shape = Int32Array.of(M,1),
        F = new NDArray( F_shape, solver.F0.slice() );

      solver.computeNewton();

      // check that J,F is unmodified by computeNewton
      expect(J).toBeAllCloseTo( tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ), {rtol:0, atol:0} );
      expect(F).toBeAllCloseTo(  new NDArray( F_shape, solver.F0.slice() ) );

      const J1 = tabulate( [M,    MX*NX], 'float64', (i,j) => solver.__DEBUG_J(i,j) ),
            R1 = tabulate( [MX*NX,MX*NX], 'float64', (i,j) => solver.__DEBUG_R(i,j) ),
            QF = new NDArray( Int32Array.of(MX*NX,1), solver.QF.slice(0,MX*NX) );

      const [Q,r1] = qr_decomp(J1);

      expect(R1).toBeAllCloseTo(r1);
      expect(QF).toBeAllCloseTo( matmul2(Q.T, F) ); // <- TODO only first (MX*NX) rows should be checked
    });


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 1024; )
        {
          const MX = _rand_int(NP,64), // <- TODO: replace with `_rand_int(1,64)` after rank-deficient implementation is finished
             shape = [MX,NX],
                 p = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
                p0 = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
               dx0 = tabulate(shape, 'float64', () => Math.random()*4 - 2);

          let x = tabulate(shape, 'float64', () => Math.random()*8 - 4),
              y = stack( [...x].map(f(p)) );

          x = x.mapElems('float64', x => x + rand_normal() / 8 );
          y = y.mapElems('float64', y => y + rand_normal() / 8 );

          yield [x,y, p0,dx0]
        }
      }()
    ).it(`initial computeNewton() QR decomposes sparse part of J correctly given random examples`, ([x, y, p0, dx0]) => {
      const solver = new TrustRegionSolverODR_NY(x,y, fgg, p0,dx0),
          {M,N,MX} = solver;

      const L = Math.min(NX,NY);

      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = Math.random()*1.8 + 0.1;
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = Math.random()*4 - 2;
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = Math.random()*4 - 2;
      for( let i=M       ; i-- > 0; ) solver. F0[i] = Math.random()*4 - 2;

      const J = tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) );

      const              F_shape = Int32Array.of(M,1),
        F = new NDArray( F_shape, solver.F0.slice() );

      solver.computeNewton();

      // check that J,F is unmodified by computeNewton
      expect(J).toBeAllCloseTo( tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ), {rtol:0, atol:0} );
      expect(F).toBeAllCloseTo(  new NDArray( F_shape, solver.F0.slice() ) );

      expect(solver.rank).toBe( Math.min(M,N) );

      if( J.data.some(isNaN) ) throw new Error('Assertion failed.');

      const [r,V] = solver.__DEBUG_RV,
            [Q,R] = qr_decomp_full( matmul2(J,V.T) );

      if( r.data.some(isNaN) ) throw new Error('Assertion failed.');
      if( V.data.some(isNaN) ) throw new Error('Assertion failed.');

      expect(r).toBeAllCloseTo(R);

      const                    {rank} = solver;
      expect( solver.QF.slice(0,rank) ).toBeAllCloseTo( matmul2(Q.T,F).data.slice(0,rank) );
    });


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 3*512; )
        {
          const MX = _rand_int(NP,32),
             shape = [MX,NX],
                 p = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
                p0 = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
               dx0 = tabulate(shape, 'float64', () => Math.random()*4 - 2);

          let x = tabulate(shape, 'float64', () => Math.random()*8 - 4),
              y = stack( [...x].map(f(p)) );

          x = x.mapElems('float64', x => x + rand_normal() / 8 );
          y = y.mapElems('float64', y => y + rand_normal() / 8 );

          yield [x,y, p0,dx0]
        }
      }()
    ).it(`initial computeNewton() solves generated over-determined examples`, ([x, y, p0, dx0]) => {
      const solver = new TrustRegionSolverODR_NY(x,y, fgg, p0,dx0),
          {M,N,MX} = solver;

      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = Math.random()*1.8 + 0.1;
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = Math.random()*4 - 2;
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = Math.random()*4 - 2;
      for( let i=M       ; i-- > 0; ) solver. F0[i] = Math.random()*4 - 2;

      const J = tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) );

      const              F_shape = Int32Array.of(M,1),
        F = new NDArray( F_shape, solver.F0.slice() );

      solver.computeNewton();

      // expect(solver.rank).toBe( Math.min(M,N) );

      // check that J,F is unmodified by computeNewton
      expect(J).toBeAllCloseTo( tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) ), {rtol:0, atol:0} );
      expect(F).toBeAllCloseTo(  new NDArray( F_shape, solver.F0.slice() ) );

      const  X = new NDArray( Int32Array.of(N,1), solver.newton_dX.slice() );
      expect(X).toBeAllCloseTo( svd_lstsq(...svd_jac_2sided(J),F).mapElems('float64', x => -x) );
    });


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 1024; )
        {
          const MX = _rand_int(1,64),
             shape = [MX,NX],
                 p = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
                p0 = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
               dx0 = tabulate(shape, 'float64', () => Math.random()*4 - 2);

          let x = tabulate(shape, 'float64', () => Math.random()*8 - 4),
              y = stack( [...x].map(f(p)) );

          x = x.mapElems('float64', x => x + rand_normal() / 8 );
          y = y.mapElems('float64', y => y + rand_normal() / 8 );

          yield [x,y, p0,dx0]
        }
      }()
    ).it(`initial computeNewton() solves generated random examples with J21=0`, ([x, y, p0, dx0]) => {
      const solver = new TrustRegionSolverODR_NY(x,y, fgg, p0,dx0),
          {M,N,MX} = solver;

      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = Math.random()*1.8 + 0.1;
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = 0;
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = Math.random()*4 - 2;
      for( let i=M       ; i-- > 0; ) solver. F0[i] = Math.random()*4 - 2;
      for( let i=N       ; i-- > 0; ) solver.  D[i] =(Math.random() < 0.99 || i < MX*NX) * (0.01 + Math.random()*2);

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
      expect(DX).toBeAllCloseTo( svd_lstsq(...svd_jac_2sided(JD),F).mapElems('float64', x => -x) );
    });


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 1024; )
        {
          const MX = _rand_int(1,64),
             shape = [MX,NX],
                 p = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
                p0 = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
               dx0 = tabulate(shape, 'float64', () => Math.random()*4 - 2);

          let x = tabulate(shape, 'float64', () => Math.random()*8 - 4),
              y = stack( [...x].map(f(p)) );

          x = x.mapElems('float64', x => x + rand_normal() / 8 );
          y = y.mapElems('float64', y => y + rand_normal() / 8 );

          yield [x,y, p0,dx0]
        }
      }()
    ).it(`initial computeNewton() solves generated random rank-deficient examples with J21=0`, ([x, y, p0, dx0]) => {
      const solver = new TrustRegionSolverODR_NY(x,y, fgg, p0,dx0),
          {M,N,MX} = solver;

      const [{data: J22},rnk] = _rand_rankdef(MX*NY,NP);
      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = Math.random()*1.8 + 0.1;
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = 0;
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = J22[i];
      for( let i=M       ; i-- > 0; ) solver. F0[i] = Math.random()*4 - 2;
      for( let i=N       ; i-- > 0; ) solver.  D[i] = (Math.random() < 0.99 || i < MX*NX) * (0.01 + Math.random()*2);

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
      function*(){
        for( let run=0; run++ < 3*512; )
        {
          const MX = _rand_int(1,32),
             shape = [MX,NX],
                 p = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
                p0 = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
               dx0 = tabulate(shape, 'float64', () => Math.random()*4 - 2);

          let x = tabulate(shape, 'float64', () => Math.random()*8 - 4),
              y = stack( [...x].map(f(p)) );

          x = x.mapElems('float64', x => x + rand_normal() / 8 );
          y = y.mapElems('float64', y => y + rand_normal() / 8 );

          yield [x,y, p0,dx0]
        }
      }()
    ).it(`initial computeNewton() solves generated random examples`, ([x, y, p0, dx0]) => {
      const solver = new TrustRegionSolverODR_NY(x,y, fgg, p0,dx0),
          {M,N,MX} = solver;

      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = Math.random()*1.8 + 0.1;
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = Math.random()*4 - 2;
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = Math.random()*4 - 2;
      for( let i=M       ; i-- > 0; ) solver. F0[i] = Math.random()*4 - 2;
      for( let i=N       ; i-- > 0; ) solver.  D[i] =(Math.random() < 0.99 || i < MX*NX) * (0.01 + Math.random()*2);

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
      expect(DX).toBeAllCloseTo( svd_lstsq(...svd_jac_2sided(JD),F).mapElems('float64', x => -x) );
    });


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 2*1024; )
        {
          const MX = _rand_int(1,32),
             shape = [MX,NX],
                 p = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
                p0 = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
               dx0 = tabulate(shape, 'float64', () => Math.random()*4 - 2);

          let x = tabulate(shape, 'float64', () => Math.random()*8 - 4),
              y = stack( [...x].map(f(p)) );

          x = x.mapElems('float64', x => x + rand_normal() / 8 );
          y = y.mapElems('float64', y => y + rand_normal() / 8 );

          yield [x,y, p0,dx0]
        }
      }()
    ).it(`initial computeNewton() solves generated random rank-deficient examples`, ([x, y, p0, dx0]) => {
      const solver = new TrustRegionSolverODR_NY(x,y, fgg, p0,dx0),
          {M,N,MX} = solver;

      const [{data: A},rnk] = _rand_rankdef(MX*NY,NP);

      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = Math.random()*1.8 + 0.1;
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = Math.random()*4 - 2;
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = A[i];
      for( let i=M       ; i-- > 0; ) solver. F0[i] = Math.random()*4 - 2;
      for( let i=N       ; i-- > 0; ) solver.  D[i] =(Math.random() < 0.99 || i < MX*NX) * (0.01 + Math.random()*2);

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
      function*(){
        for( let run=0; run++ < 3*512; )
        {
          const MX = _rand_int(1,32),
             shape = [MX,NX],
                 p = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
                p0 = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
               dx0 = tabulate(shape, 'float64', () => Math.random()*4 - 2);

          let x = tabulate(shape, 'float64', () => Math.random()*8 - 4),
              y = stack( [...x].map(f(p)) );

          x = x.mapElems('float64', x => x + rand_normal() / 8 );
          y = y.mapElems('float64', y => y + rand_normal() / 8 );

          yield [x,y, p0,dx0]
        }
      }()
    ).it(`initial computeNewton() solves generated random sparse examples`, ([x, y, p0, dx0]) => {
      const solver = new TrustRegionSolverODR_NY(x,y, fgg, p0,dx0),
          {M,N,MX} = solver;

      for( let i=MX*N    ; i-- > 0; ) solver.J11[i] = Math.random()*1.8 + 0.1;
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] =(Math.random()*4 - 2) * (Math.random() < 0.95);
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] =(Math.random()*4 - 2) * (Math.random() < 0.95);
      for( let i=M       ; i-- > 0; ) solver. F0[i] = Math.random()*4 - 2;
      for( let i=N       ; i-- > 0; ) solver.  D[i] =(Math.random() < 0.99 || i < MX*NX) * (0.01 + Math.random()*2);

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


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 1024; )
        {
          const MX = _rand_int(1,16),
             shape = [MX,NX],
                 p = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
                p0 = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
               dx0 = tabulate(shape, 'float64', () => Math.random()*4 - 2);

          let x = tabulate(shape, 'float64', () => Math.random()*8 - 4),
              y = stack( [...x].map(f(p)) );

          x = x.mapElems('float64', x => x + rand_normal() / 8 );
          y = y.mapElems('float64', y => y + rand_normal() / 8 );

          const lambdas = Array.from({length: 6}, () => Math.random()*2);
          lambdas[_rand_int(0,lambdas.length)] = 0;
          lambdas[_rand_int(0,lambdas.length)] = 0;
          Object.freeze(lambdas);

          yield [x,y, p0,dx0, lambdas];
        }
      }()
    ).it(`computeNewtonRegularized(λ) solves generated random examples`, ([x, y, p0, dx0, lambdas]) => {
      const solver = new TrustRegionSolverODR_NY(x,y, fgg, p0,dx0),
          {M,N,MX} = solver;

      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = Math.random()*1.8 + 0.1;
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = Math.random()*4 - 2;
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = Math.random()*4 - 2;
      for( let i=M       ; i-- > 0; ) solver. F0[i] = Math.random()*4 - 2;
      for( let i=N       ; i-- > 0; ) solver.  D[i] =(Math.random() < 0.99 || i < MX*NX) * (0.01 + Math.random()*2);

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
          expect(DX).toBeAllCloseTo( svd_lstsq(...svd_jac_2sided(JD),F1).mapElems('float64', x => -x) );
        }
        else
        {
          const J2 = tabulate( [N,N], 'float64', (i,j) => i!==j ? 0 : (D(i,0)*λSqrt || 1) ),
                F2 = tabulate( [N,1], 'float64',    () => 0 ),
                 J = concat([J1,J2]),
                 F = concat([F1,F2]);

          const  X = new NDArray( Int32Array.of(N,1), solver.regularized_dX.slice() );
          expect(X).toBeAllCloseTo( svd_lstsq(...svd_jac_2sided(J),F).mapElems('float64', x => -x) );
        }
      }
    });


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 1024; )
        {
          const MX = _rand_int(1,16),
             shape = [MX,NX],
                 p = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
                p0 = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
               dx0 = tabulate(shape, 'float64', () => Math.random()*4 - 2);

          let x = tabulate(shape, 'float64', () => Math.random()*8 - 4),
              y = stack( [...x].map(f(p)) );

          x = x.mapElems('float64', x => x + rand_normal() / 8 );
          y = y.mapElems('float64', y => y + rand_normal() / 8 );

          const lambdas = Array.from({length: 6}, () => Math.random()*2);
          lambdas[_rand_int(0,lambdas.length)] = 0;
          lambdas[_rand_int(0,lambdas.length)] = 0;
          Object.freeze(lambdas);

          yield [x,y, p0,dx0, lambdas];
        }
      }()
    ).it(`computeNewtonRegularized(λ) solves generated random rank-deficient examples`, ([x, y, p0, dx0, lambdas]) => {
      const solver = new TrustRegionSolverODR_NY(x,y, fgg, p0,dx0),
          {M,N,MX} = solver;

      const [{data: J22}] = _rand_rankdef(MX*NY,NP);
      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = Math.random()*1.8 + 0.1;
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = Math.random()*4 - 2;
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = J22[i];
      for( let i=M       ; i-- > 0; ) solver. F0[i] = Math.random()*4 - 2;
      for( let i=N       ; i-- > 0; ) solver.  D[i] =(Math.random() < 0.99 || i < MX*NX) * (0.01 + Math.random()*2);

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
          expect(DX).toBeAllCloseTo( svd_lstsq(...svd_jac_2sided(JD),F1).mapElems('float64', x => -x) );
        }
        else
        {
          const J2 = tabulate( [N,N], 'float64', (i,j) => i!==j ? 0 : (D(i,0)*λSqrt || 1) ),
                F2 = tabulate( [N,1], 'float64',    () => 0 ),
                J = concat([J1,J2]),
                F = concat([F1,F2]);

          const  X = new NDArray( Int32Array.of(N,1), solver.regularized_dX.slice() );
          expect(X).toBeAllCloseTo( svd_lstsq(...svd_jac_2sided(J),F).mapElems('float64', x => -x) );
        }
      }
    });


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 4096; )
        {
          const MX = _rand_int(NP,32),
             shape = [MX,NX],
                 p = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
                p0 = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
               dx0 = tabulate(shape, 'float64', () => Math.random()*4 - 2);

          let x = tabulate(shape, 'float64', () => Math.random()*8 - 4),
              y = stack( [...x].map(f(p)) );

          x = x.mapElems('float64', x => x + rand_normal() / 8 );
          y = y.mapElems('float64', y => y + rand_normal() / 8 );

          yield [x,y, p0,dx0];
        }
      }()
    ).it(`computeNewtonRegularized(0) returns correct [r,dr] given random over-determined examples`, ([x, y, p0, dx0]) => {
      const solver = new TrustRegionSolverODR_NY(x,y, fgg, p0,dx0),
          {M,N,MX} = solver;

      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = Math.random()*1.8 + 0.1;
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = Math.random()*4.0 - 2.0;
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = Math.random()*4.0 - 2.0;
      for( let i=M       ; i-- > 0; ) solver. F0[i] = Math.random()*4.0 - 2.0;
      for( let i=N       ; i-- > 0; ) solver.  D[i] = Math.random()*1.8 + 0.1;

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


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 4096; )
        {
          const MX = _rand_int(1,NP),
             shape = [MX,NX],
                 p = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
                p0 = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
               dx0 = tabulate(shape, 'float64', () => Math.random()*4 - 2);

          let x = tabulate(shape, 'float64', () => Math.random()*8 - 4),
              y = stack( [...x].map(f(p)) );

          x = x.mapElems('float64', x => x + rand_normal() / 8 );
          y = y.mapElems('float64', y => y + rand_normal() / 8 );

          yield [x,y, p0,dx0];
        }
      }()
    ).it(`computeNewtonRegularized(0) returns correct [r,dr] given random under-determined examples`, ([x, y, p0, dx0]) => {
      const solver = new TrustRegionSolverODR_NY(x,y, fgg, p0,dx0),
          {M,N,MX} = solver;

      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = Math.random()*1.8 + 0.1;
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = Math.random()*4.0 - 2.0;
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = Math.random()*4.0 - 2.0;
      for( let i=M       ; i-- > 0; ) solver.F0 [i] = Math.random()*4.0 - 2.0;
      for( let i=N       ; i-- > 0; ) solver.D  [i] = Math.random()*1.8 + 0.1;

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


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 1024; )
        {
          const MX = _rand_int(1,32),
             shape = [MX,NX],
                 p = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
                p0 = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
               dx0 = tabulate(shape, 'float64', () => Math.random()*4 - 2);

          let x = tabulate(shape, 'float64', () => Math.random()*8 - 4),
              y = stack( [...x].map(f(p)) );

          x = x.mapElems('float64', x => x + rand_normal() / 8 );
          y = y.mapElems('float64', y => y + rand_normal() / 8 );

          yield [x,y, p0,dx0];
        }
      }()
    ).it(`computeNewtonRegularized(0) returns correct [r,dr] given random rank-deficient examples`, ([x, y, p0, dx0]) => {
      const solver = new TrustRegionSolverODR_NY(x,y, fgg, p0,dx0);

      const {M,N,MX} = solver;

      const [{data: J22},rnk] = _rand_rankdef(MX*NY,NP);
      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = Math.random()*1.8 + 0.1;
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = Math.random()*4.0 - 2.0;
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = J22[i];
      for( let i=M       ; i-- > 0; ) solver. F0[i] = Math.random()*4.0 - 2.0;
      for( let i=N       ; i-- > 0; ) solver.  D[i] = Math.random()*1.8 + 0.1;

      const fJ = x => [
        new NDArray(Int32Array.of(M), solver.F0.slice()),
        tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) )
      ];

      const x0 = new NDArray(Int32Array.of(N), solver.X0.slice());

      const reference = new TrustRegionSolverLSQ(fJ,x0);
      reference.D.set(solver.D);

      const [r,dr] =    solver.computeNewtonRegularized(0);
      const [R,DR] = reference.computeNewtonRegularized(0);

      expect(   solver.rank).toBe(MX*NX + rnk);
      expect(reference.rank).toBe(MX*NX + rnk);

      expect(solver.newton_dX).toBeAllCloseTo(reference.newton_dX);
      expect( r).toBeAllCloseTo( R);
      expect(dr).toBeAllCloseTo(DR);
    });


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 4096; )
        {
          const MX = _rand_int(NP,32),
             shape = [MX,NX],
                 p = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
                p0 = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
               dx0 = tabulate(shape, 'float64', () => Math.random()*4 - 2);

          let x = tabulate(shape, 'float64', () => Math.random()*8 - 4),
              y = stack( [...x].map(f(p)) );

          x = x.mapElems('float64', x => x + rand_normal() / 8 );
          y = y.mapElems('float64', y => y + rand_normal() / 8 );

          const lambdas = Array.from({length: 6}, () => Math.random()*2);
          lambdas[_rand_int(0,lambdas.length)] = 0;
          lambdas[_rand_int(0,lambdas.length)] = 0;
          Object.freeze(lambdas);

          yield [x,y, p0,dx0, lambdas];
        }
      }()
    ).it(`computeNewtonRegularized(λ) returns correct [r,dr] given random over-determined examples`, ([x, y, p0, dx0, lambdas]) => {
      const solver = new TrustRegionSolverODR_NY(x,y, fgg, p0,dx0),
          {M,N,MX} = solver;

      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = Math.random()*1.8 + 0.1;
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = Math.random()*4.0 - 2.0;
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = Math.random()*4.0 - 2.0;
      for( let i=M       ; i-- > 0; ) solver.F0 [i] = Math.random()*4.0 - 2.0;
      for( let i=N       ; i-- > 0; ) solver.D  [i] = Math.random()*1.8 + 0.1;

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


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 4096; )
        {
          const MX = _rand_int(1,NP),
             shape = [MX,NX],
                 p = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
                p0 = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
               dx0 = tabulate(shape, 'float64', () => Math.random()*4 - 2);

          let x = tabulate(shape, 'float64', () => Math.random()*8 - 4),
              y = stack( [...x].map(f(p)) );

          x = x.mapElems('float64', x => x + rand_normal() / 8 );
          y = y.mapElems('float64', y => y + rand_normal() / 8 );

          const lambdas = Array.from({length: 6}, () => Math.random()*2);
          lambdas[_rand_int(0,lambdas.length)] = 0;
          lambdas[_rand_int(0,lambdas.length)] = 0;
          Object.freeze(lambdas);

          yield [x,y, p0,dx0, lambdas];
        }
      }()
    ).it(`computeNewtonRegularized(λ) returns correct [r,dr] given random under-determined examples`, ([x, y, p0, dx0, lambdas]) => {
      const solver = new TrustRegionSolverODR_NY(x,y, fgg, p0,dx0),
          {M,N,MX} = solver;

      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = Math.random()*1.8 + 0.1;
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = Math.random()*4.0 - 2.0;
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = Math.random()*4.0 - 2.0;
      for( let i=M       ; i-- > 0; ) solver. F0[i] = Math.random()*4.0 - 2.0;
      for( let i=N       ; i-- > 0; ) solver.  D[i] = Math.random()*1.8 + 0.1;
   
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


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 1024; )
        {
          const MX = _rand_int(1,32),
             shape = [MX,NX],
                 p = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
                p0 = tabulate( [NP], 'float64', () => Math.random()*4 - 2),
               dx0 = tabulate(shape, 'float64', () => Math.random()*4 - 2);

          let x = tabulate(shape, 'float64', () => Math.random()*8 - 4),
              y = stack( [...x].map(f(p)) );

          x = x.mapElems('float64', x => x + rand_normal() / 8 );
          y = y.mapElems('float64', y => y + rand_normal() / 8 );

          const lambdas = Array.from({length: 6}, () => Math.random()*2);
          lambdas[_rand_int(0,lambdas.length)] = 0;
          lambdas[_rand_int(0,lambdas.length)] = 0;
          Object.freeze(lambdas);

          yield [x,y, p0,dx0, lambdas];
        }
      }()
    ).it(`computeNewtonRegularized(λ) returns correct [r,dr] given random rank-deficient examples`, ([x, y, p0, dx0, lambdas]) => {
      const solver = new TrustRegionSolverODR_NY(x,y, fgg, p0,dx0),
          {M,N,MX} = solver;

      const [{data: J22},rnk] = _rand_rankdef(MX*NY,NP);
      for( let i=MX*NX   ; i-- > 0; ) solver.J11[i] = Math.random()*1.8 + 0.1;
      for( let i=MX*NY*NX; i-- > 0; ) solver.J21[i] = Math.random()*4.0 - 2.0;
      for( let i=MX*NY*NP; i-- > 0; ) solver.J22[i] = J22[i];
      for( let i=M       ; i-- > 0; ) solver. F0[i] = Math.random()*4.0 - 2.0;
      for( let i=N       ; i-- > 0; ) solver.  D[i] = Math.random()*1.8 + 0.1;

      const fJ = x => [
        new NDArray(Int32Array.of(M), solver.F0.slice()),
        tabulate( [M,N], 'float64', (i,j) => solver.__DEBUG_J(i,j) )
      ];

      const                                   x0 = new NDArray(Int32Array.of(N), solver.X0.slice()),
      reference = new TrustRegionSolverLSQ(fJ,x0);
      reference.D.set(solver.D);

      for( const λ of lambdas )
      {
        const [r,dr] =    solver.computeNewtonRegularized(λ),
              [R,DR] = reference.computeNewtonRegularized(λ);

        if( 0 === λ )
        expect(   solver.rank).toBe(MX*NX + rnk);
        expect(reference.rank).toBe(MX*NX + rnk);

        if( 0 === λ )
        expect(solver.     newton_dX).toBeAllCloseTo(reference.     newton_dX);
        expect(solver.regularized_dX).toBeAllCloseTo(reference.regularized_dX);
        expect( r).toBeAllCloseTo( R);
        expect(dr).toBeAllCloseTo(DR);
      }
    });
  });
