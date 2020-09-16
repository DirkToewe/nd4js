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
import {array, NDArray} from '../nd_array';
import {stack} from "../stack";
import {tabulate} from "../tabulate";

import {norm} from "../la/norm";

import {num_grad} from "./num_grad";
import {OptimizationNoProgressError} from './optimization_error';


const rand_poly_lin_init = rng => function rand_poly_lin( NX, NY )
{
  if( 0 !== NX%1 ) throw new Error('Assertion failed.');
  if( !(0 < NX)  ) throw new Error('Assertion failed.');

  const NP = 1 + NX;

  if( NY !== undefined )
  {
    if( 0 !== NY%1 ) throw new Error('Assertion failed.');
    if( !(0 < NY)  ) throw new Error('Assertion failed.');

    const f = Array.from({length: NY}, () => rand_poly_lin(NX));

    const result = p => {
      const fp = f.map( f => f(p) );

      return x => {
        const y   = new Float64Array(NY),
             dydp = new Float64Array(NY*NP),
             dydx = new Float64Array(NY*NX)
        for( let i=NY; i-- > 0; )
        {
          const [f,dfdp,dfdx] = fp[i](x);
                                      y  [   i  ] =  f;
          for( let j=NP; j-- > 0; )  dydp[NP*i+j] = dfdp[j];
          for( let j=NX; j-- > 0; )  dydx[NX*i+j] = dfdx[j];
        }
        return [
          new NDArray(Int32Array.of(NY),     y  ),
          new NDArray(Int32Array.of(NY,NP), dydp),
          new NDArray(Int32Array.of(NY,NX), dydx)
        ];
      };
    };

    Object.        assign(result, {NX,NY,NP});
    Object.defineProperty(result, 'name', {value: `poly_lin${NX}${NY}`, writable: false});
    Object.        freeze(result);
    return                result;
  }

  const c = new Float64Array(3*NP);
  for( let i=c.length; (i-=3) >= 0; )
  {
    c[i+0] = rng.uniform(-2  ,+2  );
    c[i+1] = rng.uniform( 0.1, 1.9);
    c[i+2] = rng.normal(0, 1e-5);
    if( rng.bool() ) {
      c[i+1] *= -1;
      c[i+2] *= -1;
    }
  }

  const fgh = p =>
  {
    if( p.ndim     !== 1  ) throw new Error('Assertion failed.');
    if( p.shape[0] !== NP ) throw new Error('Assertion failed.');
    p = p.data;

    return x => {
      if( x.ndim     !== 1  ) throw new Error('Assertion failed.');
      if( x.shape[0] !== NX ) throw new Error('Assertion failed.');
      x = x.data;

      const dfdp = new Float64Array(NP),
            dfdx = new Float64Array(NX);
      
      let  f  = c[0] + (c[1] + c[2] * p[0]*p[0])*p[0]; // <- make the coefficients of the polynomial non-linear in p by using a monotonous polynomials
      dfdp[0] =         c[1] + c[2] * p[0]*p[0]*3;

      for( let i=0; i < NX; i++ )
      {
        const pi = p[(1+i)],
              c0 = c[(1+i)*3 + 0],
              c1 = c[(1+i)*3 + 1],
              c2 = c[(1+i)*3 + 2];
               f  += x[i] * ( c0 + (c1 + c2*pi*pi)*pi );
        dfdp[1+i] += x[i] * (      (c1 + c2*pi*pi*3)  );
        dfdx[  i] +=        ( c0 + (c1 + c2*pi*pi)*pi );
      }

      return [f,dfdp,dfdx];
    };
  };
  Object.        assign(fgh,{NX,NP});
  Object.defineProperty(fgh, 'name', {value: `poly_lin${NX}`, writable: false});
  Object.        freeze(fgh);
  return                fgh;
}


/* Creates a quadratic polynomial with NX-dimensional input `x`.
 * To make the polynomial non-linear in its parameters `p`,
 * each polynomial coeffcient is itself a cubic, monotonous polynomial
 * of some entry of `p`.
 */
const rand_poly_quad_init = rng => function rand_poly_quad( NX, NY )
{
  if( 0 !== NX%1 ) throw new Error('Assertion failed.');
  if( !(0 < NX)  ) throw new Error('Assertion failed.');

  const NP = 1 + NX + NX*(NX+1)/2;

  if( NY !== undefined )
  {
    if( 0 !== NY%1 ) throw new Error('Assertion failed.');
    if( !(0 < NY)  ) throw new Error('Assertion failed.');

    const f = Array.from({length: NY}, () => rand_poly_quad(NX));

    const result = p => {
      const fp = f.map( f => f(p) );

      return x => {
        const y   = new Float64Array(NY),
             dydp = new Float64Array(NY*NP),
             dydx = new Float64Array(NY*NX)
        for( let i=NY; i-- > 0; )
        {
          const [f,dfdp,dfdx] = fp[i](x);
                                      y  [   i  ] =  f;
          for( let j=NP; j-- > 0; )  dydp[NP*i+j] = dfdp[j];
          for( let j=NX; j-- > 0; )  dydx[NX*i+j] = dfdx[j];
        }
        return [
          new NDArray(Int32Array.of(NY),     y  ),
          new NDArray(Int32Array.of(NY,NP), dydp),
          new NDArray(Int32Array.of(NY,NX), dydx)
        ];
      };
    };

    Object.        assign(result, {NX,NY,NP});
    Object.defineProperty(result, 'name', {value: `poly_lin${NX}${NY}`, writable: false});
    Object.        freeze(result);
    return                result;
  }

  const c = new Float64Array(3*NP);
  for( let i=c.length; (i-=3) >= 0; )
  {
    c[i+0] = rng.uniform(-2  ,+2  );
    c[i+1] = rng.uniform( 0.1, 1.9);
    c[i+2] = rng.normal(0,1e-5);
    if( rng.bool() ) {
      c[i+1] *= -1;
      c[i+2] *= -1;
    }
  }

  const fgh = p =>
  {
    if( p.ndim     !== 1  ) throw new Error('Assertion failed.');
    if( p.shape[0] !== NP ) throw new Error('Assertion failed.');
    p = p.data;

    return x => {
      if( x.ndim     !== 1  ) throw new Error('Assertion failed.');
      if( x.shape[0] !== NX ) throw new Error('Assertion failed.');
      x = x.data;

      const dfdp = new Float64Array(NP),
            dfdx = new Float64Array(NX);
      
      let  f  = c[0] + (c[1] + c[2] * p[0]*p[0])*p[0]; // <- make the coefficients of the polynomial non-linear in p by using a monotonous polynomials
      dfdp[0] =         c[1] + c[2] * p[0]*p[0]*3;

      for( let i=0; i < NX; i++ )
      {
        const pi = p[(1+i)],
              c0 = c[(1+i)*3 + 0],
              c1 = c[(1+i)*3 + 1],
              c2 = c[(1+i)*3 + 2];
               f  += x[i] * ( c0 + (c1 + c2*pi*pi)*pi );
        dfdp[1+i] += x[i] * (      (c1 + c2*pi*pi*3)  );
        dfdx[  i] +=        ( c0 + (c1 + c2*pi*pi)*pi );
      }

      let I = 1+NX;
      for( let i=0; i < NX; i++ )
      for( let j=0; j <= i; j++, I++ )
      {
        const pi = p[I],
              c0 = c[I*3 + 0],
              c1 = c[I*3 + 1],
              c2 = c[I*3 + 2];
             f  += x[i]*x[j] * ( c0 + (c1 + c2*pi*pi)*pi );
        dfdp[I] += x[i]*x[j] * (      (c1 + c2*pi*pi*3)  );
        dfdx[i] +=      x[j] * ( c0 + (c1 + c2*pi*pi)*pi );
        dfdx[j] += x[i]      * ( c0 + (c1 + c2*pi*pi)*pi );
      }
      expect(I).toBe(NP);

      return [f,dfdp,dfdx];
    };
  };
  Object.        assign(fgh,{NX,NP});
  Object.defineProperty(fgh, 'name', {value: `poly_quad${NX}`, writable: false});
  Object.        freeze(fgh);
  return                fgh;
}


export function generic_test_odr_gen( odr_gen )
{
  describe(`${odr_gen.name} [generic tests]`, () => {
    beforeEach( () => {
      jasmine.addMatchers(CUSTOM_MATCHERS)
    })

    for( const rand_func_init of [
      rand_poly_lin_init,
      rand_poly_quad_init
    ])
    {
      for( const NX of [1,2,3] )
        forEachItemIn(
          function*(rng){
            const rand_func = rand_func_init(rng);

            for( let run=0; run++ < 73; )
            {
              const fgh = rand_func(NX);

              for( let step=0; step++ < 73/NX; )
              {
                const x = tabulate([fgh.NX], 'float64', () => rng.uniform(-2,+2) ),
                      p = tabulate([fgh.NP], 'float64', () => rng.uniform(-2,+2) );
                yield [fgh,p,x];
              }
            }
          }
        ).it(`${rand_func_init.name}(rng)(${NX}) gradients are correct`, ([fgh, p, x]) => {
          const  F  = (p,x) => fgh(p)(x)[0],
            [f,g,h] = fgh(p)(x);

          expect(f).toBe( F(p,x) );
          expect(g).toBeAllCloseTo( num_grad( p => F(p,x) )(p) );
          expect(h).toBeAllCloseTo( num_grad( x => F(p,x) )(x) );
        });

      for( const NX of [1,2,3] )
      for( const NY of [1,2,3] )
        forEachItemIn(
          function*(rng){
            const rand_func = rand_func_init(rng);

            for( let run=0; run++ < 73; )
            {
              const fgh = rand_func(NX,NY);

              for( let step=0; step++ < 73/NX; )
              {
                const x = tabulate([fgh.NX], 'float64', () => rng.uniform(-2,+2) ),
                      p = tabulate([fgh.NP], 'float64', () => rng.uniform(-2,+2) );
                yield [fgh,p,x];
              }
            }
          }
        ).it(`${rand_func_init.name}(rng)(${NX},${NY}) gradients are correct`, ([fgh, p, x]) => {
          const  F  = (p,x) => fgh(p)(x)[0],
            [f,g,h] = fgh(p)(x);

          expect(f).toBeAllCloseTo( F(p,x) );
          expect(g).toBeAllCloseTo( num_grad( p => F(p,x) )(p) );
          expect(h).toBeAllCloseTo( num_grad( x => F(p,x) )(x) );
        });
    }

    for( const NX of [1,2,3] )
      forEachItemIn(
        function*(rng){
          const rand_poly_quad = rand_poly_quad_init(rng);

          const                   N_RUNS = 5 / NX**1.5; expect(N_RUNS).toBeGreaterThan(0);
          for( let run=0; run++ < N_RUNS; )
          {
            const fgh = rand_poly_quad(NX),
                  f   = (p,x) => fgh(p)(x)[0];

            const {NP} = fgh,
                   MX  = rng.int(1024*NP, 3072*NP);

            const p  = tabulate([NP], 'float64', () => rng.uniform(-1.5,+1.5) ),
                  p0 = tabulate([NP], 'float64', () => rng.uniform(-1.5,+1.5) );

            const samples_x = tabulate([MX,NX], 'float64', () => rng.uniform(-8,+8) ),     noise = 1e-6,
                  samples_y = array('float64', [...samples_x].map( x => rng.normal(f(p,x), noise) ) );

            for( let i=samples_x.data.length; i-- > 0; )
              samples_x.data[i] += rng.normal(0,noise);

            const dx0 = samples_x.mapElems('float64', rng.bool() ? () => 0
                                                                 : () => rng.normal(0,noise));
            yield [samples_x,samples_y, fgh,p, p0,dx0];
          }
        }
      ).it(`fits rand_poly_quad(${NX}) correctly`, ([samples_x, samples_y, fgh, P, P0, dx0]) => {

        let     p,dx, mse,dmse_dp,dmse_dx, dy, nIter=0;
        try {
          for( [p,dx, mse,dmse_dp,dmse_dx, dy] of odr_gen(samples_x, samples_y, fgh, P0, {dx0}) )
          {
            if( ++nIter >= 256 )
              throw new Error('Too many iterations.');

            if(   norm(dmse_dp) < 1e-8
               && norm(dmse_dx) < 1e-8 )
              break;

            // TODO test mse, dmse_dx, dymse_dp, dy, ...
          }
        }
        catch( err ) {
          if( ! (err instanceof OptimizationNoProgressError) )
            throw err;
        }

        expect(p).toBeAllCloseTo(P, {rtol: 1e-4, atol: 1e-4});
      });

    for( const NX of [1,2,3] )
    for( const NY of [1,2,3] )
      forEachItemIn(
        function*(rng){
          const rand_poly_lin = rand_poly_lin_init(rng);

          const                   N_RUNS = 13 / NX**1.5; expect(N_RUNS).toBeGreaterThan(0);
          for( let run=0; run++ < N_RUNS; )
          {
            const fgh = rand_poly_lin(NX,NY),
                  {NP} = fgh,
                   MX  = rng.int(1024*NP/NY | 0, 2048*NP/NY | 0);

            const p  = tabulate([NP], 'float64', () => rng.uniform(-2,+2) ),
                  p0 = tabulate([NP], 'float64', () => rng.uniform(-2,+2) );

            const noise = 1e-4,
              samples_x = tabulate([MX,NX], 'float64', () => rng.uniform(-4,+4) ),
              samples_y = stack(
                [...samples_x].map(
                  x => fgh(p)(x)[0].mapElems('float64', y => rng.normal(y,noise) )
                )
              );

            for( let i=samples_x.data.length; i-- > 0; )
              samples_x.data[i] += rng.normal(0,noise);

            const dx0 = samples_x.mapElems('float64', rng.bool() ? () => 0
                                                                 : () => rng.normal(0,noise));
            yield [samples_x,samples_y, fgh,p, p0,dx0];
          }
        }
      ).it(`fits rand_poly_lin(${NX},${NY}) correctly`, ([samples_x, samples_y, fgh, P, P0, dx0]) => {

        let     p,dx, mse,dmse_dp,dmse_dx, dy, nIter=0;
        try {
          for( [p,dx, mse,dmse_dp,dmse_dx, dy] of odr_gen(samples_x, samples_y, fgh, P0, {dx0}) )
          {
            expect(nIter++).toBeLessThan(32);

            if(   norm(dmse_dp) < 1e-5
               && norm(dmse_dx) < 1e-5 )
              break;
          }
        }
        catch( err ) {
          if( ! (err instanceof OptimizationNoProgressError) )
            throw err;
        }

        expect(p).toBeAllCloseTo(P, {rtol: 1e-4, atol: 1e-4});
      });
  });
}
