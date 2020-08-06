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
import {_rand_int} from "../_test_data_generators";

import {norm} from "../la/norm";


import {num_grad} from "./num_grad";
import { OptimizationNoProgressError } from './optimization_error';
import { rand_normal } from '../rand_normal';


/* Creates a quadratic polynomial with NX-dimensional input `x`.
 * To make the polynomial non-linear in its parameters `p`,
 * each polynomial coeffcient is itself a cubic, monotonous polynomial
 * of some entry of `p`.
 */
function rand_quadratic_poly( NX )
{
  const NP = 1 + NX + NX*(NX+1)/2;

  const c = new Float64Array(3*NP);
  for( let i=c.length; (i-=3) >= 0; )
  {
    c[i+0] = Math.random()*4 - 2;
    c[i+1] = Math.random()*2 + 0.1;
    c[i+2] = rand_normal()*1e-5;
    if( Math.random() < 0.5 ) {
      c[i+0] *= -1;
      c[i+1] *= -1;
      c[i+2] *= -1;
    }
  }

  const fgh = p =>
  {
    expect(p.shape).toEqual( Int32Array.of(NP) );
    p = p.data;

    return x => {
      expect(x.shape).toEqual( Int32Array.of(NX) );
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
  Object.assign(fgh,{NX,NP});
  Object.defineProperty(fgh, 'name', {value: `cubic_poly(${NX})`, writable: false});
  Object.freeze(fgh);

  return fgh;
}


export function generic_test_odr_gen( odr_gen )
{
  describe(`${odr_gen.name} [generic tests]`, () => {
    beforeEach( () => {
      jasmine.addMatchers(CUSTOM_MATCHERS)
    })

    for( const NX of [1,2,3,4] )
      forEachItemIn(
        function*(){
          for( let run=0; run++ < 128; )
          {
            const fgh = new rand_quadratic_poly(NX);

            for( let step=0; step++ < 128/NX; )
            {
              const x = tabulate([fgh.NX], () => Math.random()*4 - 2),
                    p = tabulate([fgh.NP], () => Math.random()*4 - 2);
              yield [fgh,p,x];
            }
          }
        }()
      ).it(`quadratic_poly(${NX}) gradients are correct`, ([fgh, p, x]) => {
        const  F  = (p,x) => fgh(p)(x)[0],
          [f,g,h] = fgh(p)(x);

        expect(f).toBe( F(p,x) );
        expect(g).toBeAllCloseTo( num_grad( p => F(p,x) )(p) );
        expect(h).toBeAllCloseTo( num_grad( x => F(p,x) )(x) );
      });

    for( const NX of [1,2,3,4] )
      forEachItemIn(
        function*(){
          for( let run=0; run++ < 48 / NX**1.5; )
          {
            const fgh = new rand_quadratic_poly(NX);

            const {NP} = fgh,
                   MX  = _rand_int(800*NP, 1600*NP);

            const p  = tabulate([NP], () => Math.random()*4 - 2),
                  p0 = tabulate([NP], () => Math.random()*4 - 2);

            const noise = 1e-6;

            const samples_x = tabulate([MX,NX], () => Math.random()*8 - 4),
                  samples_y = array('float64', [...samples_x].map( x => fgh(p)(x)[0] + rand_normal()*noise ) );

            for( let i=samples_x.data.length; i-- > 0; )
              samples_x.data[i] += rand_normal()*noise;

            const dx0 = samples_x.mapElems('float64', Math.random < 0.5 ? () => 0
                                                                        : () => rand_normal()*1e-3);
            yield [samples_x,samples_y, fgh,p, p0,dx0];
          }
        }()
      ).it(`fits random quadratic_poly(${NX}) correctly`, ([samples_x, samples_y, fgh, P, P0, dx0]) => {
        let     p,dx, mse,dmse_dp,dmse_dx, dy;
        try {
          for( [p,dx, mse,dmse_dp,dmse_dx, dy] of odr_gen(samples_x, samples_y, fgh, P0, dx0) )
          {
            if(   norm(dmse_dp) < 1e-8
               && norm(dmse_dx) < 1e-8 )
              break;
          }
        }
        catch( err ) {
          if( ! (err instanceof OptimizationNoProgressError) )
            throw err;
        }

        expect(p).toBeAllCloseTo(P, {atol: 1e-4});
      });
  });
}
