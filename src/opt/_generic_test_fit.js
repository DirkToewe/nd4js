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
import {tabulate} from "../tabulate";

import {norm} from "../la/norm";

import {num_grad} from "./num_grad";
import {OptimizationNoProgressError} from './optimization_error';


export function generic_test_fit_gen( fit_gen )
{
  describe(`${fit_gen.name} [generic tests]`, () => {
    beforeEach( () => {
      jasmine.addMatchers(CUSTOM_MATCHERS)
    })
  
  
    forEachItemIn(
      function*(rng){
        for( let run=7; run-- > 0; ) {
          const                            N      =  rng.int(1,4),
                            x0 = tabulate([N], () => rng.uniform(-2,+2) ),
                        coeffs = tabulate([N], () => rng.uniform(-2,+2) );
          Object.freeze(coeffs)
          Object.freeze(coeffs.data.buffer)
  
          const fg = p => {
            expect(p).toEqual( jasmine.any(NDArray) );
            expect(p.shape).toEqual(coeffs.shape);
            p = p.data;
    
            return x => {
              expect(x).toEqual( jasmine.any(NDArray) );
              expect(x.shape).toEqual( Int32Array.of(1) );
              x = x.data[0];
    
              return [
                p.reduce( (sum,p,i) => sum + p * x**i ),
                p.map( (p,i) => x**i )
              ];
            };
          };
          const fg_c = fg(coeffs);
      
          const           M = rng.int(128,256),
            x = tabulate([M,1], 'float64', () => rng.uniform(-2,+2) ),
            y = new NDArray(
              Int32Array.of(M),
              x.data.map(
                x => fg_c( array([x]) )[0]
              )
            );    
          Object.freeze(x);
          Object.freeze(y);
          Object.freeze(x.data.buffer);
          Object.freeze(y.data.buffer);

          yield [coeffs, fg,x,y, x0];
        }
      }
    ).it(`fits coefficients of polynomial.`, ([coeffs, fg,x,y, x0]) => {
      const [N] = coeffs.shape,
            [M] =      y.shape;

      const f = p => { const fg_p = fg(p); return x => fg_p(x)[0]; },
            g = p => { const fg_p = fg(p); return x => fg_p(x)[1]; };

      // check gradients
      ;{
        const h = p => x => num_grad( p => f(p)(x) )(p);
  
        for( let repeat=128; repeat-- > 0; )
        {
          const x = array([ Math.random()*4 - 2 ]),
                p = tabulate([N], 'float64', () => Math.random()*4 - 2);
          expect( g(p)(x) ).toBeAllCloseTo( h(p)(x) );
        }
      };
  
      const computeRes = p => x.data.map(        (x,i) =>        f(p)( array([x]) ) - y(i) ),
            computeErr = p => x.data.reduce( (sum,x,i) => sum + (f(p)( array([x]) ) - y(i))**2 / M, 0 ),
            computeGrad= num_grad(computeErr);
  
      let nIter = 0,
          mse,
          grad,
          param,
          res;
  
      try
      {
        for( [param, mse, grad, res] of fit_gen(
          x,y, fg, x0
        ))
        {
          expect(++nIter).toBeLessThan(64);
    
          expect(mse ).toBeAllCloseTo( computeErr (param) );
          expect(grad).toBeAllCloseTo( computeGrad(param), {rtol: 1e-3} );
          expect(res ).toBeAllCloseTo( computeRes (param) );
    
          if( norm(grad) <= Math.sqrt(M)*1e-12 )
            break;
        }
      }
      catch(    onpe ) {
        if( ! ( onpe instanceof OptimizationNoProgressError ) )
          throw onpe;
      }
  
      expect(param).toBeAllCloseTo(coeffs);
    })
  

    for( const [N_RUNS,N] of [
      [8,1],
      [4,2],
      [1,3],
      [1,4]
    ])
      forEachItemIn(
        function*(rng){
          for( let run=0; run++ < N_RUNS; )
          {
            const         coeffs = tabulate([N], 'float64', () => rng.uniform(-2,+2) ),
                     x0 = coeffs.mapElems('float64', x => x + rng.normal() );
            Object.freeze(coeffs);
            Object.freeze(coeffs.data.buffer);

            const fg = p => {
              expect(p).toEqual( jasmine.any(NDArray) );
              expect(p.shape).toEqual(coeffs.shape);
              p = p.data;
    
              return x => {
                expect(x).toEqual( jasmine.any(NDArray) );
                expect(x.shape).toEqual( Int32Array.of(1) );
                x = x.data;
    
                return [
                  p.reduce( (prod,p) => prod*(p-x), 1 ),
                  p.map(
                    (_,i) => p.reduce( (prod,p,j) => i===j ? prod : prod*(p-x), 1 )
                  )
                ];
              };
            };

            const           M = rng.int(128,256),
              x = tabulate([M,1], 'float64', () => rng.uniform(-2,+2) ),
              y = new NDArray(
                Int32Array.of(M),
                x.data.map(
                  x => fg(coeffs)(array([x]))[0]
                )
              );
            Object.freeze(x);
            Object.freeze(y);
            Object.freeze(x.data.buffer);
            Object.freeze(y.data.buffer);

            yield [coeffs, fg,x,y, x0];
          }
        }
      ).it(`fits coefficients of polynomial in root form with ${N} roots.`, ([coeffs, fg,x,y, x0]) => {
        const [N] = coeffs.shape,
              [M] =      y.shape;

        const f = p => { const fg_p = fg(p); return x => fg_p(x)[0]; },
              g = p => { const fg_p = fg(p); return x => fg_p(x)[1]; };

        // check gradient
        ;{      
          const h = p => x => num_grad( p => f(p)(x) )(p);

          for( let repeat=16; repeat-- > 0; )
          {
            const x = array([ Math.random()*4 - 2 ]),
                  p = tabulate([N], 'float64', () => Math.random()*4 - 2);
            expect( g(p)(x) ).toBeAllCloseTo( h(p)(x) );
          }
        };

        const computeRes = p => x.data.map(        (x,i) =>        f(p)( array([x]) ) - y(i) ),
              computeErr = p => x.data.reduce( (sum,x,i) => sum + (f(p)( array([x]) ) - y(i))**2 / M, 0 ),
              computeGrad= num_grad(computeErr);

        let nIter = 0,
            mse,
            grad,
            param,
            res;

        try
        {
          for( [param, mse, grad, res] of fit_gen(x,y, fg, x0) )
          {
            expect(++nIter).toBeLessThan(512);

            expect(mse ).toBeAllCloseTo( computeErr (param) );
            expect(grad).toBeAllCloseTo( computeGrad(param), {rtol: 1e-3} );
            expect(res ).toBeAllCloseTo( computeRes (param) );

            if( norm(grad) <= Math.sqrt(M)*1e-12 )
              break;
          } 
        }
        catch(    onpe ) {
          if( ! ( onpe instanceof OptimizationNoProgressError ) )
            throw onpe;
        }

        const par =  param.data.slice().sort( (x,y) => x-y ),
              coe = coeffs.data.slice().sort( (x,y) => x-y );

        expect(par).toBeAllCloseTo(coe);
      });


    forEachItemIn(
      function*(rng){
        for( let run=4; run-- > 0; )
        {
          const         coeffs = tabulate([2], 'float64', () => rng.uniform(-2,+2) ),
                   x0 = coeffs.mapElems('float64', x => x + rng.normal() );
          Object.freeze(coeffs);
          Object.freeze(coeffs.data.buffer);
          yield        [coeffs,x0];
        }
      }
    ).it('fits scaled exponential function.', ([coeffs,x0]) => {
  
      expect( coeffs.shape ).toEqual( Int32Array.of(2) );
  
      const [f,g] = [
        ([a,b]) => x => a * Math.exp(b*x),
        ([a,b]) => x => [Math.exp(b*x), a*x * Math.exp(b*x)]
      ].map( f => {
        expect(f).toEqual( jasmine.any(Function) );
  
        return p => {
          expect(p).toEqual( jasmine.any(NDArray) );
          expect(p.shape).toEqual( Int32Array.of(2) );
  
          const fp = f(p.data);
  
          return x => {
            expect(x).toEqual( jasmine.any(NDArray) );
            expect(x.shape).toEqual( Int32Array.of(1) );
  
            return fp(x.data[0]);
          };
        };
      });
  
      const fg = p => {
        const fp = f(p),
              gp = g(p);
        return x => [ fp(x), gp(x) ];
      };
  
      ;{      
        const h = p => x => num_grad( p => f(p)(x) )(p);
  
        for( let repeat=128; repeat-- > 0; )
        {
          const x = array([ Math.random()*4 - 2 ]),
                p = tabulate([2], 'float64', () => Math.random()*4 - 2);
          expect( g(p)(x) ).toBeAllCloseTo( h(p)(x) );
        }
      };
  
      const M = 128;
  
      const x = tabulate([M,1], 'float64', () => Math.random()*6 - 3),
            y = new NDArray(
              Int32Array.of(M),
              x.data.map(
                x => f(coeffs)( array([x]) )
              )
            );
      Object.freeze(x);
      Object.freeze(y);
      Object.freeze(x.data.buffer);
      Object.freeze(y.data.buffer);
  
      const computeRes = p => x.data.map(        (x,i) =>        f(p)( array([x]) ) - y(i) ),
            computeErr = p => x.data.reduce( (sum,x,i) => sum + (f(p)( array([x]) ) - y(i))**2 / M, 0 ),
            computeGrad= num_grad(computeErr);
  
      let nIter = 0,
          mse,
          grad,
          param,
          res;
  
      try
      {
        for( [param, mse, grad, res] of fit_gen(x,y, fg, x0) )
        {
          expect(++nIter).toBeLessThan(1024);
    
          if( ! isFinite(mse) )
            console.log({coeffs, param, mse, grad})
    
          expect(mse ).toBeAllCloseTo( computeErr (param) );
          expect(grad).toBeAllCloseTo( computeGrad(param), {rtol: 1e-3, atol:1e-5} );
          expect(res ).toBeAllCloseTo( computeRes (param) );
    
          if( norm(grad) <= Math.sqrt(M)*1e-12 )
            break;
        }
      }
      catch(    onpe ) {
        if( ! ( onpe instanceof OptimizationNoProgressError ) )
          throw onpe;
      }

      expect(param).toBeAllCloseTo(coeffs);
    })
  })
}
