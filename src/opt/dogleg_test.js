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
import {array, NDArray} from '../nd_array'
import {tabulate} from '../tabulate'

import {norm} from '../la/norm'

import {lsq_dogleg_gen,
        fit_dogleg_gen} from './dogleg'
import {num_grad} from './num_grad'
import {rosenbrock,
        rosenbrock_grad,
        rosenbrock_lsq,
        rosenbrock_lsq_jac} from './test_fn/rosenbrock'


describe('dogleg', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

//*DEBUG*/  const samples = [];

  forEachItemIn(
    function*(){                     const n = 16;
      function*       range() { for( let i=n+1; i-- > 0; ) yield Math.PI*(1-2*(i/n)); }
      for( const x of range() )
      for( const y of range() ) { yield [x,y];
      for( const z of range() ) { yield [x,y,z]; }}

 //*DEBUG*/     const avg = samples.reduce((x,y) => x+y) / samples.length,
 //*DEBUG*/           std = Math.hypot( ...samples.map( x => (x-avg) / Math.sqrt(samples.length) ) );
 //*DEBUG*/
 //*DEBUG*/     console.log('Dogleg')
 //*DEBUG*/     console.log('------')
 //*DEBUG*/     console.log('MIN:', samples.reduce((x,y) => Math.min(x,y)) );
 //*DEBUG*/     console.log('MAX:', samples.reduce((x,y) => Math.max(x,y)) );
 //*DEBUG*/     console.log('AVG:', avg );
 //*DEBUG*/     console.log('STD:', std );
    }()
  ).it('lsq_dogleg_gen works on rosenbrock_lsq', x0 => {

    let nCalls = 0
    const fJ = x => {
      expect(++nCalls).toBeLessThan(32);
      return [
        rosenbrock_lsq(x),
        rosenbrock_lsq_jac(x)
      ];
    };

    const N = x0.length,
          M = (N-1)*2;

    for( const [x, mse, mse_grad] of lsq_dogleg_gen(fJ, x0) )
    {
      expect(x       ).toEqual( jasmine.any(NDArray) )
      expect(mse     ).toEqual( jasmine.any(NDArray) )
      expect(mse_grad).toEqual( jasmine.any(NDArray) )

      expect(x       .ndim).toBe(1)
      expect(mse     .ndim).toBe(0)
      expect(mse_grad.ndim).toBe(1)

      expect(x       .shape).toEqual( Int32Array.of(N) )
      expect(mse_grad.shape).toEqual( Int32Array.of(N) )

      expect(mse     ).toBeAllCloseTo( rosenbrock(x) / M )
      expect(mse_grad).toBeAllCloseTo( rosenbrock_grad(x).mapElems(x => x/M) )

      const gNorm = norm(mse_grad);
      if(   gNorm <= 1e-8 ) {
        expect(x       ).toBeAllCloseTo(1)
        expect(mse     ).toBeAllCloseTo(0)
        expect(mse_grad).toBeAllCloseTo(0)
        break;
      }
    }

//*DEBUG*/    samples.push(nCalls);
  })


  forEachItemIn(
    function*(){
      for( let run=6; run-- > 0; ) {
        const                            N      =  Math.random()*4 + 1 | 0,
                      coeffs = tabulate([N], () => Math.random()*4 - 2);
        Object.freeze(coeffs)
        Object.freeze(coeffs.data.buffer)
        yield         coeffs;
      }
    }()
  ).it(`fit_dogleg_gen fits coefficients of polynomial.`, coeffs => {

    const [N] = coeffs.shape;

    const [f,g] = [
      p => x => p.reduce( (sum,p,i) => sum + p * x**i ),
      p => x => p.map( (p,i) => x**i )
    ].map( f => {
      expect(f).toEqual( jasmine.any(Function) );

      return p => {
        expect(p).toEqual( jasmine.any(NDArray) );
        expect(p.shape).toEqual(coeffs.shape);

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
              p = tabulate([N], 'float64', () => Math.random()*4 - 2);
        expect( g(p)(x) ).toBeAllCloseTo( h(p)(x) );
      }
    };

    const M = 256;

    const x = tabulate([M,1], 'float64', () => Math.random()*8 - 4),
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

    const computeRes = p => x.data.map( (x,i) => f(p)( array([x]) ) - y(i) ),
          computeErr = p => x.data.reduce( (sum,x,i) => sum + (f(p)( array([x]) ) - y(i))**2 / M, 0 ),
          computeGrad= num_grad(computeErr);

    let nIter = 0,
        mse,
        grad,
        param,
        res;

    for( [param, mse, grad, res] of fit_dogleg_gen(
      x,y, fg, /*p0=*/tabulate([N], () => Math.random()*4 - 2)
    ))
    {
      expect(++nIter).toBeLessThan(64);

      expect(mse ).toBeAllCloseTo( computeErr (param) );
      expect(grad).toBeAllCloseTo( computeGrad(param), {rtol: 1e-3} );
      expect(res ).toBeAllCloseTo( computeRes (param) );

      if( norm(grad) <= Math.sqrt(M)*1e-12 )
        break;
    }

    expect(param).toBeAllCloseTo(coeffs);
  })


  forEachItemIn(
    function*(){
      for( let run=6; run-- > 0; )
      { const                            N      =  Math.random()*4 + 1 | 0,
                      coeffs = tabulate([N], () => Math.random()*4 - 2);
        Object.freeze(coeffs)
        Object.freeze(coeffs.data.buffer)
        yield         coeffs;
      }
    }()
  ).it('fit_dogleg_gen fits coefficients of polynomial in root form.', coeffs => {

    const [N] = coeffs.shape;

    const [f,g] = [
      p => x => p.reduce( (prod,p) => prod*(p-x), 1 ),
      p => x => p.map(
        (_,i) => p.reduce( (prod,p,j) => i===j ? prod : prod*(p-x), 1 )
      )
    ].map( f => {
      expect(f).toEqual( jasmine.any(Function) );

      return p => {
        expect(p).toEqual( jasmine.any(NDArray) );
        expect(p.shape).toEqual(coeffs.shape);

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

    const M = 128;

    const x = tabulate([M,1], 'float64', () => Math.random()*8 - 4),
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

    for( [param, mse, grad, res] of fit_dogleg_gen(
      x,y, fg, /*p0=*/tabulate([N], () => Math.random()*4 - 2)
    ))
    {
      expect(++nIter).toBeLessThan(512);

      expect(mse ).toBeAllCloseTo( computeErr (param) );
      expect(grad).toBeAllCloseTo( computeGrad(param), {rtol: 1e-3} );
      expect(res ).toBeAllCloseTo( computeRes (param) );

      if( norm(grad) <= Math.sqrt(M)*1e-12 )
        break;
    }

    const par =  param.data.slice().sort( (x,y) => x-y ),
          coe = coeffs.data.slice().sort( (x,y) => x-y );

    expect(par).toBeAllCloseTo(coe);
  })


  forEachItemIn(
    function*(){
      for( let run=6; run-- > 0; )
      {
        const         coeffs = tabulate([2], () => Math.random()*4 - 2);
        Object.freeze(coeffs);
        Object.freeze(coeffs.data.buffer);
        yield         coeffs;
      }
    }()
  ).it('fit_dogleg_gen fits scaled exponential function.', coeffs => {

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

    for( [param, mse, grad, res] of fit_dogleg_gen(
      x,y, fg, /*p0=*/coeffs.mapElems( x => x + Math.random()*4 - 2 )
    ))
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
    expect(param).toBeAllCloseTo(coeffs);
  })
})
