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


import {cartesian_prod, linspace} from '../../iter'
import {forEachItemIn, CUSTOM_MATCHERS} from '../../jasmine_utils'
import {asarray} from '../../nd_array'
import {stack} from '../../stack'
import {tabulate} from '../../tabulate'
import {_rand_int} from '../../_test_data_generators'

import {matmul} from '../../la/matmul'

import {num_grad} from '../num_grad'


export function generic_test_test_fn( test_fn, x_range )
{
  if( x_range.length !== test_fn.nIn )
    throw new Error('Assertion failed.');
  for( const lo_hi of x_range )
  {
    if( lo_hi.length !== 2 )
      throw new Error('Assertion failed.');

    const [lo,hi] = lo_hi;
    if( ! (lo < hi))
      throw new Error('Assertion failed.');
  }


  describe(`test_fn.${test_fn.name} [generic tests]`, () => {
    beforeEach( () => {
      jasmine.addMatchers(CUSTOM_MATCHERS);
    });


    const test_fn_num_grad = num_grad(test_fn),
          test_fn_num_hess = x => {
            x = asarray(x);

            const N = x.shape[x.ndim-1];

            return stack(-1, Array.from({length: N}, (_,i) =>
              num_grad( x => test_fn.grad(x)(i) )(x)
            ));
          };


    for( const x of test_fn.roots )
      it(`has a root at x=[${x}]`, () => {
        const f = test_fn(x);

        expect(f.shape).toEqual( Int32Array.of() );
        expect(f).toBeAllCloseTo(0)
      })


    for( const x of test_fn.roots )
      it(`lsq has a root at x=[${x}]`, () => {
        const f = test_fn.lsq(x);

        expect(f.shape).toEqual( Int32Array.of(test_fn.nOut) );
        expect(f).toBeAllCloseTo(0)
      })


    it('minima_global has no duplicates', () => {
      const minima = test_fn.minima_global;

      for( let i=minima.length; i-- > 0; )
      for( let j=i            ; j-- > 0; )
        expect(minima[i]).not.toEqual(minima[j]);
    })


    it('minima starts with minima_global', () => {
      const eq = (a,b) => {
        expect(a.length).toBe(test_fn.nIn);
        expect(b.length).toBe(test_fn.nIn);

        for( let i=test_fn.nIn; i-- > 0; )
          if( a[i] !== b[i] )
            return false;

        return true;
      };

      const minima = test_fn.minima_global;
      let i=0;

      for( const x of test_fn.minima )
      {
        expect(x).toEqual(minima[i]);

        if( ++i >= minima.length )
          break;
      }
    })


    if( test_fn.minima_global.length > 1 )
      it('minima_global all produce same function value', () => {
        const [first, ...rest] = test_fn.minima_global.map(test_fn);
        expect(first).toBeAllCloseTo([...rest]);
      })


    forEachItemIn(
      function(){
        const N = Math.round( 2 ** (17/test_fn.nIn) );

        return cartesian_prod(
          ...x_range.map( r => linspace(...r,N) )
        )
      }()
    ).it('minima_global checks out given generated examples', x => {
      const f0 = test_fn(test_fn.minima_global[0]),
            f  = test_fn(x);
      expect(f0).toBeLessThanOrEqual(f);
    })


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 3*1337; )
          yield x_range.map( ([lo,hi]) => {
            const s = Math.random();
            return lo*(1-s) + s*hi;
          })
      }()
    ).it('minima_global checks out given random examples', x => {
      const f0 = test_fn(test_fn.minima_global[0]),
            f  = test_fn(x);
      expect(f0).toBeLessThanOrEqual(f);
    })


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 7*1337; )
        {
          const shape = Int32Array.from({ length: _rand_int(1,4) }, () => _rand_int(1,8) );
          shape[shape.length-1] = test_fn.nIn;

          const X = tabulate(shape, (...idx) => {
            const  i = idx.pop(),
              [lo,hi]= x_range[i],
                   s = Math.random();
            return lo*(1-s) + s*hi;
          })
          Object.freeze(X)
          Object.freeze(X.data.buffer)

          yield X;
        }
      }()
    ).it('minima_global checks out given random broadcast examples', x => {
      const f0 = test_fn(test_fn.minima_global[0]),
            f  = test_fn(x);
      expect(f0).toBeAllLessOrClose(f);
    })


    forEachItemIn(
      function*(){
        const {minima} = test_fn;
        let len = minima.length || 2**16;
  
        for( const x of minima )
        {
          if( --len < 0 )
            break;
          yield x;
        }
      }()
    ).it(`minima have zero grad`, x => {
      const g = test_fn.grad(x);

      expect(g.shape).toEqual( Int32Array.of(test_fn.nIn) );
      expect(g).toBeAllCloseTo(0, {atol: 7e-5})
    });


    for( const x of test_fn.minima_global )
      it(`has zero grad at global min. x=[${x}]`, () => {
        const g = test_fn.grad(x);

        expect(g.shape).toEqual(Int32Array.of(test_fn.nIn));
        expect(g).toBeAllCloseTo(0, {atol: 7e-5})
      })


    //
    // GRAD
    //
    forEachItemIn(
      function(){
        const N = Math.round( 2 ** (17/test_fn.nIn) );

        return cartesian_prod(
          ...x_range.map( r => linspace(...r,N) )
        )
      }()
    ).it('grad works given generated examples', x => {
      const g = test_fn.grad(x),
            G = test_fn_num_grad(x)
  
      expect(g).toBeAllCloseTo(G, {rtol:1e-3, atol:1e-5})
    })


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 13*1337; )
          yield x_range.map( ([lo,hi]) => {
            const s = Math.random();
            return lo*(1-s) + s*hi;
          })
      }()
    ).it('grad works given random examples', x => {
      const g = test_fn.grad(x),
            G = test_fn_num_grad(x)
  
      expect(g).toBeAllCloseTo(G, {rtol:1e-3, atol:1e-5})
    })


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 13*1337; )
        {
          const shape = Int32Array.from({ length: _rand_int(1,4) }, () => _rand_int(1,8) );
          shape[shape.length-1] = test_fn.nIn;

          const X = tabulate(shape, (...idx) => {
            const  i = idx.pop(),
              [lo,hi]= x_range[i],
                   s = Math.random();
            return lo*(1-s) + s*hi;
          })
          Object.freeze(X)
          Object.freeze(X.data.buffer)

          yield X;
        }
      }()
    ).it('grad works given random broadcast examples', x => {
      const g = test_fn.grad(x),
            G = stack(
              [...x.reshape(-1, test_fn.nIn) ].map(test_fn_num_grad)
            ).reshape(...x.shape)
  
      expect(g).toBeAllCloseTo(G, {rtol:1e-3, atol:1e-5})
    })


    //
    // HESS
    //
    forEachItemIn(
      function(){
        const N = Math.round( 2 ** (16/test_fn.nIn) );

        return cartesian_prod(
          ...x_range.map( r => linspace(...r,N) )
        )
      }()
    ).it('hess works given generated examples', x => {
      const h = test_fn.hess(x),
            H = test_fn_num_hess(x)
  
      expect(h).toBeAllCloseTo(H, {rtol:1e-3, atol:1e-5})
    })


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 7*1337; )
          yield x_range.map( ([lo,hi]) => {
            const s = Math.random();
            return lo*(1-s) + s*hi;
          })
      }()
    ).it('hess works given random examples', x => {
      const h = test_fn.hess(x),
            H = test_fn_num_hess(x)
  
      expect(h).toBeAllCloseTo(H, {rtol:1e-3, atol:1e-5})
    })


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 3*1337; )
        {
          const shape = Int32Array.from({ length: _rand_int(1,4) }, () => _rand_int(1,8) );
          shape[shape.length-1] = test_fn.nIn;

          const X = tabulate(shape, (...idx) => {
            const  i = idx.pop(),
              [lo,hi]= x_range[i],
                   s = Math.random();
            return lo*(1-s) + s*hi;
          })
          Object.freeze(X)
          Object.freeze(X.data.buffer)

          yield X;
        }
      }()
    ).it('hess works given random broadcast examples', x => {
      const h = test_fn.hess(x),
            H = stack(
              [...x.reshape(-1, test_fn.nIn) ].map(test_fn_num_hess)
            ).reshape(...x.shape.slice(0,-1), test_fn.nIn, test_fn.nIn)
  
      expect(h).toBeAllCloseTo(H, {rtol:1e-3, atol:1e-5})
    })


    //
    // LSQ
    //
    forEachItemIn(
      function(){
        const N = Math.round( 2 ** (18/test_fn.nIn) );

        return cartesian_prod(
          ...x_range.map( r => linspace(...r,N) )
        )
      }()
    ).it(`lsq is consistent with ${test_fn.name} given generated examples`, x => {
      const Y = test_fn(x),
            y = test_fn.lsq(x)
                  .mapElems(x.dtype, x => x*x)
                  .reduceElems(-1, x.dtype, (x,y) => x+y);
  
      expect(y).toBeAllCloseTo(Y);
    })


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 7*1337; )
          yield x_range.map( ([lo,hi]) => {
            const s = Math.random();
            return lo*(1-s) + s*hi;
          })
      }()
    ).it(`lsq is consistent with ${test_fn.name} given random examples`, x => {
      const Y = test_fn(x),
            y = test_fn.lsq(x)
                  .mapElems(x.dtype, x => x*x)
                  .reduceElems(-1, x.dtype, (x,y) => x+y);
  
      expect(y).toBeAllCloseTo(Y);
    })


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 3*1337; )
        {
          const shape = Int32Array.from({ length: _rand_int(1,4) }, () => _rand_int(1,8) );
          shape[shape.length-1] = test_fn.nIn;

          const X = tabulate(shape, (...idx) => {
            const  i = idx.pop(),
              [lo,hi]= x_range[i],
                  s = Math.random();
            return lo*(1-s) + s*hi;
          })
          Object.freeze(X)
          Object.freeze(X.data.buffer)

          yield X;
        }
      }()
    ).it(`lsq is consistent with ${test_fn.name} given random broadcast examples`, x => {
      const Y = test_fn(x),
            y = test_fn.lsq(x)
                  .mapElems(x.dtype, x => x*x)
                  .reduceElems(-1, x.dtype, (x,y) => x+y);
  
      expect(y).toBeAllCloseTo(Y);
    })


    //
    // LSQ_JAC
    //
    forEachItemIn(
      function(){
        const N = Math.round( 2 ** (17/test_fn.nIn) );

        return cartesian_prod(
          ...x_range.map( r => linspace(...r,N) )
        )
      }()
    ).it('lsq_jac is consistent with grad given generated examples', x => {
      let G = test_fn.grad(x),
          J = test_fn.lsq_jac(x),
          f = test_fn.lsq (x);
      G = G.reshape(...G.shape,1);
      f = f.reshape(...f.shape,1);
      const  g = matmul(J.T, f, [[2]]);
      expect(g).toBeAllCloseTo(G);
    })


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 73*1337; )
          yield x_range.map( ([lo,hi]) => {
            const s = Math.random();
            return lo*(1-s) + s*hi;
          })
      }()
    ).it('lsq_jac is consistent with grad given random examples', x => {
      let G = test_fn.grad(x),
          J = test_fn.lsq_jac(x),
          f = test_fn.lsq (x);
      G = G.reshape(...G.shape,1);
      f = f.reshape(...f.shape,1);
      const  g = matmul(J.T, f, [[2]]);
      expect(g).toBeAllCloseTo(G);
    })


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 13*1337; )
        {
          const shape = Int32Array.from({ length: _rand_int(1,4) }, () => _rand_int(1,8) );
          shape[shape.length-1] = test_fn.nIn;

          const X = tabulate(shape, (...idx) => {
            const  i = idx.pop(),
              [lo,hi]= x_range[i],
                  s = Math.random();
            return lo*(1-s) + s*hi;
          })
          Object.freeze(X)
          Object.freeze(X.data.buffer)

          yield X;
        }
      }()
    ).it('lsq_jac is consistent with grad given random broadcast examples', x => {
      let G = test_fn.grad(x),
          J = test_fn.lsq_jac(x),
          f = test_fn.lsq (x);
      G = G.reshape(...G.shape,1);
      f = f.reshape(...f.shape,1);
      const  g = matmul(J.T, f, [[2]]);
      expect(g).toBeAllCloseTo(G);
    })
  });
}
