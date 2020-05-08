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

import {cartesian_prod} from '../../iter'
import {asarray, NDArray} from '../../nd_array'
import {root1d_bisect} from '../root1d_bisect'


// REFERENCES
// ----------
// .. [1] https://en.wikipedia.org/wiki/Rastrigin_function


const NONGLOBAL_MINIMA = function(){
  const g = x => x + Math.PI*10 * Math.sin(Math.PI*2*x);

  const minima = [];

  for( let x=1; x < 32; x++ )
  {
    if( !( g(x-0.25) < 0 ) ) throw new Error('Assertion failed.');
    if( !( g(x+0.25) > 0 ) ) throw new Error('Assertion failed.');

    minima.push(
      root1d_bisect(g, x-0.25, x+0.25)
    );
  }

  Object.freeze(minima);
  return minima;
}();


export class Rastrigin extends Function
{
  constructor( n )
  {
    super();
    if( 0 !== n%1 ) throw new Error(`new Rastrigin(n): n not a valid integer.`);
    if( ! (0 < n) ) throw new Error(`new Rastrigin(n): n must be positive.`);

    const rastrigin_fn = x =>
    {
      x = asarray(x)

      const N = rastrigin_fn.nIn;

      if(         x.ndim     <  1 ) throw new Error(`${rastrigin_fn.name}(x): x.ndim must be at least 1.`);
      if( x.shape[x.ndim-1] !== N ) throw new Error(`${rastrigin_fn.name}(x): x.shape[-1] must be ${N}.`);

      const F_shape =      x.shape.slice(0,-1),
            F       = new (x.dtype==='float32' ? Float32Array
                                               : Float64Array)(x.data.length/N);
      x = x.data;

      for( let x_off=x.length,
               F_off=F.length;
               F_off-- > 0; )
      {        x_off -= N;
        for( let i=N; i-- > 0; )
        {
          const xi = x[x_off + i];
          F[F_off] += xi*xi + 10*( 1 - Math.cos(Math.PI*2*xi) );
        }
      }

      return new NDArray(F_shape, F);
    };

    rastrigin_fn.nIn = n;

    rastrigin_fn.roots =
    rastrigin_fn.minima_global = [ Array.from({length: n}, () => 0) ];
    rastrigin_fn.minima =
    {
      [Symbol.iterator]: function*()
      {
        yield Array.from({length: n}, () => 0);

        let inner,
            outer = [0];

        // grow the minima in a hypercube from the inside out
        for( let layer=0; layer < NONGLOBAL_MINIMA.length; layer++ )
        {
          const min = NONGLOBAL_MINIMA[layer];

          inner =           outer;
          outer = [-min, ...inner, +min];

          const ranges = Array.from({length: n}, () => outer);

          for( let i=0; i < n; i++ )
          {
            ranges[i] = [-min,+min]; yield* cartesian_prod(...ranges);
            ranges[i] = inner;
          }
        }
      }
    }

    Object.defineProperty(rastrigin_fn, 'name', {value: `rastrigin${n}d`, writable: false})
    Object.setPrototypeOf(rastrigin_fn, Rastrigin.prototype);
    return Object.freeze( rastrigin_fn );
  }


  get nOut() {
    return this.nIn;
  }


  grad( x )
  {
    x = asarray(x);

    const N = this.nIn;

    if(         x.ndim     <  1 ) throw new Error(`${this.name}.grad(x): x.ndim must be at least 1.`);
    if( x.shape[x.ndim-1] !== N ) throw new Error(`${this.name}.grad(x): x.shape[-1] must be ${N}.`);

    const G_shape = x.shape,
          G       = new (x.dtype==='float32' ? Float32Array
                                             : Float64Array)(x.data.length);
    x = x.data;

    for( let off=x.length;
             off >= 0; )
    {        off -= N;
      for( let i=N; i-- > 0; )
      {
        const xi = x[off+i];
        G[off+i] = 2*xi  +  Math.PI*20 * Math.sin(Math.PI*2*xi);
      }
    }

    return new NDArray(G_shape, G);
  }


  hess( x )
  {
    x = asarray(x);

    const N = this.nIn;

    if(         x.ndim     <  1 ) throw new Error(`${this.name}.hess(x): x.ndim must be at least 1.`);
    if( x.shape[x.ndim-1] !== N ) throw new Error(`${this.name}.hess(x): x.shape[-1] must be ${N}.`);

    const H_shape = Int32Array.of(...x.shape, N),
          H = new (x.dtype==='float32' ? Float32Array
                                       : Float64Array)(x.data.length*N);
    x = x.data;

    for( let H_off=H.length,
             x_off=x.length;
            (x_off -= N )>= 0; )
    {        H_off -= N*N;
      for( let i=N; i-- > 0; )
      {
        const xi = x[x_off + i];
        H[H_off + N*i+i] += 2  +  Math.PI*Math.PI*40 * Math.cos(Math.PI*2*xi);
      }
    }

    return new NDArray(H_shape, H);
  }


  lsq( x )
  {
    x = asarray(x);

    const N = this.nIn;

    if(         x.ndim     <  1 ) throw new Error(`${this.name}.lsq(x): x.ndim must be at least 1.`);
    if( x.shape[x.ndim-1] !== N ) throw new Error(`${this.name}.lsq(x): x.shape[-1] must be ${N}.`);

    const F_shape = x.shape,
          F       = new (x.dtype==='float32' ? Float32Array
                                             : Float64Array)(x.data.length);
    x = x.data;

    for( let x_off=x.length,
             F_off=F.length;
            (F_off -= N) >= 0; )
    {        x_off -= N;
      for( let i=N; i-- > 0; )
      {
        const xi = x[x_off + i];
        F[F_off + i] = Math.sqrt( xi*xi + 10*( 1 - Math.cos(Math.PI*2*xi) ) );
      }
    }

    return new NDArray(F_shape, F);
  }


  lsq_jac( x )
  {
    x = asarray(x);

    const N = this.nIn;

    if(         x.ndim     <  1 ) throw new Error(`${this.name}.lsq_jac(x): x.ndim must be at least 1.`);
    if( x.shape[x.ndim-1] !== N ) throw new Error(`${this.name}.lsq_jac(x): x.shape[-1] must be ${N}.`);
  
    const J_shape = Int32Array.of(...x.shape, N),
          J       = new (x.dtype==='float32' ? Float32Array
                                             : Float64Array)(x.data.length*N);
    x = x.data;
  
    for( let x_off=x.length,
             J_off=J.length;
            (J_off -= N*N) >= 0; )
    {        x_off -= N;
      for( let i=N; i-- > 0; )
      {
        const xi = x[x_off + i];
        const    denom  = Math.sqrt( xi*xi + 10*( 1 - Math.cos(Math.PI*2*xi) ) );
        if(0 !== denom)
          J[J_off + N*i+i] = ( xi + Math.PI*10*Math.sin(Math.PI*2*xi) ) / denom;
      }
    }
  
    return new NDArray(J_shape, J);
  }
}
