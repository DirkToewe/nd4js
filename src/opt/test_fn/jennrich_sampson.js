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

import {asarray, NDArray} from '../../nd_array'


// REFERENCES
// ----------
// .. [1] "Testing Unconstrained Optimization Software"
//         Jorge J. Moré, Burton S. Garbow, Kenneth E. Hillstrom
//         ACM Transactions on Mathematical Software, Vol 7, No. 1, March 1982, pp. 17-41
// .. [2] "Application of stepwise regression to nonlinear estimation"
//         R. I. Jennrich, P. F. Sampson
//         Technometrics 10 (1968), pp. 64-72


export class JennrichSampson extends Function
{
  constructor( m )
  {
    super();
    if( 0 !== m%1 ) throw new Error(`new JennrichSampson(m): m not a valid integer.`);
    if( ! (1 < m) ) throw new Error(`new JennrichSampson(m): m must be greater than 1.`);

    const jennrich_sampson_fn = x =>
    {
      x = asarray(x)

      if(         x.ndim     <  1 ) throw new Error(`${jennrich_sampson_fn.name}(x): x.ndim must be at least 1.`);
      if( x.shape[x.ndim-1] !== 2 ) throw new Error(`${jennrich_sampson_fn.name}(x): x.shape[-1] must be 2.`);

      const F_shape =      x.shape.slice(0,-1),
            F       = new (x.dtype==='float32' ? Float32Array
                                               : Float64Array)(x.data.length/2);
      x = x.data;

      const M = jennrich_sampson_fn.nOut;

      for( let x_off=x.length,
               F_off=F.length;
               F_off-- > 0; )
      {        x_off -= 2;
        const x1 = x[x_off+0],
              x2 = x[x_off+1];

        let sum = 0;

        for( let  i=M; i-- > 0; ) {
          const  fi = 4 + 2*i - Math.exp( (i+1)*x1 ) - Math.exp( (i+1)*x2 );
          sum += fi*fi;
        }

        F[F_off] = sum;
      }

      return new NDArray(F_shape, F);
    };

    if( m === 10 ) {
      jennrich_sampson_fn.minima_global=
      jennrich_sampson_fn.minima       = [
        [0.25782521321500883,
         0.25782521321500883]
      ];
    }
    jennrich_sampson_fn.roots = [];
    jennrich_sampson_fn.nOut = m;

    Object.defineProperty(jennrich_sampson_fn, 'name', {value: `jennrich_sampson${m}d`, writable: false})
    Object.setPrototypeOf(jennrich_sampson_fn, JennrichSampson.prototype);
    return Object.freeze( jennrich_sampson_fn );
  }


  get nIn() {
    return 2;
  }


  grad( x )
  {
    x = asarray(x);

    if(         x.ndim     <  1 ) throw new Error(`${this.name}.grad(x): x.ndim must be at least 1.`);
    if( x.shape[x.ndim-1] !== 2 ) throw new Error(`${this.name}.grad(x): x.shape[-1] must be 2.`);

    const G_shape =      x.shape,
          G       = new (x.dtype==='float32' ? Float32Array
                                             : Float64Array)(x.data.length);
    x = x.data;

    const M = this.nOut;

    for( let off=x.length;
             off >= 0; )
    {        off -= 2;
      const x1 = x[off+0],
            x2 = x[off+1];

      let sum = 0;

      for( let  i=M; i-- > 0; ) {
        const  fi = 4 + 2*i - Math.exp( (i+1)*x1 ) - Math.exp( (i+1)*x2 );
        G[off + 0] -= 2*fi*(i+1) * Math.exp( (i+1)*x1 );
        G[off + 1] -= 2*fi*(i+1) * Math.exp( (i+1)*x2 );
      }
    }

    return new NDArray(G_shape, G);
  }


  hess( x )
  {
    x = asarray(x);

    if(         x.ndim     <  1 ) throw new Error(`${this.name}.hess(x): x.ndim must be at least 1.`);
    if( x.shape[x.ndim-1] !== 2 ) throw new Error(`${this.name}.hess(x): x.shape[-1] must be 2.`);

    const H_shape = Int32Array.of(...x.shape, 2),
          H = new (x.dtype==='float32' ? Float32Array
                                      : Float64Array)(x.data.length*2);
    x = x.data;

    const M = this.nOut;

    for( let H_off=H.length,
             x_off=x.length;
            (x_off -= 2) >= 0; )
    {        H_off -= 2*2;
      const x1 = x[x_off+0],
            x2 = x[x_off+1];

      let sum = 0;

      for( let i=M; i-- > 0; ) {
        const  fi = 4 + 2*i - Math.exp( (i+1)*x1 ) - Math.exp( (i+1)*x2 );

        H[H_off + 0] += 2*(i+1)*(i+1) * ( Math.exp( (i+1)* x1*2   ) - fi * Math.exp( (i+1)*x1 ) );
        H[H_off + 3] += 2*(i+1)*(i+1) * ( Math.exp( (i+1)* x2*2   ) - fi * Math.exp( (i+1)*x2 ) );
        H[H_off + 1] += 2*(i+1)*(i+1) *   Math.exp( (i+1)*(x1+x2) );
      }
      H[H_off + 2] = H[H_off + 1];
    }

    return new NDArray(H_shape, H);
  }


  lsq( x )
  {
    x = asarray(x);

    if(         x.ndim     <  1 ) throw new Error(`${this.name}.lsq(x): x.ndim must be at least 1.`);
    if( x.shape[x.ndim-1] !== 2 ) throw new Error(`${this.name}.lsq(x): x.shape[-1] must be 2.`);

    const M = this.nOut;

    const F_shape = x.shape.slice(),
          F       = new (x.dtype==='float32' ? Float32Array
                                             : Float64Array)(x.data.length/2*M);
    F_shape[F_shape.length-1] = M;
    x = x.data;

    for( let x_off=x.length,
             F_off=F.length;
            (F_off -= M) >= 0; )
    {        x_off -= 2;
      const x1 = x[x_off+0],
            x2 = x[x_off+1];

      for( let  i=M; i-- > 0; )
        F[F_off + i] = 4 + 2*i - Math.exp( (i+1)*x1 ) - Math.exp( (i+1)*x2 );
    }

    return new NDArray(F_shape, F);
  }


  lsq_jac( x )
  {
    x = asarray(x);

    const M = this.nOut;

    if(         x.ndim     <  1 ) throw new Error(`${this.name}.lsq_jac(x): x.ndim must be at least 1.`);
    if( x.shape[x.ndim-1] !== 2 ) throw new Error(`${this.name}.lsq_jac(x): x.shape[-1] must be 2.`);
  
    const J_shape = Int32Array.of(...x.shape.slice(0,-1),M,2),
          J       = new (x.dtype==='float32' ? Float32Array
                                             : Float64Array)(x.data.length*M);
    x = x.data;
  
    for( let x_off=x.length,
             J_off=J.length;
            (J_off -= M*2) >= 0; )
    {        x_off -= 2;
      const x1 = x[x_off+0],
            x2 = x[x_off+1];

      for( let  i=M; i-- > 0; ) {
        J[J_off + 2*i+0] = - (i+1) * Math.exp( (i+1)*x1 );
        J[J_off + 2*i+1] = - (i+1) * Math.exp( (i+1)*x2 );
      }
    }
  
    return new NDArray(J_shape, J);
  }
}
