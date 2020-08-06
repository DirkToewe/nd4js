'use strict';

/* This file is part of ND.JS.
 *
 * ND.JS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ND.JS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ND.JS. If not, see <http://www.gnu.org/licenses/>.
 */

import {array, NDArray} from '../nd_array'
import {ARRAY_TYPES, eps} from '../dt'


// https://en.wikipedia.org/wiki/Finite_difference#Forward,_backward,_and_central_differences
// https://www.geometrictools.com/Documentation/FiniteDifferences.pdf
const NUM_GRAD_D = Float64Array.of(+2,+1,-1,-2),
      NUM_GRAD_W = Float64Array.of(-1, 8,-8,+1),
      NUM_GRAD_S = 12;


export const num_grad = ( f, {h_rel=undefined, h_abs=undefined}={} ) => x => {
  x = array('float', x);

  const dtype = x.dtype==='float32' ? 'float32' : 'float64',
        DTypeArray = ARRAY_TYPES[dtype],
        x_shape = x.shape;

  x = x.data;

  let g = null,
      g_shape = null;

  const 
    epsRel = null != h_rel ? h_rel : eps(dtype)**(1/3),
    epsAbs = null != h_abs ? h_abs : eps(dtype)**(1/3);

  let h = Math.max(Math.abs(x[0])*epsRel, epsAbs);

  // first function evaluation to determine output shape
  let        y = function() {
    const    X = new NDArray( x_shape, DTypeArray.from(x) );
             X.data[0] += h * NUM_GRAD_D[0];
    return f(X);
  }();

  if( 'number' === typeof y )
  { // SCALAR OUTPUT
    // -------------
    g_shape = x_shape;
    g = new Float64Array(x.length);
    g[0] = y * NUM_GRAD_W[0];

    outer_loop:for( let i=0,
                        j=1;; j=0 )
    {
      for( ; j < NUM_GRAD_D.length; j++ )
      {
        const X = new NDArray( x_shape, DTypeArray.from(x) );
              X.data[i] += h * NUM_GRAD_D[j];
        y = f(X);
        g[i] += y * NUM_GRAD_W[j];
      }
      g[i] /= h * NUM_GRAD_S;

      if( ! (++i < x.length) )
        break outer_loop;

      h = Math.max(Math.abs(x[i])*epsRel, epsAbs);
    }
  }
  else
  { // ND-ARRAY OUTPUT
    // ---------------
    y = array('float', y);

    const     y_shape = y.shape;
          y = y.data;
    const M = y.length,
          N = x.length;

    g_shape = Int32Array.of(...y_shape, ...x_shape);
    g = new Float64Array(y.length * x.length);

    for( let k=0; k < M; k++ )
      g[N*k] = y[k] * NUM_GRAD_W[0];

    outer_loop:for( let i=0,
                        j=1;; j=0 )
    {
      for( ; j < NUM_GRAD_D.length; j++ )
      {
        const X = new NDArray( x_shape, DTypeArray.from(x) );
              X.data[i] += h * NUM_GRAD_D[j];
        y = array('float', f(X));

        if( y.ndim !== y_shape.length )
          throw new Error('num_grad(f)(x): return shape of f must be consistent.');
        for( let i=y_shape.length; i-- > 0; )
          if( y.shape[i] !== y_shape[i] )
            throw new Error('num_grad(f)(x): return shape of f must be consistent.');

        y = y.data;

        for( let k=0; k < M; k++ )
          g[N*k+i] += y[k] * NUM_GRAD_W[j];
      }

      for( let k=0; k < M; k++ )
        g[N*k+i] /= h * NUM_GRAD_S;

      if( ! (++i < N) )
        break outer_loop;

      h = Math.max(Math.abs(x[i])*epsRel, epsAbs);
    }
  }

  return new NDArray(g_shape, g)
}


export const num_grad_forward = ( f, {h_rel=undefined, h_abs=undefined}={} ) => x => {
  x = array(x)

  const dtype = x.dtype==='float32' ? 'float32' : 'float64',
        DTypeArray = ARRAY_TYPES[dtype],
        shape = x.shape
  x = x.data;

  const  g = new DTypeArray(x.length),
    epsRel = null != h_rel ? h_rel : eps(dtype)**(1/3),
    epsAbs = null != h_abs ? h_abs : eps(dtype)**(1/3);

  const f0 = f( new NDArray(shape, x.slice()) );
  g.fill(-1.5*f0);

  // http://macs.citadel.edu/chenm/344.dir/temp.dir/lect4_1.pdf
  for( let i=x.length; i-- > 0; )
  {
    const h = Math.max(Math.abs(x[i])*epsRel, epsAbs) // TODO maybe a check for h===0 might be sensible
    for( const [d,w] of [[1*h,+2.0],
                         [2*h,-0.5]] )
    {
      const     X = new NDArray( shape, DTypeArray.from(x) );
                X.data[i] += d;
      g[i] += f(X) * w;
    }
    g[i] /= h;
  }

  return new NDArray(shape, g)
}
