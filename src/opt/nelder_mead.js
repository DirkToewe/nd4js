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

import {array,
      NDArray} from '../nd_array'

import {regular_simplex} from "../geom/simplex";

import {matmul2} from "../la/matmul";
import {FrobeniusNorm} from '../la/norm';

import {argmax,
        argmin} from '../iter/min_max'

import {AleaRNG} from "../rand/alea_rng";

import {OptimizationNoProgressError} from "./optimization_error";


// REFERENCES
// ----------
// .. [1] "A  simplex  method  for  function  minimization."
//         J. A. Nelder and R. Mead
const RNG = new AleaRNG('opt/nelder_mead.js');


export function* min_nelder_mead_gen(
  f,
  x0,
  {
    scale   = 0.001,
    reflect = 1,
    expand  = Math.E/8,
    worstRatio    = 1e-10,
    contractOuter = 2 / (1 + Math.sqrt(5)),
    contractInner = 2 / (1 + Math.sqrt(5)),
    shrink        = 2.5/Math.PI,
    rng = RNG
  } = {}
)
{
  if( !(f instanceof Function) ) throw new Error('min_nelder_mead_gen(f,x0): f has to be a Function.');

  if( ! (0 < scale) ) throw new Error('min_nelder_mead_gen(f,x0,opt={}): opt.scale must be greater than 0.');

  if( ! (0 < reflect) ) throw new Error('min_nelder_mead_gen(f,x0,opt={}): opt.reflect must be greater than 0.');

  if( ! (0 <= expand ) ) throw new Error('min_nelder_mead_gen(f,x0,opt={}): opt.expand must be non-negative.');
  if( ! (0 <  expand ) )    console.warn('min_nelder_mead_gen(f,x0,opt={}): opt.expand should be greater than 0.');

  if( ! (0 <  contractOuter) ) throw new Error('min_nelder_mead_gen(f,x0,opt={}): opt.contractOuter must be greater than 0.');
  if( ! (1 >= contractOuter) ) throw new Error('min_nelder_mead_gen(f,x0,opt={}): opt.contractOuter must not be greater than 1.');
  if( ! (1 >  contractOuter) )    console.warn('min_nelder_mead_gen(f,x0,opt={}): opt.contractOuter must should be less than 1.');

  if( ! (0 <  contractInner) ) throw new Error('min_nelder_mead_gen(f,x0,opt={}): opt.contractInner must be greater than 0.');
  if( ! (1 >= contractInner) ) throw new Error('min_nelder_mead_gen(f,x0,opt={}): opt.contractInner must not be greater than 1.');
  if( ! (1 >  contractInner) )    console.warn('min_nelder_mead_gen(f,x0,opt={}): opt.contractInner must should be less than 1.');

  if( ! (0 <  shrink) ) throw new Error('min_nelder_mead_gen(f,x0,opt={}): opt.shrink must be greater than 0.');
  if( ! (1 >= shrink) ) throw new Error('min_nelder_mead_gen(f,x0,opt={}): opt.shrink must not be greater than 1.');
  if( ! (1 >  shrink) )    console.warn('min_nelder_mead_gen(f,x0,opt={}): opt.shrink must should be less than 1.');

  x0 = array('float64', x0);

  if( x0.ndim !== 1 ) throw new Error('min_nelder_mead_gen(f,x0): x0.ndim has to be 1.');

  const   N = x0.shape[0],
    shape_x = x0.shape;

  const verts_f = new Float64Array(N+1),
        centroid= new Float64Array(N  ),
              x = x0.data;
                  x0 = undefined;

  f = function(){
    const F = f;
    return (x,off=0) => {
      if( !(off >= 0         ) ) throw new Error('Assertion failed.');
      if( !(off <= x.length-N) ) throw new Error('Assertion failed.');
      if(   off%N !== 0        ) throw new Error('Assertion failed.');

      const    f = 1 * F( new NDArray(shape_x, x.slice(off,off+N)) );
      if(isNaN(f)) throw new Error('Assertion failed.');
      return   f;
    };
  }();

  function distance( x, y, y_off )
  {
    if( !(y_off >= 0         ) ) throw new Error('Assertion failed.');
    if( !(y_off <= y.length-N) ) throw new Error('Assertion failed.');
    if(   y_off%N !== 0        ) throw new Error('Assertion failed.');

    const norm = new FrobeniusNorm();
    for( let i=N; i-- > 0; )
      norm.include(x[i] - y[i+y_off]);
    return norm.result;
  }

  // INIT SIMPLEX
  const {data: verts_x} = matmul2( regular_simplex(N), rng.ortho(N) );

  for( let i=N+1; i-- > 0; )
  for( let j=N  ; j-- > 0; )
    verts_x[N*i+j] = verts_x[N*i+j]*scale + x[j];

  for( let i=N+1; i-- > 0; )
    verts_f[i] = f(verts_x, N*i);

  ;{
    const b = argmin(verts_f);
    // YIELD STARTING POINT (and simplex)
    yield [
      new NDArray( shape_x, verts_x.slice(N*b,N*(b+1)) ),
      verts_f[b]
    ];
  }

  loop: for(;;)
  {
    // FIND WORST VERTEX (and move it to row 0)
    ;{
      const i = argmax(verts_f);

      for( let j=N; j-- > 0; ) {
        const x_ij = verts_x[N*i+j];
                     verts_x[N*i+j] = verts_x[j];
                                      verts_x[j] = x_ij;
      }
      const f0 = verts_f[0];
                 verts_f[0] = verts_f[i];
                              verts_f[i] = f0;
    }

    const  b = argmin(verts_f.subarray(1)) + 1,
      f_best =        verts_f[b],
      f_worst=        verts_f[0];

    // COMPUTE CENTROID OF REMAINING VERTICES (write to x2)
    centroid.fill(0.0);
    for( let i=N+1; i-- > 1; )
    for( let j=N  ; j-- > 0; )
      centroid[j] += verts_x[N*i+j];

    for( let j=N; j-- > 0; )
      centroid[j] /= N;

    let dist_min = +Infinity,
        dist_max = -Infinity;

    for( let i=N+1; i-- > 1; )
    { const               dist = distance(centroid, verts_x,N*i);
      dist_min = Math.min(dist, dist_min);
      dist_max = Math.max(dist, dist_max);
    }

    const dist = 0.5 * distance(centroid, verts_x,0);
    // if( !(0 < dist) )
    //   throw new OptimizationNoProgressError();

    if( 0 < dist )
    {
      // REFLECT
      for( let i=N; i-- > 0; )
        x[i] = centroid[i] + reflect*(centroid[i] - verts_x[i]);

      const f_reflect = f(x);
      if(   f_reflect < f_worst )
      {
        for( let i=N; i-- > 0; )
          verts_x[i] = x[i];
        verts_f[0] = f_reflect;

        if( f_reflect < f_best )
        {
          // EXPAND (if reflection new best)
          if( dist*worstRatio <= dist_min )
          {
            for( let i=N; i-- > 0; )
              x[i] = centroid[i] + expand*(x[i] - centroid[i]);

            const f_expand = f(x);
            if(   f_expand < f_reflect )
            {
              for( let i=N; i-- > 0; )
                verts_x[i] = x[i];
              verts_f[0] = f_expand;
            }
          }

          // YIELD NEW BEST POINT (and simplex)
          yield [
            new NDArray( shape_x, verts_x.slice(0,N) ),
            verts_f[0]
          ];
        }

        continue loop;
      }
    }

    if( dist >= dist_max*worstRatio )
    {
      // OUTER CONTRACTION
      ;{
        for( let i=N; i-- > 0; )
          x[i] = centroid[i] + contractOuter*(centroid[i] - verts_x[i]);

        const f_contract = f(x);
        if(   f_contract < f_worst )
        {
          for( let i=N; i-- > 0; )
            verts_x[i] = x[i];
          verts_f[0] = f_contract;

          if( f_contract < f_best )
            yield [
              new NDArray( shape_x, verts_x.slice(0,N) ),
              verts_f[0]
            ];

          continue loop;
        }
      }

      // INNER CONTRACTION
      ;{
        for( let i=N; i-- > 0; )
          x[i] = centroid[i] - contractInner*(centroid[i] - verts_x[i]);

        const f_contract = f(x);
        if(   f_contract < f_worst )
        {
          for( let i=N; i-- > 0; )
            verts_x[i] = x[i];
          verts_f[0] = f_contract;

          if( f_contract < f_best )
            yield [
              new NDArray( shape_x, verts_x.slice(0,N) ),
              verts_f[0]
            ];

          continue loop;
        }
      }
    }

    let stuck = true;

    // SHRINK
    for( let i=N+1; i-- > 0; )
    if( i !== b )
    {
      for( let j=N  ; j-- > 0; ) { const  x_ij = verts_x[N*i+j] + (1-shrink)*(verts_x[N*b+j] - verts_x[N*i+j]);
        stuck = stuck && verts_x[N*i+j]===x_ij;
                         verts_x[N*i+j] = x_ij;
        if( ! isFinite(x_ij) )
          throw new Error('Assertion failed.');
      }
      verts_f[i] = f(verts_x, N*i);
    }

    if(stuck)
      throw new OptimizationNoProgressError();

    const       c = argmin(verts_f);
    if( verts_f[c] < f_best )
      yield [
        new NDArray( shape_x, verts_x.slice(N*c,N*(c+1)) ),
        verts_f[c]
      ];
  }
}
