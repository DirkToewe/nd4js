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

import {array, NDArray} from '../nd_array'

import {LineSearchError} from './line_search/line_search_error'
import {more_thuente_u123} from './line_search/more_thuente_u123'


// REFERENCES
// ----------
// .. [1] "On the Limited Memory BFGS Method for Large Scale Optimization."
//         Dong C. Liu and Jorge Nocedal
// .. [2] "Updating Quasi-Newton Matrices with Limited Storage"
//         Jorge Nocedal 
//         MATHEMATICS OF COMPUTATION, VOL. 35, NO. 151, JULY 1980, PAGES 773-78
//         https://courses.engr.illinois.edu/ece544na/fa2014/nocedal80.pdf


// TODO: Implement scaling as described in [1].


function dot( u, v )
{
  if(  u.length !== v.length ) throw new Error()
  let result = 0
  for( let i=u.length; i-- > 0; )
    result += u[i] * v[i]
  return result
}


export function* min_lbfgs_gen(
  fg,
  x0,
  {
    historySize=8,
    lineSearch=more_thuente_u123(),
    negDir0 = g => {
      const G = g.data;
      for( let i=G.length; i-- > 0; )
        G[i] /= 1024;
      return g;
    },
    // scaling = function(){
    //   let scale = 1e-3;
    //   return {
    //     update( dx, dg ) {
    //       if( dx.ndim !== 1 ) throw new Error('Assertion failed.');
    //       if( dg.ndim !== 1 ) throw new Error('Assertion failed.');
    //       dx = dx.data;
    //       dg = dg.data;
    //       const len = dx.length;
    //       if(   len!==dg.length )
    //         throw new Error('Assertion failed.');

    //       let dot=0,
    //           div=0;
    //       for( let i=len; i-- > 0; )
    //       {
    //         dot += dx[i]*dg[i];
    //         div += dg[i]*dg[i];
    //       }

    //       scale = dot/div;
    //     },
    //     apply( negDir ) {
    //       if( negDir.ndim !== 1 )
    //         throw new Error('Assertion failed.');
    //       const d = negDir.data;

    //       for( let i=d.length; i-- > 0; )
    //         d[i] *= scale;

    //       return negDir;
    //     }
    //   };
    // }()
  } = {}
)
{
  if(  historySize%1 !== 0 ) throw new Error('min_lbfgs_gen(fg, x0, opt): opt.historySize must be an integer.');
  if(!(historySize    >  0)) throw new Error('min_lbfgs_gen(fg, x0, opt): opt.historySize must be positive.');

  let x = array('float64', x0);
                           x0 = undefined;
  let [f,g] = fg(x);
       f *= 1;
  if(                      isNaN(f) ) throw new Error('min_lbfgs_gen(fg, x0, opt): fg must return [float,float[]].');
  if( ! g.dtype.startsWith('float') ) throw new Error('min_lbfgs_gen(fg, x0, opt): fg must return [float,float[]].');

  const  [L] = x.shape,
     shape   = x.shape,
     shape_f = Int32Array.of();
  lineSearch = lineSearch(fg);

  let dX = [],
    dGdX = [],
    dG   = [];

  // computes product H⋅g (H: inverse Hessian, g: gradient vector)
  const compute_Hg = g => {
    // SEE [2] FOR COMPUTATION OF INVERSE HESSIAN
          g =     Float64Array.from(g.data)
    const α = new Float64Array(dX.length)

    for( let i=dGdX.length; i-- > 0; )
    {
      α[i] = dot(dX[i],g) / dGdX[i]
      for( let j=L; j-- > 0; )
        g[j] -= α[i] * dG[i][j]
    }

    // the product H0⋅g
    g = negDir0( new NDArray(shape,g) );
    g = Float64Array.from(g.data);
    // negDir = scaling.apply(negDir)
  
    for( let i=0; i < dGdX.length; i++ )
    {
      const βi = dot(dG[i],g) / dGdX[i]
      for( let j=L; j-- > 0; )
        g[j] += (α[i] - βi) * dX[i][j]
    }
  
    return new NDArray(shape, g)
  }

  const step = () =>
  {
    let negDir = compute_Hg(g)

    let [X,F,G] = lineSearch(
      new NDArray(shape, x.data.slice()), f,
      new NDArray(shape, g.data.slice()), negDir
    );

    X = array('float64', X);
    G = array('float64', G);
    f*= 1;

    if( X.ndim !== 1 ) throw new Error('Assertion failed.');
    if( G.ndim !== 1 ) throw new Error('Assertion failed.');
    if( isNaN(f)     ) throw new Error('Assertion failed.');

    const _x = x.data,
          _g = g.data,
          dg = Float64Array.from(G.data, (Gi,i) => Gi - _g[i]),
          dx = Float64Array.from(X.data, (Xi,i) => Xi - _x[i]);

    if( !(0 < dot(dx,dg)) ) throw new LineSearchError()

    // scaling.update(
    //   new NDArray(shape, dx.slice()),
    //   new NDArray(shape, dg.slice())
    // )

    dG  .push(    dg    )
    dGdX.push(dot(dg,dx))
      dX.push(       dx )

    // LIMIT THE NUMBER OF MEMOIZED GRADIENTS
    // (while loop in case historySize was changed)
    if( dX.length > historySize ) {
        dX.shift()
      dGdX.shift()
      dG  .shift()
    }
    if( dX.length  >  historySize ) throw new Error('Assertion#1 failed.')
    if( dX.length !== dGdX.length ) throw new Error('Assertion#2 failed.')
    if( dX.length !== dG  .length ) throw new Error('Assertion#3 failed.')

    x = X;
    f = F*1;
    g = G;
  }

  for(;;)
  {
//    if( ! (x instanceof NDArray) ) throw new Error('Assertion failed.')
//    if( ! (f instanceof NDArray) ) throw new Error('Assertion failed.')
//    if( ! (g instanceof NDArray) ) throw new Error('Assertion failed.')
    if( x.ndim !== 1 ) throw new Error('Assertion#4 failed.')
    if( isNaN(f*1)   ) throw new Error('Assertion#5 failed.')
    if( g.ndim !== 1 ) throw new Error('Assertion#6 failed.')
    if( x.shape[0] !== L ) throw new Error('Assertion#7 failed.')
    if( g.shape[0] !== L ) throw new Error('Assertion#8 failed.')

    yield [
      new NDArray(shape, x.data.slice()), f*1,
      new NDArray(shape, g.data.slice())
    ];

    loop:for(;;)
      try{
        step()
        break loop;
      }
      catch(err) {
        if( err instanceof LineSearchError && dX.length > 0 ) {
          // forget half and retry
          const             mid = dX.length+1 >>> 1;
            dX =   dX.slice(mid);
          dGdX = dGdX.slice(mid);
          dG   = dG  .slice(mid);
        }
        else throw err
      }
  }
}
