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

import {strong_wolfe} from './line_search/strong_wolfe'
import {array, NDArray} from '../nd_array'
import {LineSearchNoProgressError, LineSearchError} from './line_search/line_search_error'


// TODO: add scaling as described in:
//  "On the Limited Memory BFGS Method for Large Scale Optimization."
//   by Dong C. Liu and Jorge Nocedal


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
    lineSearch=strong_wolfe(),
    negDir0 = G => G.mapElems('float64',g=>g/1024)
  } = {}
)
{
  let x = array('float64', x0);
                           x0 = undefined;
  let [f,g] = fg(x);
  if( ! f.dtype.startsWith('float') ) throw new Error('min_lbfgs_gen(fg, x0, opt): fg must return [float,float[]].');
  if( ! g.dtype.startsWith('float') ) throw new Error('min_lbfgs_gen(fg, x0, opt): fg must return [float,float[]].');

  const [L] = x.shape,
      shape = x.shape
  lineSearch = lineSearch(fg)

  let dX = [],
    dGdX = [],
    dG   = []

  const negDir = g => {
    //SEE:
    // Jorge Nocedal "Updating Quasi-Newton Matrices with Limited Storage"
    // MATHEMATICS OF COMPUTATION, VOL. 35, NO. 151, JULY 1980, PAGES 773-78
    // https://courses.engr.illinois.edu/ece544na/fa2014/nocedal80.pdf

          g =     Float64Array.from(g.data)
    const α = new Float64Array(dX.length)

    for( let i=dGdX.length; i-- > 0; )
    {
      α[i] = dot(dX[i],g) / dGdX[i]
      for( let j=L; j-- > 0; )
        g[j] -= α[i] * dG[i][j]
    }

    g = negDir0(new NDArray(shape,g))
    g = Float64Array.from(g.data)
  
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
    const [X,F,G] = lineSearch( x,f,g, negDir(g) )

    x = x.data
    g = g.data
    const dg = Float64Array.from(G.data, (Gi,i) => Gi - g[i]),
          dx = Float64Array.from(X.data, (Xi,i) => Xi - x[i])

    if( !(0 < dot(dx,dg)) ) throw new LineSearchNoProgressError();

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
    if( dX.length >   historySize ) throw new Error('Assertion#1 failed.')
    if( dX.length !== dGdX.length ) throw new Error('Assertion#2 failed.')
    if( dX.length !== dG  .length ) throw new Error('Assertion#3 failed.')

    x = X
    f = F
    g = G
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

    try{
      step()
    }
    catch(err) {
      if( err instanceof LineSearchError ) {
        // single attempt to restart
          dX = [],
        dGdX = [],
        dG   = []
        step()
      }
      else throw err
    }
  }
}
