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

import {strong_wolfe} from './line_search/strong_wolfe'
import {asarray, NDArray} from '../nd_array'
import {LineSearchNoProgressError} from './line_search/line_search_error'


function dot( u, v )
{
  if(  u.length !== v.length ) throw new Error()
  let result = 0
  for( let i=u.length; i-- > 0; )
    result += u[i] * v[i]
  return result
}

export function* min_lbfgs_gen( fg, x0, {historySize=8, lineSearch=strong_wolfe(), negDir0 = g=>g} = {} )
{
  let x = asarray(x0)//; x0 = undefined
  let [f,g] = fg(x)

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
    if( dX.length >   historySize ) throw new Error('Assertion failed.')
    if( dX.length !== dGdX.length ) throw new Error('Assertion failed.')
    if( dX.length !== dG  .length ) throw new Error('Assertion failed.')

    x = X
    f = F
    g = G
  }

  for(;;)
  {
//    if( ! (x instanceof NDArray) ) throw new Error('Assertion failed.')
//    if( ! (f instanceof NDArray) ) throw new Error('Assertion failed.')
//    if( ! (g instanceof NDArray) ) throw new Error('Assertion failed.')
    if( x.ndim !== 1 ) throw new Error('Assertion failed.')
    if( f.ndim !== 0 ) throw new Error('Assertion failed.')
    if( g.ndim !== 1 ) throw new Error('Assertion failed.')
    if( x.shape[0] !== L ) throw new Error('Assertion failed.')
    if( g.shape[0] !== L ) throw new Error('Assertion failed.')

    yield [x,f,g].map( a => new NDArray(a.shape, a.data.slice()) )

    try{
      step()
    }
    catch(err) {
      if( err instanceof LineSearchNoProgressError ) {
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

//export function min_fbgs( fg, x0, opt={} )
//{
//  const {gtol=1e-4, maxIter=128} = opt
//
//  
//  for( const [x,f,g] of min_lbfgs_gen(fg,x0,opt) )
//  {
//    
//  }
//}
