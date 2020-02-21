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

import {zip_elems} from '../../zip_elems'
import {LineSearchError, LineSearchNoProgressError} from './line_search_error'
import {asarray} from '../../nd_array'


export const strong_wolfe = ({c1=0.4, c2=0.8, c3=1.6}={}) => {
  // SEE:
  //   "Numerical Optimization" 2n Edition,
  //   Jorge Nocedal Stephen J. Wright,
  //   Chapter 3. Line Search Methods, page 60.

  // CHECK 0 < c1 < c2 < 1 < c3
  if( c1 <=  0 ) throw new Error('strong_wolfe({c1,c2,c3}): c1 must be positive.'  );
  if( c1 >= c2 ) throw new Error('strong_wolfe({c1,c2,c3}): c1 must less than c2.' );
  if(  1 <= c2 ) throw new Error('strong_wolfe({c1,c2,c3}): c2 must less than 1.'  );
  if(  1 >= c3 ) throw new Error('strong_wolfe({c1,c2,c3}): c3 must larger than 1.');

  return fg => (X0,f0,G0, negDir) => {
    X0 = asarray(X0)
    G0 = asarray(G0)
    if(X0.ndim !== 1) throw new Error('strong_wolfe()(fg)(X0,f0,G0): X0.ndim must be 1.')
    if(G0.ndim !== 1) throw new Error('strong_wolfe()(fg)(X0,f0,G0): G0.ndim must be 1.')
    if(X0.shape[0] !== G0.shape[0]) throw new Error('strong_wolfe()(fg)(X0,f0,G0): X0 and G0 must have the same shape.')
    if( isNaN(f0) ) throw new Error('strong_wolfe()(fg)(X0,f0,G0): f0 is NaN.')

    const dtype = [X0,G0].every(x => x.dtype==='float32') ? 'float32' : 'float64'

    const projGrad = g => { // <- projected gradient
      let pg = 0
      for( let [i]=negDir.shape; i-- > 0; )
        pg -= negDir.data[i] * g.data[i]
      return pg
    }

    const p0 = projGrad(G0); // <- projected gradient

    if( p0 >= 0 ) throw new Error('strong_wolfe: Initial projected gradient not negative.');

    let αMin =  0, α = 1, αMax = Infinity,
        fMin = f0;

    // STEP 1: BRACKETING PHASE
    //   Find a range guaranteed to contain an α satisfying strong Wolfe.
    bracketing: while(true)
    {
      const X = zip_elems([X0,negDir], dtype, (x,nDir) => x - α*nDir),
           [f,G] = fg(X),
            p = projGrad(G);

      if( f - f0 > c1*α*p0 || (0 < αMin && f >= fMin) )
      {
        αMax = α;
        break bracketing;
      }

      if( Math.abs(p) <= -c2*p0 )
        return [X,f,G];

      if( p >= 0 )
      {
        αMax = αMin;
        αMin = α;
        fMin = f;
        break bracketing;
      }

      if( ! (α < αMax) )
        throw new LineSearchError('strong_wolfe: Strong Wolfe condition not satisfiable in range.');

      αMin = α; α *= c3;
      fMin = f;
    }

    if( αMin === αMax ) throw new LineSearchError('strong_wolfe: bracketing failed.');

    const noProgress = () => {
      if( αMin === 0 ) throw new LineSearchNoProgressError('strongWolfeLineSearch(): no progress.');
      throw new LineSearchError(`strong_wolfe: bisection failed ${JSON.stringify({αMin, α, αMax})}.`);
    }

    // STEP 2: BISECTION PHASE
    //   Given a range that is guaranteed to contain a valid
    //   strong Wolfe α values, this method finds such a value.
    while(true)
    {
      α = (αMin + αMax) / 2; // <- TODO: use quadratic polynomial to find new point

      const X = zip_elems([X0,negDir], dtype, (x,nDir) => x - α*nDir),
           [f,G] = fg(X),
            p = projGrad(G);

      if( f - f0 > c1*α*p0 || f >= fMin ) {
        if( αMax === α )
          noProgress()
        αMax = α;
      }
      else {
        if( Math.abs(p)  <= -c2*p0 )
          return [X,f,G];

        if( p * (αMax - αMin)  >=  0 )
          αMax = αMin;

        if( αMin === α )
          noProgress()
        αMin = α;
        fMin = f;
      }
    }
  }
}
