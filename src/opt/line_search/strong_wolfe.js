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
import {asarray} from '../../nd_array'

import {_min1d_interp_quad} from '../_opt_utils'
//*DEBUG*/import {num_grad} from '../num_grad'

import {LineSearchError, LineSearchNoProgressError} from './line_search_error'


export const strong_wolfe = ({c1=0.01, c2=0.9, c3=2, minRed=0.2}={}) => {
  // SEE:
  //   "Numerical Optimization" 2n Edition,
  //   Jorge Nocedal Stephen J. Wright,
  //   Chapter 3. Line Search Methods, page 60.

  c1    *= 1; // (1) Wolfe cond. 1: Sufficient decrease of objective.
  c2    *= 1; // (2) Wolfe cond. 2: Sufficient decrease of projected gradient.
  c3    *= 1; // (3) Growth factor for 1st phase of line search (bracketing).
  minRed*= 1; // (4) Minimum reduction of search range during 2nd phase of line search (zoom).

  // CHECK 0 < c1 < c2 < 1 < c3
  if( !(c1 >  0) ) throw new Error('strong_wolfe(opt): opt.c1 must be positive.'  );
  if( !(c1 < c2) ) throw new Error('strong_wolfe(opt): opt.c1 must less than opt.c2.' );
  if( !(1  > c2) ) throw new Error('strong_wolfe(opt): opt.c2 must less than 1.'  );
  if( !(1  < c3) ) throw new Error('strong_wolfe(opt): opt.c3 must larger than 1.');
  if( !(minRed >= 0.0) ) throw new Error('Assertion failed.');
  if( !(minRed <= 0.5) ) throw new Error('Assertion failed.');

  const Λ = 1 / minRed,
        λ = Λ - 1;

  return fg => (X0,f0,G0, negDir) => {
    X0 = asarray(X0)
    G0 = asarray(G0)
    if    (X0.ndim !== 1) throw new Error('strong_wolfe()(fg)(X0,f0,G0): X0.ndim must be 1.')
    if(    G0.ndim !== 1) throw new Error('strong_wolfe()(fg)(X0,f0,G0): G0.ndim must be 1.')
    if(negDir.ndim !== 1) throw new Error('strong_wolfe()(fg)(X0,f0,G0): negDir.ndim must be 1.')
    if(X0.shape[0] !== G0.shape[0]) throw new Error('strong_wolfe()(fg)(X0,f0,G0): X0 and G0 must have the same shape.')
    if( isNaN(f0) ) throw new Error('strong_wolfe()(fg)(X0,f0,G0): f0 is NaN.')

    const dtype = [X0,G0].every(x => x.dtype==='float32') ? 'float32' : 'float64'

//*DEBUG*/    const projGradTest = num_grad( α => {
//*DEBUG*/      const X = zip_elems([X0,negDir], dtype, (x,nDir) => x - α*nDir);
//*DEBUG*/      return fg(X)[0];
//*DEBUG*/    });

    const projGrad = g => { // <- projected gradient
      let pg = 0
      for( let [i]=negDir.shape; i-- > 0; )
        pg -= negDir.data[i] * g.data[i]
      return pg
    }

    const p0 = projGrad(G0); // <- projected gradient

    if( p0 >= 0 ) throw new Error('strong_wolfe: Initial projected gradient not negative.');

    let αMin =  0, α = 1,
        fMin = f0,
        pMin = p0,

        αMax = Infinity,
        fMax = NaN;

    // STEP 1: BRACKETING PHASE
    //   Find a range guaranteed to contain an α satisfying strong Wolfe.
    bracketing: for(;;)
    {
      const X = zip_elems([X0,negDir], dtype, (x,nDir) => x - α*nDir),
           [f,G] = fg(X),
            p = projGrad(G);

//*DEBUG*/      const q = projGradTest(α);
//*DEBUG*/      if( !(Math.abs(p-q) <= Math.max(Math.abs(p),Math.abs(q),1)*1e-5) ) throw Error(`Assertion failed: ${p} !== ${q}`);

      if( f - f0 > c1*α*p0 || (0 < αMin && f >= fMin) )
      {
        αMax = α;
        fMax = f;
        break bracketing;
      }

      if( Math.abs(p) <= -c2*p0 )
        return [X,f,G];

      if( p >= 0 )
      {
        αMax = αMin;
        fMax = fMin;

        αMin = α;
        fMin = f;
        pMin = p;
        break bracketing;
      }

      if( ! (α < αMax) )
        throw new LineSearchError('strong_wolfe: Strong Wolfe condition not satisfiable in range.');

      αMin = α; α *= c3;
      fMin = f;
      pMin = p;
    }

    if( αMin === αMax ) throw new LineSearchError('strong_wolfe: bracketing failed.');

    const noProgress = () => {
      if( αMin === 0 ) throw new LineSearchNoProgressError('strongWolfeLineSearch(): no progress.');
      throw new LineSearchError(`strong_wolfe: bisection failed ${JSON.stringify({αMin, α, αMax})}.`);
    }

    // STEP 2: ZOOM PHASE
    //   Given a range that is guaranteed to contain a valid
    //   strong Wolfe α values, this method finds such a value.
    for( let run=1;; run++ )
    {
//*DEBUG*/      const q = projGradTest(αMin);
//*DEBUG*/      if( !(Math.abs(pMin-q) <= Math.max(Math.abs(pMin),Math.abs(q),1)*1e-5) ) throw Error(`Assertion failed: ${pMin} !== ${q}`);

      if( 2===Λ )
        α = (αMin + αMax) / 2;
      else {
        α = _min1d_interp_quad(
          αMin,αMax,
          fMin,fMax,
          pMin
        );

        const αLo = (αMax + λ*αMin) / Λ,
              αHi = (αMin + λ*αMax) / Λ;
              
//*DEBUG*/        if( !(Math.min(αMin,αMax) <= α) ) throw new Error('Assertion failed.');
//*DEBUG*/        if( !(Math.max(αMin,αMax) >= α) ) throw new Error('Assertion failed.');

             if( !(αLo <= α) ) α = αLo; // < handles NaN
        else if( !(αHi >= α) ) α = αHi; // < handles NaN
      }

      const X = zip_elems([X0,negDir], dtype, (x,nDir) => x - α*nDir),
           [f,G] = fg(X),
            p = projGrad(G);

      if( f - f0 > c1*α*p0 || f >= fMin ) {
        if( αMax === α )
          noProgress()
        αMax = α;
        fMax = f;
      }
      else {
        if( Math.abs(p)  <= -c2*p0 )
          return [X,f,G];

        if( p * (αMax - αMin) >= 0 ) {
          αMax = αMin;
          fMax = fMin;
        }

        if( αMin === α )
          noProgress()
        αMin = α;
        fMin = f;
        pMin = p;
      }
    }
  }
}
