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

import {array, NDArray} from '../../nd_array'

import {nextUp} from "../../dt/float64_utils";

import {_min1d_interp_quad} from '../_opt_utils'

import {LineSearchError,
        LineSearchNoProgressError,
        LineSearchBoundReachedError,
        LineSearchBisectionError} from './line_search_error'


// References
// ----------
// .. [1] "Numerical Optimization", 2nd Edition,
//         Jorge Nocedal & Stephen J. Wright,
//         Chapter 3. Line Search Methods, page 60f
// .. [2] "Line Search Algorithms with Guaranteed Sufficient Decrease"
//         Jorge J. Moré and David J. Thuente
//         ACM Transactions on Mathematical Software,
//         Vol 20, No. 3, Septermber 1994, Pages 286-307


// TODO: study dcsrch.f and dcstep.f from MINPACK2 for possible improvements
//       https://github.com/scipy/scipy/tree/master/scipy/optimize/minpack2


// implements the pseudocode lines (U1,U2,U3) from [2].
export const more_thuente_u123 = (opt={}) => {
  let   {fRed=0.01, gRed=0.9, grow=Math.PI/3, shrinkLeast=0.2} = opt;
  opt = {fRed,      gRed,     grow,           shrinkLeast    };

  fRed       *= 1; // (1) Wolfe cond. 1: Sufficient decrease of objective.
  gRed       *= 1; // (2) Wolfe cond. 2: Sufficient decrease of projected gradient.
  grow       *= 1; // (3) Growth factor for 1st phase of line search (bracketing).
  shrinkLeast*= 1; // (4) Minimum reduction of search range during 2nd phase of line search (zoom).

  // CHECK 0 < fRed < gRed < 1 < c3
  if( !(fRed   >  0  ) ) throw new Error('more_thuente_u123(opt): opt.fRed must be positive.'  );
  if( !(fRed   < gRed) ) throw new Error('more_thuente_u123(opt): opt.fRed must less than opt.gRed.' );
  if( !(1      > gRed) ) throw new Error('more_thuente_u123(opt): opt.gRed must less than 1.'  );
  if( !(1      < grow) ) throw new Error('more_thuente_u123(opt): opt.grow must larger than 1.');
  if( !(fRed   <  0.5) )    console.warn('more_thuente_u123(opt): opt.fRed should be less than 0.5 to work properly with (quasi-)Newton methods.');

  if( !(shrinkLeast >= 0.0) ) throw new Error('more_thuente_u123(opt): opt.shrinkLeast must be non-negative.');
  if( !(shrinkLeast <= 0.5) ) throw new Error('more_thuente_u123(opt): opt.shrinkLeast must not be greater than 0.5.');

  const Λ = 1 / shrinkLeast,
        λ = Λ - 1;

  const line_search = fg_raw => (X0,f0,G0, negDir, αMin=0, α0=undefined, αMax=Infinity) => {
    X0 = array('float64', X0);
    G0 = array('float64', G0);
    f0*= 1;

    if( !(negDir instanceof NDArray) ) throw new Error('Assertion failed.');

    const   shape = negDir.shape,
      [L] = shape;

    if    (X0.ndim !== 1 ) throw new Error('more_thuente_u123()(fg)(X0,f0,G0): X0.ndim must be 1.')
    if(    G0.ndim !== 1 ) throw new Error('more_thuente_u123()(fg)(X0,f0,G0): G0.ndim must be 1.')
    if(negDir.ndim !== 1 ) throw new Error('more_thuente_u123()(fg)(X0,f0,G0): negDir.ndim must be 1.')
    if(X0.shape[0] !== L ) throw new Error('more_thuente_u123()(fg)(X0,f0,G0): negDir and X0 must have the same shape.')
    if(G0.shape[0] !== L ) throw new Error('more_thuente_u123()(fg)(X0,f0,G0): negDir and G0 must have the same shape.')
    if(  isNaN(f0)       ) throw new Error('more_thuente_u123()(fg)(X0,f0,G0): f0 is NaN.')
    if(       αMin !== 0 ) throw new Error('more_thuente_u123()(fg)(X0,f0,G0, αMin,α0,αMax): αMin !== 0 not (yet) supported.')

    negDir = negDir.data;

    // TODO: Instead of a specified αMin, we could compute
    //       the smallest αMin that still results in an
    //       X different from X0, i.e. an αMin that
    //       does not result in complete underflow.

    if( null == α0 )
      α0 = Math.min(1, αMax/2);

    if( !(αMin <= α0) ) throw new Error('more_thuente_u123()(fg)(X0,f0,G0, αMin,α0,αMax): αMin must not be greater than α0.')
    if( !(αMax >= α0) ) throw new Error('more_thuente_u123()(fg)(X0,f0,G0, αMin,α0,αMax): αMax must not be less than α0.')

    if( 0 === αMax )
      throw new LineSearchNoProgressError();

    const projGrad = g => g.data.reduce( (p, gi,i) => p - negDir[i]*gi, 0 );

    // wrap some checks aroung fg_raw
    const xfgp = α => {
      const  X =         new NDArray(shape, X0.data.map((x,i) => x - α*negDir[i]) );
      let [f,G]= fg_raw( new NDArray(shape, X .data.slice()) );

      f *= 1;
      G = array('float64', G);

      if( isNaN(f)        ) throw new Error('more_thuente_u123()(fg)(X0,f0,G0): fg returned NaN.');
      if(G.ndim     !== 1 ) throw new Error('more_thuente_u123()(fg)(X0,f0,G0): fg returned gradient with invalid shape.');
      if(G.shape[0] !== L ) throw new Error('more_thuente_u123()(fg)(X0,f0,G0): fg returned gradient with invalid shape.');

      const         p = projGrad(G);
      return [X,f,G,p];
    };

    const  p0 = projGrad(G0); // <- projected gradient
    if(0===p0) throw new LineSearchNoProgressError();
    if(0 < p0) throw new                     Error('more_thuente_u123: Initial projected gradient not negative.');

    let αLo =  0, αHi = Infinity, α = α0,
        fLo = f0, fHi = NaN,
        pLo = p0, pHi = NaN;

    // PHASE 1: BRACKETING PHASE
    // ------------------------
    //   Find a range guaranteed to contain an α satisfying strong Wolfe.
    bracketing: for(;;)
    {
      if( ! isFinite(α) )
        throw new Error('Assertion failed.');

      const [X,f,G,p] = xfgp(α);

      // check convergence
      if(  f-f0 <= fRed*α*p0  &&  Math.abs(p) <= -gRed*p0  )
        return [X,f,G];

      if( f - fLo > fRed*(α-αLo)*p0 )
      {
        αHi = α;
        fHi = f;
        pHi = p;
        break bracketing;
      }

      if( p > 0 )
      {
        αHi = αLo; αLo = α;
        fHi = fLo; fLo = f;
        pHi = pLo; pLo = p;
        break bracketing;
      }

      if( ! (α < αMax) ) throw new LineSearchBoundReachedError(X,f,G);
      if( ! (α < αHi ) ) throw new LineSearchError            ('more_thuente_u123: Strong Wolfe condition not satisfiable in range.');

      αLo = α; α = Math.min(Math.max(nextUp(α), α*grow), αMax);
      fLo = f;
      pLo = p;
    }

    if( αLo === αHi ) throw new LineSearchError('more_thuente_u123: bracketing failed.');

    const noProgress = (X,f,G) => {
      if( αLo === 0 ) // <- TODO FIXME
        throw new LineSearchNoProgressError('strongWolfeLineSearch(): no progress.');
      throw new LineSearchBisectionError(X,f,G);
    };

    // PHASE 2: ZOOM PHASE
    // ------------------
    //   Given a range that is guaranteed to contain a valid
    //   strong Wolfe α values, this method finds such a value.
    for( let run=1;; run++ )
    {
      // SELECT A NEW TRIAL VALUE FOR α
      if( 2===Λ || αLo === αHi || ! isFinite(fHi) )
        α = αLo + (αHi-αLo) / 2;
      else {
        // TODO: the minimum will satisfy gRed but might not satisfy fRed.
        //       It might be smarter to use the interpolation to find
        //       a point that potentially satisfies both.
        α = _min1d_interp_quad(
          αLo,αHi,
          fLo,fHi,
          pLo
        );
        const αLst = (αHi + λ*αLo) / Λ, // least
              αMst = (αLo + λ*αHi) / Λ; // most
             if( !(αLst <= α) ) α = αLst; // < handles NaN
        else if( !(αMst >= α) ) α = αMst; // < handles NaN
      }

      const [X,f,G,p] = xfgp(α);

      // check convergence
      if(  f-f0 <= fRed*α*p0  &&  Math.abs(p) <= -gRed*p0  )
        return [X,f,G];

      if( f - fLo > fRed*(α-αLo)*p0 ) {
        if( αHi === α )
          noProgress(X,f,G);
        αHi = α;
        fHi = f;
        pHi = p;
      }
      else {
        if( Math.sign(αHi - αLo) * p > 0 ) {
          αHi = αLo;
          fHi = fLo;
          pHi = pLo;
        }

        if( αLo === α )
          noProgress(X,f,G);

        αLo = α;
        fLo = f;
        pLo = p;
      }
    }
  };

  Object.defineProperty(line_search, 'name', {value: `more_thuente_u123(${JSON.stringify(opt)})`, writable: false});
  Object.assign(line_search, opt);
  Object.freeze(line_search);
  return line_search;
}
