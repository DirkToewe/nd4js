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

import {nextUp,
        nextDown} from "../../dt/float64_utils";

import {LineSearchError,
        LineSearchNoProgressError,
        LineSearchBoundReachedError,
        LineSearchBisectionError} from './line_search_error'
import {_min1d_interp_ffg,
        _min1d_interp_ffgg,
        _min1d_interp_gg} from './_line_search_utils'

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


const DEFAULT_OPTIONS = Object.freeze({
  fRed: 1e-2,
  gRed: 0.9,
  growMin: Math.PI/3,
  growMax: Math.E - 1.5,
  shrinkLeast: 0.1
});


// implements the pseudocode lines (a,b,c) from [2].
export const more_thuente_abc = (opt={}) =>
{
  for( const key of Object.keys(opt) )
    if( ! (key in DEFAULT_OPTIONS) )
      console.warn(`more_thuente_abc(opt): unknown parameter "${key}" in opt.`);

  let {fRed, gRed, growMin, growMax, shrinkLeast} = opt = Object.assign({}, DEFAULT_OPTIONS, opt);

  fRed       *= 1; // (1) Wolfe cond. 1: Sufficient decrease of objective.
  gRed       *= 1; // (2) Wolfe cond. 2: Sufficient decrease of projected gradient.
  growMin    *= 1; // (3) Lower growth factor bound for 1st phase of line search (bracketing).
  growMax    *= 1; // (4) Upper growth factor bound for 1st phase of line search (bracketing).
  shrinkLeast*= 1; // (4) Minimum reduction of search range during 2nd phase of line search (zoom).

  // CHECK 0 < fRed < gRed < 1 < c3
  if( !(fRed    >   0     ) ) throw new Error('more_thuente_abc(opt): opt.fRed must be positive.'  );
  if( !(fRed    <  gRed   ) ) throw new Error('more_thuente_abc(opt): opt.fRed must less than opt.gRed.' );
  if( !(1       >  gRed   ) ) throw new Error('more_thuente_abc(opt): opt.gRed must less than 1.'  );
  if( !(1       <  growMin) ) throw new Error('more_thuente_abc(opt): opt.growMin must larger than 1.');
  if( !(growMax >= growMin) ) throw new Error('more_thuente_abc(opt): opt.growMax must not be less than opt.growMin.');
  if( !(fRed    <  0.5    ) )    console.warn('more_thuente_abc(opt): opt.fRed should be less than 0.5 to work properly with (quasi-)Newton methods.');

  if( !(shrinkLeast >= 0.0) ) throw new Error('more_thuente_abc(opt): opt.shrinkLeast must be non-negative.');
  if( !(shrinkLeast <= 0.5) ) throw new Error('more_thuente_abc(opt): opt.shrinkLeast must not be greater than 0.5.');

  const Λ = 1 / shrinkLeast,
        λ = Λ - 1;

  const line_search = fg_raw => (X0,f0,G0, negDir, αMin=0, α0=undefined, αMax=Infinity) => {
    X0 = array('float64', X0);
    G0 = array('float64', G0);
    f0*= 1;

    if( !(negDir instanceof NDArray) ) throw new Error('Assertion failed.');

    const   shape = negDir.shape,
      [L] = shape;

    if    (X0.ndim !== 1 ) throw new Error('more_thuente_abc()(fg)(X0,f0,G0): X0.ndim must be 1.')
    if(    G0.ndim !== 1 ) throw new Error('more_thuente_abc()(fg)(X0,f0,G0): G0.ndim must be 1.')
    if(negDir.ndim !== 1 ) throw new Error('more_thuente_abc()(fg)(X0,f0,G0): negDir.ndim must be 1.')
    if(X0.shape[0] !== L ) throw new Error('more_thuente_abc()(fg)(X0,f0,G0): negDir and X0 must have the same shape.')
    if(G0.shape[0] !== L ) throw new Error('more_thuente_abc()(fg)(X0,f0,G0): negDir and G0 must have the same shape.')
    if(  isNaN(f0)       ) throw new Error('more_thuente_abc()(fg)(X0,f0,G0): f0 is NaN.')
    if(       αMin !== 0 ) throw new Error('more_thuente_abc()(fg)(X0,f0,G0, αMin,α0,αMax): αMin !== 0 not (yet) supported.')

    negDir = negDir.data;

    // TODO: Instead of a specified αMin, we could compute
    //       the smallest αMin that still results in an
    //       X different from X0, i.e. an αMin that
    //       does not result in complete underflow.

    if( null == α0 )
      α0 = Math.min(1, αMax/2);

    if( !(αMin <= α0) ) throw new Error('more_thuente_abc()(fg)(X0,f0,G0, αMin,α0,αMax): αMin must not be greater than α0.')
    if( !(αMax >= α0) ) throw new Error('more_thuente_abc()(fg)(X0,f0,G0, αMin,α0,αMax): αMax must not be less than α0.')

    if( 0 === αMax )
      throw new LineSearchNoProgressError();

    const projGrad = g => g.data.reduce( (p, gi,i) => p - negDir[i]*gi, 0 );

    // wrap some checks aroung fg_raw
    const xfgp = α => {
      const  X =         new NDArray(shape, X0.data.map((x,i) => x - α*negDir[i]) );
      let [f,G]= fg_raw( new NDArray(shape, X .data.slice()) );

      f *= 1;
      G = array('float64', G);

      if( isNaN(f)        ) throw new Error('more_thuente_abc()(fg)(X0,f0,G0): fg returned NaN.');
      if(G.ndim     !== 1 ) throw new Error('more_thuente_abc()(fg)(X0,f0,G0): fg returned gradient with invalid shape.');
      if(G.shape[0] !== L ) throw new Error('more_thuente_abc()(fg)(X0,f0,G0): fg returned gradient with invalid shape.');

      const         p = projGrad(G);
      return [X,f,G,p];
    };

    const  p0 = projGrad(G0); // <- projected gradient
    if(0===p0) throw new LineSearchNoProgressError();
    if(0 < p0) throw new                     Error('more_thuente_abc: Initial projected gradient not negative.');

    let αLo =  0, αHi = Infinity, α = α0,
        fLo = f0, fHi = NaN,
        pLo = p0, pHi = NaN;

    // PHASE 1: BRACKETING PHASE
    // ------------------------
    //   Find a range guaranteed to contain an α satisfying strong Wolfe.
    bracketing: for(;;)
    {
      if( ! isFinite(α) )
        throw new Error('Assertion failed: ' + α);

      const [X,f,G,p] = xfgp(α);

      // check convergence
      if(  f-f0 <= fRed*α*p0  &&  Math.abs(p) <= -gRed*p0  )
        return [X,f,G];

      if( f > fLo )
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

      if( ! (α < αMax            ) ) throw new LineSearchBoundReachedError(X,f,G);
      if( ! (α < Number.MAX_VALUE) ) throw new Error('more_thuente_abc: Infinity reached.');

      let αTrial;

      if( pLo < p )
        αTrial = _min1d_interp_gg(αLo,α, pLo,p);
      else {
        αTrial = α*growMin;
      }

      αTrial = Math.min(αTrial, α*growMax);
      αTrial = Math.max(αTrial, α*growMin);

      if( ! (α    <  αTrial) ) αTrial = nextUp(α); // <- handles NaN
      if( ! (αMax >= αTrial) ) αTrial = αMax;      // <- handles NaN

      αLo = α; α = αTrial;
      fLo = f;
      pLo = p;
    }

    if( αLo === αHi ) throw new LineSearchError('more_thuente_abc: bracketing failed.');

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
      {
        const αLil = Math.min(αLo,αHi),
              αBig = Math.max(αLo,αHi);

        const αLst = Math.max( nextUp  (αLil), shrinkLeast*αBig + (1-shrinkLeast)*αLil ), // <- Least
              αMst = Math.min( nextDown(αBig), shrinkLeast*αLil + (1-shrinkLeast)*αBig ); // <- Most

        const αMid = αLo + (αHi-αLo) / 2;

        if( shrinkLeast===0.5 || !(αLst < αMst) || ! isFinite(fHi) || ! isFinite(pHi) )
          α = αMid;
        else {
          // TODO: the minimum will satisfy gRed but might not satisfy fRed.
          //       It might be smarter to use the interpolation to find
          //       a point that potentially satisfies both.

          if( fLo < fHi )
          { // See [2], section "Trial Value Selection", Case 1
            const αc = _min1d_interp_ffgg(αLo, αHi, fLo, fHi, pLo, pHi),
                  αq = _min1d_interp_ffg (αLo, αHi, fLo, fHi, pLo);
            α = Math.abs(αc - αLo) < Math.abs(αq - αLo) ? αc : (αc+αq)/2; // <- TODO: Looking for a minimum might not be smart, we should consider trying to find a strong wolfe range instead
          }
          else if( Math.sign(pLo)*pHi < 0 )
          { // See [2], section "Trial Value Selection", Case 1
            const αc = _min1d_interp_ffgg(αLo, αHi, fLo, fHi, pLo, pHi),
                  αs = _min1d_interp_gg  (αLo, αHi,           pLo, pHi);
            α = Math.abs(αs - αHi) <= Math.abs(αc - αHi) ? αc : αs; // <- TODO: Looking for a minimum might not be smart, we should consider trying to find a strong wolfe range instead
          }
          else
          { // α = αMid;
            α = _min1d_interp_ffg(
              αLo,αHi,
              fLo,fHi,
              pLo
            );
          }

               if( !(αLst <= α) ) α = αLst; // < handles NaN
          else if( !(αMst >= α) ) α = αMst; // < handles NaN
        }
      }

      const [X,f,G,p] = xfgp(α);

      // check convergence
      if(  f-f0 <= fRed*α*p0  &&  Math.abs(p) <= -gRed*p0  )
        return [X,f,G];

      if( f > fLo ) {
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

  Object.defineProperty(line_search, 'name', {value: `more_thuente_abc(${JSON.stringify(opt)})`, writable: false});
  Object.assign(line_search, opt);
  Object.freeze(line_search);
  return line_search;
}
