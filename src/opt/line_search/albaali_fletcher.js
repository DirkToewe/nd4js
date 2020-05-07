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

import {asarray, NDArray} from '../../nd_array'

import {nextUp} from "../../dt/float64_utils";

import {_min1d_interp_quad} from '../_opt_utils'
//*DEBUG*/import {num_grad} from '../num_grad'
//*DEBUG*/import {zip_elems} from '../../zip_elems'

import {LineSearchError,
        LineSearchNoProgressError,
        LineSearchBoundReachedError} from './line_search_error'

// References
// ----------
// .. [1] "Numerical Optimization", 2nd Edition,
//         Jorge Nocedal & Stephen J. Wright,
//         Chapter 3. Line Search Methods, page 60f
// .. [2] "An Efficient Line Search for Nonlinear Least Squares",
//         M. Al-Baali & R. Fletcher
//         Journal Of Optimization Theory and Applications: Vol. 48, No. 3, MARCH 1986
// .. [3] "Line Search Algorithms with Guaranteed Sufficient Decrease"
//         Jorge J. Moré and David J. Thuente
//         ACM Transactions on Mathematical Software,
//         Vol 20, No. 3, Septermber 1994, Pages 286-307


// Implementation of Scheme S2 from [2].
export const albaali_fletcher = ({c1=0.1, c2=0.9, c3=Math.PI/3, minRed=0.2}={}) => {
  c1    *= 1; // (1) Wolfe cond. 1: Sufficient decrease of objective.
  c2    *= 1; // (2) Wolfe cond. 2: Sufficient decrease of projected gradient.
  c3    *= 1; // (3) Growth factor for 1st phase of line search (bracketing).
  minRed*= 1; // (4) Minimum reduction of search range during 2nd phase of line search (zoom).

  // CHECK 0 < c1 < c2 < 1 < c3
  if( !(c1 >  0) ) throw new Error('albaali_fletcher(opt): opt.c1 must be positive.'  );
  if( !(c1 < c2) ) throw new Error('albaali_fletcher(opt): opt.c1 must less than opt.c2.' );
  if( !(1  > c2) ) throw new Error('albaali_fletcher(opt): opt.c2 must less than 1.'  );
  if( !(1  < c3) ) throw new Error('albaali_fletcher(opt): opt.c3 must larger than 1.');
  if( !(c1 < 0.5)) console.warn('albaali_fletcher(opt): opt.c1 should be less than 0.5 to work properly with (quasi-)Newton methods.');
  if( !(minRed >= 0.0) ) throw new Error('Assertion failed.');
  if( !(minRed <= 0.5) ) throw new Error('Assertion failed.');

  const Λ = 1 / minRed,
        λ = Λ - 1;

  return fg_raw => (X0,f0,G0, negDir, αMin=0, α0=undefined, αMax=Infinity) => {
    X0 = asarray(X0);
    G0 = asarray(G0);
    f0*= 1;

    if( !(negDir instanceof NDArray) ) throw new Error('Assertion failed.');

    const shape = negDir.shape;

    if    (X0.ndim !== 1       ) throw new Error('albaali_fletcher()(fg)(X0,f0,G0): X0.ndim must be 1.')
    if(    G0.ndim !== 1       ) throw new Error('albaali_fletcher()(fg)(X0,f0,G0): G0.ndim must be 1.')
    if(negDir.ndim !== 1       ) throw new Error('albaali_fletcher()(fg)(X0,f0,G0): negDir.ndim must be 1.')
    if(X0.shape[0] !== shape[0]) throw new Error('albaali_fletcher()(fg)(X0,f0,G0): negDir and X0 must have the same shape.')
    if(G0.shape[0] !== shape[0]) throw new Error('albaali_fletcher()(fg)(X0,f0,G0): negDir and G0 must have the same shape.')
    if(  isNaN(f0)             ) throw new Error('albaali_fletcher()(fg)(X0,f0,G0): f0 is NaN.')

    if(   αMin !== 0               ) throw new Error('albaali_fletcher()(fg)(X0,f0,G0, αMin,α0,αMax): αMin !== 0 not (yet) supported.')
    if( !(αMax > Number.MAX_VALUE) ) throw new Error('albaali_fletcher()(fg)(X0,f0,G0, αMin,α0,αMax): αMin !== 0 not (yet) supported.')

    if( null == α0 )
      α0 = Math.min(1, αMax/2);

    if( !(αMin <  α0) ) throw new Error('albaali_fletcher()(fg)(X0,f0,G0, αMin,α0,αMax): αMin !== 0 not (yet) supported.')
    if( !(αMax >= α0) ) throw new Error('albaali_fletcher()(fg)(X0,f0,G0, αMin,α0,αMax): αMin !== 0 not (yet) supported.')

    const dtype = [X0,G0].every(x => x.dtype==='float32') ? 'float32' : 'float64'

    // wrap some checks aroung fg_raw
    const fg = x => {
      let [f,g] = fg_raw(x);

      f *= 1;
      g = asarray(g);

      if( isNaN(f)              ) throw new Error('albaali_fletcher()(fg)(X0,f0,G0): fg returned NaN.');
      if(g.ndim     !== 1       ) throw new Error('albaali_fletcher()(fg)(X0,f0,G0): fg returned gradient with invalid shape.');
      if(g.shape[0] !== shape[0]) throw new Error('albaali_fletcher()(fg)(X0,f0,G0): fg returned gradient with invalid shape.');

      return [f,g];
    };

    negDir = negDir.data;

//*DEBUG*/    const projGradTest = num_grad( α => {
//*DEBUG*/      const X = zip_elems([X0,negDir], dtype, (x,nDir) => x - α*nDir);
//*DEBUG*/      return fg(X)[0];
//*DEBUG*/    });

    // computes the projected gradient
    const projGrad = g => { // <- projected gradient
      let pg = 0
      for( let i=negDir.length; i-- > 0; )
        pg -= negDir[i] * g.data[i]
      return pg
    }

    const p0 = projGrad(G0); // <- projected gradient

    if( p0 >= 0 ) throw new Error('albaali_fletcher: Initial projected gradient not negative.');

    let αLo =  0, αHi = Infinity, α = 1,
        fLo = f0, fHi = NaN,
        pLo = p0, pHi = NaN;

    // STEP 1: BRACKETING PHASE
    // ------------------------
    //   Find a range guaranteed to contain an α satisfying strong Wolfe.
    bracketing: for(;;)
    {
      if( ! isFinite(α) )
        throw new Error('Assertion failed.');

      const X = new NDArray( shape, X0.data.map((x,i) => x - α*negDir[i]) ),
         [f,G]= fg(X),
            p = projGrad(G);

//*DEBUG*/      const q = projGradTest(α);
//*DEBUG*/      if( !(Math.abs(p-q) <= Math.max(Math.abs(p),Math.abs(q),1)*1e-5) ) throw Error(`Assertion failed: ${p} !== ${q}`);

      if( f - f0 > c1*α*p0 || (0 < αLo && f >= fLo) )
      {
        αHi = α;
        fHi = f;
        pHi = p;
        break bracketing;
      }

      if( Math.abs(p) <= -c2*p0 )
        return [X,f,G];

      if( p >= 0 )
      {
        αHi = αLo; αLo = α;
        fHi = fLo; fLo = f;
        pHi = pLo; pLo = p;
        break bracketing;
      }

      if( α === αMax )
        throw new LineSearchBoundReachedError('albaali_fletcher: αMax reached without finding solution.');

      if( ! (α < αHi) )
        throw new LineSearchError('albaali_fletcher: Strong Wolfe condition not satisfiable in range.');

      αLo = α; α = Math.min(Math.max(nextUp(α), α*c3), αMax);
      fLo = f;
      pLo = p;
    }

    if( αLo === αHi ) throw new LineSearchError('albaali_fletcher: bracketing failed.');

    const noProgress = () => {
      if( true || αLo === 0 ) // <- TODO FIXME
        throw new LineSearchNoProgressError('strongWolfeLineSearch(): no progress.');
      const str = Object.entries({p0, αLo, α, αHi, fLo, fHi, pLo, pHi}).map(([k,v]) => k + ':\t' + v).join('\n  ');
      throw new LineSearchError(`albaali_fletcher: bisection failed {\n  ${str}\n}.`);
    }

    // STEP 2: ZOOM PHASE
    // ------------------
    //   Given a range that is guaranteed to contain a valid
    //   strong Wolfe α values, this method finds such a value.
    for( let run=1;; run++ )
    {
//*DEBUG*/    ;{
//*DEBUG*/      const q = projGradTest(αLo);
//*DEBUG*/      if( !(Math.abs(pLo-q) <= Math.max(Math.abs(pLo),Math.abs(q),1)*1e-5) ) throw Error(`Assertion failed: ${pLo} !== ${q}`);
//*DEBUG*/    };
//*DEBUG*/    ;{
//*DEBUG*/      const q = projGradTest(αHi);
//*DEBUG*/      if( !(Math.abs(pHi-q) <= Math.max(Math.abs(pHi),Math.abs(q),1)*1e-5) ) throw Error(`Assertion failed: ${pHi} !== ${q}`);
//*DEBUG*/    };

      if( 2===Λ || αLo === αHi || ! isFinite(fHi) )
        α = (αLo + αHi) / 2;
      else {
        α = _min1d_interp_quad(
          αLo,αHi,
          fLo,fHi,
          pLo
        );

        const αLst = (αHi + λ*αLo) / Λ, // least
              αMst = (αLo + λ*αHi) / Λ; // most
            
//*DEBUG*/        if( !(Math.min(αLo,αHi) <= α) ) throw new Error('Assertion failed.');
//*DEBUG*/        if( !(Math.max(αLo,αHi) >= α) ) throw new Error('Assertion failed.');

             if( !(αLst <= α) ) α = αLst; // < handles NaN
        else if( !(αMst >= α) ) α = αMst; // < handles NaN
      }

      const X = new NDArray( shape, X0.data.map((x,i) => x - α*negDir[i]) ),
         [f,G]= fg(X),
            p = projGrad(G);

      if( f - f0 > c1*α*p0 || f >= fLo ) {
        if( αHi === α )
          noProgress()
        αHi = α;
        fHi = f;
        pHi = p;
      }
      else {
        if( Math.abs(p)  <= -c2*p0 )
          return [X,f,G];

        if( Math.sign(αHi - αLo) * p >= 0 ) {
          αHi = αLo;
          fHi = fLo;
          pHi = pLo;
        }

        if( αLo === α )
          noProgress()

        αLo = α;
        fLo = f;
        pLo = p;
      }
    }
  }
}
