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

import {asarray, NDArray} from '../nd_array'

import {nextUp} from "../dt/float64_utils";

import {FrobeniusNorm} from '../la/norm'

import {OptimizationNoProgressError} from "./optimization_error";
import {TrustRegionSolverLSQ} from './_trust_region_solver_lsq'


// References
// ----------
// .. [1] "THE LEVENBERG-MARQUARDT ALGORITHM: IMPLEMENTATION AND THEORY"
//         Jorge J. Moré.


/** Solves a nonlinear least-squares problem using a
 *  trust-region variant Levenberg-Marquard algorithm
 *  as described in:
 *
 *    "THE LEVENBERG-MARQUARDT ALGORITHM: IMPLEMENTATION AND THEORY"
 *    by Jorge J. Moré.
 */
export const lsq_lm_gen = (fJ, x0, opt) => _lm(new TrustRegionSolverLSQ(fJ,x0), opt);

export function* _lm(
  solver,
  /*opt=*/{
    r0   = 1.1,               // <-   START TRUST REGION RADIUS. r0=1 MEANS THE TRUST REGION CONTAINS THE CAUCHY POINT ON ITS RADIUS.
    rMin = 0,                 // <- MINIMUM TRUST REGION RADIUS (meaning radius AFTER scaling) TODO should be scaled, i.e. relative to norm(D) or something // TODO rMin could be computed to be the smallest radius that still allows some progress, i.e. x[i] + rMin = 
    rMax = Infinity,          // <- MAXIMUM TRUST REGION RADIUS (meaning radius AFTER scaling) TODO should be scaled, i.e. relative to norm(D) or something
    rNewton  = 2,             // <- IF GAUSS-NEWTON INSIDE TRUST RADIUS AND IMPROVEMENT OKAY (according to expectGainMin), RADIUS IS SET TO rGN*distance(x,gaussNewton)
    rTol = 0.05,              // <-         TRUST REGION RADIUS TOLERANCE
    lmLower = 0.001,          // <- LOWER RELATIVE PER-STEP BOUND OF THE LEVENBERG-MARQUARDT PARAMETER ITERATION (USED TO AVOID EXCEEDINGLY SMALL LAMBDA PARAMETERS)
    shrinkLower= 0.05,        // <-  SHRINK TRUST REGION RADIUS BY THIS FACTOR AT MOST
    shrinkUpper= 0.95,        // <-  SHRINK TRUST REGION RADIUS BY THIS FACTOR AT LEAST
    grow = 1.4142135623730951,// <-    GROW TRUST REGION RADIUS BY THIS FACTOR
    expectGainMin = 0.25,     // <- IF ACTUAL IMPROVEMENT RELATIVE TO PREDICTION IS WORSE  THAN THIS FACTOR, SHRINK TRUST REGION
    expectGainMax = 0.75,     // <- IF ACTUAL IMPROVEMENT RELATIVE TO PREDICTION IS BETTER THAN THIS FACTOR, GROW   TRUST REGION
    stuckLimit = 64           // <- MAX. NUMBER OF CONSECUTIVE ITERATIONS WITHOUT PROGRESS
  } = {}
)
{
  if( !( 0 <= lmLower                ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.lmLower must be within [0,1).');
  if( !(      lmLower < 1            ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.lmLower must be within [0,1).');
  if( !(            1 >  rTol        ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.rTol must be less than 1.');
  if( !(            0 <  rTol        ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.rTol must be positive number.');
  if( !(            0 <  r0          ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.r0 must be positive number.');
  if( !(            0 <= rMin        ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.rMin must be non-negative.');
  if( !(         rMin <= rMax        ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.rMin must not be greater than opt.rMax.');
  if( !(         rMin <= r0          ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.r0 must not be less than opt.rMin.');
  if( !(         rMax >= r0          ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.r0 must not be greater than opt.rMax.');
  if( !(   1/(1+rTol) <  rNewton     ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.rGN must be greater than 1/(1+rTol).');
  if( !(            0 <  shrinkLower ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.shrinkLower must be positive number.');
  if( !(  shrinkLower <= shrinkUpper ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.shrinkLower must not be greater than opt.shrinkUpper.');
  if( !(  shrinkUpper < 1            ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.shrinkUpper must be less than 1.');
  if( !(            1 < grow         ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.grow must be number greater than 1.');
  if( !(            0 < expectGainMin) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.expectGainMin must be positive number.');
  if( !(expectGainMin < expectGainMax) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.expectGainMin must be less than expectGainMax.');
  if( !(expectGainMax < 1            ) )    console.warn('lsq_lm_gen(fJ, x0, opt): opt.expectGainMax should be less than 1.');
  if(   0!== stuckLimit%1              ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.stuckLimit must be an integer.');
  if( !(0 <= stuckLimit              ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.stuckLimit must not be negative.');

  const {N, D, regularized_dX: dX, G0: G} = solver,
                                  NORM = new FrobeniusNorm();

//  let R = r0; // <- TODO maybe there is a better starting value for R0 like r0*solver.scaledNorm(G) order something like that
  let         R = -r0 * solver.cauchyTravel() * solver.scaledNorm(G),
           loss = solver.loss,
    stuckometer = 0;

  for(;;)
  {
    if( 0 === stuckometer )
      yield solver.report();

    let [r,dr] = solver.computeNewtonRegularized(0);

    if( ! isFinite(r) ) throw new Error('Assertion failed.');

    const gnInRadius = r <= R*(1+rTol);
    if( ! gnInRadius )
    {
      NORM.reset();
      for( let i=N; i-- > 0; ) {
        const d = D[i];
        let   g = G[i];
        if(0 !== d)  g /= d;
        NORM.include(g);
      }
      r -= R;

      if( ! isFinite(r) ) throw new Error('Assertion failed.');
      if( !(0 <  r) ) throw new Error('Assertion failed.');
      if( !(0 > dr) ) throw new Error('Assertion failed.');

      let λMin = -r / dr,
          λMax = NORM.result / R,
          λ = 0;

/*DEBUG*/      let nIter = -1;
      for(;;)
      {
/*DEBUG*/        if( 64 < ++nIter ) throw new Error('Assertion failed.');
//*DEBUG*/        if( ! (solver.computeNewtonRegularized(λMin)[0] > R) ) throw new Error('Assertion failed.');
//*DEBUG*/        if( ! (solver.computeNewtonRegularized(λMax)[0] < R) ) throw new Error('Assertion failed.');

/*DEBUG*/        if( ! (λMin >= 0   ) ) throw new Error('Assertion failed.');
/*DEBUG*/        if( ! (   0 <= λMax) ) throw new Error('Assertion failed.');
/*DEBUG*/        if( ! (λMin <  λMax) ) throw new Error('Assertion failed.');

        // Algorithm (5.5) (a)
        if( ! (λMin < λ && λ < λMax) )
          λ = Math.max(lmLower*λMax, Math.sqrt(λMin*λMax));

/*DEBUG*/        if( ! (λMin < λ && λ < λMax) ) throw new Error('Assertion failed.');

        // Algorithm (5.5) (b)
        [r,dr] = solver.computeNewtonRegularized(λ);

        if( ! isFinite(r) ) throw new Error('Assertion failed.');

        if( Math.abs(R-r) <= R*rTol )
          break;

        if( r < R )
          λMax = λ;

        λMin = Math.max(λMin, λ + (R-r) / dr);

        // Algorithm (5.5) (c)
        λ += (R-r)/dr * (r/R);
      }

/*DEBUG*/      if( !(Math.abs(solver.scaledNorm(dX)-R) <= R*rTol) ) // <- check that solution is in radius
/*DEBUG*/        throw new Error('Assertion failed.');
    }

/*DEBUG*/    if( !(solver.scaledNorm(dX) <= R*(1+rTol)) )
/*DEBUG*/      throw new Error('Assertion failed.');

    const [
      loss_predict,
      loss_consider
    ] = solver.considerMove(dX);

    /*DEBUG*/    if( !(loss_predict <= loss+1e-6) ) throw new Error('Assertion failed.');

    const rScale = () => 1;
    // const rScale = () => solver.scaledNorm(G); // <- TODO: examine scaling alternatives for the radius limits

    const predict = loss - loss_predict,
           actual = loss - loss_consider;
    if(    actual < predict*expectGainMin )
    {
      // PREDICTION BAD -> REDUCE TRUST RADIUS
      // (by fitting a quadratic polynomial through err_now, grad_now and err_next)
      let grad = 0;
      for( let i=N; i-- > 0; )
        grad += dX[i] * G[i];
      grad *= 2;

      let   shrink = grad / (2*grad + loss - loss_consider);
      if( !(shrink >= shrinkLower) ) shrink = shrinkLower;
      if( !(shrink <= shrinkUpper) ) shrink = shrinkUpper;
      const r = Math.max(rMin*rScale(), R*shrink);
      if(   r === R && !(loss > loss_consider) )
        throw new OptimizationNoProgressError();
      R = r;
    }
    else if( gnInRadius )
    {
      // See [1] Algorithm (7.1) (d). If Gauss-Newton point is inside the trust radius and
      // prediction is not bad, set the radius to rGN times the distance to the Gauss-Newton
      // point. This also avoids that the trust radius grows uncontrollably while the
      // Gauss-Newton point is inside the trust region for a few iterations.
      const        rs = rScale();
      R = Math.max(rs*rMin,
          Math.min(rs*rMax, rNewton*r));
    }
    else if( actual > predict*expectGainMax || actual===0 )
    {
      // PREDICTION GREAT -> INCREASE TRUST RADIUS
      R = Math.min(rMax*rScale(),
          Math.max(nextUp(R), R*grow));
    }

    // ONLY ACCEPT NEW X IF BETTER
    if( loss > loss_consider ) {
        loss = loss_consider;
      solver.makeConsideredMove();
      stuckometer = 0;
    }
    else if( ++stuckometer > stuckLimit )
      throw new OptimizationNoProgressError('Too many unsuccessfull iterations.');
  }
}


export function* fit_lm_gen(
  x, //: float[nSamples,nDim]
  y, //: float[nSamples]
  fg,//: (params: float[nParam]) => (x: float[nDim]) => [f: float, g: float[nParam]],
  p0,//: float[nParam]
  opt
)
{
  if( ! (fg instanceof Function) )
    throw new Error('fit_lm_gen(x,y, fg, p0): fg must be a function.');

  x = asarray(          x ); // <- TODO allow non-ndarray x ?
  y = asarray('float64',y );
  p0= asarray('float64',p0);

  const FloatArray = Float64Array;

  if(x .ndim !== 2 ) throw new Error('fit_lm_gen(x,y, fg, p0): x.ndim must be 2.');
  if(y .ndim !== 1 ) throw new Error('fit_lm_gen(x,y, fg, p0): y.ndim must be 1.');
  if(p0.ndim !== 1 ) throw new Error('fit_lm_gen(x,y, fg, p0): p0.ndim must be 1.');

  const [P] = p0.shape,
      [M,N] =  x.shape,
    x_shape = Int32Array.of(N);

  if( M != y.shape[0] )
    throw new Error('fit_lm_gen(x,y, fg, p0): x.shape[0] must be equal to y.shape[0].');

  x = x.data.slice(); // <- TODO: TypedArray.subarray could be used here
  y = y.data;

  // if x is frozen, no protection copies are necessary while passing it to fg
  Object.freeze(x.buffer);
  x = Array.from(
    {length: M},
    (_,i) => Object.freeze(
      new NDArray(x_shape, x.subarray(N*i,N*i+1))
    )
  );
  Object.freeze(x);

  const R = new FloatArray(M  ),
        J = new FloatArray(M*P),
       RJ =[new NDArray(Int32Array.of(M  ), R),
            new NDArray(Int32Array.of(M,P), J)];

  const fJ = p =>
  {
    const fgp = fg(p);
    if( !(fgp instanceof Function) )
      throw new Error();

    for( let i=M; i-- > 0; )
    {
      let [f,g] = fgp(x[i]);

      f = asarray(f);
      g = asarray(g);

      if( f.ndim     !== 0 ) throw new Error('fit_lm_gen(x,y, fg, p0): fg must have signature float[nParams] => float[nDim] => [float, float[nParams]].');
      if( g.ndim     !== 1 ) throw new Error('fit_lm_gen(x,y, fg, p0): fg must have signature float[nParams] => float[nDim] => [float, float[nParams]].');
      if( g.shape[0] !== P ) throw new Error('fit_lm_gen(x,y, fg, p0): fg must have signature float[nParams] => float[nDim] => [float, float[nParams]].');

      f = f.data[0];
      g = g.data;

      if(       isNaN(f)) throw new Error('fit_lm_gen(x,y, fg, p0): NaN encountered.');
      if(g.some(isNaN)  ) throw new Error('fit_lm_gen(x,y, fg, p0): NaN encountered.');

      R[i] = f - y[i];

      for( let j=P; j-- > 0; )
        J[P*i+j] = g[j];
    }
    return RJ;
  };

  yield* lsq_lm_gen(fJ, p0, opt);
}
