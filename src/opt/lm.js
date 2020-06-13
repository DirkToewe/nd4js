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

import {array, asarray, NDArray} from '../nd_array'

import {nextUp} from "../dt/float64_utils";

import {FrobeniusNorm} from '../la/norm'

import {OptimizationNoProgressError} from "./optimization_error";
import {TrustRegionSolverLSQ} from './_trust_region_solver'

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
export function* lsq_lm_gen(
  fJ,//: (x0: float[nInputs]) => [f: float[nSamples], j: float[nSamples,nInputs]] - THE OPTIMIZATION FUNCTION AND ITS JACOBIAN
  x0,//: float[nInputs] - THE START POINT OF THE OPTIMIZATION
  /*opt=*/{
    r0   = 1.1,               // <-   START TRUST REGION RADIUS. r0=1 MEANS THE TRUST REGION CONTAINS THE CAUCHY POINT ON ITS RADIUS.
    rMin = 0,                 // <- MINIMUM TRUST REGION RADIUS (meaning radius AFTER scaling) TODO should be scaled, i.e. relative to norm(D) or something // TODO rMin could be computed to be the smallest radius that still allows some progress, i.e. x[i] + rMin = 
    rMax = Infinity,          // <- MAXIMUM TRUST REGION RADIUS (meaning radius AFTER scaling) TODO should be scaled, i.e. relative to norm(D) or something
    rGN  = 2,                 // <- IF GAUSS-NEWTON INSIDE TRUST RADIUS AND IMPROVEMENT OKAY (according to expectGainMin), RADIUS IS SET TO rGN*distance(x,gaussNewton)
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
  x0 = array('float64', x0);
  if( x0.ndim !== 1 )
    throw new Error('lsq_dogleg_gen(fJ, x0): x0.ndim must be 1.');

  if( !(fJ instanceof Function) )
    throw new Error('lsq_dogleg_gen(fJ, x0): fJ must be Function.');

  if( !( 0 <= lmLower                ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.lmLower must be within [0,1).');
  if( !(      lmLower < 1            ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.lmLower must be within [0,1).');
  if( !(            1 >  rTol        ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.rTol must be less than 1.');
  if( !(            0 <  rTol        ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.rTol must be positive number.');
  if( !(            0 <  r0          ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.r0 must be positive number.');
  if( !(            0 <= rMin        ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.rMin must be non-negative.');
  if( !(         rMin <= rMax        ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.rMin must not be greater than opt.rMax.');
  if( !(         rMin <= r0          ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.r0 must not be less than opt.rMin.');
  if( !(         rMax >= r0          ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.r0 must not be greater than opt.rMax.');
  if( !(   1/(1+rTol) <  rGN         ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.rGN must be greater than 1/(1+rTol).');
  if( !(            0 <  shrinkLower ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.shrinkLower must be positive number.');
  if( !(  shrinkLower <= shrinkUpper ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.shrinkLower must not be greater than opt.shrinkUpper.');
  if( !(  shrinkUpper < 1            ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.shrinkUpper must be less than 1.');
  if( !(            1 < grow         ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.grow must be number greater than 1.');
  if( !(            0 < expectGainMin) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.expectGainMin must be positive number.');
  if( !(expectGainMin < expectGainMax) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.expectGainMin must be less than expectGainMax.');
  if( !(expectGainMax < 1            ) )    console.warn('lsq_lm_gen(fJ, x0, opt): opt.expectGainMax should be less than 1.');
  if(   0!== stuckLimit%1              ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.stuckLimit must be an integer.');
  if( !(0 <= stuckLimit              ) ) throw new Error('lsq_lm_gen(fJ, x0, opt): opt.stuckLimit must not be negative.');

  // TODO: Instead of a specified rMin, we could compute
  //       the smallest rMin that still results in an
  //       X different from X0, i.e. an rMin that
  //       does not result in complete underflow.

  let X = x0.data,
      W = new Float64Array(X.length);

  const NORM = new FrobeniusNorm();

  const X_shape = x0.shape,
      mse_shape = new Int32Array(0);

  let [f,J] = fJ( new NDArray(X_shape, X.slice()) );
  f = array('float64', f);
  J = array('float64', J);
  if( f.ndim !== 1 ) throw new Error('lsq_dogleg_gen(fJ, x0): fJ(x) returned [f,J] with f.ndim !== 1.');
  if( J.ndim !== 2 ) throw new Error('lsq_dogleg_gen(fJ, x0): fJ(x) returned [f,J] with J.ndim !== 2.');

  const f_shape = f.shape,
        J_shape = J.shape;

  const                  solver = new TrustRegionSolverLSQ(...J.shape),
    {M,N, D, X: dX, G} = solver;

  if( x0.shape[0] !== N )
    throw new Error('lsq_dogleg_gen(fJ, x0): fJ(x) returned [f,J] where J.shape[1] did not match x0.shape[0].');

  solver.update(f,J);

  f = f.data,
  J = J.data;

//  let R = r0; // <- TODO maybe there is a better starting value for R0 like r0*solver.scaledNorm(G) order something like that
  let R = -r0 * solver.cauchyPointTravel() * solver.scaledNorm(G);

  let err_now = 0;
  for( let i=M; i-- > 0; ) {
    const      s = f[i];
    err_now += s*s;
  }

  let stuckometer = 0; // <- indicates how stuck we are, i.e. counts the number of consecutive operations without progress

  for(;;)
  {
    if( stuckometer === 0 )
      yield [
        new NDArray(  X_shape, X.slice()),
        new NDArray(mse_shape, Float64Array.of(err_now / M) ),
        new NDArray(  X_shape, G.map(g => g*2/M) ),
        new NDArray( f_shape, f.slice() ),
        new NDArray( J_shape, J.slice() )
      ];

    solver.computeMinGlobal();

    const              distanceToGN = solver.scaledNorm(dX), // <- distance to Gauss-Newton point
          gnInRadius = distanceToGN <= R*(1+rTol);
    if( ! gnInRadius )
    {
      NORM.reset();
      for( let i=N; i-- > 0; ) {
        const d = D[i];
        let   g = G[i];
        if(0 !== d)  g /= d;
        NORM.include(g);
      }

      // GLOBAL MIN OUTSIDE TRUST REGION -> FIND LAMBDA
      let [r,dr] = solver.computeMinRegularized(0);
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
//*DEBUG*/        if( ! (solver.computeMinRegularized(λMin)[0] > R) ) throw new Error('Assertion failed.');
//*DEBUG*/        if( ! (solver.computeMinRegularized(λMax)[0] < R) ) throw new Error('Assertion failed.');

/*DEBUG*/        if( ! (λMin >= 0   ) ) throw new Error('Assertion failed.');
/*DEBUG*/        if( ! (   0 <= λMax) ) throw new Error('Assertion failed.');
/*DEBUG*/        if( ! (λMin <  λMax) ) throw new Error('Assertion failed.');

        // Algorithm (5.5) (a)
        if( ! (λMin < λ && λ < λMax) )
          λ = Math.max(lmLower*λMax, Math.sqrt(λMin*λMax));

/*DEBUG*/        if( ! (λMin < λ && λ < λMax) ) throw new Error('Assertion failed.');

        // Algorithm (5.5) (b)
        [r,dr] = solver.computeMinRegularized(λ);

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

    let same_x = true; // <- TRUE IF X REMAINS COMPLETELY UNCHANGED (THIS COULD BE AVOIDED BY COMPUTING AN APPROPRIATE rMin)
    for( let i=N; i-- > 0; ) {
      W[i] = dX[i] + X[i];
      same_x &= (W[i] === X[i]);
    }

    let err_predict = 0;
    for( let i=M; i-- > 0; )
    {
      let s = 0;
      for( let j=N; j-- > 0; )
        s += J[N*i+j] * (W[j] - X[j]);
      s += f[i];
      err_predict += s*s;
    }

//*DEBUG*/    if( !(err_predict <= err_now) ) throw new Error('Assertion failed.');

    let [f_next,
         J_next] = fJ( new NDArray(X_shape, W.slice()) );
    f_next = array('float64', f_next);
    J_next = array('float64', J_next);
    if( f_next.ndim !== 1 ) throw new Error('lsq_dogleg_gen(fJ, x0): fJ(x) returned [f,J] with f.ndim !== 1.');
    if( J_next.ndim !== 2 ) throw new Error('lsq_dogleg_gen(fJ, x0): fJ(x) returned [f,J] with J.ndim !== 2.');

    if( M !== f_next.shape[0] ) throw new Error('lsq_dogleg_gen(fJ, x0): fJ must return [f,J] with consistent shapes');
    if( M !== J_next.shape[0] ) throw new Error('lsq_dogleg_gen(fJ, x0): fJ must return [f,J] with consistent shapes');
    if( N !== J_next.shape[1] ) throw new Error('lsq_dogleg_gen(fJ, x0): fJ must return [f,J] with consistent shapes');

    const f_now = f,
          J_now = J;
    f = f_next.data;
    J = J_next.data;

    let err_next = 0;
    for( let i=M; i-- > 0; ) {
      const       s = f[i];
      err_next += s*s;
    }

    const rScale = () => 1;
    // const rScale = () => solver.scaledNorm(G); // <- TODO: examine scaling alternatives for the radius limits

    const predict = err_now - err_predict,
           actual = err_now - err_next;
    if(    actual < predict*expectGainMin )
    {
      // PREDICTION BAD -> REDUCE TRUST RADIUS
      // (by fitting a quadratic polynomial through err_now, grad_now and err_next)
      let grad_now = 0;
      for( let i=N; i-- > 0; )
        grad_now += dX[i] * G[i];
      grad_now *= 2;

      let   shrink = 0.5 * grad_now / (grad_now + err_now - err_next);
      if( !(shrink >= shrinkLower) ) shrink = shrinkLower;
      if( !(shrink <= shrinkUpper) ) shrink = shrinkUpper;
      const r = Math.max(rMin*rScale(), R*shrink);
      if(   r === R && !(err_now > err_next) )
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
          Math.min(rs*rMax, rGN*distanceToGN));
    }
    else if( actual > predict*expectGainMax || same_x )
    {
      // PREDICTION GREAT -> INCREASE TRUST RADIUS
      R = Math.min(rMax*rScale(),
          Math.max(nextUp(R), R*grow));
    }

    // ONLY ACCEPT NEW X IF BETTER
    if( err_now > err_next ) {
        err_now = err_next;
      [W,X] = [X,W];
      solver.update(f_next, J_next);
      stuckometer = 0;
    }
    else
    {
      if( ++stuckometer > stuckLimit )
        throw new OptimizationNoProgressError();
      // Rollback
      f = f_now;
      J = J_now;
    }
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
