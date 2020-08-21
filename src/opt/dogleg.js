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

import {array, asarray, NDArray} from '../nd_array';

import {nextUp} from "../dt/float64_utils";

import {OptimizationNoProgressError} from "./optimization_error";
import {roots1d_polyquad} from './polyquad';
import {TrustRegionSolverLBFGS} from './_trust_region_solver_lbfgs';
import {TrustRegionSolverLSQ  } from './_trust_region_solver_lsq';
import {fit_odr_gen} from "./_trust_region_solver_odr_ny";


export function min_dogleg_gen(fg, x0, opt={})
{
  const { updateTol = 1e-14,
        historySize = 8,
          scaleInit = 1e-2 } = opt;
  return _dogleg(
    new TrustRegionSolverLBFGS(fg, x0, {updateTol, historySize, scaleInit}),
    opt
  );
}


// References
// ----------
// .. [1] "An Improved Optimization Method for iSAM2"
//         Rana Talaei Shahir and Hamid D. Taghirad Senior Member, IEEE
//         Proceeding of the 2nd RSI/ISM International Conference on Robotics and Mechatronics
//         October 15-17, 2014, Tehran, Iran
// .. [2] "THE LEVENBERG-MARQUARDT ALGORITHM: IMPLEMENTATION AND THEORY"
//         Jorge J. Moré.

// TODO implemented double-dogleg method as well as described in [1]
// TODO implement dogbox method as well which allows for box-constrained optimization

/** Solves a nonlinear least-squares problem using the dogleg method.
 */
export const lsq_dogleg_gen = (fJ, x0, opt) => _dogleg(new TrustRegionSolverLSQ(fJ,x0), opt);

export function* _dogleg(
  solver,
  /*opt=*/{
    r0   = 1.1,               // <-   START TRUST REGION RADIUS (meaning radius AFTER scaling)
    rMin = 0,                 // <- MINIMUM TRUST REGION RADIUS (meaning radius AFTER scaling)
    rMax = Infinity,          // <- MAXIMUM TRUST REGION RADIUS (meaning radius AFTER scaling)
    rNewton = 2,              // <- IF NEWTON-POINT INSIDE TRUST RADIUS AND IMPROVEMENT OKAY (according to expectGainMin), RADIUS IS SET TO rNewton*distance(x,xNewton)
    shrinkLower= 0.05,        // <-  SHRINK TRUST REGION RADIUS BY THIS FACTOR AT MOST
    shrinkUpper= 0.95,        // <-  SHRINK TRUST REGION RADIUS BY THIS FACTOR AT LEAST
    grow = 1.4142135623730951,// <-    GROW TRUST REGION RADIUS BY THIS FACTOR
    expectGainMin = 0.25,     // <- IF ACTUAL IMPROVEMENT RELATIVE TO PREDICTION IS WORSE  THAN THIS FACTOR, SHRINK TRUST REGION
    expectGainMax = 0.75,     // <- IF ACTUAL IMPROVEMENT RELATIVE TO PREDICTION IS BETTER THAN THIS FACTOR, GROW   TRUST REGION
    stuckLimit = 64           // <- MAX. NUMBER OF CONSECUTIVE ITERATIONS WITHOUT PROGRESS
  } = {}
)
{
  if( !(            0 <  r0          ) ) throw new Error('dogleg_gen(fJ, x0, opt): opt.r0 must be positive number.');
  if( !(            0 <= rMin        ) ) throw new Error('dogleg_gen(fJ, x0, opt): opt.rMin must be non-negative.');
  if( !(         rMin <= rMax        ) ) throw new Error('dogleg_gen(fJ, x0, opt): opt.rMin must not be greater than opt.rMax.');
  if( !(         rMin <= r0          ) ) throw new Error('dogleg_gen(fJ, x0, opt): opt.r0 must not be less than opt.rMin.');
  if( !(         rMax >= r0          ) ) throw new Error('dogleg_gen(fJ, x0, opt): opt.r0 must not be greater than opt.rMax.');
  if( !(            1 <  rNewton     ) ) throw new Error('dogleg_gen(fJ, x0, opt): opt.rGN must be greater than 1.');
  if( !(            0 <  shrinkLower ) ) throw new Error('dogleg_gen(fJ, x0, opt): opt.shrinkLower must be positive number.');
  if( !(  shrinkLower <= shrinkUpper ) ) throw new Error('dogleg_gen(fJ, x0, opt): opt.shrinkLower must not be greater than opt.shrinkUpper.');
  if( !(  shrinkUpper < 1            ) ) throw new Error('dogleg_gen(fJ, x0, opt): opt.shrinkUpper must be less than 1.');
  if( !(            1 < grow         ) ) throw new Error('dogleg_gen(fJ, x0, opt): opt.grow must be number greater than 1.');
  if( !(            0 < expectGainMin) ) throw new Error('dogleg_gen(fJ, x0, opt): opt.expectGainMin must be positive number.');
  if( !(expectGainMin < expectGainMax) ) throw new Error('dogleg_gen(fJ, x0, opt): opt.expectGainMin must be less than expectGainMax.');
  if( !(expectGainMax < 1            ) )    console.warn('dogleg_gen(fJ, x0, opt): opt.expectGainMax should be less than 1.');
  if(   0!== stuckLimit%1              ) throw new Error('dogleg_gen(fJ, x0, opt): opt.stuckLimit must be an integer.');
  if( !(0 <= stuckLimit              ) ) throw new Error('dogleg_gen(fJ, x0, opt): opt.stuckLimit must not be negative.');

  const {N, D, newton_dX, G0: G} = solver,
                      dX = new Float64Array(N);
  let     gNorm = solver.scaledNorm(G),
             cp = solver.cauchyTravel(),
              R = r0 * -gNorm*cp,
           loss = solver.loss,
    stuckometer = 0; // <- indicates how stuck we are, i.e. counts the number of consecutive operations without progress

  for(;;)
  {
/*DEBUG*/    if( !(cp <= 0) ) throw new Error('Assertion failed.');

    if( 0 === stuckometer )
      yield solver.report();

    { const t = Math.max(-R/gNorm, cp); // <- don't go beyond trust region radius
      for( let i=N; i-- > 0; )
        dX[i] = t*G[i]; // <- FIXME: Shouldn't the scaling be used here?!?!
    };

    let gnInRadius = false;

    if( -gNorm*cp < R ) // <- IF CAUCHY POINT INSIDE TRUST RADIUS
    {
      solver.computeNewton(); // <- computes Gauss-Newton point

      // TODO In [2], a method is described to adjust trust radius when global min is
      //      inside trust region. Maybe it could be used in the Dogleg method as well.

      // let's find the travel t from cauchy point (cp) to gauss-newton point (gp) which intersects
      // with the trust region boundary. t=0 means cauchy point t=1 mean gauss-newton point.
      // The travel is the solution to following the quadratic equation:
      //   R² = a + 2bt + ct²
      let a = 0,
          b = 0,
          c = 0;
      for( let i=N; i-- > 0; ) {
        const d =  D[i],
          u = d * dX[i],
          v = d * (newton_dX[i] - dX[i]);
        a += u*u;
        b += u*v;
        c += v*v;
      } a -= R*R;
        b *= 2;

      let t = roots1d_polyquad(a,b,c)[1];

/*DEBUG*/      if( !(0 <= t) ) throw new Error('Assertion failed.');

      if( t > 1 ) { // <- don't go beyond gauss-newton point (i.e. when it lies inside rust region)
          t = 1;
          gnInRadius = true;
      }

      for( let i=N; i-- > 0; ) {
        const    v = newton_dX[i] - dX[i];
        dX[i] += v*t;
      }
    }

/*DEBUG*/    if( !gnInRadius ) { if( !(Math.abs(solver.scaledNorm(dX) - R) <= R*   1e-5 + 1e-8 ) ) throw new Error( 'Assertion failed: ' + JSON.stringify({cp, err: Math.abs(solver.scaledNorm(dX) - R), D: [...D], G: [...G], dX: [...dX]}) ); } // <- check that solution is on ellipsoid boundary
/*DEBUG*/    else                if( !(         solver.scaledNorm(dX)      <= R*(1+1e-5) ) ) throw new Error('Assertion failed.');   // <- check that solution is in ellipsoid

    const [
      loss_predict,
      loss_consider
    ] = solver.considerMove(dX);

//*DEBUG*/    if( !(loss_predict <= loss+1e-12) ) throw new Error('Assertion failed: ' + JSON.stringify({loss, loss_consider, loss_predict, gnInRadius}));

    const rScale = () => 1;
    // const rScale = () => gNorm; // <- TODO: examine scaling alternatives

    const predict = loss - loss_predict,
           actual = loss - loss_consider;
    if(    actual < predict*expectGainMin )
    {
      // PREDICTION BAD -> REDUCE TRUST RADIUS
      // (by fitting a quadratic polynomial through err_now, grad_now and err_next)
      let grad = 0;
      for( let i=N; i-- > 0; )
        grad += dX[i] * G[i];

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
      // See [2] Algorithm (7.1) (d). If Gauss-Newton point is inside the trust radius and
      // prediction is not bad, set the radius to rGN times the distance to the Gauss-Newton
      // point. This also avoids that the trust radius grows uncontrollably while the
      // Gauss-Newton point is inside the trust region for a few iterations.
      const distanceToNewton = solver.scaledNorm(dX),
                        rs   = rScale();
      R = Math.max(rMin*rs,
          Math.min(rMax*rs, rNewton*distanceToNewton));
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
      solver.wiggle();

    gNorm = solver.scaledNorm(G);
    cp = solver.cauchyTravel();
  }
}


// export const odr_dogleg_gen = (fgg, p0, dx0, opt) => 
//   ???


export const fit_odr_dogleg_gen = fit_odr_gen(_dogleg);


export function* fit_dogleg_gen(
  x, //: float[nSamples,nDim]
  y, //: float[nSamples]
  fg,//: (params: float[nParam]) => (x: float[nDim]) => [f: float, g: float[nParam]],
  p0,//: float[nParam]
  opt
)
{
  if( ! (fg instanceof Function) )
    throw new Error('fit_dogleg_gen(x,y, fg, p0): fg must be a function.');

  x = array(           x ); // <- TODO allow non-ndarray x ?
  y = array('float64', y );
  p0= array('float64', p0);

  const FloatArray = Float64Array;

  if(x .ndim !== 2 ) throw new Error('fit_dogleg_gen(x,y, fg, p0): x.ndim must be 2.');
  if(y .ndim !== 1 ) throw new Error('fit_dogleg_gen(x,y, fg, p0): y.ndim must be 1.');
  if(p0.ndim !== 1 ) throw new Error('fit_dogleg_gen(x,y, fg, p0): p0.ndim must be 1.');

  const [P] = p0.shape,
      [M,N] =  x.shape,
    x_shape = Int32Array.of(N);

  if( M != y.shape[0] )
    throw new Error('fit_dogleg_gen(x,y, fg, p0): x.shape[0] must be equal to y.shape[0].');

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

  const R  = new FloatArray(M  ),
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

      if( f.ndim     !== 0 ) throw new Error('fit_dogleg_gen(x,y, fg, p0): fg must have signature float[nParams] => float[nDim] => [float, float[nParams]].');
      if( g.ndim     !== 1 ) throw new Error('fit_dogleg_gen(x,y, fg, p0): fg must have signature float[nParams] => float[nDim] => [float, float[nParams]].');
      if( g.shape[0] !== P ) throw new Error('fit_dogleg_gen(x,y, fg, p0): fg must have signature float[nParams] => float[nDim] => [float, float[nParams]].');

      f = f.data[0];
      g = g.data;

      if(       isNaN(f)) throw new Error('fit_dogleg_gen(x,y, fg, p0): NaN encountered.');
      if(g.some(isNaN)  ) throw new Error('fit_dogleg_gen(x,y, fg, p0): NaN encountered.');

      R[i] = f - y[i];

      for( let j=P; j-- > 0; )
        J[P*i+j] = g[j];
    }
    return RJ;
  };

  yield* lsq_dogleg_gen(fJ, p0, opt);
}
