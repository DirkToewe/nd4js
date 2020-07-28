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

import {array, NDArray} from '../nd_array'

import {FrobeniusNorm} from '../la/norm'

import {LBFGSB_Solver} from "./_lbfgsb_solver";
import {OptimizationNoProgressError} from './optimization_error';


const REPORT_STATE_READY    = 1,
      REPORT_STATE_NA       = 2,
      REPORT_STATE_CONSIDER = 3;


/** Computes the dot product of two arrays.
 */
const dot = (u,v) => {
  if( u.length%1 !== 0       ) throw new Error('Assertion failed.');
  if( v.length%1 !== 0       ) throw new Error('Assertion failed.');
  if( u.length   !== v.length) throw new Error('Assertion failed.');

  let uv = 0;
  for( let i=u.length; i-- > 0; )
    uv += u[i]*v[i];

  return uv;
};


export class TrustRegionSolverLBFGS
{
  constructor( fg, x0, {updateTol, historySize, scaleInit} )
  {
    if( !(fg instanceof Function) )
      throw new Error('new TrustRegionSolverLBFGS(fg, x0, updateTol): fg must be Function.');

    scaleInit *= 1;
    if( ! (0   <   scaleInit) ) throw new Error('TrustRegionSolverLBFGS(fg,x0,opt): opt.scaleInit must be greater than 0.');
    if( ! isFinite(scaleInit) ) throw new Error('TrustRegionSolverLBFGS(fg,x0,opt): opt.scaleInit must be greater than 0.');

    this.fg = fg;

    this.report_x = x0 = array('float64', x0);
    if( x0.ndim !== 1 )
      throw new Error('new TrustRegionSolverLBFGS(fg, x0, updateTol): x0.ndim must be 1.');

    let [f,g] = fg( new NDArray(x0.shape, x0.data.slice()) );

    this.report_loss      = f *= 1;
    this.report_loss_grad = g  = array('float64', g);

    if( isNaN(f)     ) throw new Error('new TrustRegionSolverLBFGS(fg, x0, updateTol): fg must have signature float[N] => [float,float[N]].');
    if( g.ndim !== 1 ) throw new Error('new TrustRegionSolverLBFGS(fg, x0, updateTol): fg must have signature float[N] => [float,float[N]].');
  
    const             [N] = x0.shape;
    if( g.shape[0] !== N ) throw new Error('new TrustRegionSolverLBFGS(fg, x0, updateTol): fg must have signature float[N] => [float,float[N]].');

    this.M = historySize;
    this.N = N;

    this.scaleInit = scaleInit = 1 / scaleInit;

    this.lbfgs = new LBFGSB_Solver(historySize, N);
    this.lbfgs.scale = scaleInit;

    this.D = new Float64Array(N); this.D.fill( Math.sqrt(scaleInit) );

    this.X0 = Float64Array.from(x0);
    this.G0 = Float64Array.from( g);

    // the last X,G that was used as LBFGS update
    this.last_X = Float64Array.from(x0);
    this.last_G = Float64Array.from( g);

    this.dg = new Float64Array(N);
    this.loss = f;
    this.bv = new Float64Array(2*historySize);

    this.newton_dX = new Float64Array(N);
    // this.regularized_dX = new Float64Array(N);

    this.updateTol = updateTol;
    this._report_state = REPORT_STATE_READY;

    Object.seal(this);
  }


  wiggle()
  {
    const {lbfgs, D, scaleInit} = this,
                            {m} = lbfgs;
    if( 0===m )
      throw new OptimizationNoProgressError('Too many unsuccessfull iterations.'); 

    if( !(0 < m) )
      throw new Error('Assertion failed.');

    lbfgs.forget(lbfgs.m+1 >>> 1);
    // lbfgs.scale = scaleInit;
    // D.fill( Math.sqrt(scaleInit) );
  }


  report()
  {
    if( this._report_state !== REPORT_STATE_READY )
      throw new Error('TrustRegionSolverLSQ::report: can only be called once after each makeConsideredMove() but not directly after considerMove(dX).');
    this._report_state = REPORT_STATE_NA;

    const result = [
      this.report_x,
      this.report_loss,
      this.report_loss_grad
    ];

    this.report_x         =
    this.report_loss_grad = null;
    this.report_loss      = NaN;

    return result;
  }


  considerMove( dX )
  {
    const {N, D, updateTol, lbfgs, loss, bv, X0, G0, dg, last_X, last_G, scaleSharpness} = this;

    this._report_state = REPORT_STATE_CONSIDER;

    if( N !== dX.length )
      throw new Error('Assertion failed.');

    const X_shape = Int32Array.of(N),
          X = new Float64Array(N);

    // EVALUATE FUNCTION AT CONSIDERED X

    // compute considered x
    for( let i=N; i-- > 0; )
      X[i] = X0[i] + dX[i];

    // evaluate function at considered x
    let [f,g] = this.fg(
      new NDArray( X_shape, X.slice() )
    );

    f *= 1;
    g  = array('float64', g);

    if( isNaN(f) ) throw new Error('Assertion failed.');
    if( g.ndim !== 1 ) throw new Error('Assertion failed.');
    if( g.shape[0] !== N ) throw new Error('Assertion failed.');

    this.report_x         = new NDArray(X_shape, X);
    this.report_loss      = f;
    this.report_loss_grad = g;

    g = g.data;

    // compute actual deltaX
    for( let i=N; i-- > 0; )
      X[i] -= X0[i];

    // COMPUTE PREDICTION
    let predict_loss = loss;
    for( let i=N; i-- > 0; )
      predict_loss += G0[i] * X[i];
    lbfgs.compute_bv(X, bv);
    predict_loss += 0.5 * lbfgs.compute_ubbv(bv, dot(X,X), bv);

    // UPDATE LBFGSB MODEL
    ;{
      const dx = X;

      for( let i=N; i-- > 0; ) dg[i] = g[i] - last_G[i];
      for( let i=N; i-- > 0; ) {
        dx[i]  = X0[i] + dX[i];
        dx[i] -= last_X[i];
      }

      const   mx = dg.reduce((mx,x) => Math.max(mx,Math.abs(x)), 0); // <- avoid underflow
      if( 0 < mx )
      {
        let gg = 0,
            gx = 0;

        for( let i=N; i-- > 0; )
        { const xi = dx[i]/mx,
                gi = dg[i]/mx;
          gg += gi*gi;
          gx += xi*gi;
        }

        if( 1 <= gg ) // <- checks NaN
        {
          if( ! isFinite(gx) )
            throw new Error('Assertion failed.');

          if( updateTol*gg < gx )
          {
            lbfgs.update(dx,dg);
            lbfgs.scale = gg/gx;
            // lbfgs.scale = Math.max(lbfgs.scale, gg/gx);

            // UPDATE SCALING
            for( let i=N; i-- > 0; )
            {                  lbfgs.compute_be(i,bv);
              const     B_ii = lbfgs.compute_ubbv(bv,1,bv);
              if( !(0 < B_ii) )
                throw new Error('Assertion failed: LBFGS model not positive definite.');
              D[i] = Math.max(D[i], Math.sqrt(B_ii));
            }

            for( let i=N; i-- > 0; ) last_X[i] = X0[i] + dX[i];
            for( let i=N; i-- > 0; ) last_G[i] =  g[i];
          }
        }
      }
    }

    // compute considered x (again)
    for( let i=N; i-- > 0; )
      X[i] = X0[i] + dX[i];

    return [ predict_loss,
         this.report_loss ];
  }


  makeConsideredMove()
  {
    if( this._report_state !== REPORT_STATE_CONSIDER )
      throw new Error('Assertion failed.');

    this._report_state = REPORT_STATE_READY;

    // discard consideration
    this.loss = this.report_loss;

    const {N, X0, G0} = this;

    const {data: x0} = this.report_x,
          {data: g0} = this.report_loss_grad;

    for( let i=N; i-- > 0; ) X0[i] = x0[i];
    for( let i=N; i-- > 0; ) G0[i] = g0[i];
  }


  scaledNorm( X )
  {
    const {N,D} = this;

    if( X.length !== N ) throw new Error('Assertion failed.');

    const norm = new FrobeniusNorm();

    for( let i=N; i-- > 0; )
      norm.include(D[i]*X[i]);

    return norm.result;
  }


  cauchyTravel()
  {
    const {N, lbfgs, G0, bv} = this;

    // polynomial along the gradient const1 + 2ax + bx² = const2 + (a/b + x)² * b
    let a=0;
    for( let i=N; i-- > 0; )
      a += G0[i]*G0[i];

    lbfgs.compute_bv(G0,bv);
    const b = lbfgs.compute_ubbv(bv, dot(G0,G0), bv);

    if( isNaN(b) )
      throw new Error('Assertion failed: ' + b);

    return 0===a ? 0 :
           0===b ? -Infinity : -a/b;
  }


  computeNewton()
  {
    const {N, lbfgs, G0, newton_dX} = this;

    lbfgs.compute_Hv(G0, newton_dX);

    if( newton_dX.some(isNaN) ) 
      throw new Error('Assertion failed.');

    for( let i=N; i-- > 0; )
      newton_dX[i] *= -1;
  }
}
