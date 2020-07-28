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

import {_giv_rot_qr,
        _giv_rot_rows} from "../la/_giv_rot";
import {FrobeniusNorm} from "../la/norm";
import {_rrqr_decomp_inplace,
        _rrqr_rank} from "../la/rrqr";
import {_triu_solve} from "../la/tri";

import {OptimizationNoProgressError} from "./optimization_error";


const REPORT_STATE_READY    = 1,
      REPORT_STATE_NA       = 2,
      REPORT_STATE_CONSIDER = 3;


// This class is used by DogLeg (and possibly LevenbergMarquardt) to solve the following minimization problem:
//
// minimize ∑ {f(p, X[i,:] + Δ[i,:]) - y[i]}²  +  ∑ (Δ[i,j])²
//          i                                   (i,j)
//
//     Param              | Description
//  ----------------------|-------------
//  f: float[NX] -> float | The curve to be fit
//  X: float[MX,NX]       | X-values of the curve sampling points
//  y: float[MX]          | y-values of the curve sampling points
//  p: float[NP]          | The curve parameters
//  Δ: float[MX,NX]       | X-correction summand for the sampling point X (is called dx in the code below)
//
// This minimization problem is also known as orthogonal distance regression (ODR). In other words ODR takes
// noisy sampling points (X,y) as input and fits a curve function `f` through those sampling point by finding
// the parameters `Δ` and `p`, such that `f` has minimal orthogonal distance to the sampling points. `p` are the
// actual curve parameters of `f`, while `Δ` are summands accounting for the noise in `X`.


export class TrustRegionSolverODR
{
  constructor( samples_x, samples_y, fgg, p0, dx0 )
  {
    if( ! (fgg instanceof Function) )
      throw new Error('Assertion failed.');

    samples_x = array('float64', samples_x);
    samples_y = array('float64', samples_y);
           p0 = array('float64',  p0);
          dx0 = array('float64', dx0);

    if(samples_x.ndim !== 1 &&
       samples_x.ndim !== 2      ) throw new Error('Assertion failed.');
    if(samples_y.ndim !== 1      ) throw new Error('Assertion failed.');
    if(samples_x.ndim !==dx0.ndim) throw new Error('Assertion failed.');
    if(             1 !== p0.ndim) throw new Error('Assertion failed.');
    

    const [MX,NX=1] = dx0.shape,
             [NP  ] =  p0.shape;

    if(   1 !==       dx0.ndim
      && NX !==       dx0.shape[1] ) throw new Error('Assertion failed.');
    if(  MX !==       dx0.shape[0] ) throw new Error('Assertion failed.');
    if(  MX !== samples_y.shape[0] ) throw new Error('Assertion failed.');

    const M = MX*NX + MX,
          N = MX*NX + NP,
        J11 = new Float64Array(MX*NX),
        J21 = new Float64Array(MX*NX),
        J22 = new Float64Array(MX*NP),
          D = new Float64Array(N),
      _consider_J11 = new Float64Array(J11.length),
      _consider_J21 = new Float64Array(J21.length),
      _consider_J22 = new Float64Array(J22.length),

      R11 = new Float64Array( MX * (NX*(NX+1)>>>1) ),
      R22 = new Float64Array( MX * NP ),

      QF = new Float64Array(M),

      X0 = new Float64Array(N),
      F0 = new Float64Array(M),
      G0 = new Float64Array(N),

      newton_dX = new Float64Array(N),

      P = new Int32Array(NP),

      tmp = new Float64Array(MX*NX),
      norm= new Float64Array(2*NP);

    Object.assign(this, {
      MX,NX, NP, M,N,
      loss: 0.0,
      rank: -1,
      _report_state: REPORT_STATE_NA,
      fgg,
      report_p        : null,
      report_dx       : null,
      report_loss     : NaN,
      report_dloss_dp : null,
      report_dloss_ddx: null,
      report_dy       : null,
      samples_x: samples_x.data,
      samples_y: samples_y.data,
      x_shape: samples_x.shape,
      y_shape: samples_y.shape,
      p_shape:        p0.shape,
      J11,     R11,
      J21,J22, R22,
      _consider_J11,
      _consider_J21,
      _consider_J22,
      D, tmp, norm,
      X0,F0,G0, QF, P,
      newton_dX
    });
    Object.seal(this);

           p0 =        p0.data;
          dx0 =       dx0.data;
    samples_x = samples_x.data;

    for( let i=NP; i-- > 0; )
      newton_dX[MX*NX + i] = p0[i];

    for( let j=NX; j-- > 0; )
    for( let i=MX; i-- > 0; )
      newton_dX[MX*j+i] = dx0[NX*i+j];

    this._considerMove(newton_dX);
                       newton_dX.fill(0.0);
    this.makeConsideredMove();
  }


  _considerMove( dX )
  {
    const {
      M,N, MX,NX,NP,
      fgg,
      x_shape, samples_x,
      y_shape, samples_y,
      p_shape,
      _consider_J11,
      _consider_J21,
      _consider_J22,
      X0
    } = this;

    if( dX.length !== N  ) throw new Error('Assertion failed.');

    this._report_state = REPORT_STATE_CONSIDER;

    const report_p        = new Float64Array(NP),
          report_dx       = new Float64Array(MX*NX),
          report_dy       = new Float64Array(MX),
          report_dloss_dp = new Float64Array(NP),
          report_dloss_ddx= new Float64Array(MX*NX);
     this.report_p        = new NDArray( p_shape, report_p );
     this.report_dx       = new NDArray( x_shape, report_dx);
     this.report_dy       = new NDArray( y_shape, report_dy);
     this.report_dloss_dp = new NDArray( p_shape, report_dloss_dp );
     this.report_dloss_ddx= new NDArray( x_shape, report_dloss_ddx);

    for( let i=NP; i-- > 0; ) {
      const                    I = MX*NX + i;
      report_p[i] = X0[I] + dX[I];
    }

    for( let i=MX; i-- > 0; )
    for( let j=NX; j-- > 0; ) {
      const                           ji = MX*j+i;
      report_dx[NX*i+j] = X0[ji] + dX[ji];
    }

    const fgg_p = fgg(
      new NDArray(p_shape, report_p.slice())
    );

    for( let i=MX; i-- > 0; )
    {
      let [f,gp,gx] = function(){
        const xi = new Float64Array(NX);

        for( let j=NX; j-- > 0; ) { const   ij = NX*i+j;
          xi[j] = samples_x[ij] + report_dx[ij];
        }

        return fgg_p(
          new NDArray(x_shape.slice(1), xi)
        );
      }();

      f *= 1;
      gp = asarray('float64',gp);
      gx = asarray('float64',gx);

      if( isNaN(f)                   ) throw new Error('Assertion failed.');
      if(gp.ndim !== 1               ) throw new Error('Assertion failed.');
      if(gx.ndim !== x_shape.length-1) throw new Error('Assertion failed.');

      report_dy[i] = f - samples_y[i];

      gp = gp.data;
      gx = gx.data;

      for( let j=NP; j-- > 0; ) _consider_J22[NP*i+j] = gp[j];
      for( let j=NX; j-- > 0; ) _consider_J21[MX*j+i] = gx[j]; // <- TODO: consider "transpose" memory organization of J21, so that it's same memory organization as R21 in `computeNewton`
    }

    _consider_J11.fill(1.0); // <- TODO: use weights instead;

    // COMPUTE LOSS GRADIENT w.r.t. P
    for( let i=MX; i-- > 0; )
    for( let j=NP; j-- > 0; )
      report_dloss_dp[j] += _consider_J22[NP*i+j] * report_dy[i] / M * 2;

    // COMPUTE LOSS GRADIENT w.r.t. ΔX
    for( let i=MX; i-- > 0; )
    for( let j=NX; j-- > 0; )
      report_dloss_ddx[NX*i+j] = (
        _consider_J11[MX*j+i] * report_dx[NX*i+j] +
        _consider_J21[MX*j+i] * report_dy[   i  ]
      ) / M * 2;


    // COMPUTE LOSS (mean squared error)
    let report_loss = 0.0;
    for( let i=MX*NX; i-- > 0; ) {
      const          s = report_dx[i];
      report_loss += s*s / M;
    }
    for( let i=MX; i-- > 0; ) {
      const          s = report_dy[i];
      report_loss += s*s / M;
    }

    this.report_dloss_dp = new NDArray(p_shape, report_dloss_dp );
    this.report_dloss_ddx= new NDArray(x_shape, report_dloss_ddx);
    this.report_loss     = report_loss;
  }


  considerMove( dX )
  {
    this._considerMove(dX);

    const {
      M,N, MX,NX,NP,

      J11,
      J21,J22,

      F0,X0
    } = this;

    const report_p  = this.report_p .data,
          report_dx = this.report_dx.data;

    let predict_loss = 0.0;
    for( let i=MX; i-- > 0; )
    for( let j=NX; j-- > 0; ) {
      const ji = MX*j+i,
            ij = NX*i+j,
                      f = F0[ji] + J11[ji] * (report_dx[ij] - X0[ji]);
      predict_loss += f*f/M;
    }

    for( let i=MX; i-- > 0; )
    {
      let f = F0[MX*NX + i];
      for( let j=NX; j-- > 0; ) f += J21[MX*j+i] * (report_dx[NX*i+j] - X0[MX*j+i]);
      for( let j=NP; j-- > 0; ) f += J22[NP*i+j] * (report_p[j]       - X0[MX*NX + j]);
      predict_loss += f*f/M;
    }

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
    this.rank = -1;

    [this._consider_J11, this.J11] = [this.J11, this._consider_J11];
    [this._consider_J21, this.J21] = [this.J21, this._consider_J21];
    [this._consider_J22, this.J22] = [this.J22, this._consider_J22];

    const {
      MX,NX,NP,
      J11,
      J21,J22,
      D,
      X0,F0,G0
    } = this;
  
    ;{
      const p = this.report_p .data,
           dx = this.report_dx.data,
           dy = this.report_dy.data;

      for( let i=NP   ; i-- > 0; ) X0[MX*NX + i] = p[i];
      for( let i=MX   ; i-- > 0; )
      for( let j=NX   ; j-- > 0; ) X0[MX*j+i] = dx[NX*i+j];

      for( let i=MX   ; i-- > 0; ) F0[MX*NX + i] = dy[i];
      for( let i=MX*NX; i-- > 0; ) F0[        i] = X0[i];
    }

    // COMPUTE GRADIENT OF (HALF) SQUARED ERROR (MIND THE HALF :P)
    for( let i=MX*NX; i-- > 0; ) G0[i] = J11[i] * F0[i];

    for( let j=NX; j-- > 0; )
    for( let i=MX; i-- > 0; )
      G0[MX*j+i] += J21[MX*j+i] * F0[MX*NX + i];

    G0.fill(0.0, MX*NX, MX*NX+NP);
    for( let i=MX; i-- > 0; )
    for( let j=NP; j-- > 0; )
      G0[MX*NX + j] += J22[NP*i+j] * F0[MX*NX + i];

    for( let i=MX*NX; i-- > 0; )
      D[i] = Math.hypot(J11[i],J21[i]);

    const norm = new FrobeniusNorm();
    for( let j=NP; j-- > 0; ) { norm.reset();
      for( let i=MX; i-- > 0; ) norm.include(J22[NP*i+j]);
      D[MX*NX + j] =            norm.result;
    }
  }


  cauchyTravel()
  {
    const {
      N,MX,NX,NP,
      J11,
      J21,J22,
      G0: G
    } = this;

    // polynomial along the gradient const1 + 2ax + bx² = const2 + (a/b + x)² * b
    let a=0,
        b=0;

    for( let i=N; i-- > 0; )
      a += G[i]*G[i];

    for( let i=MX; i-- > 0; ) {
      let Jg = 0;

      for( let j=NP; j-- > 0; )
        Jg += J22[NP*i+j] * G[MX*NX + j];

      for( let j=NX; j-- > 0; ) {
        const     ji = MX*j+i;
        Jg += J21[ji] * G[ji];
      }

      b += Jg * Jg;
    }

    for( let i=MX*NX; i-- > 0; ) {
      const     Jg = J11[i] * G[i];
      b += Jg * Jg;
    }

    return 0===a ? 0 :
           0===b ? -Infinity : -a/b;
  }


  report()
  {
    if( this._report_state !== REPORT_STATE_READY )
      throw new Error('TrustRegionSolverLSQ::report: can only be called once after each makeConsideredMove() but not directly after considerMove(dX).');
    this._report_state = REPORT_STATE_NA;

    const result = [
      this.report_p,
      this.report_dx,
      this.report_loss,
      this.report_dloss_dp,
      this.report_dloss_ddx,
      this.report_dy
    ];

    this.report_p        =
    this.report_dx       =
    this.report_dy       =
    this.report_dloss_dp =
    this.report_dloss_ddx= null;
    this.report_loss     = NaN;

    return result;
  }


  wiggle()
  {
    throw new OptimizationNoProgressError('Too many unsuccessfull iterations.');
  }


  __DEBUG_J( i, j ) // <- meant for debugging only!
  {
    if( i%1 !== 0 ) throw new Error('Assertion failed.');
    if( j%1 !== 0 ) throw new Error('Assertion failed.');

    i |= 0;
    j |= 0;

    if( !( 0 <= i ) ) throw new Error('Assertion failed.');
    if( !( 0 <= j ) ) throw new Error('Assertion failed.');

    const {
      M,N, MX,NX,NP,
      J11,
      J21,J22
    } = this;

    if( !( i < M ) ) throw new Error('Assertion failed.');
    if( !( j < N ) ) throw new Error('Assertion failed.');

    if( i < MX*NX ) return i===j    ? J11[i] : 0;
        i-= MX*NX;
    if( j < MX*NX ) return i===j%MX ? J21[j] : 0;
        j-= MX*NX;

    return J22[NP*i+j];
  }


  scaledNorm( X )
  {
    const {N,D} = this;

    if( X.length !== N ) throw new Error('Assertion failed.');

    const                    norm = new FrobeniusNorm();
    for( let i=N; i-- > 0; ) norm.include(D[i]*X[i]);
    return                   norm.result;
  }


  __DEBUG_R( i, j ) // <- meant for debugging only!
  {
    if( i%1 !== 0 ) throw new Error('Assertion failed.');
    if( j%1 !== 0 ) throw new Error('Assertion failed.');

    i |= 0;
    j |= 0;

    if( !( 0 <= i ) ) throw new Error('Assertion failed.');
    if( !( 0 <= j ) ) throw new Error('Assertion failed.');

    const {
      M,N, MX,NX,NP,
           R11,
      tmp: R21,R22
    } = this;

    if( !( i < M ) ) throw new Error('Assertion failed.');
    if( !( j < N ) ) throw new Error('Assertion failed.');

    if( j < MX*NX )
    {
      if( i%MX !== j%MX ) return 0;

      if( i < MX*NX )
      {
        if( j < i ) return 0;

        const I = Math.floor(i/MX);

        const off = MX*( NX*I - (I*(I-1) >>> 1) );

        j = Math.floor(j/MX) - I;
        i = i%MX;

        return R11[off + (NX-I)*i+j];
      }

      i -= MX*NX;
      j = Math.floor(j/MX);
      return R21[NX*i+j];
    }

    if( i < MX*NX )
      return NaN;

    i -= MX*NX;
    j -= MX*NX;

    return R22[NP*i+j];
  }


  // The Jacobian of the orthogonal least squares problem is sparse with the following structure:
  //
  //     ┏                 ╷    ┓   ┏                 ╷    ┓   
  //     ┃  ╲              ┊    ┃   ┃                 ┊    ┃
  //     ┃    ╲            ┊    ┃   ┃                 ┊    ┃
  //     ┃      ╲          ┊    ┃   ┃                 ┊    ┃
  //     ┃        ╲        ┊ 0  ┃   ┃       J11       ┊    ┃
  // J = ┃          ╲      ┊    ┃ = ┃                 ┊    ┃
  //     ┃            ╲    ┊    ┃   ┃                 ┊    ┃
  //     ┃              ╲  ┊    ┃   ┃                 ┊    ┃
  //     ┃┄┄┄┄┄┬┄┄┄┄┄┬┄┄┄┄┄┼┄┄┄┄┃   ┃┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┃
  //     ┃ ╲   ┊ ... ┊ ╲   ┊ ██ ┃   ┃                 ┊    ┃
  //     ┃   ╲ ┊     ┊   ╲ ┊ ██ ┃   ┃       J21       ┊ J22┃
  //     ┗     ╵     ╵     ╵    ┛   ┗                 ╵    ┛
  //
  //  J  : float[MX*(NX+1),MX*NX+NP];
  //  J11: float[MX*NX    ,MX*NX   ]; J11 is a Diagonal Matrix; Diag. represent the weights on Δ; Assumed to NOT be rank-deficient;
  //  J21: float[MX       ,MX*NX   ]; J21 is a row of diagonal matrix blocks
  //  J22: float[MX       ,        ]; J22 is a dense matrix
  //
  // If there are no rank deficiencies, J can be sparsely QR decomposed to solve ODR problem:
  // 
  //       ┏     ╷     ╷     ╷    ┓     ┏                 ╷    ┓  
  //       ┃ ╲   ┊ ... ┊ ╲   ┊ ██ ┃     ┃                 ┊    ┃  
  //       ┃   ╲ ┊     ┊   ╲ ┊ ██ ┃     ┃                 ┊    ┃  
  //       ┃┄┄┄┄┄┼┄┄┄┄┄┼┄┄┄┄┄┤ ██ ┃     ┃                 ┊    ┃  
  //       ┃     ┊ ⋱   ┊  ⋮  ┊ ██ ┃     ┃       R1        ┊    ┃  
  // J = Q·┃     ┊   ⋱ ┊  ⋮  ┊ ██ ┃ = Q·┃                 ┊ R2 ┃  
  //       ┃     └┄┄┄┄┄┼┄┄┄┄┄┤ ██ ┃     ┃                 ┊    ┃  
  //       ┃           ┊ ╲   ┊ ██ ┃     ┃                 ┊    ┃  
  //       ┃           ┊   ╲ ┊ ██ ┃     ┃                 ┊    ┃  
  //       ┃     0     └┄┄┄┄┄┤ ██ ┃     ┃┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┤    ┃  
  //       ┃                 ┊ ▜█ ┃     ┃        0        ┊    ┃  
  //       ┃                 ┊  ▜ ┃     ┃                 ┊    ┃  
  //       ┗                 ╵    ┛     ┗                 ╵    ┛  
  //
  // R1: float[MX*NX    ,MX*NX]; Upper (block) triangular matrix consisting of diagonal matrix blocks.
  // R2: float[MX*(NX+1),   NP]; Dense Matrix. R2[-MX:, :] is upper triangular.
  //
  computeNewton()
  {
    const {
      M, MX,NX,NP,
      J11,
      J21,J22,
           R11,
      tmp: R21,R22,
      F0, QF,
      newton_dX: X,
      P, norm
    } = this;

    if( !(R21.length >= MX*NX) ) throw new Error('Assertion failed.');

    if( ! (0 > this.rank) )
      throw new Error('Assertion failed.');

    //
    // INIT R
    //
    // In J11, memory order is like (example with MX=3,NX=2):
    //   0
    //     4
    //       8
    //        12
    //          15
    //            18
    //              21
    //                23
    //                  25
    //                    27
    //                      28
    //                        29
    //
    // In R1, memory order is like:
    //
    //   0     1     2     3
    //     4     5     6     7
    //       8     9    10    11
    //        12    13    14
    //          15    16    17
    //            18    19    20
    //              21    22
    //                23    24
    //                  25    26
    //                    27
    //                      28
    //                        29
    //
    // Where the new indices are initialized to zero.

    //
    // STEP 1: MEMORY INITIALIZATION
    //
    R11.fill(0.0);
    for( let j=NX; j-- > 0; )
    {
      const off = MX*( NX*j - (j*(j-1) >>> 1) );

      for( let i=MX; i-- > 0; )
        R11[off + (NX-j)*i] = J11[MX*j+i];
    }

    // J21 memory is orderd left to right, basically going down the diagonals one after another
    // R21 on the other hand is ordered row by row top to bottom
    //
    // In other words: Memory order in J21 is like (example with MX=3,NX=4):
    //   0     3     6     9
    //     1     4     7    10
    //       2     5     8    11
    //
    // And memory order in R21 is like:
    //   0     1     2     3
    //     4     5     6     7    
    //       8     9    10    11
    for( let i=MX; i-- > 0; )
    for( let j=NX; j-- > 0; )
      R21[NX*i+j] = J21[MX*j+i];

    for( let i=MX*NP; i-- > 0; )
      R22[i] = J22[i];

    for( let i=M; i-- > 0; )
      QF[i] = F0[i];

    //
    // STEP 2.1: ELIMINATE R21 USING GIVENS ROTATIONS
    //           O( MX*NX² ) operations
    //
    X.fill(1.0, 0,MX); // <- temp. stores the (accumulated) cosines of the givens rotations

    for( let I=0; I < NX; I++ ) // <- for each block top to bottom
    {
      const R11_off = MX*( NX*I - (I*(I-1) >>> 1) );

      for( let i=0; i < MX; i++ ) // <- for each diagonal entry in block
      {
        // COMPUTE GIVENS ROTATION
        const    i1 = R11_off + (NX-I)* i,
                 i2 =            NX   * i+I,
                 r1 = R11[i1],
                 r2 = R21[i2],
         [c,s,norm] =_giv_rot_qr(r1,r2);

        R21[i2] = s * X[i];
                      X[i] *= c; if( s === 0 ) continue;
        R11[i1] = norm;

        // ROTATE R1 & R21
        for( let j=0; ++j < NX-I; ) // <- for earch non-zero entry in row
        { const         rj = R21[i2+j];
          R11[i1+j] = s*rj;
          R21[i2+j] = c*rj;
        }

        // ROTATE QF
        _giv_rot_rows(QF,1, MX*I +i,
                            MX*NX+i, c,s);
      }
    }

    //
    // STEP 2.2: APPLY TO R22 THE SCALE WHICH RESULTS FROM THE GIVENS ROTATIONS IN STEP 2.1
    //           O( MX*NP ) operations
    //
    for( let i=MX; i-- > 0; )
    for( let j=NP; j-- > 0; )
      R22[NP*i+j] *= X[i];

    //
    // STEP 2.3: RRQR-DECOMPOSE R22
    //           O( min(MX,NP)*NP*MX ) operations
    //
    for( let i=NP; i-- > 0; ) P[i] = i;

     _rrqr_decomp_inplace(MX,NP,1, R22,0, QF,MX*NX, P,0, norm);
    const rnk =_rrqr_rank(MX,NP,   R22,0,                norm);
    this.rank = rnk + MX*NX;

    for( let i=MX*NX+rnk; i-- > 0; )
      X[i] = QF[i];

    if( rnk !== NP )
    { //
      // STEP 2.4: ELIMINATE RANK-DEFICIENT COLUMS
      //           O( rnk² * (NP-rnk) ) operations
      X.fill(0.0, MX*NX+rnk);

      // eliminate lower part of linear dependent columns of R22
      for( let i=rnk; i-- > 0; ) {     const ii = NP*i+i;
        for( let j= NP; j-- > rnk; ) { const ij = NP*i+j,
                                  R_ii = R22[ii],
                                  R_ij = R22[ij];
          // compute Givens rot.
          let [c,s,nrm] = _giv_rot_qr(R_ii,R_ij);
          if( s !== 0 )
          { if( c < 0 ) {
                c *= -1;
                s *= -1;
              nrm *= -1;
            }
            // apply Givens rot.
            for( let k=i; k-- > 0; )
            { const         R_ki = R22[NP*k+i],
                            R_kj = R22[NP*k+j];
              R22[NP*k+i] = R_kj*s + c*R_ki;
              R22[NP*k+j] = R_kj*c - s*R_ki;
            }
            R22[ii] = nrm;
          } R22[ij] = s;
        }
      }
    }

    //
    // STEP 3.1: BACKWARDS SUBSITUTION OF R2-PART
    //           O( rnk² ) operations
    //
    _triu_solve(rnk,NP,1, R22,0, X,MX*NX);

    if( rnk != NP )
    { //
      // STEP 3.2: APPLY GIVENS ROTATIONS TO X (UNDO 2.4)
      //           O( rnk * (NP-rnk) ) operations
      //
      for( let i= 0 ; i < rnk; i++ )
      for( let j=rnk; j < NP ; j++ )
      {
        const s = R22[NP*i+j]; if( 0 === s ) continue;
        const c = Math.sqrt(1-s*s);
        _giv_rot_rows(X,1, MX*NX + j,
                           MX*NX + i, c,s);
      }
    }

    //
    // STEP 3.3: APPLY COLUMN PERMUTATIONS P TO X (UNDO 2.3)
    //           O( NP ) operations
    //
    for( let i=NP; i-- > 0; )                norm[P[i]] = X[MX*NX + i];
    for( let i=NP; i-- > 0; ) X[MX*NX + i] = norm[  i ];

    //
    // STEP 3.4: MOVE SOLVED R22-PART TO THE RIGHT
    //           O( max(NP,NX)*MX ) operations
    // -------------------------------------------
    // As a result of steps 2.1 and 2.2, each row in R12 can be described as
    // a scaled row of J22, i.e.:
    //
    // R12[MX*j+i,:] = s[i,j] * J22[i,:]
    //
    // Where:
    //
    // s[i,j] = R21[i,j] * c[i,j-1]
    // s[i,0] = R21[i,0]
    // c[i,j] = sqrt( 1 - (s[i,j])² )
    //
    for( let i=MX; i-- > 0; )
    {
      let xj = 0.0;
      for( let j=NP; j-- > 0; )
        xj += J22[NP*i+j] * X[MX*NX + j];

      for( let j=0; j < NX; j++ )
      { const        s = R21[NX*i+j];
        X[MX*j+i] -= s*xj;
      }
    }

    //
    // STEP 3.5: BACKWARD SUBSTITUTION OF R11
    //           O( MX*NX² ) operations
    for( let I=NX; I-- > 0; ) // <- for each block bottom to top
    {
      const R11_off = MX*( NX*I - (I*(I-1) >>> 1) );

      for( let i=MX; i-- > 0; ) // <- for each diagonal entry in block bottom to top
      {
        const iX = MX*I+i,
              iR = R11_off + (NX-I)*i;

        for( let j=NX-I; --j > 0; )
          X[iX] -= X[iX + MX*j] * R11[iR + j];

        X[iX] /= R11[iR];
      }
    }
  }


  // computeNewtonRegularized( λ )
  // {
  //   throw new Error('Not yet implemented.');
  // }
}
