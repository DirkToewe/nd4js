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

import {lstsq  } from '../la/lstsq'
import {matmul2} from '../la/matmul'
import {solve  } from '../la/solve'


// Regularized least squares (RLS) solver. Solves the follwing equation system:
//
//   (JᵀJ + λI)x = Jᵀy
//
const RLS = (J,y, λ) => {
  if( J.ndim !== 2 ) throw new Error('Assertion failed.');
  if( y.ndim !== 2 ) throw new Error('Assertion failed.');
  if( y.shape[0] !== J.shape[0] ) throw new Error('Assertion failed.');
  if( y.shape[1] !== 1          ) throw new Error('Assertion failed.');

  if( ! (0 <= λ && λ < Infinity) )
    throw new Error(`fit_lm_gen(x,y, fg, p0): Illegal damping factor λ=${λ} during iteration.`);

  if( 0 === λ )
    return lstsq(J,y);

  const [M,N] = J.shape;

  J = J.data;
  y = y.data;

  // The equation above can be solved as the following least squares problem:
  //
  //   min ║Ax - z║₂
  //    x
  //
  // where
  //       ┏       ┓        ┏   ┓
  //       ┃       ┃        ┃   ┃
  //       ┃       ┃        ┃   ┃
  //       ┃   J   ┃        ┃ y ┃
  //       ┃       ┃        ┃   ┃
  //       ┃       ┃        ┃   ┃
  //   A = ┃┄┄┄┄┄┄┄┃    z = ┃┄┄┄┃
  //       ┃√λ     ┃        ┃ 0 ┃
  //       ┃   ⋱   ┃        ┃ ⋮ ┃
  //       ┃    √λ ┃        ┃ 0 ┃
  //       ┗       ┛        ┗   ┛
  let A = new Float64Array((M+N)*N),
      z = new Float64Array( M+N );

  λ = Math.sqrt(λ);

  // Init A
  for( let i=N; i-- > 0; )
    A[M*N + N*i+i] = λ;

  for( let i=M*N; i-- > 0; )
    A[i] = J[i];

  // Init z
  for( let i=M; i-- > 0; )
    z[i] = y[i];

  A = new NDArray(Int32Array.of(M+N,N), A);
  z = new NDArray(Int32Array.of(M+N,1), z);

  return lstsq(A,z);
};


/** Solves a nonlinear least-squares problem using the
 *  Levenberg-Marquard algorithm.
 */
export function* lsq_lm_gen(
  fj,//: (x0: float[nInputs]) => [f: float[nSamples], j: float[nSamples,nInputs]]
  x0,//: float[nInputs]
  {
    lambda0 = 10,
    lambdaFactor = [1 / (Math.E-1), Math.PI-2],
    improvement: {
      accept    = 0,
      expectLow = 0.25,
      expectHigh= 0.75
    } = {},
    rls = RLS
  } = {}
)
{
  const FloatArray = Float64Array;

  let x = asarray(x0, 'float64');
                  x0 = undefined;

  if( x.ndim !== 1 ) throw new Error('lsq_lm_gen(fj, x0): x0.ndim must be 1.');

  if( accept     < 0         )    console.warn('lsq_lm_gen(fj, x0, opt): opt.improvement.accept should not be negative.');
  if( accept     > expectLow ) throw new Error('lsq_lm_gen(fj, x0, opt): opt.improvement.accept must be less than opt.improvement.expectLow.');
  if( expectHigh < expectLow ) throw new Error('lsq_lm_gen(fj, x0, opt): opt.improvement.expectLow must be less than opt.improvement.expectHigh.');
  if( expectHigh > 1         )    console.warn('lsq_lm_gen(fj, x0, opt): opt.improvement.expectHigh should not be greater than 1.');

  if( ! (0 <= lambda0 && lambda0 < Infinity) )
    throw new Error(`lsq_lm_gen(fj, x0, opt): Invalid opt.lambda0=${lambda0}.`);

  if( lambdaFactor instanceof Array )
  {
    if(  lambdaFactor.length !== 2  ) throw new Error(`lsq_lm_gen(fj, x0, opt): opt.lambdaFactor must contain two numbers if Array.`);
    const [lo,hi] = lambdaFactor;
    if( ! (0 < lo && lo < 1       ) ) throw new Error(`lsq_lm_gen(fj, x0, opt): opt.lambdaFactor[0] must be in (0,1).`);
    if( ! (1 < hi && hi < Infinity) ) throw new Error(`lsq_lm_gen(fj, x0, opt): opt.lambdaFactor[1] must be in (1,∞)`);
  }
  else if( ! (1 <= lambdaFactor && lambdaFactor < Infinity) )
    throw new Error(`lsq_lm_gen(fj, x0, opt): Invalid opt.lambdaFactor=${lambda0}.`);

  const [N] = x.shape,
    x_shape = x.shape;

  let     M = -1,
    r_shape = undefined,
    J_shape = undefined,
    R_shape = undefined;

  x = x.data;

  let R = undefined,
      G = undefined,
      J = undefined,
    mse = NaN;

  const recompute = () =>
  {
    let [r,j] = fj( new NDArray(x_shape, x.slice()) );
    r = asarray(r);
    j = asarray(j);

    if(r.ndim     !== 1         ) throw new Error('lsq_lm_gen(fj, x0): fj must have signature float[nInputs] => [float[nSamples], float[nSamples,nInputs]].');
    if(j.ndim     !== 2         ) throw new Error('lsq_lm_gen(fj, x0): fj must have signature float[nInputs] => [float[nSamples], float[nSamples,nInputs]].');
    if(j.shape[0] !== r.shape[0]) throw new Error('lsq_lm_gen(fj, x0): fj must have signature float[nInputs] => [float[nSamples], float[nSamples,nInputs]].');
    if(j.shape[1] !== N         ) throw new Error('lsq_lm_gen(fj, x0): fj must have signature float[nInputs] => [float[nSamples], float[nSamples,nInputs]].');

    if( M < 0 ) { // <- 1st time around we have to determine M
      [M] = r.shape;
      r_shape = r.shape;
      J_shape = j.shape;
      R_shape = Int32Array.of(M,1);
      R = new FloatArray(M  );
      G = new FloatArray(  N);
      J = new FloatArray(M*N);
    }
    else if( j.shape[0] !== M )
      throw new Error('lsq_lm_gen(fj, x0): fj must be consistent.');

    r = r.data;
    j = j.data;

    // copy j -> J
    for( let i=M*N; i-- > 0; )
      J[i] = j[i];

    // copy r -> R
    for( let i=M; i-- > 0; )
      R[i] = r[i];

    // compute G
    G.fill(0);
    for( let i=M; i-- > 0; )
    for( let j=N; j-- > 0; )
      G[j] += J[N*i+j] * R[i] / M; // <- TODO check underflow?

    // compute mean square error
    mse = 0;
    for( let i=M; i-- > 0; ) {
      const  r = R[i];
      mse += r*r / M; // <- TODO check underflow?
    }
    mse *= 0.5;
  };

  const compute_dx = (J,R,lambda) => {
    let dx = rls(
      new NDArray(J_shape, J.slice()),
      new NDArray(R_shape, R.slice()),
      lambda
    );  dx = asarray(dx);
    if( dx.ndim     !== 2 ) throw new Error('Assertion failed.');
    if( dx.shape[0] !== N ) throw new Error('Assertion failed.');
    if( dx.shape[1] !== 1 ) throw new Error('Assertion failed.');
    return dx.data;
  };

  if( lambdaFactor === 1 )
    for(;;)
    {
      recompute();

      yield [
        /*meanSquareErr=*/mse,
        /*gradient =*/new NDArray(x_shape, G.slice()),
        /*params   =*/new NDArray(x_shape, x.slice()),
        /*residuals=*/new NDArray(r_shape, R.slice())
      ];

      const dx = compute_dx(J,R, lambda0);

      for( let i=N; i-- > 0; )
        x[i] -= dx[i];
    }
  else
  {
    const [shrink,grow] = lambdaFactor instanceof Array
      ? lambdaFactor
      : [1/lambdaFactor, lambdaFactor];

    recompute();

    let   _mse,
          _R = R.slice(),
          _J = J.slice(),
          _x = x.slice(),
      lambda = lambda0;

    for(;;)
    {
      yield [
        mse,
        /*gradient =*/new NDArray(x_shape, G.slice()),
        /*params   =*/new NDArray(x_shape, x.slice()),
        /*residuals=*/new NDArray(r_shape, R.slice())
      ];

      [_R,R] = [R,_R];
      [_J,J] = [J,_J];
      [_x,x] = [x,_x];
      _mse = mse;

      for(;;)
      {
        const dx = compute_dx(_J,_R, lambda);

        for( let i=N; i-- > 0; )
          x[i] = _x[i] - dx[i];

        recompute();

        let est = 0; // <- estimated mean square error
        for( let i=0; i < M; i++ ) { let r = 0;
        for( let j=0; j < N; j++ )       r-= _J[N*i+j] * dx[j];
            r += _R[i]
          est += r*r/M;
        }
        est *= 0.5;

        if( est > _mse ) throw new Error('Assertion failed.'); // <- TODO comment out after testing

        const improvement = (_mse - mse)
                          / (_mse - est);

        if( improvement < expectLow  ) // <- if worse  than expected, shift towards gradient descent
        {
          let l = lambda*grow;
          if( l===lambda ) l += Number.MIN_VALUE;
          if( l===lambda ) throw new Error('Assertion failed.');
          if(!isFinite(l)) throw new Error('Assertion failed.');
          lambda = l;
        }
        else if( improvement > expectHigh ) // <- if better than expected, shift towards gauss-newton
        {
          if( 0===lambda ) break;
          let l = lambda*shrink;
          if( l===lambda ) l -= Number.MIN_VALUE;
          if( l===lambda ) throw new Error('Assertion failed.');
          lambda = l;
        }

        if( improvement >= accept )
          break;
      }
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
  x = asarray(x, 'float64');
  y = asarray(y, 'float64');
  p0= asarray(p0,'float64');

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

  const fj = p =>
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

  yield* lsq_lm_gen(fj, p0, opt);
}
