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

import {LineSearchBisectionError,
        LineSearchError,
        LineSearchNoProgressError} from './line_search/line_search_error'
import {LBFGS_Solver} from './_lbfgs_solver';
import {more_thuente_abc} from './line_search/more_thuente_abc';


// REFERENCES
// ----------
// .. [1] "On the Limited Memory BFGS Method for Large Scale Optimization."
//         Dong C. Liu and Jorge Nocedal
// .. [2] "Updating Quasi-Newton Matrices with Limited Storage"
//         Jorge Nocedal 
//         MATHEMATICS OF COMPUTATION, VOL. 35, NO. 151, JULY 1980, PAGES 773-78
//         https://courses.engr.illinois.edu/ece544na/fa2014/nocedal80.pdf


// TODO: Implement scaling as described in [1].


function dot( u, v )
{
  if(  u.length !== v.length ) throw new Error()
  let result = 0
  for( let i=u.length; i-- > 0; )
    result += u[i] * v[i]
  return result
}


export function* min_lbfgs_gen(
  fg,
  x0,
  {
    historySize=8,
    lineSearch=more_thuente_abc(),
    updateTol = 1e-13,
    negDir0 = g => g,
    scaling = function(){
      let scale = 1/1024;
      return {
        update( dx, dg ) {
          if( dx.ndim !== 1 ) throw new Error('Assertion failed.');
          if( dg.ndim !== 1 ) throw new Error('Assertion failed.');
          dx = dx.data;
          dg = dg.data;
          const len = dx.length;
          if(   len!==dg.length )
            throw new Error('Assertion failed.');

          const mx = dg.reduce((mx,x) => Math.max(mx,Math.abs(x)), 0); // <- avoid underflow
    
          let gg = 0,
              gx = 0;
    
          for( let i=dg.length; i-- > 0; )
          { const       gi = dg[i]/mx,
                        xi = dx[i]/mx;
            gg += gi*gi;
            gx += xi*gi;
          }

          scale = gx/gg;

          if( ! isFinite(scale) )
            throw new Error('Assertion failed.');
        },
        apply( negDir ) {
          if( negDir.ndim !== 1 )
            throw new Error('Assertion failed.');
          const d = negDir.data;

          for( let i=d.length; i-- > 0; )
            d[i] *= scale;

          return negDir;
        }
      };
    }()
  } = {}
)
{
  if(  historySize%1 !== 0 ) throw new Error('min_lbfgs_gen(fg, x0, opt): opt.historySize must be an integer.');
  if(!(historySize    >  0)) throw new Error('min_lbfgs_gen(fg, x0, opt): opt.historySize must be positive.');

  let x = array('float64', x0);
                           x0 = undefined;

  if( x.ndim !== 1 ) throw new Error('min_lbfgs_gen(fg, x0, opt): x0.ndim must be 1.');

  const [L] = x.shape,
      shape = x.shape;

  let [f,g] = fg(x);

  g = array('float64', g);
  f*= 1;

  if( isNaN(f)         ) throw new Error('min_lbfgs_gen(fg, x0, opt): fg must return [float,float[]].');
  if( g.ndim     !== 1 ) throw new Error('min_lbfgs_gen(fg, x0, opt): fg must return [float,float[]].');
  if( g.shape[0] !== L ) throw new Error('min_lbfgs_gen(fg, x0, opt): fg must return [float,float[]].');

  lineSearch = lineSearch(fg);

  const solver = new LBFGS_Solver(historySize, L);

  const dx = new Float64Array(L),
        dg = new Float64Array(L);

  const step = () =>
  {
    let                      negDir = Float64Array.from(g.data);
    solver.compute_Hv_phase1(negDir);
                             negDir = negDir0( new NDArray(shape, negDir) );
                             negDir = array('float64', negDir);
                             negDir = scaling.apply(negDir);
                             negDir = array('float64', negDir);
    solver.compute_Hv_phase2(negDir.data);

    let X,F,G;

    let update = true;
    try {
      [X,F,G] = lineSearch(
        new NDArray(shape, x.data.slice()), f,
        new NDArray(shape, g.data.slice()), negDir,
      );
    }
    catch(err) {
      if( !(err instanceof LineSearchBisectionError) )
        throw err;
      ({x: X,
        f: F,
        g: G} = err);
      update = false;
    }

    X = array('float64', X);
    G = array('float64', G);
    F*= 1;

    if( X.ndim     !== 1 ) throw new Error('Assertion failed.');
    if( G.ndim     !== 1 ) throw new Error('Assertion failed.');
    if( X.shape[0] !== L ) throw new Error('Assertion failed.');
    if( G.shape[0] !== L ) throw new Error('Assertion failed.');
    if( isNaN(F)         ) throw new Error('Assertion failed.');

    for( let i= L; i-- > 0; ) dg[i] = G.data[i] - g.data[i];
    for( let i= L; i-- > 0; ) dx[i] = X.data[i] - x.data[i];

    if( dx.every(x => x===0) )
      throw new LineSearchNoProgressError();

    update = update || function(){
      const mx = dg.reduce((mx,x) => Math.max(mx,Math.abs(x)), 0); // <- avoid underflow

      let gg = 0,
          gx = 0;

      for( let i=dg.length; i-- > 0; )
      { const       gi = dg[i]/mx,
                    xi = dx[i]/mx;
        gg += gi*gi;
        gx += xi*gi;
      }

      return updateTol*gg < gx;
    }();

    if( update ) {
      solver .update(dx,dg);
      scaling.update(
        new NDArray(shape, dx.slice()),
        new NDArray(shape, dg.slice())
      );
    }

    x = X;
    f = F;
    g = G;
  }

  for(;;)
  {
    if( ! (x instanceof NDArray) ) throw new Error('Assertion failed.')
    if( ! (g instanceof NDArray) ) throw new Error('Assertion failed.')
    if( isNaN(f) ) throw new Error('Assertion#5 failed.')
    if( x.ndim     !== 1 ) throw new Error('Assertion#4 failed.')
    if( g.ndim     !== 1 ) throw new Error('Assertion#6 failed.')
    if( x.shape[0] !== L ) throw new Error('Assertion#7 failed.')
    if( g.shape[0] !== L ) throw new Error('Assertion#8 failed.')

    yield [
      new NDArray(shape, x.data.slice()), f,
      new NDArray(shape, g.data.slice())
    ];

    loop:for(;;)
      try{
        step()
        break loop;
      }
      catch(err) {
        if( err instanceof LineSearchError && solver.m > 0 )
          solver.forget(solver.m+1 >>> 1);
        else
          throw err
      }
  }
}


export function* lsq_lbfgs_gen(
  fJ,//: (x0: float[nInputs]) => [f: float[nSamples], j: float[nSamples,nInputs]] - THE OPTIMIZATION FUNCTION AND ITS JACOBIAN
  x0,//: float[nInputs] - THE START POINT OF THE OPTIMIZATION
  {
    historySize=8,
    lineSearch=more_thuente_abc(),
    updateTol = 1e-13
  } = {}
)
{
  if(  historySize%1 !== 0 ) throw new Error('min_lbfgs_gen(fg, x0, opt): opt.historySize must be an integer.');
  if(!(historySize    >  0)) throw new Error('min_lbfgs_gen(fg, x0, opt): opt.historySize must be positive.');

  // param | meaning
  // ------|-------------------------
  //  x    | function input
  //  f    | function output f(x)
  //  j    | jacobian f'(x)
  //  e    | mean squared error (MSE)
  //  g    | gradient of MSE

  let x = array('float64', x0);
                           x0 = undefined;

  if( x.ndim !== 1 ) throw new Error('min_lbfgs_gen(fg, x0, opt): x0.ndim must be 1.');

  let [f,j] = fJ(x);

  f = array('float64', f);
  j = array('float64', j);

  if( f.ndim !== 1 ) throw new Error('lsq_lbfgs_gen(fg, x0: float[L], opt): fg must return [float[K],float[K,L]].');
  if( j.ndim !== 2 ) throw new Error('lsq_lbfgs_gen(fg, x0: float[L], opt): fg must return [float[K],float[K,L]].');

  const [K,L]   = j.shape,
        shape_x = x.shape,
        shape_f = f.shape,
        shape_j = j.shape,
        shape_e = Int32Array.of();

  if( f.shape[0] !== K ) throw new Error('lsq_lbfgs_gen(fg, x0: float[L], opt): fg must return [float[K],float[K,L]].');
  if( x.shape[0] !== L ) throw new Error('lsq_lbfgs_gen(fg, x0: float[L], opt): fg must return [float[K],float[K,L]].');

  fJ = function(){
    const FJ = fJ;

    return x => {
      x = new NDArray(shape_x, x.data);

      let [f,j] = FJ(x);
      f = array('float64', f);
      j = array('float64', j);

      if( f.data.some(isNaN) ) throw new Error('Assertion failed.');
      if( j.data.some(isNaN) ) throw new Error('Assertion failed.');

      if( f.ndim     !== 1 ) throw new Error('lsq_lbfgs_gen(fg, x0: float[L], opt): fJ(x: float[L]) must return [float[K],float[K,L]].');
      if( j.ndim     !== 2 ) throw new Error('lsq_lbfgs_gen(fg, x0: float[L], opt): fJ(x: float[L]) must return [float[K],float[K,L]].');
      if( f.shape[0] !== K ) throw new Error('lsq_lbfgs_gen(fg, x0: float[L], opt): fJ(x: float[L]) must return [float[K],float[K,L]].');
      if( j.shape[0] !== K ) throw new Error('lsq_lbfgs_gen(fg, x0: float[L], opt): fJ(x: float[L]) must return [float[K],float[K,L]].');
      if( j.shape[1] !== L ) throw new Error('lsq_lbfgs_gen(fg, x0: float[L], opt): fJ(x: float[L]) must return [float[K],float[K,L]].');

      return [f,j];
    };
  }();

  const eg = (f,j) => {
    const F = f.data,
          J = j.data;

    let e = 0,
        G = new Float64Array(L);

    for( let i=K; i-- > 0; )
      e += F[i]*F[i] / K;

    for( let i=K; i-- > 0; )
    for( let j=L; j-- > 0; )
      G[j] += 2*F[i]*J[L*i+j] / K;

    return [e, new NDArray(shape_x, G)];
  };

  let last_x,
      last_f,
      last_j;

  lineSearch = lineSearch( x => {
    last_x = array('float64', x);
                              x = undefined;
    if( last_x.ndim     !== 1 ) throw new Error('Assertion failed.');
    if( last_x.shape[0] !== L ) throw new Error('Assertion failed.');
    [last_f,
     last_j] = fJ(last_x);

    return eg(last_f, last_j);
  });

  const negDir0 = negDir => {
    // OPTION 1: USE CAUCHY POINT OF THE GAUSS NEWTON APPROXIMATION
    const G = g.data,
          J = j.data;

    // compute the minimum along negDir (similar to cauchy point)
    let a=0;
    for( let i=L; i-- > 0; )
      a += G[i]*negDir[i];

    if( 0 < a )
    {
      let b=0;
      for( let i=K; i-- > 0; ) {
        let Jd = 0;
        for( let j=L; j-- > 0; )
          Jd += J[L*i+j] * negDir[j];
        b += Jd*Jd * 2/K;
      }

      if( ! (0 < b) )
        throw new Error('Assertion failed.');

      a /= b;
      for( let i=L; i-- > 0; )
        negDir[i] *= a;
    }

    // // OPTION 1: USE GAUSS NEWTON APPROXIMATION
    // const nd = lstsq(
    //   matmul2(j.T, j),
    //   new NDArray(Int32Array.of(L,1), negDir)
    // );

    // negDir.set(nd.data);

    // for( let i=L; i-- > 0; )
    //   negDir[i] *= K/2;
  };

  const dx = new Float64Array(L),
        dg = new Float64Array(L),
    solver = new LBFGS_Solver(historySize, L);
  let[e,g] = eg(f,j);

  const step = () =>
  {
    let                      negDir = Float64Array.from(g.data);
    solver.compute_Hv_phase1(negDir);
                     negDir0(negDir);
    solver.compute_Hv_phase2(negDir);

    let X,E,G,update = true;
    try {
      [X,E,G] = lineSearch(
        new NDArray(shape_x, x.data.slice()), e,
        new NDArray(shape_x, g.data.slice()),
        new NDArray(shape_x, negDir)
      );
    }
    catch(err) {
      if( ! (err instanceof LineSearchBisectionError) )
        throw err;
      ({x: X,
        f: E,
        g: G} = err);
      update = false;
    }

    x = array('float64', x); e *= 1;
    g = array('float64', g);

    if( x.ndim     !== 1 ) throw new Error('Assertion failed.');
    if( g.ndim     !== 1 ) throw new Error('Assertion failed.');
    if( x.shape[0] !== L ) throw new Error('Assertion failed.');
    if( g.shape[0] !== L ) throw new Error('Assertion failed.');
    if( isNaN(e)         ) throw new Error('Assertion failed.');

    let F = last_f,
        J = last_j;

    // only recompute f,j is necessary
    if( last_x.data.some( (xi,i) => xi !== X.data[i] ) )
      [F,J] = fJ(X);

    for( let i= L; i-- > 0; ) dx[i] = X.data[i] - x.data[i];
    for( let i= L; i-- > 0; ) dg[i] = G.data[i] - g.data[i];

    if( dx.every(x => x===0) )
      throw new LineSearchNoProgressError();

    update = update || function(){
      const mx = dg.reduce((mx,x) => Math.max(mx,Math.abs(x)), 0); // <- avoid underflow

      let gg = 0,
          gx = 0;

      for( let i=dg.length; i-- > 0; )
      { const       gi = dg[i]/mx,
                    xi = dx[i]/mx;
        gg += gi*gi;
        gx += xi*gi;
      }

      return updateTol*gg < gx;
    }();

    if(update)
      solver.update(dx,dg);

    x = X;
    e = E;
    g = G;
    f = F;
    j = J;
  }

  for(;;)
  {
    yield [
      new NDArray(shape_x, x.data.slice()), e,
      new NDArray(shape_x, g.data.slice()),
      new NDArray(shape_f, f.data.slice()),
      new NDArray(shape_j, j.data.slice())
    ];

    loop:for(;;)
      try{
        step()
        break loop;
      }
      catch(err) {
        if( err instanceof LineSearchError && solver.m > 0 )
          solver.forget(solver.m+1 >>> 1);
        else
          throw err
      }
  }
}
