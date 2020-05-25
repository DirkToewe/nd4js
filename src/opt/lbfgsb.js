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

import {array,
      asarray,
      NDArray} from '../nd_array'

import {nextDown} from '../dt/float64_utils';

import {LineSearchError, LineSearchNoProgressError, LineSearchBoundReachedError, LineSearchBisectionError} from './line_search/line_search_error'
import {more_thuente_u123} from './line_search/more_thuente_u123'

import {LBFGSB_Solver} from "./_lbfgsb_solver";


// REFERENCES
// ----------
// .. [1] "On the Limited Memory BFGS Method for Large Scale Optimization."
//         Dong C. Liu and Jorge Nocedal


function dot( u, v )
{
  if(  u.length !== v.length ) throw new Error()
  let result = 0
  for( let i=u.length; i-- > 0; )
    result += u[i] * v[i]
  return result
}


export function* min_lbfgsb_gen(
  fg,
  x0,
  bounds,
  {
    historySize = 8,
    lineSearch = more_thuente_u123({fRed: 1e-3, gRed: 0.8}), // <- there's not always a minimum along search direction (due to bounds) -> use _u123 instead of _abc
    updateTol = 1e-14
  } = {}
)
{

  if(  historySize%1 !== 0 ) throw new Error('min_lbfgsb_gen(fg, x0, opt): opt.historySize must be an integer.');
  if(!(historySize    >  0)) throw new Error('min_lbfgsb_gen(fg, x0, opt): opt.historySize must be positive.');

  if( ! (fg instanceof Function) )
    throw new Error('min_lbfgsb_gen(fg,x0,bounds): fg must be a function.');

  x0    =  array('float64', x0);
  bounds=asarray(bounds);

  const [L] = x0.shape,
      shape = x0.shape;

  if(    x0.ndim     !== 1 ) throw new Error('min_lbfgsb_gen(fg,x0,bounds): x0.ndim must be 1.');
  if(bounds.ndim     !== 2 ) throw new Error('min_lbfgsb_gen(fg,x0,bounds): bounds.ndim must be 2.');
  if(bounds.shape[0] !== L ) throw new Error('min_lbfgsb_gen(fg,x0,bounds): bounds.shape[0] must be x0.shape[0].');
  if(bounds.shape[1] !== 2 ) throw new Error('min_lbfgsb_gen(fg,x0,bounds): bounds.shape[1] must be 2.');
  if(  ! (updateTol >= 0)  ) throw new Error('min_lbfgsb_gen(fg,x0,bounds,opt): opt.updateTol most be non-negative.');

  bounds = bounds.data;
  updateTol *= 1;

  for( let i=bounds.length; (i-=2) >= 0; )
    if( ! (bounds[i+0] <= bounds[i+1]) )
      throw new Error('Assertion failed.');

  let x = x0;
          x0 = undefined;
  let [f,g] = fg(x);
  f *= 1;
  g = array('float64', g);

  if( ! isFinite(f) ) throw new Error('min_lbfgsb_gen(fg, x0, bounds): fg must return [float,float[]].');
  if( ! g.dtype.startsWith('float') ) throw new Error('min_lbfgsb_gen(fg, x0, bounds): fg must return [float,float[]].');

  lineSearch = lineSearch(fg);

  const solver = new LBFGSB_Solver(historySize, L);
  solver.scale = 1024; // TODO: Implement scaling as described in [1]

  const indices =       Int32Array.from({length: L}, (_,i) => i),
             dx = new Float64Array(L),
             dg = new Float64Array(L);

  const step = () =>
  {
     const g_arg = g.data.slice(),
           x_arg = x.data.slice(),
           n_dir = new Float64Array(L);

     const [n_fix, df] = solver.compute_cauchyGeneralized(g_arg, x_arg, bounds, indices, /*complete=*/false);
                         solver.compute_subspace_Hv(g_arg, n_dir, indices, n_fix);

//*DEBUG*/    if( 0 !== n_fix ) throw new Error('Not yet tested.');                        
//*DEBUG*/    if( 0 !== df ) throw new Error('Not yet tested.');
//*DEBUG*/    if( ! g_arg.every( (gi,i) => gi === g.data[i] ) ) throw new Error('Not yet tested.');

/*DEBUG*/    if( ! indices.subarray(0,n_fix).every(i => n_dir[i]===0) )
/*DEBUG*/      throw new Error('Assertion failed.');

    const αMax = function(){
      let αMax = Infinity;

      for( const i of indices.subarray(n_fix) )
        if( 0 !== n_dir[i] ) { const      end = bounds[2*i + (n_dir[i] < 0)];
          let        travel = (x_arg[i] - end) / n_dir[i];
          if( !(0 <= travel) )
            throw new Error('Assertion failed: ' + travel);

          let n=0;
          while( Math.sign(n_dir[i]) * (x_arg[i] - travel*n_dir[i] - end) < 0 ) {
            if( ++n > 2 )
              throw new Error('Assertion failed.');
            travel = nextDown(travel);
          }

          if( n_dir[i] < 0 && x_arg[i] - travel*n_dir[i] > end ) throw new Error('Assertion failed.');
          if( n_dir[i] > 0 && x_arg[i] - travel*n_dir[i] < end ) throw new Error('Assertion failed.');

          αMax = Math.min(αMax, travel);
        }

      return αMax;
    }();

    // TODO: g_arg is only the quasi newton approximation of the gradient at x_arg.
    //       Alternatively we could compute the gradient at x_arg insteal which might
    //       improve robustness.

    let    X,F,G;
    try { [X,F,G] = lineSearch(
            new NDArray(shape, x_arg.slice()), f+df,
            new NDArray(shape, g_arg),
            new NDArray(shape, n_dir),
            /*αMin=*/undefined,
            /*α0  =*/undefined,
            /*αMax=*/αMax
          );
    }
    catch(err) {
      if(  err instanceof LineSearchBoundReachedError
        || err instanceof LineSearchBisectionError )
      {
        ({x: X,
          f: F,
          g: G} = err);
      }
      else if( err instanceof LineSearchNoProgressError && df < 0 )
      {
          X  =      new NDArray(shape, x_arg);
        [F,G] = fg( new NDArray(shape, x_arg.slice()) );
      }
      else
        throw err;
    }

    F *= 1;
    X = array('float64', X);
    G = array('float64', G);

    for( let i=L; i-- > 0; ) dx[i] = X.data[i] - x.data[i];
    for( let i=L; i-- > 0; ) dg[i] = G.data[i] - g.data[i];

    if( dx.every(x => x===0) )
      throw new LineSearchNoProgressError();

    const mx = dg.reduce((mx,x) => Math.max(mx,Math.abs(x)), 0); // <- avoid underflow

    let gg = 0,
        gx = 0;

    for( let i=dg.length; i-- > 0; )
    { const xi = dx[i]/mx,
            gi = dg[i]/mx;
      gg += gi*gi;
      gx += xi*gi;
    }

    if( updateTol*gg < gx )
    {
      solver.update(dx,dg);
      solver.scale = gg/gx;
    }

    x = X;
    f = F;
    g = G;
  }

  for(;;)
  {
    if( ! (x instanceof NDArray) ) throw new Error('Assertion failed.')
    if( ! (g instanceof NDArray) ) throw new Error('Assertion failed.')
    if( x.shape[0] !== L ) throw new Error('Assertion#7 failed.');
    if( g.shape[0] !== L ) throw new Error('Assertion#8 failed.');
    if( x.ndim !== 1 ) throw new Error('Assertion#4 failed.');
    if( g.ndim !== 1 ) throw new Error('Assertion#6 failed.');
    if( isNaN(f) ) throw new Error('Assertion#5 failed.');

    const γ = g.data.slice();
    x.data.forEach( (xi,i) => {
      const [x_min,x_max] = bounds.subarray(2*i,2+2*i),
                       gi = γ[i];
           if( 0 < gi && xi===x_min ) γ[i] = 0;
      else if( 0 > gi && xi===x_max ) γ[i] = 0;
    });
    
    yield [
      new NDArray(shape, x.data.slice()), f,
      new NDArray(shape, γ),
      new NDArray(shape, g.data.slice())
    ];

    loop:for(;;)
      try{
        step();
        break loop;
      }
      catch(err) {
        if( err instanceof LineSearchError && solver.m > 0 ) {
          solver.forget(solver.m+1 >>> 1);
          // if( solver.m === 0 )
          //   solver.scale = 1024;
        }
        else
          throw err;
      }
  }
}
