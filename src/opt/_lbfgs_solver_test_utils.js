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

import {NDArray} from '../nd_array'
import {_rand_int,
        _shuffle} from '../_test_data_generators'
import {zip_elems} from '../zip_elems'

import {eye} from '../la/eye'
import {matmul,
        matmul2} from '../la/matmul'


/** Produces an infinite series of LBFGS update input pairs.
 *  The are filtered to result in a sufficiently well-conditioned
 *  (positive definite) Hessian approximation.
 */
export function* _rand_updates_v2( rng, N )
{
  if( N%1 !== 0 ) throw new Error('Assertion failed.');
  if( ! (0 < N) ) throw new Error('Assertion failed.');

  const shape = Int32Array.of(N,1);

  const updateTol = 1/16;

  loop:for(;;)
  {
    const dx = Float64Array.from({length: N}, () => rng.normal() ),
          dg = Float64Array.from({length: N}, () => rng.normal() ); // <- avoid underflow

    const sx = rng.uniform(0.5, 2) / Math.hypot(...dx),
          sg = rng.uniform(0.5, 2) / Math.hypot(...dg);

    for( let i=N; i-- > 0; ) {
      dx[i] *= sx;
      dg[i] *= sg;
    }

    const mx = dg.reduce((mx,x) => Math.max(mx,Math.abs(x)), 0);

    let gg = 0,
        gx = 0;

    for( let i=dg.length; i-- > 0; )
    { const xi = dx[i]/mx,
            gi = dg[i]/mx;
      gg += gi*gi;
      gx += xi*gi;
    }

    if( ! (updateTol*gg < gx) )
      continue loop;

    Object.freeze(dx.buffer);
    Object.freeze(dg.buffer);
    yield Object.freeze([
      new NDArray(shape,dx),
      new NDArray(shape,dg)
    ]);
  }
}


/** A slow reference implementation of the LBFGS Hessian approximation
 *  used to test the much more efficient `LBFGSSolver`.
 */
export class LBFGS_SolverRefImpl
{
  constructor( M, N )
  {
    if( M%1 !== 0 ) throw new Error('Assertion failed.');
    if( N%1 !== 0 ) throw new Error('Assertion failed.');
    if( !(0 < M) ) throw new Error('Assertion failed.');
    if( !(0 < N) ) throw new Error('Assertion failed.');
    this.M = M|0;
    this.N = N|0;
    this.dx = [];
    this.dg = [];
    this.scale = 1;
  }

  update( dx, dg )
  {
    if( dx.length%1 !== 0 ) throw new Error('Assertion failed.');
    if( dg.length%1 !== 0 ) throw new Error('Assertion failed.');

    const {dx: dX, dg: dG, N} = this;

    dx = Float64Array.from(dx);
    dg = Float64Array.from(dg);

    if( dx.length !== N ) throw new Error('Assertion failed: ' + dx.length);
    if( dg.length !== N ) throw new Error('Assertion failed: ' + dg.length);

    const shape = Int32Array.of(N,1);

    dx = new NDArray(shape, dx);
    dg = new NDArray(shape, dg);
    Object.freeze(dx.buffer);
    Object.freeze(dg.buffer);

    if( dX.length != dG.length ) throw new Error('Assertion failed.');
    dX.push(dx);
    dG.push(dg);
    if( dX.length > this.M ) {
        dX.shift();
        dG.shift();
    }
    if( dX.length > this.M ) throw new Error('Assertion failed.');
  }

  *_history()
  {
    const {dx,dg} = this;
    if( dx.length != dg.length ) throw new Error('Assertion failed.');
    for( let i=0; i < dx.length; i++ )
      yield [dx[i], dg[i]];
  }

  forget( k )
  {
    const {dx,dg} = this;
    if( 0 !== k%1 ) throw new Error('Assertion failed.')
    if( !(0 <= k) ) throw new Error('Assertion failed.');
    if( !(k <= dx.length) ) throw new Error('Assertion failed.');
    dx.splice(0,k);
    dg.splice(0,k);
  }

  get B()
  {
    const {N} = this;

    let B = eye(N).mapElems(x => x * this.scale);

    for( const [dx,dg] of this._history() )
      B = zip_elems(
        [ B,
          matmul2(dg, dg.T),
          matmul2(dg.T, dx),
          matmul(B,dx, dx.T, B.T),
          matmul(dx.T, B, dx) ],
        (B, gg, gx, BxxB, xBx) => B + gg/gx - BxxB/xBx
      );

    Object.freeze(B.data.buffer);
    Object.freeze(B);
    return B;
  }

  get H()
  {
    const {N} = this;

    const I = eye(N);

    let H = I.mapElems(x => x / this.scale);

    for( const [dx,dg] of this._history() )
    {
      const gx = matmul2(dg.T, dx),
             W = zip_elems(
               [I, matmul2(dx, dg.T), gx],
               (I,xg,gx) => I - xg/gx
             );

      H = zip_elems(
        [ matmul(W, H, W.T),
          matmul2(dx, dx.T),
          gx ],
        (WHW,xx,gx) => WHW + xx/gx
      );
    }

    Object.freeze(H.data.buffer);
    Object.freeze(H);
    return        H;
  }
}
