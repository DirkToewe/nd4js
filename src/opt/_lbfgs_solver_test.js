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

import {forEachItemIn, CUSTOM_MATCHERS} from '../jasmine_utils'
import {NDArray} from '../nd_array'
import {zip_elems} from '../zip_elems'

import {eye} from '../la/eye'
import {matmul,
        matmul2} from '../la/matmul'
import {norm} from '../la/norm'

import {LBFGSSolver} from './_lbfgs_solver'


class LBFGSSolverRef
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


const dot = (u,v) => {
  if( u.length%1 !== 0       ) throw new Error('Assertion failed.');
  if( v.length%1 !== 0       ) throw new Error('Assertion failed.');
  if( u.length   !== v.length) throw new Error('Assertion failed.');

  let uv = 0;
  for( let i=u.length; i-- > 0; )
    uv += u[i]*v[i];

  return uv;
};


function* _rand_updates( N )
{
  if( N%1 !== 0 ) throw new Error('Assertion failed.');
  if( ! (0 < N) ) throw new Error('Assertion failed.');

  const shape = Int32Array.of(N,1);

  loop:for(;;)
  {
    const dx = Float64Array.from({length: N}, () => Math.random()*2-1),
          dg = Float64Array.from({length: N}, () => Math.random()*2-1);

    const nx = Math.hypot(...dx),
          ng = Math.hypot(...dg);

    for( let i=N; i-- > 0; ) {
      dx[i] /= nx;
      dg[i] /= ng;
    }

    let dot_xg = dot(dx,dg);
    if( dot_xg < 0 ) {
      for( let i=N; i-- > 0; )
        dx[i] *= -1;
      dot_xg *= -1;
    }
    if( dot_xg < 0.05 )
      continue loop;

    const sx = Math.random()*1 + 0.5,
          sg = Math.random()*1 + 0.5;

    for( let i=N; i-- > 0; ) {
      dx[i] *= sx;
      dg[i] *= sg;
    }

    Object.freeze(dx.buffer);
    Object.freeze(dg.buffer);
    yield Object.freeze([
      new NDArray(shape,dx),
      new NDArray(shape,dg)
    ]);
  }
}


describe('LBFGSSolver', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  for( const [suffix, rand_scale] of [
    [''           , () => 1                  ],
    [' and scales', () => Math.random() + 0.5]
  ])
  {
    forEachItemIn(
      function*(){
        for( let M=0; M++ <  8; )
        for( let N=0; N++ < 16; )
          yield [M,N];
      }()
    ).it('compute_bv and compute_ubbv work given random updates' + suffix, ([M,N]) => {
      const ref = new LBFGSSolverRef(M,N),
            tst = new LBFGSSolver   (M,N);

      expect(M).toBeGreaterThan(0);
      expect(N).toBeGreaterThan(0);

      const shape = Int32Array.of(N,1);

      let i=0;
      for( let [dx,dg] of _rand_updates(N) )
      {
        ref.update(dx.data, dg.data);
        tst.update(dx.data, dg.data);
        for( let outmost=0; outmost++ < 2; )
        {
          const       scale = rand_scale();
          ref.scale = scale;
          tst.scale = scale;

          const {B} = ref;

          const u = Float64Array.from({length: N}, () => Math.random()*8-4),
                v = Float64Array.from({length: N}, () => Math.random()*8-4),
               bu = new Float64Array(2*M),
               bv = new Float64Array(2*M);
          Object.freeze(u.buffer);
          Object.freeze(v.buffer);

          const U = new NDArray(shape, u),
                V = new NDArray(shape, v);

          const UBU = matmul(U.T, B, U),
                UBV = matmul(U.T, B, V),
                VBV = matmul(V.T, B, V);
          Object.freeze(UBU); Object.freeze(UBU.data.buffer);
          Object.freeze(UBV); Object.freeze(UBV.data.buffer);
          Object.freeze(VBV); Object.freeze(VBV.data.buffer);
          
          for( let outer=0; outer++ < 2; )
          {
            tst.compute_bv(u, bu);

            const    uBu = tst.compute_ubbv(bu, dot(u,u), bu);
            expect([[uBu]]).toBeAllCloseTo(UBU);

            for( let inner=0; inner++ < 2; )
            {
              tst.compute_bv(v, bv);

              const    uBv = tst.compute_ubbv(bu, dot(u,v), bv);
              expect([[uBv]]).toBeAllCloseTo(UBV);

              const    vBu = tst.compute_ubbv(bv, dot(v,u), bu);
              expect([[vBu]]).toBeAllCloseTo(UBV);

              const    vBv = tst.compute_ubbv(bv, dot(v,v), bv);
              expect([[vBv]]).toBeAllCloseTo(VBV);
            }
          }
        }

        if( ++i >= 64 ) break;
      }
    })


    forEachItemIn(
      function*(){
        for( let M=0; M++ < 12; )
        for( let N=0; N++ < 24; )
          yield [M,N];
      }()
    ).it('compute_be works given random updates' + suffix, ([M,N]) => {
      const tst = new LBFGSSolver   (M,N);

      expect(M).toBeGreaterThan(0);
      expect(N).toBeGreaterThan(0);

      const shape = Int32Array.of(N,1);

      let i=0;
      for( let [dx,dg] of _rand_updates(N) )
      {
        tst.update(dx.data, dg.data);

        const m = Math.min(M,i+1);

        const BEI = new Float64Array(2*M),
              bei = new Float64Array(2*M);

        for( let i=0; i < N; i++ )
        {
          const EI = Float64Array.from({length: N}, (_,j) => (i===j)*1);
          Object.freeze(EI.buffer);

          for( let run=2; run-- > 0; )
          {
            for( let k=2*M; k-- > 0; )
            {
              BEI[k] = Math.random()*8-4;
              bei[k] = Math.random()*8-4;
            }
            tst.scale = rand_scale();
            tst.compute_bv(EI, BEI);
            tst.compute_be( i, bei);

              expect( bei.subarray(0,m*2) )
            .toEqual( BEI.subarray(0,m*2) );
          }
        }

        if( ++i >= 128 ) break;
      }
    })
  }
})


describe('LBFGSSolverRef', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  for( const [suffix, rand_scale] of [
    [''           , () => 1                  ],
    [' and scales', () => Math.random() + 0.5]
  ])
  {
    forEachItemIn(
      function*(){
        for( let M=0; M++ <  8; )
        for( let N=0; N++ < 16; )
          yield [M,N];
      }()
    ).it('computes B correctly given random updates' + suffix, ([M,N]) => {
      const lbfgs = new LBFGSSolverRef(M,N);
    
      expect(M).toBeGreaterThan(0);
      expect(N).toBeGreaterThan(0);
    
      let i=0;
      for( let [dx,dg] of _rand_updates(N) )
      {
        lbfgs.update(dx.data, dg.data);
        lbfgs.scale = rand_scale();
    
        const {B} = lbfgs;
        expect( matmul2(B,dx) ).toBeAllCloseTo(dg);
    
        if( ++i >= 128 ) break;
      }
    })
    
    
    forEachItemIn(
      function*(){
        for( let M=0; M++ <  8; )
        for( let N=0; N++ < 16; )
          yield [M,N];
      }()
    ).it('computes H correctly given random updates' + suffix, ([M,N]) => {
      const lbfgs = new LBFGSSolverRef(M,N);
    
      expect(M).toBeGreaterThan(0);
      expect(N).toBeGreaterThan(0);
    
      const shape = Int32Array.of(N,1);
      Object.freeze(shape.buffer);
    
      let i=0;
      for( let [dx,dg] of _rand_updates(N) )
      {
        lbfgs.update(dx.data, dg.data);
        lbfgs.scale = rand_scale();
    
        const {H} = lbfgs;
    
        const atol = Number.EPSILON * norm(H) * norm(dx);
        expect( matmul2(H,dg) ).toBeAllCloseTo(dx, {atol});
    
        if( ++i >= 128 ) break;
      }
    })


    forEachItemIn(
      function*(){
        for( let M=0; M++ <  8; )
        for( let N=0; N++ < 16; )
          yield [M,N];
      }()
    ).it('H=inv(B) given random updates' + suffix, ([M,N]) => {
      const lbfgs = new LBFGSSolverRef(M,N);
    
      expect(M).toBeGreaterThan(0);
      expect(N).toBeGreaterThan(0);
    
      const shape = Int32Array.of(N,1);
      Object.freeze(shape.buffer);
    
      const I = eye(N);
    
      let i=0;
      for( let [dx,dg] of _rand_updates(N) )
      {
        lbfgs.update(dx.data, dg.data);
        lbfgs.scale = rand_scale();
    
        const {H,B} = lbfgs;
    
        const atol = Number.EPSILON * (1<<15) * norm(H) * norm(B);
        expect( matmul2(B,H) ).toBeAllCloseTo(I, {atol});
        expect( matmul2(H,B) ).toBeAllCloseTo(I, {atol});
    
        if( ++i >= 64 ) break;
      }
    })
  }
})
