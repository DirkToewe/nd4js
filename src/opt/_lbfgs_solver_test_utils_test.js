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
import {_rand_int,
        _shuffle} from '../_test_data_generators'

import {diag_mat} from '../la/diag'
import {eye} from '../la/eye'
import {matmul,
        matmul2} from '../la/matmul'
import {norm} from '../la/norm'
import {rand_ortho} from '../la/rand_ortho'

import {LBFGS_SolverRefImpl, _rand_updates} from './_lbfgs_solver_test_utils'


describe('LBFGS_SolverRefImpl', () => {
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
      const lbfgs = new LBFGS_SolverRefImpl(M,N);
    
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
      const lbfgs = new LBFGS_SolverRefImpl(M,N);
    
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
      const lbfgs = new LBFGS_SolverRefImpl(M,N);
    
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
    
        const tol = {
          rtol: 0,
          atol: Number.EPSILON * (1<<18) * norm(H) * norm(B)
        };
        expect( matmul2(B,H) ).toBeAllCloseTo(I, tol);
        expect( matmul2(H,B) ).toBeAllCloseTo(I, tol);
    
        if( ++i >= 64 ) break;
      }
    })


    forEachItemIn(
      function*(){
        for( let run=0; run++ < 512; )
        for( let N=0; N++ < 16; )
        {
          const   d = Float64Array.from({length: N}, () => Math.random()*1.8 + 0.2),
                  D = diag_mat(d),
                  Q = rand_ortho(N),
                  B = matmul(Q.T, D, Q),
            updates = [];

          let i=0;
          for( const {data: row} of Q )
          {
            const dxScale = Math.random()*1.8 + 0.2,
                  dgScale = d[i];
            i++;
            const dx = row.map(dx => dx*dxScale),
                  dg =  dx.map(dg => dg*dgScale);
            Object.freeze(dx.buffer);
            Object.freeze(dg.buffer);
            updates.push([dx,dg]);
          }

          Object.freeze(updates);
          Object.freeze(B.data.buffer);
          Object.freeze(B);

          yield [updates,B];
        }
      }()
    ).it('approximates Hessian given random multivariate polynomials' + suffix, ([updates,B]) => {
      const N = B.shape[0];

      const lbfgs = new LBFGS_SolverRefImpl(N,N);

      for( const [dx,dg] of updates )
        lbfgs.update(dx,dg);
      lbfgs.scale = rand_scale();

      const  b = lbfgs.B;
      expect(b).toBeAllCloseTo(B);
    });
  }
})
