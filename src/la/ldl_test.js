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

import {forEachItemIn, CUSTOM_MATCHERS} from '../jasmine_utils'
import {NDArray} from '../nd_array'
import {tabulate} from '../tabulate'
import {zip_elems} from '../zip_elems'

import {ldl_decomp,
        ldl_solve} from './ldl'
import {matmul} from './matmul'


describe('LDLᵀ Decomposition', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  forEachItemIn(
    function*(rng){
      for( let run=2048; run-- > 0; )
      {
        let ndim = rng.int(0,4),
          shapes = [ Array.from({length: ndim}, () => rng.int(1,4)) ]
        shapes.splice( rng.int(0,2), 0, shapes[0].slice( rng.int(0,ndim+1) ) )

        for( let d=ndim; d > 0; d-- )
        for( let i=rng.int(0,2); i-- > 0; ) {
          const    shape = shapes[rng.int(0,2)],
               j = shape.length - d
          if(0<=j) shape[j] = 1
        }

        const M = rng.int(1,32); shapes[0].push(M,M)
        const N = rng.int(1,32); shapes[1].push(M,N)

        const y = tabulate(shapes[1],'float64', () => rng.uniform(-4,+4) ),
              LD= tabulate(shapes[0],'float64', (...indices) => {
                const [i,j] = indices.slice(-2);
                return i===j ? rng.uniform( 0.4, 1.6) * (rng.bool() ? +1 : -1)
                     : i < j ? 0
                     :         rng.uniform(-1.2,+1.2);
              })

        yield [LD,y]
      }
    }
  ).it('ldl_solve works on random examples', ([LD,y]) => {
    const x = ldl_solve(LD,y),
          L = LD.mapElems( 'float64', (L,...idx) => { const [i,j] = idx.slice(-2); return i===j ? 1 : L } ),
          D = LD.mapElems( 'float64', (L,...idx) => { const [i,j] = idx.slice(-2); return i===j ? L : 0 } ),
          Y = matmul(L, D, L.T, x);

    expect(Y).toBeAllCloseTo(y)
  });


  forEachItemIn(
    function*(rng){
      for( let run=733; run-- > 0; )
      {
        const shape = Int32Array.from({ length: rng.int(2,6) }, () => rng.int(1,4) );
              shape[shape.length-2] =
              shape[shape.length-1] = rng.int(1,32);

        const L = tabulate(shape,'float64',(...indices) => {
          const [i,j] = indices.slice(-2);
          if( i < j ) return 0;
          if( i===j ) return 1;
                      return rng.uniform(-1.2,+1.2);
        });

        const D = tabulate(shape,'float64',(...indices) => {
          const [i,j] = indices.slice(-2);
          if( i===j ) return rng.uniform(0.4,1.6) * (rng.bool() ? +1 : -1);
                      return  0;
        });

        Object.freeze(L); Object.freeze(L.data.buffer);
        Object.freeze(D); Object.freeze(D.data.buffer);

        yield [L,D];
      }
    }
  ).it('ldl_decomp works on random examples', ([L,D]) => {
    expect(L).toBeLowerTriangular();
    expect(D).toBeDiagonal();
    expect(L.shape).toEqual(D.shape);
    expect(L.shape[L.ndim-2]).toBe(L.shape[L.ndim-1]);

    const LDLT= matmul(L,D,L.T),
          ld  = ldl_decomp(LDLT),
          l   = ld.mapElems( 'float64', (L,...idx) => { const [i,j] = idx.slice(-2); return i===j ? 1 : L } ),
           d  = ld.mapElems( 'float64', (L,...idx) => { const [i,j] = idx.slice(-2); return i===j ? L : 0 } );

    expect(ld.shape).toEqual(L.shape);
    expect(ld).toBeLowerTriangular();

    expect(d.shape).toEqual(D.shape);
    expect(d).toBeDiagonal();
    expect(d).toBeAllCloseTo(D);

    expect(l.shape).toEqual(L.shape);
    expect(l).toBeLowerTriangular();
    expect(l).toBeAllCloseTo(L);
  });


  forEachItemIn(
    function*(rng){
      for( let N=3; N < 768; N = N*1.4 | 0 )
      {
        const A = new NDArray(
            Int32Array.of(N,N),
          Float64Array.from({length: N*N}, () => rng.uniform(-4,+4) )
        );
        const         LDLT = zip_elems([A,A.T], 'float64', (a,at) => (a+at)/2);
        Object.freeze(LDLT.data.buffer)
        Object.freeze(LDLT)
        yield         LDLT
      }
    }
  ).it('ldl_decomp works on random large examples', LDLT => {
    const ld  = ldl_decomp(LDLT),
          l   = ld.mapElems('float64', (L,i,j) => i===j ? 1 : L),
           d  = ld.mapElems('float64', (L,i,j) => i===j ? L : 0),
          ldlt= matmul(l, d, l.T);

    expect(ld.shape).toEqual(LDLT.shape);
    expect(d .shape).toEqual(LDLT.shape);
    expect( l.shape).toEqual(LDLT.shape);

    expect(ld).toBeLowerTriangular();
    expect( d).toBeDiagonal();
    expect(l ).toBeLowerTriangular();

    expect(ldlt).toBeAllCloseTo(LDLT);  
  });
})
