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
import math from '../math'
import {NDArray} from '../nd_array'
import {tabulate} from '../tabulate'
import {zip_elems} from '../zip_elems'

import {eye} from './eye'
import {matmul,
        matmul2} from './matmul'
import {permute_cols} from './permute'
import {qr_decomp} from './qr'
import {rand_ortho} from './rand_ortho'
import { rrqr_decomp,
         rrqr_decomp_full,
         rrqr_rank,
         rrqr_solve,
         rrqr_lstsq}       from  './rrqr'
import {srrqr_decomp_full} from './srrqr'


const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from


function _rand_rank_deficient(...shape)
{
  if( shape.length < 2 )
    throw new Error('Assertion failed.');

  const [M,N]= shape.slice(-2),
           L = Math.min(M,N);

  const ndim = randInt(0,7);

  let r_shape = shape.slice(0,-2);

  // use random, mostly rank-deficient SVD to generate test matrix A
  const ranks = tabulate(r_shape, 'int32', () => randInt(0,L+1) ), // <- ranks
            U = rand_ortho('float64', ...r_shape,M,L),
            V = rand_ortho('float64', ...r_shape,L,N),
            S = tabulate([...r_shape,L,L], 'float64', (...idx) => {
              const j = idx.pop(),
                    i = idx.pop(), rank = ranks(...idx);

              if( i !== j || rank <= i ) return 0; 

              return Math.random()*8 + 0.1
            }),
            A = matmul(U,S,V);

  // add random scaling
  for( let S_off=S.data.length; (S_off -= L*L) >= 0; )
  {
    const scale = Math.random()*1e300 + 1;

    for( let i=L; i-- > 0; )
      S.data[S_off + L*i+i] *= scale;
  }

  Object.freeze(ranks.data.buffer);
  Object.freeze(A.data.buffer);
  return [ranks,A];
}


describe('srrqr', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  const SRRQR_METHODS = Object.entries({
    srrqr_decomp_full
  });
  Object.freeze(SRRQR_METHODS);

  for( const [srrqr_name,srrqr_deco] of SRRQR_METHODS )

  forEachItemIn(
    function*(){
      for( let run=1024; run-- > 0; )
      {
        const ndim = randInt(0,5);

        yield _rand_rank_deficient(
          ...Array.from({ length: ndim-2 }, () => randInt(1,8) ),
          randInt(1,64),
          randInt(1,64)
        );
      }
    }()
  ).it(`${srrqr_name} computes rank correctly for random examples`, ([RANKS,A]) => {
//    console.log(`A.shape: [${A.shape}].`);

    const [Q,R,P, ranks] = srrqr_deco(A);

    expect(ranks.shape).toEqual(A.shape.slice(0,-2));
    expect(RANKS.shape).toEqual(A.shape.slice(0,-2));
    expect(ranks.dtype).toBe('int32');
    expect(RANKS.dtype).toBe('int32');
    expect(ranks).toBeAllCloseTo(RANKS, {rtol:0, atol:0});
  })
})


describe('(s)rrqr', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  const RRQR_METHODS = Object.entries({
    srrqr_decomp_full,
     rrqr_decomp,
     rrqr_decomp_full,
  });
  Object.freeze(RRQR_METHODS);


  for( const [rrqr_name,rrqr_deco] of RRQR_METHODS )
    forEachItemIn(
      function*(){
        for( let run=1024; run-- > 0; )
        {
          const ndim = randInt(2,5),
              shape = Int32Array.from({ length: ndim }, () => randInt(1,24) )
          const A = tabulate(shape, 'float64', () => Math.random() < 0.1 ? 0 : Math.random()*2-1)
          Object.freeze(A.data.buffer)
          yield A
        }
      }()
    ).it(`${rrqr_name} works on random examples`, A => {
      const [M,N] = A.shape.slice(-2),
          [Q,R,P] = rrqr_deco(A),
               L  = rrqr_name.endsWith('_full') ? M : Math.min(M,N),
               a  = matmul2(Q,R)
      Object.freeze(Q.data.buffer); expect(Q.dtype).toBe('float64')
      Object.freeze(R.data.buffer); expect(R.dtype).toBe('float64')
      Object.freeze(P.data.buffer); expect(P.dtype).toBe(  'int32')
      Object.freeze(a.data.buffer)

      expect( Q.shape ).toEqual( Int32Array.of(...A.shape.slice(0,-2), M, L) )
      expect( R.shape ).toEqual( Int32Array.of(...A.shape.slice(0,-2), L, N) )
      expect( P.shape ).toEqual( Int32Array.of(...A.shape.slice(0,-2),    N) )

      A = permute_cols(A,P);

      expect(R).toBeUpperTriangular()
      expect(a).toBeAllCloseTo(A)

      const I = eye(L)
      Object.freeze(I.data.buffer)
                  expect( matmul2(Q.T,Q) ).toBeAllCloseTo(I)
      if( L===M ) expect( matmul2(Q,Q.T) ).toBeAllCloseTo(I)

      if( rrqr_name.startsWith('rrqr_') )
      {
        const RR = R.data;
        for( let off=RR.length; (off -= L*N) >= 0; )
          for( let i=Math.min(M,N); --i > 0; )
          {
            const  R_II = math.abs(RR[off + N* i   + i   ]),
                   R_ii = math.abs(RR[off + N*(i-1)+(i-1)])
            expect(R_ii).not.toBeLessThan(R_II)
          }
      }
    })


  for( const [rrqr_name,rrqr_deco] of RRQR_METHODS )
    forEachItemIn(
      function*(){
        for( let run=1024; run-- > 0; )
        {
          let ndim = randInt(0,5),
            shapes = [ Array.from({length: ndim}, () => randInt(1,8)) ]
          shapes.splice( randInt(0,2), 0, shapes[0].slice( randInt(0,ndim) ) )

          for( let d=ndim; d > 0; d-- )
          for( let i=randInt(0,2); i-- > 0; ) {
            const    shape = shapes[randInt(0,2)],
                j = shape.length - d
            if(0<=j) shape[j] = 1
          }

          const M = randInt(1,32); shapes[0].push(M,M)
          const N = randInt(1,32); shapes[1].push(M,N)

          yield shapes.map(
            s => tabulate(s,'float64', () => Math.random()*2-1)
          )
        }
      }()
    ).it(`${rrqr_name}+rrqr_solve solves random square examples`, ([QR,y]) => {
      const [Q,R,P] = rrqr_deco(QR),
               x    = rrqr_solve(Q,R,P, y)
    
      expect( matmul2(QR,x) ).toBeAllCloseTo(y)
    })


  for( const [rrqr_name,rrqr_deco] of RRQR_METHODS )
    forEachItemIn(
      function*(){
        for( let run=1024; run-- > 0; )
        {
          let ndim = randInt(0,5),
            shapes = [ Array.from({length: ndim}, () => randInt(1,8)) ]
          shapes.splice( randInt(0,2), 0, shapes[0].slice( randInt(0,ndim) ) )

          for( let d=ndim; d > 0; d-- )
          for( let i=randInt(0,2); i-- > 0; ) {
            const    shape = shapes[randInt(0,2)],
                j = shape.length - d
            if(0<=j) shape[j] = 1
          }

          const M = randInt(1,32); shapes[0].push(M,M)
          const N = randInt(1,32); shapes[1].push(M,N)

          yield shapes.map(
            s => tabulate(s,'float64', () => Math.random()*2-1)
          )
        }
      }()
    ).it(`${rrqr_name}+rrqr_lstsq solves random square examples`, ([QR,y]) => {
      const [Q,R,P] = rrqr_deco(QR),
               x    = rrqr_lstsq(Q,R,P, y)
    
      expect( matmul2(QR,x) ).toBeAllCloseTo(y)
    })


  for( const [rrqr_name,rrqr_deco] of RRQR_METHODS )
    forEachItemIn(
      function*(){
        for( let run=1024; run-- > 0; )
        {
          let ndim = randInt(0,5),
            shapes = [ Array.from({length: ndim}, () => randInt(1,8)) ]
          shapes.splice( randInt(0,2), 0, shapes[0].slice( randInt(0,ndim) ) )

          for( let d=ndim; d > 0; d-- )
          for( let i=randInt(0,2); i-- > 0; ) { // <- 50/50 chance to collapse dimension to test broadcasting
            const    shape = shapes[randInt(0,2)],
                j = shape.length - d
            if(0<=j) shape[j] = 1
          }

          let M = randInt(1,32),
              N = randInt(1,32)
          if( M < N )
              [M,N] = [N,M];

          if( M < N )
            throw new Error('Assertion failed.');

          const J = randInt(1,32)
          shapes[0].push(M,N)
          shapes[1].push(M,J)

          yield shapes.map(
            s => tabulate(s,'float64', () => Math.random()*4-2)
          )
        }
      }()
    ).it(`${rrqr_name}+rrqr_lstsq solves least squares of random over-determined examples`, ([A,y]) => {
      const [Q,R,P]= rrqr_deco(A),
              x    = rrqr_lstsq(Q,R,P, y),
              Ax   = matmul2(A,x)

      // every least square solution satisfies the normal equaltion Aᵀ(Ax - y) = 0
      expect( matmul2(A.T, zip_elems([Ax,y], math.sub) ) ).toBeAllCloseTo(0)
    })


  for( const [rrqr_name,rrqr_deco] of RRQR_METHODS )
    forEachItemIn(
      function*(){
        for( let run=1024; run-- > 0; )
        {
          let ndim = randInt(0,5),
            shapes = [ Array.from({length: ndim}, () => randInt(1,8)) ]
          shapes.splice( randInt(0,2), 0, shapes[0].slice( randInt(0,ndim) ) )

          for( let d=ndim; d > 0; d-- )
          for( let i=randInt(0,2); i-- > 0; ) {
            const   shape = shapes[randInt(0,2)],
                  j=shape.length - d
            if(0<=j)shape[j] = 1
          }

          let M = randInt(1,32),
              N = randInt(1,32)
          if( M > N )
            [M,N] = [N,M];

          if( M > N )
            throw new Error('Assertion failed.');

          const J = randInt(1,32)
          shapes[0].push(M,N)
          shapes[1].push(M,J)

          yield shapes.map(
            s => tabulate(s,'float64', () => Math.random()*4-2)
          )
        }
      }()
    ).it(`${rrqr_name}+rrqr_lstsq computes one solution of random under-determined examples`, ([A,y]) => {
      const [Q,R,P]= rrqr_deco(A),
               x   = rrqr_lstsq(Q,R,P, y)

      expect( matmul2(A,x) ).toBeAllCloseTo(y)
    })


  for( const [rrqr_name,rrqr_deco] of RRQR_METHODS )
    forEachItemIn(
      function*(){
        for( let run=1024; run-- > 0; )
        {
          const ndim = randInt(0,5),
              shapes = [ Array.from({length: ndim}, () => randInt(1,8)) ]

          // add some broadcasting test cases (A.ndim != y.ndim)
          shapes.splice( randInt(0,2), 0, shapes[0].slice( randInt(0,ndim) ) )

          // add some broadcasting test cases (A.shape[i] = 1 || y.shape[i] = 1)
          for( let d=ndim; d > 0; d-- )
          for( let i=randInt(0,2); i-- > 0; ) {
            const   shape = shapes[randInt(0,2)],
                  j=shape.length - d
            if(0<=j)shape[j] = 1
          }

          const M = randInt(1,32),
                N = randInt(1,32);

          const J = randInt(1,32);
          shapes[0].push(M,N);
          shapes[1].push(M,J);

          yield [
            _rand_rank_deficient(...shapes[0])[1],
            tabulate(shapes[1], 'float64', () => Math.random()*4-2)
          ];
        }
      }()
    ).it(`${rrqr_name}+rrqr_lstsq solves least squares for random rank-deficient examples`, ([A,y]) => {
  //    console.log('SHAPES:', [...A.shape], [...y.shape])
      const [Q,R,P]= rrqr_deco(A),
               x   = rrqr_lstsq(Q,R,P, y),
              Ax   = matmul2(A,x)

      // every least square solution satisfies the normal equaltion Aᵀ(Ax - y) = 0
      expect( matmul2(A.T, zip_elems([Ax,y], (Ax,y) => Ax-y) ) ).toBeAllCloseTo(0)
    })


  for( const [rrqr_name,rrqr_deco] of RRQR_METHODS )
    forEachItemIn(
      function*(){
        for( let run=1024; run-- > 0; )
        {
          const ndim = randInt(0,5);

          yield _rand_rank_deficient(
            ...Array.from({ length: ndim-2 }, () => randInt(1,8) ),
            randInt(1,32),
            randInt(1,32)
          );
        }
      }()
    ).it(`${rrqr_name}+rrqr_rank works on random examples`, ([RANKS,A]) => {
      const [Q,R,P] = rrqr_deco(A),
             ranks  = rrqr_rank(R)

      expect(ranks.shape).toEqual(A.shape.slice(0,-2))
      expect(RANKS.shape).toEqual(A.shape.slice(0,-2))
      expect(ranks.dtype).toBe('int32')
      expect(RANKS.dtype).toBe('int32')
      expect(ranks).toBeAllCloseTo(RANKS, {rtol:0, atol:0})
    })
})
