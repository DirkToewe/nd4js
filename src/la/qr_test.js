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
import {tabulate} from '../tabulate'
import {_rand_rankdef,
        _rand_rows0,
        _rand_cols0} from '../_test_data_generators'
import {zip_elems} from '../zip_elems'

import {eye} from './eye'
import {matmul2} from './matmul'
import {qr_decomp,
        qr_decomp_full,
       _qr_decomp_inplace,
        qr_lstsq} from './qr'


describe('qr', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from;


  for( const [rng,suffix] of [
    [() =>                           Math.random()*8 - 4, ''                      ],
    [() => Math.random() < 0.1 ? 0 : Math.random()*8 - 4, ' with occasional zeros']
  ])
  forEachItemIn(
    function*(){
      for( let run=1024; run-- > 0; )
      {
        const                       shape = Int32Array.from({ length: randInt(2,5) }, () => randInt(1,24)),
                      QR = tabulate(shape, 'float64', rng);
        Object.freeze(QR.data.buffer);
        yield         QR;
      }
    }()
  ).it('qr_decomp_full works on random examples' + suffix, QR => {
    const [M,N] = QR.shape.slice(-2),
          [Q,R] = qr_decomp_full(QR),
             L  = Math.min(M,N),
            qr  = matmul2(Q,R)
    Object.freeze(Q .data.buffer)
    Object.freeze( R.data.buffer)
    Object.freeze(qr.data.buffer)

    expect( Q.shape ).toEqual( Int32Array.of(...QR.shape.slice(0,-1), M) )
    expect( R.shape ).toEqual( QR.shape )

    expect(R).toBeUpperTriangular()
    expect(qr).toBeAllCloseTo(QR)

    const I = eye(M)
    Object.freeze(I.data.buffer)
    expect( matmul2(Q.T,Q) ).toBeAllCloseTo(I)
    expect( matmul2(Q,Q.T) ).toBeAllCloseTo(I)
  })


  for( const [rng,suffix] of [
    [() =>                           Math.random()*8 - 4, ''                      ],
    [() => Math.random() < 0.1 ? 0 : Math.random()*8 - 4, ' with occasional zeros']
  ])
  forEachItemIn(
    function*(){
      for( let run=1024; run-- > 0; )
      {
        const                       shape = Int32Array.from({ length: randInt(2,5) }, () => randInt(1,32)),
                      QR = tabulate(shape, 'float64', rng);
        Object.freeze(QR.data.buffer);
        yield         QR;
      }
    }()
  ).it('qr_decomp works on random examples' + suffix, QR => {
    const [M,N] = QR.shape.slice(-2),
          [Q,R] = qr_decomp(QR),
             L  = Math.min(M,N),
            qr  = matmul2(Q,R)
    Object.freeze(Q .data.buffer)
    Object.freeze( R.data.buffer)
    Object.freeze(qr.data.buffer)

    expect( Q.shape ).toEqual( Int32Array.of(...QR.shape.slice(0,-2), M, L) )
    expect( R.shape ).toEqual( Int32Array.of(...QR.shape.slice(0,-2), L, N) )

    expect(R).toBeUpperTriangular()
    expect(qr).toBeAllCloseTo(QR, {atol: 1e-7})

    const I = eye(L)
    Object.freeze(I.data.buffer)
                 expect( matmul2(Q.T,Q) ).toBeAllCloseTo(I)
    if( M <= N ) expect( matmul2(Q,Q.T) ).toBeAllCloseTo(I)
  })


  forEachItemIn(
    function*(){
      function* shapes() {
        const steps_per_binade = 4;

        for( let M=8; M > 1; M-- )
        for( let N=8; N > 1; N-- )
          yield [M,N];

        for( let run=1024; run-- > 0; )
        {
          const M = randInt(1,128),
                N = randInt(1,128);
          yield [M,N];
        }
      }

      for( const [M,N] of shapes() )
        yield _rand_rows0(M,N);
    }()
  ).it('qr_decomp works on random matrices with zero rows', QR => {
    const [M,N] = QR.shape.slice(-2),
          [Q,R] = qr_decomp(QR),
             L  = Math.min(M,N),
            qr  = matmul2(Q,R)
    Object.freeze(Q .data.buffer)
    Object.freeze( R.data.buffer)
    Object.freeze(qr.data.buffer)

    expect( Q.shape ).toEqual( Int32Array.of(...QR.shape.slice(0,-2), M, L) )
    expect( R.shape ).toEqual( Int32Array.of(...QR.shape.slice(0,-2), L, N) )

    expect(R).toBeUpperTriangular()
    expect(qr).toBeAllCloseTo(QR)

    const I = eye(L)
    Object.freeze(I.data.buffer)
                 expect( matmul2(Q.T,Q) ).toBeAllCloseTo(I)
    if( M <= N ) expect( matmul2(Q,Q.T) ).toBeAllCloseTo(I)
  })


  forEachItemIn(
    function*(){
      function* shapes() {
        const steps_per_binade = 4;

        for( let M=8; M > 1; M-- )
        for( let N=8; N > 1; N-- )
          yield [M,N];

        for( let run=1024; run-- > 0; )
        {
          const M = randInt(1,128),
                N = randInt(1,128);
          yield [M,N];
        }
      }

      for( const [M,N] of shapes() )
        yield _rand_cols0(M,N);
    }()
  ).it('qr_decomp works on random matrices with zero columns', QR => {
    const [M,N] = QR.shape.slice(-2),
          [Q,R] = qr_decomp(QR),
             L  = Math.min(M,N),
            qr  = matmul2(Q,R)
    Object.freeze(Q .data.buffer)
    Object.freeze( R.data.buffer)
    Object.freeze(qr.data.buffer)

    expect( Q.shape ).toEqual( Int32Array.of(...QR.shape.slice(0,-2), M, L) )
    expect( R.shape ).toEqual( Int32Array.of(...QR.shape.slice(0,-2), L, N) )

    expect(R).toBeUpperTriangular()
    expect(qr).toBeAllCloseTo(QR)

    const I = eye(L)
    Object.freeze(I.data.buffer)
                 expect( matmul2(Q.T,Q) ).toBeAllCloseTo(I)
    if( M <= N ) expect( matmul2(Q,Q.T) ).toBeAllCloseTo(I)
  })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

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
  ).it('qr_lstsq solves random square examples', ([QR,y]) => {
    const [Q,R]= qr_decomp(QR),
             x = qr_lstsq(Q,R, y)
  
    expect( matmul2(QR,x) ).toBeAllCloseTo(y)
  })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

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
        if( M < N ) {
          const L=M; M=N; N=L
        }

        const J = randInt(1,32)
        shapes[0].push(M,N)
        shapes[1].push(M,J)

        yield shapes.map(
          s => tabulate(s,'float64', () => Math.random()*2-1)
        )
      }
    }()
  ).it('qr_lstsq+qr_decomp solves least squares of random over-determined examples', ([A,y]) => {
    const [Q,R]= qr_decomp(A),
             x = qr_lstsq(Q,R, y),
            Ax = matmul2(A,x)

    // every least square solution satisfies the normal equaltion Aᵀ(Ax - y) = 0
    expect( matmul2(A.T, zip_elems([Ax,y], math.sub) ) ).toBeAllCloseTo(0)
  })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

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
        if( M < N ) {
          const L=M; M=N; N=L
        }

        const J = randInt(1,32)
        shapes[0].push(M,N)
        shapes[1].push(M,J)

        yield shapes.map(
          s => tabulate(s,'float64', () => Math.random()*2-1)
        )
      }
    }()
  ).it('qr_lstsq+qr_decomp_full solve least squares of random over-determined examples', ([A,y]) => {
    const [Q,R]= qr_decomp_full(A),
             x = qr_lstsq(Q,R, y),
            Ax = matmul2(A,x)

    // every least square solution satisfies the normal equaltion Aᵀ(Ax - y) = 0
    expect( matmul2(A.T, zip_elems([Ax,y], math.sub) ) ).toBeAllCloseTo(0)
  })


  forEachItemIn(
    function*(){
      for( let run=4096; run-- > 0; )
      {
        const M = randInt(1,64),
              N = randInt(1,64),
              L = randInt(1,64),
             [A]= _rand_rankdef(M,N),
              Y = tabulate([M,L], 'float64', () => (Math.random() < 0.01)*(Math.random()*8-4) );
        Object.freeze(Y.data.buffer);
        Object.freeze(Y);
        yield [A,Y];
      }
    }()
  ).it(`_qr_decomp_inplace works on random rank-deficient examples`, ([A,Y]) => {
    const [M,N] = A.shape,
             L  = Y.shape[1];

    const [Q,R] = qr_decomp_full(A);

    const a = A.mapElems(),
          y = Y.mapElems();

    _qr_decomp_inplace(M,N,L, a.data,0, y.data,0); // <- TODO test with offset

    expect(a).toBeUpperTriangular();
    expect(a).toBeAllCloseTo(R);
    expect(y).toBeAllCloseTo( matmul2(Q.T,Y) );
  })
})
