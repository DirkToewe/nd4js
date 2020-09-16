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
import {generic_test_lstsq} from "./_generic_test_lstsq";
import {tabulate} from '../tabulate'
import {_rand_rankdef,
        _rand_rows0,
        _rand_cols0} from '../_test_data_generators'

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

  // TODO: test qr_lstsq without relying on qr_decomp(_full)

});


for( const lstsq of Object.values({
  [`qr_decomp_full + qr_lstsq   `]: (A,y) => qr_lstsq(   qr_decomp_full(A), y ),
  [`qr_decomp_full + qr_lstsq...`]: (A,y) => qr_lstsq(...qr_decomp_full(A), y ),
  [`qr_decomp      + qr_lstsq   `]: (A,y) => qr_lstsq(   qr_decomp     (A), y ),
  [`qr_decomp      + qr_lstsq...`]: (A,y) => qr_lstsq(...qr_decomp     (A), y )
}) )
  generic_test_lstsq(lstsq, 'no underdet.');


describe('qr', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  // TODO: test qr_lstsq without relying on qr_decomp(_full)

});


const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from;


const test_qr_decomp = (decomp_name, test_body) => describe(decomp_name, () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  for( const [rng,suffix] of [
    [() =>                           Math.random()*8 - 4, ''                      ],
    [() => Math.random() < 0.1 ? 0 : Math.random()*8 - 4, ' with occasional zeros']
  ])
    forEachItemIn(
      function*(){
        for( let run=256; run-- > 0; )
        { const                       shape = Int32Array.from({ length: randInt(2,5) }, () => randInt(1,24)),
                        QR = tabulate(shape, 'float64', rng);
          Object.freeze(QR.data.buffer);
          yield         QR;
        }
      }()
    ).it(`correctly decomposes random examples${suffix}`, test_body);


  forEachItemIn(
    function*(){
      function* shapes() {
        const steps_per_binade = 4;

        for( let M=8; M > 1; M-- )
        for( let N=8; N > 1; N-- )
          yield [M,N];

        for( let run=256; run-- > 0; )
        { const  M = randInt(1,128),
                 N = randInt(1,128);
          yield [M,N];
        }
      }

      for( const [M,N] of shapes() )
        yield _rand_rows0(M,N);
    }()
  ).it(`correctly decomposes random matrices with zero rows`, test_body);


  forEachItemIn(
    function*(){
      function* shapes() {
        const steps_per_binade = 4;

        for( let M=8; M > 1; M-- )
        for( let N=8; N > 1; N-- )
          yield [M,N];

        for( let run=256; run-- > 0; )
        {
          const M = randInt(1,128),
                N = randInt(1,128);
          yield [M,N];
        }
      }

      for( const [M,N] of shapes() )
        yield _rand_cols0(M,N);
    }()
  ).it(`correctly decomposes random matrices with zero columns`, test_body);


  forEachItemIn(
    function*(){
      for( let run=256; run-- > 0; )
      {
        const                                shape = Int32Array.from({ length: randInt(2,5) }, () => randInt(1,24)),
                     [QR] = _rand_rankdef(...shape);
        Object.freeze(QR.data.buffer);
        yield         QR;
      }
    }()
  ).it(`correctly decomposes random rank-deficient examples`, test_body);

});


test_qr_decomp('qr_decomp_full', QR => {
  const [M,N] = QR.shape.slice(-2),
        [Q,R] = qr_decomp_full(QR),
          qr  = matmul2(Q,R)
  Object.freeze(Q .data.buffer)
  Object.freeze( R.data.buffer)
  Object.freeze(qr.data.buffer)

  expect( Q.shape ).toEqual( Int32Array.of(...QR.shape.slice(0,-1), M) )
  expect( R.shape ).toEqual( QR.shape )

  expect(R).toBeUpperTriangular()
  expect(qr).toBeAllCloseTo(QR)

  const                                   I = eye(M)
  expect( matmul2(Q.T,Q) ).toBeAllCloseTo(I)
  expect( matmul2(Q,Q.T) ).toBeAllCloseTo(I)
});


test_qr_decomp('qr_decomp', QR => {
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

  const                                                I = eye(L)
               expect( matmul2(Q.T,Q) ).toBeAllCloseTo(I)
  if( M <= N ) expect( matmul2(Q,Q.T) ).toBeAllCloseTo(I)
});


describe('_qr_decomp_inplace', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS);
  });


  function* test_data( rand_matrix )
  {
    for( let run=1024; run-- > 0; )
    {
      const M = randInt(1,64),
            N = randInt(1,64),
            L = randInt(1,64),
            A = rand_matrix(M,N),
            Y = tabulate([M,L], 'float64', () => (Math.random()*8-4) * (Math.random() < 0.9) );
      Object.freeze(Y.data.buffer);
      Object.freeze(Y);
      yield [A,Y];
    }
  }


  const test_body = ([A,Y]) => {
    const [M,N] = A.shape,
             L  = Y.shape[1];

    const [Q,R] = qr_decomp_full(A);

    const a = A.mapElems(),
          y = Y.mapElems();

    _qr_decomp_inplace(M,N,L, a.data,0, y.data,0); // <- TODO test with offset

    expect(a).toBeUpperTriangular();
    expect(a).toBeAllCloseTo(R);
    expect(y).toBeAllCloseTo( matmul2(Q.T,Y) );
  };


  for( const [rng,suffix] of [
    [() =>                           Math.random()*8 - 4, ''                      ],
    [() => Math.random() < 0.1 ? 0 : Math.random()*8 - 4, ' with occasional zeros']
  ])
    forEachItemIn(
      test_data( (M,N) => tabulate([M,N], 'float64', rng ) )
    ).it(`correctly decomposes random matrices${suffix}`, test_body);


  forEachItemIn(
    test_data(_rand_rows0)
  ).it(`correctly decomposes random matrices with zero rows`, test_body);


  forEachItemIn(
    test_data(_rand_cols0)
  ).it(`correctly decomposes random matrices with zero columns`, test_body);


  forEachItemIn(
    test_data( (M,N) => _rand_rankdef(M,N)[0] )
  ).it(`correctly decomposes random rank-deficient matrices`, test_body);
});
