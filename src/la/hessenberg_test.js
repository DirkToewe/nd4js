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
import {tabulate} from '../tabulate'
import {hessenberg_decomp} from './hessenberg'
import {matmul} from './matmul'
import {eye} from './eye'


describe('hessenberg', () => {
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
      function* shapes()
      {
        for( let run = 0; run < 16; run++ )
        for( let N=1; N < 16; N++ )
          yield [N,N];

        for( let run=1024; run-- > 0; )
        {
          const N = randInt(1,64),
             ndim = randInt(2,5),
            shape = Array.from({ length: ndim-2 }, () => randInt(1,8) );
          yield [...shape, N,N];
        }
      }

      for( const shape of shapes() )
      {
        const A = tabulate(shape, 'float64', rng);
        Object.freeze(A.data.buffer);
        yield A;
      }
    }()
  ).it('hessenberg_decomp works on random examples' + suffix, A => {
    const [N]= A.shape.slice(-1),
        [U,H]= hessenberg_decomp(A),
           a = matmul(U,H,U.T)

    expect(U.shape).toEqual(A.shape)
    expect(H.shape).toEqual(A.shape)

    const I = eye(N)
    expect( matmul(U,U.T) ).toBeAllCloseTo(I)
    expect( matmul(U.T,U) ).toBeAllCloseTo(I)
    expect(a).toBeAllCloseTo(A)

    expect(H).toBeUpperHessenberg()
  })


  forEachItemIn(
    function*(){
      function* shapes()
      {
        for( let run = 0; run < 16; run++ )
        for( let N=1; N < 16; N++ )
          yield [N,N];

        for( let run=1024; run-- > 0; )
        {
          const N = randInt(1,64),
             ndim = randInt(2,5),
            shape = Array.from({ length: ndim-2 }, () => randInt(1,8) );
          yield [...shape, N,N];
        }
      }

      for( const shape of shapes() )
      {
        const                         sparseness = Math.random(),
          rng = () => Math.random() < sparseness ? 0 : Math.random()*2 - 1;

        const A = tabulate(shape, 'float64', rng);
        Object.freeze(A.data.buffer);
        yield A;
      }
    }()
  ).it('hessenberg_decomp works on random sparse examples', A => {
    const [N]= A.shape.slice(-1),
        [U,H]= hessenberg_decomp(A),
           a = matmul(U,H,U.T)

    expect(U.shape).toEqual(A.shape)
    expect(H.shape).toEqual(A.shape)

    const I = eye(N)
    expect( matmul(U,U.T) ).toBeAllCloseTo(I)
    expect( matmul(U.T,U) ).toBeAllCloseTo(I)
    expect(a).toBeAllCloseTo(A)

    expect(H).toBeUpperHessenberg()
  })
})
