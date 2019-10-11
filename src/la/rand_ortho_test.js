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
import {matmul2} from './matmul'
import {eye} from './eye'
import {rand_ortho} from './rand_ortho'
import {zip_elems} from '../zip_elems'

import math from '../math'


describe('nd.la.rand_ortho', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
    jasmine.DEFAULT_TIMEOUT_INTERVAL = 60*60*1000;
  })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

      const steps_per_binade = 3;

      for( let run=0*steps_per_binade; run <= 8*steps_per_binade; run++ )
        yield Math.round(2**(run/steps_per_binade));
    }()
  ).it('Generates orthogonal matrices', N => {
    const I = eye(N);

    const E = new Float64Array(N*N),
          s = new Float64Array(N*N);

    for( const Q of [rand_ortho(          N),
                     rand_ortho(          N,N),
                     rand_ortho('float64',N),
                     rand_ortho('float64',N,N)] )
    {
      const Q_T = Q.T
      expect( matmul2(Q,Q_T) ).toBeAllCloseTo(I);
      expect( matmul2(Q_T,Q) ).toBeAllCloseTo(I);
    }
  })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

      const steps_per_binade = 3;

      for( let repeat=0; repeat < 4; repeat++ )
      for( let run=0*steps_per_binade; run <= 3*steps_per_binade; run++ )
      {
        const  N = Math.round(2**(run/steps_per_binade)),
               N_RUNS = randInt(14*1024, 16*1024);
        yield [N,N_RUNS];
      }
    }()
  ).it('Generates evenly distributed matrices', ([N,N_RUNS]) => {
    const I = eye(N);

    const E = new Float64Array(N*N),
          s = new Float64Array(N*N);

    for( let run=0; run < N_RUNS; run++ )
    {
      for( const Q of [rand_ortho(          N),
                       rand_ortho(          N,N),
                       rand_ortho('float64',N),
                       rand_ortho('float64',N,N)] )
      {
        for( let i=N*N; i-- > 0; )
        {
          const   q = Q.data[i];
          E[i] += q;
          s[i] += q*q;
        }
      }
    }

    for( let i=N*N; i-- > 0; )
      E[i] /= (4*N_RUNS);

    for( let i=N*N; i-- > 0; )
      s[i] = s[i] / (4*N_RUNS) - E[i]*E[i];

    const s_mean = s.reduce((x,y) => x+y) / (N*N);
    expect(E).toBeAllCloseTo(0,      {atol:1e-2, rtol:0});
    expect(s).toBeAllCloseTo(s_mean, {atol:1e-2, rtol:0});
  })


  for( const [rndo_name, rndo_method] of Object.entries({
    'rand_ortho(M,N,N)'          : (M,N) => rand_ortho(          M,N,N),
    "rand_ortho('float64',M,N,N)": (M,N) => rand_ortho('float64',M,N,N)
  }))
    forEachItemIn(
      function*(){
        const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from;

        const steps_per_binade = 3;

        for( let run=0*steps_per_binade; run <= 8*steps_per_binade; run++ )
        {
          const  N = Math.round(2**(run/steps_per_binade)),
                 M = randInt(1,32);
          yield [M,N];
        }
      }()
    ).it(`${rndo_name} generates batches of orthogonal matrices`, ([M,N]) => {
      const I = eye(N);

      const Q = rndo_method(M,N),
            Q_T = Q.T;

      expect(Q.dtype).toBe('float64');

      expect( matmul2(Q,Q_T) ).toBeAllCloseTo(I);
      expect( matmul2(Q_T,Q) ).toBeAllCloseTo(I);
    })


  for( const [rndo_name, rndo_method] of Object.entries({
    'rand_ortho(M,N,N)'          : (N_RUNS,N) => rand_ortho(           N_RUNS,N,N),
    "rand_ortho('float64',M,N,N)": (N_RUNS,N) => rand_ortho('float64', N_RUNS,N,N)
  }))
    forEachItemIn(
      function*(){
        const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from;

        const steps_per_binade = 3;

        for( let run=0*steps_per_binade; run <= 3*steps_per_binade; run++ )
        {
          const  N = Math.round(2**(run/steps_per_binade)),
                 N_RUNS = randInt(56*1024, 64*1024);
          yield [N,N_RUNS];
        }
      }()
    ).it(`${rndo_name} generates evenly distributed batches of matrices`, ([N,N_RUNS]) => {
      const I = eye(N);

      const Q = rndo_method(N_RUNS,N),
            Q_T = Q.T;

      expect(Q.dtype).toBe('float64');

      const E =            Q                                     .reduceElems(0, 'float64', (x,y) => x+y).mapElems('float64', x => x/N_RUNS),
            s = zip_elems([Q,E], 'float64', (q,e) => (q-e)*(q-e)).reduceElems(0, 'float64', (x,y) => x+y).mapElems('float64', x => x/N_RUNS);

      const s_mean = s.reduceElems((x,y) => x+y) / (N*N);
      expect(E).toBeAllCloseTo(0,      {atol:1e-2, rtol:0});
      expect(s).toBeAllCloseTo(s_mean, {atol:1e-2, rtol:0});
    })
})
