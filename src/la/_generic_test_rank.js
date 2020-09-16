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
import {array,
      NDArray} from '../nd_array';

import {FrobeniusNorm} from './norm';
import { permute_rows } from "./permute";


export function generic_test_rank( rank )
{
  describe(`${rank.name} [generic RANK tests]`, () => {

    beforeEach( () => {
      jasmine.addMatchers(CUSTOM_MATCHERS)
    });


    forEachItemIn(
      function*(rng){
        for( let run=773; run-- > 0; )
        {
          const M = rng.int(1,32),
                N = rng.int(1,32),
            shape = Array.from({length: rng.int(0,4) }, () => rng.int(1,4));

          yield rng.rankDef(...shape, M,N);
        }
      }
    ).it(`correctly computes the ranks given random batches rank-deficient examples`, ([A,R]) => {
      const r = rank(A);

      expect(r.shape).toEqual(A.shape.slice(0,-2))
      expect(R.shape).toEqual(A.shape.slice(0,-2))
      expect(r.dtype).toBe('int32')
      expect(R.dtype).toBe('int32')
      expect(r).toBeAllCloseTo(R, {rtol:0, atol:0})
    });


    forEachItemIn(
      function*(rng){
        function* sizes()
        {
          for( let M=0; M++ < 16; )
          for( let N=0; N++ < 16; )
            yield [M,N];

          for( let run=512; run-- > 0; )
            yield [
              rng.int(1,64),
              rng.int(1,64)
            ];
        }

        const NORM = new FrobeniusNorm();

        for( const [M,N] of sizes() )
        {
          const K = Math.max(M,N),
                L = Math.min(M,N),
                R = rng.int(0,L+1);

          const W = new Float64Array(R),
                I =       Int32Array.from({length: K}, (_,i) => i);
          let   A = new Float64Array(K*L);

          for( let i=0; i < R; i++ )
          for( let j=0; j < L; j++ )
            A[L*i+j] = rng.uniform(-4,+4);

          // build remaining/rank-deficients rows as linear combination of the top 
          for( let i=R; i < K; i++ )
          {
            const n = rng.int(0,R+1);

            rng.shuffle( I.subarray(0,R) );
                                     NORM.reset();
            for( let j=n; j-- > 0; ) NORM.include( W[j]  = rng.uniform(-4,+4) );

            const scale = rng.uniform(0,2) / NORM.result;

            for( let j=n; j-- > 0; ) W[j] *= scale;
            
            for( let j=0; j < n; j++ )
            for( let k=0; k < L; k++ )
              A[L*i+k] += W[j] * A[L*I[j]+k];
          }

          A = new NDArray(Int32Array.of(K,L), A);
          rng.shuffle(I);
          A = permute_rows(A,I);

          yield [
            M < N ? A.T : A,
            array(R)
          ];
        }
      }
    ).it(`correctly computes the ranks given random batches rank-deficient matrices`, ([A,R]) => {
      const  r = rank(A);
      expect(r.shape).toEqual(A.shape.slice(0,-2))
      expect(R.shape).toEqual(A.shape.slice(0,-2))
      expect(r.dtype).toBe('int32')
      expect(R.dtype).toBe('int32')
      expect(r).toBeAllCloseTo(R, {rtol:0, atol:0})
    });

  });
}
