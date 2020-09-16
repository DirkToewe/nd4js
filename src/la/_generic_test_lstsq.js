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

import {forEachItemIn, CUSTOM_MATCHERS} from '../jasmine_utils';
import {array} from "../nd_array";
import {tabulate} from "../tabulate";
import {zip_elems} from "../zip_elems";

import {generic_test_solve} from "./_generic_test_solve";
import {matmul2} from "./matmul";
import {solve} from "./solve";


export function generic_test_lstsq( lstsq, check_underdet = 'least norm solution' )
{
  if(  check_underdet !== 'no underdet.'
    && check_underdet !== 'some solution'
    && check_underdet !== 'least norm solution' )
    throw new Error('Assertion failed.');

  generic_test_solve(lstsq);

  describe(`${lstsq.name} [generic LSTSQ tests]`, () => {

    beforeEach( () => {
      jasmine.addMatchers(CUSTOM_MATCHERS)
    });


    it('solves hand-crafted over-determined least-squares example', () => {
      for( const dt1 of ['int32', 'float32', 'float64', 'object'] )
      for( const dt2 of ['int32', 'float32', 'float64', 'object'] )
      {
        const A = array(dt1,
          [[1,1],
           [1,2],
           [1,3],
           [1,4],
           [1,5]]
        );
        const y = array(dt2, [[2,3,4,5,6]]).T
        Object.freeze(A.data.buffer);
        Object.freeze(y.data.buffer);
  
        const  x = lstsq(A,y);
        expect(x.dtype).toMatch(/^float/)
        expect(x).toBeAllCloseTo([[1],
                                  [1]]);
      }
    });


    forEachItemIn(
      function*(rng){
        for( let run=512; run-- > 0; )
        {
          let ndim = rng.int(0,3),
            shapes = [ Array.from({length: ndim}, () => rng.int(1,4)) ]
            shapes.splice( rng.int(0,2), 0, shapes[0].slice( rng.int(0,ndim+1) ) )
  
          for( let d=ndim; d > 0; d-- )
          for( let i=rng.int(0,2); i-- > 0; ) { // <- 50/50 chance to collapse dimension to test broadcasting
            const    shape = shapes[rng.int(0,2)],
                 j = shape.length - d
            if(0<=j) shape[j] = 1
          }
  
          let N = rng.int(1,24),
              M = rng.int(N,24);

          const J = rng.int(1,32)
          shapes[0].push(M,N)
          shapes[1].push(M,J)
  
          yield shapes.map(
            s => tabulate(s,'float64', () => rng.uniform(-4,+4))
          )
        }
      }
    ).it('solves random over-determined least-squares examples', ([A,y]) => {
      const x =   lstsq(A,y),
           Ax = matmul2(A,x)
  
      // check least squares solution via normal equaltion Aᵀ(Ax - y) = 0
      expect( matmul2(A.T, zip_elems([Ax,y], (x,y) => x-y) ) ).toBeAllCloseTo(0, {atol: 1e-7})
    });


    function* underdet(rng){
      for( let run=773; run-- > 0; )
      {
        let ndim = rng.int(0,3),
          shapes = [ Array.from({length: ndim}, () => rng.int(1,4)) ]
          shapes.splice( rng.int(0,2), 0, shapes[0].slice( rng.int(0,ndim+1) ) )

        for( let d=ndim; d > 0; d-- )
        for( let i=rng.int(0,2); i-- > 0; ) {
          const     shape = shapes[rng.int(0,2)],
                j = shape.length - d;
          if(0<=j)  shape[j] = 1;
        }

        let   M = rng.int(1,24),
              N = rng.int(M,24)+1;
        const J = rng.int(1,24);
        shapes[0].push(M,N);
        shapes[1].push(M,J);

        yield shapes.map(
          s => tabulate(s,'float64', () => rng.uniform(-4,+4))
        );
      }
    };


    if( check_underdet !== 'no underdet.' )
      forEachItemIn(
        underdet
      ).it('computes      some solution of random under-determined least-squares examples', ([A,y]) => {
        const x = lstsq(A,y);

        // check that it's one solution
        expect( matmul2(A,x) ).toBeAllCloseTo(y);
      });


    if( check_underdet !== 'no underdet.' )
      forEachItemIn(
        function*(rng){
          for( let run=512; run-- > 0; )
          {
            let ndim = rng.int(0,3),
              shapes = [ Array.from({length: ndim}, () => rng.int(1,4)) ]
              shapes.splice( rng.int(0,2), 0, shapes[0].slice( rng.int(0,ndim+1) ) )
    
            for( let d=ndim; d > 0; d-- )
            for( let i=rng.int(0,2); i-- > 0; ) {
              const    shape = shapes[rng.int(0,2)],
                  j = shape.length - d
              if(0<=j) shape[j] = 1
            }
    
            let M = rng.int(1,32),
                N = rng.int(1,32),
                L = rng.int(1,32);
    
            shapes[0].push(M,N)
            shapes[1].push(M,L)
    
            const [A]= rng.rankDef(...shapes[0]),
                  y = tabulate(shapes[1], 'float64', () => rng.uniform(-4,+4));
            Object.freeze(y.data.buffer);
            Object.freeze(y);
    
            yield [A,y];
          }
        }
      ).it('computes      some solution of random  rank-deficient  least-squares examples', ([A,y]) => {
        const x =   lstsq(A,y),
            Ax = matmul2(A,x);
        // check least squares solution via normal equaltion Aᵀ(Ax - y) = 0
        expect( matmul2(A.T, zip_elems([Ax,y], (x,y) => x-y) ) ).toBeAllCloseTo(0, {atol: 1e-6})
      });


    if( check_underdet === 'least norm solution' )
      forEachItemIn(
        underdet
      ).it('computes min.-norm solution of random under-determined least-squares examples', ([A,y]) => {
        const x = lstsq(A,y);

        // check that it's one solution
        expect( matmul2(A,x) ).toBeAllCloseTo(y);

        // check that it's the minimum norm solution
        // https://www.math.usm.edu/lambers/mat419/lecture15.pdf
        const AAT = matmul2(A,A.T);
        expect(x).toBeAllCloseTo( matmul2(A.T, solve(AAT,y)) );
      });

  });
}