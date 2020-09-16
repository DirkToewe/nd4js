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
import {matmul, matmul2} from './matmul'
import {array, NDArray} from '../nd_array'
import {zip_elems} from '../zip_elems'
import {tabulate} from '../tabulate'
import {Complex128Array} from '../dt/complex_array'
import math from '../math'


describe('matmul', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  it('matmul2 works on example factors of shape [2,1] and [1,3]', () => {
    const
      a = array([[1],
                 [2]]),
      b = array([[30,40,50]]),
      c = matmul2(a,b)
    expect(c.shape).toEqual( Int32Array.of(2,3) )
    expect(c.data ).toEqual( Int32Array.of(
      1*30, 1*40, 1*50,
      2*30, 2*40, 2*50
    ))
  })

  it('matmul2 works on example factors of shape [2,3] and [3,2]', () => {
    const
      a = array([
        [1,2,3],
        [4,5,6]
      ]),
      b = array([
        [ 70, 80],
        [ 90,100],
        [110,120]
      ]),
      c = matmul2(a,b)
    expect(c.shape).toEqual( Int32Array.of(2,2) )
    expect(c.data ).toEqual( Int32Array.of(
      (1*70 + 2*90 + 3*110), (1*80+2*100+3*120),
      (4*70 + 5*90 + 6*110), (4*80+5*100+6*120)
    ))
  })

  it('matmul works on example triple of shapes [1,4], [4,3], [3,2]', () => {
    const a = array([[1,2,3,4]]),
          b = array([[11,12,13],
                     [21,22,23],
                     [31,32,33],
                     [41,42,43]]),
          c = array([[5, 6],
                     [7, 8],
                     [9,10]]),
          abc = matmul(a,b,c),
          ABC = array([[6760, 7720]])

     expect(abc.shape).toEqual(ABC.shape)
     expect(abc).toBeAllCloseTo(ABC)
  })

  for( const dtype1 of ['int32', 'float64', 'complex128'] )
  for( const dtype2 of ['int32', 'float64', 'complex128'] )
    forEachItemIn(
      function*(){
        const randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from
  
        for( let run=373; run-- > 0; )
        {
          let cShape = Int32Array.from({ length: randInt(2,6) }, () => randInt(1,4) ),
              aShape = cShape.slice(),
              bShape = cShape.slice( randInt(0,cShape.length-1) );
      
          for( let a=aShape.length-2, b=bShape.length-2; a-- > 0 && b-- > 0; )
            switch( randInt(0,3) )
            {
              case 0: break;
              case 1: aShape[a] = 1; break;
              case 2: bShape[b] = 1; break;
            }
      
          if( Math.random() < 0.5 ) {
            const tmp = aShape;
                        aShape = bShape;
                                 bShape = tmp;
          };
          aShape[aShape.length-2] = randInt(1,16)
          aShape[aShape.length-1] =
          bShape[bShape.length-2] = randInt(1,16);
          bShape[bShape.length-1] = randInt(1,16);

          let a = ( dtype1==='int32' ? Int32Array : Float64Array ).from( {length: aShape.reduce(math.mul,1)*(1 + (dtype1==='complex128'))}, () => Math.random()*4-2 ),
              b = ( dtype2==='int32' ? Int32Array : Float64Array ).from( {length: bShape.reduce(math.mul,1)*(1 + (dtype2==='complex128'))}, () => Math.random()*4-2 );
          if( dtype1==='complex128' ) a = new Complex128Array(a.buffer);
          if( dtype2==='complex128' ) b = new Complex128Array(b.buffer);
          yield [
            new NDArray(aShape, a),
            new NDArray(bShape, b)
          ];
        }
      }()
    ).it(`matmul2 works on random examples with dtypes ${dtype1} and ${dtype2}`, ([a,b]) => {
      expect(a.dtype).toBe(dtype1)
      expect(b.dtype).toBe(dtype2)
      const dtype3 = function(){
              if( dtype1 === 'complex128' ) return 'complex128';
              if( dtype2 === 'complex128' ) return 'complex128';
              if( dtype1 !==      'int32' ) return    'float64';
              if( dtype2 !==      'int32' ) return    'float64';
              return              'int32';
            }(),
            c = matmul2(a,b),
            C = zip_elems([
              a.sliceElems('...','new'),
              b.sliceElems('...','new',[],[])
            ], dtype3, math.mul).reduceElems([-2], math.add)
  
      expect(c.shape).toEqual(C.shape)
      expect(c.dtype).toBe(dtype3)
      expect(c).toBeAllCloseTo(C)
    })
})
