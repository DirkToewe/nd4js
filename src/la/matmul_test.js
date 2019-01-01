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
import {array} from '../nd_array'
import {zip_elems} from '../zip_elems'
import {tabulate} from '../tabulate'


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

  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from

      for( let run=1024; run-- > 0; )
      {
        let cShape = Int32Array.from({ length: randInt(2,8) }, () => randInt(1,8) ),
            aShape = cShape.slice(),
            bShape = cShape.slice( randInt(0,cShape.length-2) );
    
        for( let a=aShape.length-2, b=bShape.length-2; a-- > 0 && b-- > 0; )
          switch( randInt(0,3) )
          {
            case 0: break;
            case 1: aShape[a] = 1; break;
            case 2: bShape[b] = 1; break;
          }
    
        if( Math.random() < 0.5 ) {
          const tmp = aShape; aShape = bShape; bShape = tmp;
        }
        aShape[aShape.length-2] = randInt(1,5)
        aShape[aShape.length-1] = bShape[bShape.length-2]
        bShape[bShape.length-1] = randInt(1,5)
    
        yield [
          tabulate(aShape, 'float64', () => Math.random()*2 - 1 ),
          tabulate(bShape, 'float64', () => Math.random()*2 - 1 )
        ]
      }
    }()
  ).it('matmul2 works on random examples', ([a,b]) => {
    const c = matmul2(a,b),
          C = zip_elems([
            a.sliceElems('...','new'),
            b.sliceElems('...','new',[],[])
          ], 'float64', (x,y) => x*y ).reduceElems([-2], (x,y) => x+y )

    expect(c.shape).toEqual(C.shape)
    expect(c.dtype).toBe('float64')
    expect(c).toBeAllCloseTo(C)
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
})
