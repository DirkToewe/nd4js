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
import {num_grad,
        num_grad_forward} from './num_grad'
import {array} from '../nd_array'


describe('num_grad_forward', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  forEachItemIn(
    function*(){
      for( let run=1024; run-- > 0; )
      {
        const [a,b,c] = function*(){ while(true) yield Math.random()*2-1 }(),
              F = (x,y) => a + (b + c*x)*x,
              G = (x,y) =>      b + c*x*2
        yield [F,G]
      }
    }()
  ).it('works on random 1d polynomials of degree 2', ([F,G]) => {
    const g = num_grad_forward( xy => F(...xy.data) )

    for( let x=-3.14; x <=3.14; x+=0.42 )
      expect(g([x])).toBeAllCloseTo(G(x), {atol:1e-7})
  })


  forEachItemIn(
    function*(){
      for( let run=1024; run-- > 0; )
      {
        const [a,b,c,d,e,f] = function*(){ while(true) yield Math.random()*2-1 }(),
              F = (x,y) => a + (b + d*y + e*x)*x + (c + f*y)*y,
              G = (x,y) => array([
                b + d*y + 2*e*x,
                c + d*x + 2*f*y
              ])
        yield [F,G]
      }
    }()
  ).it('works on random 2d polynomials of degree 2', ([F,G]) => {
    const g = num_grad_forward( xy => F(...xy.data) )

    for( let x=-3.14; x <=3.14; x+=0.42 )
    for( let y=-3.14; y <=3.14; y+=0.42 )
      expect(g([x,y])).toBeAllCloseTo(G(x,y), {atol:1e-7})
  })


  forEachItemIn(
    function*(){
      for( let run=128; run-- > 0; )
      {
        const a = Array.from({length: 10}, () => Math.random()*2-1),
              F = (x,y,z) => {
                let result = 0
                for( let i=0; i<=2; i++ )
                for( let j=0; j<=2; j++ )
                for( let k=0; k<=2; k++ )
                  if( i+j+k <= 2 )
                    result += a[4*i+2*j+k] * x**i * y**j * z**k
                return result
              },
              G = (x,y,z) => {
                let result = [0,0,0]
                for( let i=0; i<=2; i++ )
                for( let j=0; j<=2; j++ )
                for( let k=0; k<=2; k++ )
                  if( i+j+k <= 2 ) {
                    if(i > 0) result[0] += a[4*i+2*j+k] * i * x**(i-1) * y** j    * z** k
                    if(j > 0) result[1] += a[4*i+2*j+k] * j * x** i    * y**(j-1) * z** k
                    if(k > 0) result[2] += a[4*i+2*j+k] * k * x** i    * y** j    * z**(k-1)
                  }
                return array(result)
              }
        yield [F,G]
      }
    }()
  ).it('works on random 3d polynomials of degree 2', ([F,G]) => {
    const g = num_grad_forward( xy => F(...xy.data) )

    for( let x=-3.14; x <=3.14; x+=0.42 )
    for( let y=-3.14; y <=3.14; y+=0.42 )
    for( let z=-3.14; z <=3.14; z+=0.42 )
      expect(g([x,y,z])).toBeAllCloseTo(G(x,y,z), {atol:1e-7})
  })
})


describe('num_grad', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  forEachItemIn(
    function*(){
      for( let run=1024; run-- > 0; )
      {
        const [a,b,c] = function*(){ while(true) yield Math.random()*2-1 }(),
              F = (x,y) => a + (b + c*x)*x,
              G = (x,y) =>      b + c*x*2
        yield [F,G]
      }
    }()
  ).it('works on random 1d polynomials of degree 2', ([F,G]) => {
    const g = num_grad( xy => F(...xy.data) )

    for( let x=-3.14; x <=3.14; x+=0.42 )
      expect(g([x])).toBeAllCloseTo(G(x))//, {rtol:1e-4, atol:1e-6})
  })


  forEachItemIn(
    function*(){
      for( let run=1024; run-- > 0; )
      {
        const [a,b,c,d,e,f] = function*(){ while(true) yield Math.random()*2-1 }(),
              F = (x,y) => a + (b + d*y + e*x)*x + (c + f*y)*y,
              G = (x,y) => array([
                b + d*y + 2*e*x,
                c + d*x + 2*f*y
              ])
        yield [F,G]
      }
    }()
  ).it('works on random 2d polynomials of degree 2', ([F,G]) => {
    const g = num_grad( xy => F(...xy.data) )

    for( let x=-3.14; x <=3.14; x+=0.42 )
    for( let y=-3.14; y <=3.14; y+=0.42 )
      expect(g([x,y])).toBeAllCloseTo(G(x,y))//, {rtol:1e-5, atol:1e-6})
  })


  forEachItemIn(
    function*(){
      for( let run=128; run-- > 0; )
      {
        const a = Array.from({length: 10}, () => Math.random()*2-1),
              F = (x,y,z) => {
                let result = 0
                for( let i=0; i<=2; i++ )
                for( let j=0; j<=2; j++ )
                for( let k=0; k<=2; k++ )
                  if( i+j+k <= 2 )
                    result += a[4*i+2*j+k] * x**i * y**j * z**k
                return result
              },
              G = (x,y,z) => {
                let result = [0,0,0]
                for( let i=0; i<=2; i++ )
                for( let j=0; j<=2; j++ )
                for( let k=0; k<=2; k++ )
                  if( i+j+k <= 2 ) {
                    if(i > 0) result[0] += a[4*i+2*j+k] * i * x**(i-1) * y** j    * z** k
                    if(j > 0) result[1] += a[4*i+2*j+k] * j * x** i    * y**(j-1) * z** k
                    if(k > 0) result[2] += a[4*i+2*j+k] * k * x** i    * y** j    * z**(k-1)
                  }
                return array(result)
              }
        yield [F,G]
      }
    }()
  ).it('works on random 3d polynomials of degree 2', ([F,G]) => {
    const g = num_grad( xy => F(...xy.data) )

    for( let x=-3.14; x <=3.14; x+=0.42 )
    for( let y=-3.14; y <=3.14; y+=0.42 )
    for( let z=-3.14; z <=3.14; z+=0.42 )
      expect(g([x,y,z])).toBeAllCloseTo(G(x,y,z))//, {rtol:1e-4, atol:1e-6})
  })
})
