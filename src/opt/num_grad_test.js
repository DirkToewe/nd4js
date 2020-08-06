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
import {array, NDArray} from '../nd_array'


describe('num_grad_forward', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  forEachItemIn(
    function*(){
      for( let run=8192; run-- > 0; )
      {
        const [a,b,c] = function*(){ for(;;) yield Math.random()*2-1 }(),
              F = x => a + (b + c*x)*x,
              G = x =>      b + c*x*2
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
        const [a,b,c,d,e,f] = function*(){ for(;;) yield Math.random()*2-1 }(),
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
        const [a,b,c] = function*(){ for(;;) yield Math.random()*2-1 }(),
              F = (x,y) => a + (b + c*x)*x,
              G = (x,y) =>      b + c*x*2
        yield [F,G]
      }
    }()
  ).it('works given random float      -> float      polynomials of degree 2', ([F,G]) => {
    const g = num_grad( xy => F(...xy.data) )

    for( let x=-3.14; x <=3.14; x+=0.42 )
      expect(g([x])).toBeAllCloseTo(G(x))//, {rtol:1e-4, atol:1e-6})
  });


  forEachItemIn(
    function*(){
      for( let run=1024; run-- > 0; )
      {
        const [a,b,c,d,e,f] = function*(){ for(;;) yield Math.random()*2-1 }(),
              F = (x,y) => a + (b + d*y + e*x)*x + (c + f*y)*y,
              G = (x,y) => array([
                b + d*y + 2*e*x,
                c + d*x + 2*f*y
              ])
        yield [F,G]
      }
    }()
  ).it('works given random float[2]   -> float      polynomials of degree 2', ([F,G]) => {
    const g = num_grad( xy => F(...xy.data) )

    for( let x=-3.14; x <=3.14; x+=0.42 )
    for( let y=-3.14; y <=3.14; y+=0.42 )
      expect(g([x,y])).toBeAllCloseTo(G(x,y))//, {rtol:1e-5, atol:1e-6})
  });


  forEachItemIn(
    function*(){
      for( let run=128; run-- > 0; )
      {
        const a = Float64Array.from({length: 3*3*2 + 1}, () => Math.random()*2-1),
              F = (x,y,z) => {
                let result = 0
                for( let i=0; i    <=2; i++ )
                for( let j=0; i+j  <=2; j++ )
                for( let k=0; i+j+k<=2; k++ )
                    result += a[3*(3*i + j) + k] * x**i * y**j * z**k
                return result
              },
              G = (x,y,z) => {
                let result = [0,0,0]
                for( let i=0; i    <=2; i++ )
                for( let j=0; i+j  <=2; j++ )
                for( let k=0; i+j+k<=2; k++ )
                { const                  coeff = a[3*(3*i + j) + k];
                  if(i > 0) result[0] += coeff * i * x**(i-1) * y** j    * z** k
                  if(j > 0) result[1] += coeff * j * x** i    * y**(j-1) * z** k
                  if(k > 0) result[2] += coeff * k * x** i    * y** j    * z**(k-1)
                }
                return array(result)
              }
        yield [F,G]
      }
    }()
  ).it('works given random float[3]   -> float      polynomials of degree 2', ([F,G]) => {
    const g = num_grad( xy => F(...xy.data) )

    for( let x=-3.14; x <=3.14; x+=0.42 )
    for( let y=-3.14; y <=3.14; y+=0.42 )
    for( let z=-3.14; z <=3.14; z+=0.42 )
      expect(g([x,y,z])).toBeAllCloseTo(G(x,y,z))//, {rtol:1e-4, atol:1e-6})
  });


  forEachItemIn(
    function*(){
      for( let run=8; run-- > 0; )
      {
        const a = Float64Array.from({length: 3*3*3*2 + 1}, () => Math.random()*2-1),
              F = (x,y,z,w) => {
                let result = 0
                for( let i=0; i      <=2; i++ )
                for( let j=0; i+j    <=2; j++ )
                for( let k=0; i+j+k  <=2; k++ )
                for( let l=0; i+j+k+l<=2; l++ )
                       result += a[3*(3*(3*i + j) + k) + l] * x**i * y**j * z**k * w**l;
                return result
              },
              G = (x,y,z,w) => {
                let result = Float64Array.of(0,0,0,0);
                for( let i=0; i      <=2; i++ )
                for( let j=0; i+j    <=2; j++ )
                for( let k=0; i+j+k  <=2; k++ )
                for( let l=0; i+j+k+l<=2; l++ )
                { const                  coeff = a[3*(3*(3*i + j) + k) + l];
                  if(i > 0) result[0] += coeff * i * x**(i-1) * y** j    * z** k    * w** l;
                  if(j > 0) result[1] += coeff * j * x** i    * y**(j-1) * z** k    * w** l;
                  if(k > 0) result[2] += coeff * k * x** i    * y** j    * z**(k-1) * w** l;
                  if(l > 0) result[3] += coeff * l * x** i    * y** j    * z** k    * w**(l-1);
                }
                return new NDArray(Int32Array.of(2,2), result);
              }
        yield [F,G]
      }
    }()
  ).it('works given random float[2,2] -> float      polynomials of degree 2', ([F,G]) => {
    const g = num_grad( xy => F(...xy.data) )

    for( let x=-3.14; x <=3.14; x+=0.42 )
    for( let y=-3.14; y <=3.14; y+=0.42 )
    for( let z=-3.14; z <=3.14; z+=0.42 )
    for( let w=-3.14; w <=3.14; w+=0.42 )
      expect( g([[x,y],
                 [z,w]]) ).toBeAllCloseTo(G(x,y,z,w))//, {rtol:1e-4, atol:1e-6})
  });


  forEachItemIn(
    function*(){
      for( let run=8192; run-- > 0; )
      {
        const [a1,b1,c1,
               a2,b2,c2,
               a3,b3,c3,
               a4,b4,c4,
               a5,b5,c5,
               a6,b6,c6] = function*(){ for(;;) yield Math.random()*2-1 }(),
                       F = (x) => [
                         [ a1 + (b1 + c1*x)*x,
                           a2 + (b2 + c2*x)*x ],
                         [ a3 + (b3 + c3*x)*x,
                           a4 + (b4 + c4*x)*x ],
                         [ a5 + (b5 + c5*x)*x,
                           a6 + (b6 + c6*x)*x ],
                       ],
                       G = (x) => [
                         [ b1 + c1*x*2,
                           b2 + c2*x*2 ],
                         [ b3 + c3*x*2,
                           b4 + c4*x*2 ],
                         [ b5 + c5*x*2,
                           b6 + c6*x*2 ]
                       ];
        yield [F,G]
      }
    }()
  ).it('works given random float      -> float[3,2] polynomials of degree 2', ([F,G]) => {
    const g = num_grad( xy => F(...xy.data) )

    for( let x=-3.14; x <=3.14; x+=0.42 )
      expect(g(x)).toBeAllCloseTo(G(x))//, {rtol:1e-4, atol:1e-6})
  });


  forEachItemIn(
    function*(){
      for( let run=512; run-- > 0; )
      {
        const [a1,b1,c1,d1,e1,f1,
               a2,b2,c2,d2,e2,f2,
               a3,b3,c3,d3,e3,f3,
               a4,b4,c4,d4,e4,f4,
               a5,b5,c5,d5,e5,f5,
               a6,b6,c6,d6,e6,f6] = function*(){ for(;;) yield Math.random()*2-1 }(),
              F = (x,y) => [
                [ a1 + (b1 + d1*y + e1*x)*x + (c1 + f1*y)*y,
                  a2 + (b2 + d2*y + e2*x)*x + (c2 + f2*y)*y,
                  a3 + (b3 + d3*y + e3*x)*x + (c3 + f3*y)*y ],
                [ a4 + (b4 + d4*y + e4*x)*x + (c4 + f4*y)*y,
                  a5 + (b5 + d5*y + e5*x)*x + (c5 + f5*y)*y,
                  a6 + (b6 + d6*y + e6*x)*x + (c6 + f6*y)*y ],
              ],
              G = (x,y) => array([
                [[ b1 + d1*y + 2*e1*x,
                   c1 + d1*x + 2*f1*y ],
                 [ b2 + d2*y + 2*e2*x,
                   c2 + d2*x + 2*f2*y ],
                 [ b3 + d3*y + 2*e3*x,
                   c3 + d3*x + 2*f3*y ]],
                [[ b4 + d4*y + 2*e4*x,
                   c4 + d4*x + 2*f4*y ],
                 [ b5 + d5*y + 2*e5*x,
                   c5 + d5*x + 2*f5*y ],
                 [ b6 + d6*y + 2*e6*x,
                   c6 + d6*x + 2*f6*y ]]
              ])
        yield [F,G]
      }
    }()
  ).it('works given random float[2]   -> float[2,3] polynomials of degree 2', ([F,G]) => {
    const g = num_grad( xy => F(...xy.data) )

    for( let x=-3.14; x <=3.14; x+=0.42 )
    for( let y=-3.14; y <=3.14; y+=0.42 )
      expect(g([x,y])).toBeAllCloseTo(G(x,y), {atol:1e-7})
  })
})
