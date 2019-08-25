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

import {forEachItemIn, CUSTOM_MATCHERS} from '../../jasmine_utils'
import {array, NDArray} from '../../nd_array'
import {tabulate} from '../../tabulate'
import {zip_elems} from '../../zip_elems'
import {rosenbrock, rosenbrock_grad} from './rosenbrock'
import {num_grad} from '../num_grad'


describe('rosenbrock', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  it('rosenbrock_grad has zero grad for x=[1,1] and x=[1,1,1]', () => {
    const g11  = rosenbrock_grad([1,1]),
          g111 = rosenbrock_grad([1,1,1])

    expect(g11 .shape).toEqual(Int32Array.of(2))
    expect(g111.shape).toEqual(Int32Array.of(3))
    expect(g11 ).toBeAllCloseTo(0)
    expect(g111).toBeAllCloseTo(0)
  })


  const rosenbrock_num_grad = num_grad(rosenbrock,2**-16)


  forEachItemIn(
    function*(){
      const  S = 3.14, Δ = 0.042

      for( let x = -S; x <= +S; x+=Δ )
      for( let y = -S; y <= +S; y+=Δ )
        yield array([x,y])
    }()
  ).it('rosenbrock_grad (2d) works for generated examples', (xy) => {
    const g = rosenbrock_grad(xy),
          G = rosenbrock_num_grad(xy)

    expect(g).toBeAllCloseTo(G, {rtol:1e-4, atol:1e-6})
  })


  forEachItemIn(
    function*(){
      const  S = 1.57, Δ = 0.136

      for( let x = -S; x <= +S; x+=Δ )
      for( let y = -S; y <= +S; y+=Δ )
      for( let z = -S; z <= +S; z+=Δ )
        yield array([x,y,z])
    }()
  ).it('rosenbrock_grad (3d) works on generated examples', (xyz) => {
    const g = rosenbrock_grad(xyz),
          G = rosenbrock_num_grad(xyz)

    expect(g).toBeAllCloseTo(G, {rtol:1e-4, atol:1e-6})
  })


  forEachItemIn(
    function*(){
      const  S = 1.57, Δ = 0.272

      for( let w = -S; w <= +S; w+=Δ )
      for( let x = -S; x <= +S; x+=Δ )
      for( let y = -S; y <= +S; y+=Δ )
      for( let z = -S; z <= +S; z+=Δ )
        yield array([w,x,y,z])
    }()
  ).it('rosenbrock_grad (4d) works on generated examples', (wxyz) => {
    const g = rosenbrock_grad(wxyz),
          G = rosenbrock_num_grad(wxyz)

    expect(g).toBeAllCloseTo(G, {rtol:1e-4, atol:1e-6})
  })


  forEachItemIn(
    function*(){
      for( let d=2; d <= 8; d++ )
      for( let run=0; run < 128; run++ )
        yield tabulate([d], 'float64', () => Math.random()*2-1)
    }()
  ).it('rosenbrock_grad works on random examples', (x) => {
    const g = rosenbrock_grad(x),
          G = rosenbrock_num_grad(x)

    expect(g).toBeAllCloseTo(G, {rtol:1e-4, atol:1e-6})
  })


  forEachItemIn(
    function*(){
      const  S = 3.14, Δ = 0.042, N = 42,
        shape2 = Int32Array.of(2)

      for( let x = -S; x <= +S; x+=Δ )
      for( let y = -S; y <= +S; y+=Δ ) {
        yield new NDArray(shape2, Float32Array.of(x,y))
        yield new NDArray(shape2, Float64Array.of(x,y))
      }

      for( let z = -S; z <= +S; z+=Δ )
        for( const dtype of ['float32','float64'] ) {
          yield tabulate([2*N+1, 2], dtype, (i,j) => [  (i-N)/N*S,z][j])
          yield tabulate([2*N+1, 2], dtype, (i,j) => [z,(i-N)/N*S  ][j])
        }

      for( const dtype of ['float32','float64'] ) {
        yield tabulate([2*N+1, 2*N+1, 2], (i,j,k) => [(i-N)/N*S, (j-N)/N*S][k])
        yield tabulate([2*N+1, 2*N+1, 2], (i,j,k) => [(j-N)/N*S, (i-N)/N*S][k])
      }
    }()
  ).it('rosenbrock (2d) has a global minimum at x=[1,1]', x => {
    const y_min = rosenbrock([1,1])

    expect(y_min.shape).toEqual(Int32Array.of())
    expect(y_min).toBeAllCloseTo(0)

    const y = rosenbrock(x)

    expect( y.mapElems(x => x >= 0) ).toBeAllCloseTo(true, {rtol:0, atol:0})
  })


  forEachItemIn(
    function*(){
      const  S = 3.14, Δ = 0.271, N = 13,
        shape3 = Int32Array.of(3)

      for( let x = -S; x <= +S; x+=Δ )
      for( let y = -S; y <= +S; y+=Δ )
      for( let z = -S; z <= +S; z+=Δ ) {
        yield new NDArray(shape3, Float32Array.of(x,y,z))
        yield new NDArray(shape3, Float64Array.of(x,y,z))
      }

      for( let x = -S; x <= +S; x+=Δ )
      for( let y = -S; y <= +S; y+=Δ )
        for( const dtype of ['float32','float64'] ) {
          yield tabulate([2*N+1, 3], dtype, (i,j) => [x,y,(i-N)/N*S    ][j])
          yield tabulate([2*N+1, 3], dtype, (i,j) => [x,  (i-N)/N*S,  y][j])
          yield tabulate([2*N+1, 3], dtype, (i,j) => [    (i-N)/N*S,x,y][j])
        }

      for( let x = -S; x <= +S; x+=Δ )
        for( const dtype of ['float32','float64'] ) {
          yield tabulate([2*N+1, 2*N+1, 3], (i,j,k) => [x,(i-N)/N*S,  (j-N)/N*S  ][k])
          yield tabulate([2*N+1, 2*N+1, 3], (i,j,k) => [  (i-N)/N*S,x,(j-N)/N*S  ][k])
          yield tabulate([2*N+1, 2*N+1, 3], (i,j,k) => [  (i-N)/N*S,  (j-N)/N*S,x][k])
          yield tabulate([2*N+1, 2*N+1, 3], (i,j,k) => [x,(j-N)/N*S,  (i-N)/N*S  ][k])
          yield tabulate([2*N+1, 2*N+1, 3], (i,j,k) => [  (j-N)/N*S,x,(i-N)/N*S  ][k])
          yield tabulate([2*N+1, 2*N+1, 3], (i,j,k) => [  (j-N)/N*S,  (i-N)/N*S,x][k])
        }

      for( const dtype of ['float32','float64'] ) {
        yield tabulate([2*N+1, 2*N+1, 2*N+1, 3], (i,j,k,l) => [(i-N)/N*S, (j-N)/N*S, (k-N)/N*S][l])
        yield tabulate([2*N+1, 2*N+1, 2*N+1, 3], (i,j,k,l) => [(i-N)/N*S, (k-N)/N*S, (j-N)/N*S][l])
        yield tabulate([2*N+1, 2*N+1, 2*N+1, 3], (i,j,k,l) => [(j-N)/N*S, (i-N)/N*S, (k-N)/N*S][l])
        yield tabulate([2*N+1, 2*N+1, 2*N+1, 3], (i,j,k,l) => [(j-N)/N*S, (k-N)/N*S, (i-N)/N*S][l])
        yield tabulate([2*N+1, 2*N+1, 2*N+1, 3], (i,j,k,l) => [(k-N)/N*S, (i-N)/N*S, (j-N)/N*S][l])
        yield tabulate([2*N+1, 2*N+1, 2*N+1, 3], (i,j,k,l) => [(k-N)/N*S, (j-N)/N*S, (i-N)/N*S][l])
      }
    }()
  ).it('rosenbrock (3d) has a global minimum at x=[1,1,1]', x => {
    const y_min = rosenbrock([1,1,1])

    expect(y_min.shape).toEqual(Int32Array.of())
    expect(y_min).toBeAllCloseTo(0)

    const y = rosenbrock(x)

    expect( y.mapElems(x => x >= 0) ).toBeAllCloseTo(true, {rtol:0, atol:0})
  })
})