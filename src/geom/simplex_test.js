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

import {range} from "../iter";

import {regular_simplex} from './simplex'
import {FrobeniusNorm} from '../la/norm';
import {rank} from "../la/rank";


describe('simplex', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  for( const args of [
    [],
    ['float32'],
    ['float64']
  ])
    forEachItemIn(
      range(1,128)
    ).it(`regular_simplex(${[...args,''].join(', ')}N) generates N+1 vertices with distance 1 to each other`, N => {
      let verts = regular_simplex(...args,N);

      expect(verts.dtype).toBe( args[0] || 'float64' )
      expect(verts.ndim ).toBe(2);
      expect(verts.shape).toEqual( Int32Array.of(N+1,N) );

      verts = verts.data;

      const norm = new FrobeniusNorm();

      for( let i=0; i < N+1; i++ )
      for( let j=0; j < N+1; j++ )
        if( i !== j )
        {
          norm.reset();
          for( let k=0; k < N; k++ )
            norm.include( verts[N*i+k] - verts[N*j+k] );

          expect(norm.resultSquare).toBeAllCloseTo(1);
        }
    })


  forEachItemIn(
    range(1,128)
  ).it('regular_simplex(N) returns a matrix of rank N and 1st row filled with zeros', N => {
    const verts = regular_simplex(N);

    expect( verts.sliceElems(0) ).toBeAllCloseTo(0, {rtol:0, atol:0});
    expect( rank(verts)         ).toBeAllCloseTo(N, {rtol:0, atol:0});
  })


  for( const args of [
    [],
    ['float32'],
    ['float64']
  ])
    it('regular_simplex(1024) generates N+1 vertices with distance 1 to each other', () => {
      const N = 1024;

      let verts = regular_simplex(...args,N);

      expect(verts.dtype).toBe( args[0] || 'float64' )
      expect(verts.ndim ).toBe(2);
      expect(verts.shape).toEqual( Int32Array.of(N+1,N) );

      verts = verts.data;

      const norm = new FrobeniusNorm();

      for( let i=0; i < N+1; i++ )
      for( let j=0; j < N+1; j++ )
        if( i !== j )
        {
          norm.reset();
          for( let k=0; k < N; k++ )
            norm.include( verts[N*i+k] - verts[N*j+k] );

          expect(norm.resultSquare).toBeAllCloseTo(1);
        }
    })
})
