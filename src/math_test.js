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


import {Complex} from './dt/complex'
import math from './math'

describe('math', () => {
  beforeEach( () => {
    jasmine.addCustomEqualityTester( (a,b) => {
      if( a instanceof Complex ) return a.equals(b);
      if( b instanceof Complex ) return b.equals(a);
    })
  })

  it('does add() correctly', () => {
    expect( math.add(1,2) ).toBe(3)
    expect( math.add(1.25,2.25) ).toBe(3.5)
    expect( math.add(
      new Complex(1.25,2.25),
      new Complex(3.25,4.25)
    )).toEqual( new Complex(4.5,6.5) )
  })

  it('does sub() correctly', () => {
    expect( math.sub(1,2) ).toBe(-1)
    expect( math.sub(1,2.5) ).toBe(-1.5)
    expect( math.sub(
      new Complex(1,2),
      new Complex(3.5,4.5)
    )).toEqual( new Complex(-2.5,-2.5) )
  })

  it('does mul() correctly', () => {
    expect( math.mul(2,3) ).toBe(6)
    expect( math.mul(4.5,6) ).toBe(27)
    expect( math.mul(
      new Complex(1, 2),
      new Complex(3.5, 4.5)
    )).toEqual( new Complex(-5.5, 11.5) )
  })

  it('does div() correctly', () => {
    expect( math.div(1,2) ).toBe(0.5)
    expect( math.div(4.5,3) ).toBe(1.5)
    expect( math.div(
      new Complex(2.5, 5),
      new Complex(3, 4)
    )).toEqual( new Complex(1.1, 0.2) )    
  })
})
