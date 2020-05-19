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

import {array} from './nd_array'
import {tabulate} from './tabulate'
import {CUSTOM_MATCHERS} from './jasmine_utils'


describe('CUSTOM_MATCHERS.toBeAllCloseTo', () => {

  const expectAllClose = CUSTOM_MATCHERS.toBeAllCloseTo().compare

  it('fails for broadcast-incompatible arrays', () => {
    {
      const test = expectAllClose(
        tabulate([3,1,2], () => 1),
        tabulate([2,1,2], () => 1)
      )
      expect(test.pass).toBe(false)
      expect(test.message).toBe('Expected shape [3,1,2] to be broadcast-compatible to [2,1,2].')
    };{
      const test = expectAllClose(
        tabulate([3,1,2], () => 1),
        tabulate(    [3], () => 1)
      )
      expect(test.pass).toBe(false)
      expect(test.message).toBe('Expected shape [3,1,2] to be broadcast-compatible to [3].')
    };{
      const test = expectAllClose(
        tabulate([3,2,1], () => 1),
        tabulate(  [3,2], () => 1)
      )
      expect(test.pass).toBe(false)
      expect(test.message).toBe('Expected shape [3,2,1] to be broadcast-compatible to [3,2].')
    }
  })

  it('passes for broadcast-equal arrays', () => {
    {
      const test = expectAllClose(
        tabulate([3,1,2], (i,j,k) => 100*i+10*j+k+111),
        tabulate([3,1,2], (i,j,k) => 100*i+10*j+k+111)
      )
      expect(test.pass).withContext(test.message).toBe(true)
    };{
      const test = expectAllClose(
        tabulate([3,1,2], (i,j,k) => 10*i+k+11),
        tabulate([3,4,2], (i,j,k) => 10*i+k+11)
      )
      expect(test.pass).withContext(test.message).toBe(true)
    };{
      const test = expectAllClose(
        tabulate([3,1,2], (i,j,k) => k+1),
        tabulate(  [4,2], (i,j)   => j+1)
      )
      expect(test.pass).withContext(test.message).toBe(true)
    };{
      const test = expectAllClose(
        tabulate([3,2,1], (i,j,k) => j+1),
        tabulate(  [2,4], (i,j)   => i+1)
      )
      expect(test.pass).withContext(test.message).toBe(true)
    }
  })

  it('fails for vastly different arrays', () => {
    {
      const test = expectAllClose(
        tabulate([3,1,2], (i,j,k) => 100*i+10*j+k+111),
        tabulate([3,1,2], (i,j,k) => 100*i+10*j+k+111 + (i==1 ? 1337 : 0) )
      )
      expect(test.pass).withContext(test.message).toBe(false)
    };{
      const test = expectAllClose(
        tabulate([3,1,2], (i,j,k) => 10*i+k+11),
        tabulate([3,4,2], (i,j,k) => 10*i+k+11 + (j==2 ? 1337 : 0) )
      )
      expect(test.pass).withContext(test.message).toBe(false)
    };{
      const test = expectAllClose(
        tabulate([3,1,2], (i,j,k) => k+1),
        tabulate(  [4,2], (i,j)   => j+1 + (i==2 ? 1337 : 0) )
      )
      expect(test.pass).withContext(test.message).toBe(false)
    };{
      const test = expectAllClose(
        tabulate([3,2,1], (i,j,k) => j+1),
        tabulate(  [2,4], (i,j)   => i+1 + (j==3 ? 1337 : 0) )
      )
      expect(test.pass).withContext(test.message).toBe(false)
    }
  })

  it('handles NaN and Infinity correctly', () => {
    for( const inf of [-Infinity,+Infinity] )
    {
      const  test = expectAllClose(inf,inf)
      expect(test.pass).toBe(true)
    }

    for( const x of [-Infinity, 1337, NaN, +Infinity] )
    for( const y of [-x, NaN] )
    {
      const  test = expectAllClose(x,y)
      expect(test.pass).toBe(false)
    }
  })
})
