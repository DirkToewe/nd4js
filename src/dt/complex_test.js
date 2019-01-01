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

import {Complex} from './complex'

describe('Complex', () => {
  it('does toString() real values correctly', () => {
    expect( new Complex(1337) .toString() ).toBe('1337');
    expect( new Complex(0.125).toString() ).toBe('0.125');
  })

  it('does toString() imaginary values correctly', () => {
    expect( new Complex(0, 1337) .toString() ).toBe('1337j');
    expect( new Complex(0, 0.125).toString() ).toBe('0.125j');
  })

  it('does toString() complex values correctly', () => {
    expect( new Complex(0.125, 1337).toString() ).toBe('0.125 + 1337j');
    expect( new Complex(1337, 0.125).toString() ).toBe('1337 + 0.125j');
  })

  it('constructor() works correctly', () => {
    for( let run=4; run-- > 0; )
    {
      const re = Math.random()*2048 - 1024,
            im = Math.random()*2048 - 1024,
            c = new Complex(re,im);
      expect(c.re).toBe(re);
      expect(c.im).toBe(im);
    }
  })

  it('does conj() correctly', () => {
    for( let run=4; run-- > 0; )
    {
      const re = Math.random()*2048 - 1024,
            im = Math.random()*2048 - 1024,
            c = new Complex(re,im).conj();
      expect(c.re).toBe(+re);
      expect(c.im).toBe(-im);
    }
  })

  it('does add() correctly', () => {
    for( let run=4; run-- > 0; )
    {
      const re1 = Math.random()*2048 - 1024,
            im1 = Math.random()*2048 - 1024,
            re2 = Math.random()*2048 - 1024,
            im2 = Math.random()*2048 - 1024,
            c1 = new Complex(re1,im1),
            c2 = new Complex(re2,im2),
            s = c1.add(c2);
      expect(s.re).toBe(re1+re2);
      expect(s.im).toBe(im1+im2);
    }
  })

  it('does sub() correctly', () => {
    for( let run=4; run-- > 0; )
    {
      const re1 = Math.random()*2048 - 1024,
            im1 = Math.random()*2048 - 1024,
            re2 = Math.random()*2048 - 1024,
            im2 = Math.random()*2048 - 1024,
            c1 = new Complex(re1,im1),
            c2 = new Complex(re2,im2),
            s = c1.sub(c2);
      expect(s.re).toBe(re1-re2);
      expect(s.im).toBe(im1-im2);
    }
  })

  it('does mul() correctly', () => {
    for( let run=4; run-- > 0; )
    {
      const re1 = Math.random()*2048 - 1024,
            im1 = Math.random()*2048 - 1024,
            re2 = Math.random()*2048 - 1024,
            im2 = Math.random()*2048 - 1024,
            c1 = new Complex(re1,im1),
            c2 = new Complex(re2,im2),
            p = c1.mul(c2);
      expect(p.re).toBe(re1*re2 - im1*im2);
      expect(p.im).toBe(re1*im2 + im1*re2);
    }
  })

  it('does div() correctly', () => {
    for( let run=4; run-- > 0; )
    {
      const re1 = Math.random()*2048 - 1024,
            im1 = Math.random()*2048 - 1024,
            re2 = Math.random()*2048 - 1024,
            im2 = Math.random()*2048 - 1024,
            c1 = new Complex(re1,im1),
            c2 = new Complex(re2,im2),
            q = c1.div(c2);
      const div = re2**2 + im2**2;
      expect(q.re).toBeCloseTo( (re1*re2 + im1*im2) / div, 4 );
      expect(q.im).toBeCloseTo( (im1*re2 - re1*im2) / div, 4 );
    }
  })
})
