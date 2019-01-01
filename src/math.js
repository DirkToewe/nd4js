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


export const math = {
  nextUp(x)
  {
    // FIXME implement this completely
    // https://gist.github.com/Yaffle/4654250
    const
      EPSILON   = Number.EPSILON,
      MIN_VALUE = Number.MIN_VALUE;
    let y = x * (x < 0 ? 1 - EPSILON/2 : 1 + EPSILON);
    if( ! isFinite(x) ) throw new Error('Assertion failed!');
    if( ! isFinite(y) ) throw new Error('Assertion failed!');
    if (y === x) y += MIN_VALUE;
    let b = x + (y - x) / 2; if (x < b && b < y) { y = b; }
    let c =     (y + x) / 2; if (x < c && c < y) { y = c; }
    return y === 0 ? -0 : y;
  },

  add(x,y) {
    if( x instanceof Complex ) return x.add(y);
    if( y instanceof Complex ) return y.add(x);
    return x + y;
  },

  sub(x,y) {
    if( x instanceof Complex ) return             x .sub(y);
    if( y instanceof Complex ) return new Complex(x).sub(y);
    return x - y;
  },

  mul(x,y) {
    if( x instanceof Complex ) return x.mul(y);
    if( y instanceof Complex ) return y.mul(x);
    return x * y;
  },

  div(x,y) {
    if( x instanceof Complex ) return             x .div(y);
    if( y instanceof Complex ) return new Complex(x).div(y);
    return x / y;
  },

  cast(x,dtype)
  {
    if( dtype ===      'int32' ) return x & 0xFFFFFFFF;
    if( dtype ===    'float32' ) return Math.fround(x);
    if( dtype === 'complex128' ) return x instanceof Complex ? x : new Complex(x);
    return x;
  },

  zero: dtype => math.cast(0,dtype),
  one : dtype => math.cast(1,dtype),

  neg: x => x instanceof Complex ? x.neg() : -x,
  abs: x => x instanceof Complex ? x.abs() : Math.abs(x),

  sqrt: x => x instanceof Complex
    ?   x.sqrt()
    : ( x >= 0
      ?           Math.sqrt(x)
      : new Complex(x).sqrt()
    ),

  exp: x => x instanceof Complex ? x.exp() : Math.exp(x),

  min: (x,y) => Math.min(x,y),
  max: (x,y) => Math.max(x,y),

  hypot: (x,y) => Math.hypot(x,y),
  atan2: (x,y) => Math.atan2(x,y),

  conj: x => x instanceof Complex ? x.conj() : x,

  is_equal(x,y) {
    if( x instanceof Complex ) return x.equals(y);
    if( y instanceof Complex ) return y.equals(x);
    return x == y;
  },

  is_close(x,y) {
    const atol = 1e-8,
          rtol = 1e-5,
           tol = atol + rtol * math.max(
            math.abs(x),
            math.abs(y)
          );
    return math.abs( math.sub(x,y) ) <= tol;
  } 
}

export default math
