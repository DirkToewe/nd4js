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

import {ARRAY_TYPES} from "./dt";
import {tabulate} from './tabulate';

import {matmul2} from './la/matmul'

import {AleaRNG} from "./rand/alea_rng";


export class TestRNG extends AleaRNG
{
  rankDef(...shape)
  {
    if( !(shape.length >= 2       ) ) throw new Error('TestRNG::rankDef(...shape): shape.length must be at least 2.');

    const dtype = shape[0] in ARRAY_TYPES ? shape.shift() : 'float64';

    if( ! shape.every(s => s%1===0) ) throw new Error('TestRNG::rankDef(...shape): shape must be all integers.');
    if( ! shape.every(s => s  >  0) ) throw new Error('TestRNG::rankDef(...shape): shape must be all positive integers.');

    const N = shape.pop(),
          M = shape.pop(),
          L = Math.min(M,N);

    // use random, mostly rank-deficient SVD to generate test matrix A
    const ranks = tabulate(shape, 'int32', () => this.int(0,L+1) ), // <- ranks
      u = this.ortho(dtype, ...shape,M,L),
      v = this.ortho(dtype, ...shape,L,N),
      V =     v.data,
      R = ranks.data;

    // V = S @ V, where S is a batch of diagonal scaling matrices with given `ranks`
    for( let i=R.length; i-- > 0; )
    for( let j=L       ; j-- > 0; )
    {
      const scale = R[i] <= j ? 0 : this.uniform(1e-4, 1e+4);

      for( let k=N; k-- > 0; )
        V[N*(L*i+j)+k] *= scale;
    }

    const a = matmul2(u,v); // <- A = U @ S @ V

    Object.freeze(ranks.data.buffer);
    Object.freeze(    a.data.buffer);
    return [a, ranks];
  }
}
