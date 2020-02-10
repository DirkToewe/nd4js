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


import {tabulate} from './tabulate'
import {matmul} from './la/matmul'
import {rand_ortho} from './la/rand_ortho'


const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from


export function _rand_rankdef(...shape)
{
  if( shape.length < 2 )
    throw new Error('Assertion failed.');

  const [M,N]= shape.slice(-2),
           L = Math.min(M,N);

  const ndim = randInt(0,7);

  let r_shape = shape.slice(0,-2);

  // use random, mostly rank-deficient SVD to generate test matrix A
  const ranks = tabulate(r_shape, 'int32', () => randInt(0,L+1) ), // <- ranks
            U = rand_ortho('float64', ...r_shape,M,L),
            V = rand_ortho('float64', ...r_shape,L,N),
            S = tabulate([...r_shape,L,L], 'float64', (...idx) => {
              const j = idx.pop(),
                    i = idx.pop(), rank = ranks(...idx);

              if( i !== j || rank <= i ) return 0; 

              return Math.random()*8 + 0.1
            }),
            A = matmul(U,S,V);

  // add random scaling
  for( let S_off=S.data.length; (S_off -= L*L) >= 0; )
  {
    const scale = Math.random()*1e300 + 1;

    for( let i=L; i-- > 0; )
      S.data[S_off + L*i+i] *= scale;
  }

  Object.freeze(ranks.data.buffer);
  Object.freeze(A.data.buffer);
  return [A, ranks];
}
