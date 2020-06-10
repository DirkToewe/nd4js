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

import {ARRAY_TYPES} from '../dt'
import {NDArray} from '../nd_array'
import {rand_normal} from '../rand_normal'
import {FrobeniusNorm} from './norm'


export function rand_ortho( dtype, ...shape )
{
  // REFERENCES:
  //   - https://blogs.sas.com/content/iml/2012/03/28/generating-a-random-orthogonal-matrix.html
  if( ! (dtype in ARRAY_TYPES) ) {
    shape.unshift(dtype);
    dtype = 'float64';
  }

  if( ! dtype.startsWith('float') )
    throw new Error(`Unsupported dtype: '${dtype}'.`);

  if( shape.length === 1 )
    shape.push(shape[0]);

  for( let i=shape.length; i-- > 0; )
    if( shape[i] % 1 !== 0 || shape[i] < 1 )
      throw new Error(`shape[${i}] = ${shape[i]} not a positive integer.`);
  shape = Int32Array.from(shape);

  const [M,N] = shape.slice(-2),
           L  = Math.max(M,N);

  const DTypeArray = ARRAY_TYPES[dtype],
        U = new DTypeArray( shape.reduce((x,y) => x*y) ),
        Q = new DTypeArray(L*L),
        x = new DTypeArray(L),
        d = new DTypeArray(L);

  const NORM = new FrobeniusNorm();

  // TODO: this should be computable more efficiently for non-square matrices

  for( let off=0; off < U.length; off += M*N )
  {
    // INIT Q TO IDENTITY
    for( let i=0; i < L; i++ )
    for( let j=0; j < L; j++ )
      Q[L*i+j] = (i===j)*1;

    d[0] = Math.random() <= 0.5 ? -1 : +1;

    // apply series of random householder transformations/reflections
    for( let k=L; k-- > 1; )
    {
      NORM.reset();
      for( let i=0;; i++ ) {
        const  x_i = x[i] = rand_normal(); if(i >= k) break;
        NORM.include(x_i);
      }                                              d[k] = x[k] > 0 ? -1 : +1;
      const  norm =          NORM.resultIncl(x[k]) * d[k];
      if(0===norm) continue; NORM.   include(x[k] -= norm);
      const   max =          NORM.max,
              div =Math.sqrt(NORM.sum);
      for( let i=k; i >= 0; i-- )
        x[i] = x[i] / max / div;

      // apply householder to right of Q
      for(let i=L; i-- > 0;){
        let sum = 0; for(let j=k; j >= 0; j--) sum += Q[L*i+j] * x[j];
            sum*= 2; for(let j=k; j >= 0; j--)        Q[L*i+j]-= x[j] * sum;
      }
    }

    // flip row signs according to d
    for( let i=0; i < L; i++ )
    for( let j=0; j < L; j++ )
      Q[L*i+j] *= d[i];

    // MOVE Q -> U
    for( let i=0; i < M; i++ )
    for( let j=0; j < N; j++ )
      U[off + N*i+j] = Q[L*i+j];
  }

  return new NDArray(shape, U);
}
