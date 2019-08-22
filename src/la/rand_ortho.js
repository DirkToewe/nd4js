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
        y = new DTypeArray(L),
        d = new DTypeArray(L);

  const NORM = new FrobeniusNorm();

  for( let off=0; off < U.length; off += M*N )
  {
    // INIT Q TO IDENTITY
    for( let i=0; i < L; i++ )
    for( let j=0; j < L; j++ )
      Q[L*i+j] = i===j;

    d[0] = Math.random() <= 0.5 ? -1 : +1;

    // apply series of random householder transformations/reflections
    for( let k=1; k < L; k++ )
    {
      NORM.reset();
      for( let i=0; i <= k; i++ )
        NORM.include( x[i] = rand_normal() );

      const sgn = x[k] <= 0 ? -1 : +1,
        s = sgn * NORM.result;

      d[k]  = -sgn;
      x[k] +=  s;

      y.fill(0, 0,k+1);

      for( let j=0; j <= k; j++ )
      for( let i=0; i <= k; i++ ) y[i] += x[j]*Q[L*j+i];

      const beta = s*x[k];

      for( let i=0; i <= k; i++ )
        y[i] /= beta;

      // apply householder transformation/reflection
      for( let i=0; i <= k; i++ )
      for( let j=0; j <= k; j++ )
        Q[L*i+j] -= x[i]*y[j];
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
