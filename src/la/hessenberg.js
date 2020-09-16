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

import {asarray, NDArray} from '../nd_array'
import {ARRAY_TYPES} from '../dt'
import {FrobeniusNorm} from './norm'


// TODO: bidiag_solve


export function _hessenberg_decomp(N, U,H, off)
{
  N   |= 0;
  off |= 0;

  // INIT U TO IDENTITY
  for( let i=N-1; i-- > 0; ) U[off + N*i+i] = 1.0;

  const NORM = new FrobeniusNorm(),
     lastRow = off + N*(N-1);

  // ELIMINATE ELEMENTS LEFT OF (i,i-1) VIA HOUSEHOLDER
  for( let i=N; --i > 1; )
  {
    const  rowI = off + N*i,
      ii = rowI+(i-1);
    // COMPUTE HOUSEHOLDER VECTOR
    NORM.reset();
    for( let j=i-1; j-- > 0; )
        NORM.include(H[rowI+j]);
    if( NORM.max === 0 ) continue;
    const norm =          NORM.resultIncl(H[ii]) * (H[ii] > 0 ? -1 : +1); // <- avoid cancellation error
                          NORM.   include(H[ii] -= norm);
    const  max =          NORM.max,
           div =Math.sqrt(NORM.sum);
    for( let j=i; j-- > 0; )
      H[rowI+j] = H[rowI+j] / max * Math.SQRT2 / div;

    // APPLY HOUSEHOLDER TO RIGHT OF H
    for( let j=i; j-- > 0; ) {
      let sum = 0;
      for( let k=i; k-- > 0; ) sum += H[off + N*j+k] * H[rowI+k];
      for( let k=i; k-- > 0; )        H[off + N*j+k]-= H[rowI+k] * sum;
    }

    // APPLY HOUSEHOLDER TO LEFT OF H
    U.fill(0.0, lastRow, off + N*N);
    for( let j=i; j-- > 0; )
    for( let k=N; k-- > 0; ) U[lastRow+k] += H[off + N*j+k] * H[rowI+j];
    for( let j=i; j-- > 0; )
    for( let k=N; k-- > 0; )                 H[off + N*j+k]-= H[rowI+j]*U[lastRow+k];

    // APPLY HOUSEHOLDER TO RIGHT OF U
    for( let j=N-1; j-- > 0; ) {
      let sum = 0;
      for( let k=i; k-- > 0; ) sum += U[off + N*j+k] * H[rowI+k];
      for( let k=i; k-- > 0; )        U[off + N*j+k]-= H[rowI+k] * sum;
    }

    // FINISH ROW i
    H.fill(0.0, rowI, ii);
                    H[ii] = norm;
  }

  // LAST ROW OF U WAS USED AS TEMP MEMORY
  U.fill(0.0, lastRow, off + N*N-1);
                     U[off + N*N-1] = 1;
}


export function hessenberg_decomp(A)
{
  A = asarray(A);
  if( A.ndim < 2 ) throw new Error('hessenberg_decomp(A): A must at least be 2D.');
  const
      DType      = A.dtype==='float32' ? 'float32' : 'float64',
      DTypeArray = ARRAY_TYPES[A.dtype==='float32' ? 'float32' : 'float64'], // <- ensure at least double precision
           shape = A.shape,
      N  = shape[shape.length-1] | 0;
  if( N != shape[shape.length-2] ) throw new Error('hessenberg_decomp(A): A must be square.');

  const H =     DTypeArray.from(A.data); A = undefined;
  const U = new DTypeArray(H.length);

  for( let off=H.length; (off -= N*N) >= 0; )
    _hessenberg_decomp(N, U,H, off);

  return [
    new NDArray(shape,U),
    new NDArray(shape,H)
  ];
};
