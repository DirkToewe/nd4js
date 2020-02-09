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
import {_triu_solve} from './tri'


export function lu_decomp(A)
{
  A = asarray(A)
  const LU = ARRAY_TYPES[A.dtype==='float32' ? 'float32' : 'float64'].from(A.data),
      [N,M]= A.shape.slice(-2),
         P = new Int32Array(LU.length/M)

  if( N != M )
    throw new Error('Last two dimensions must be quadratic.')

  for(
    let LU_off=0,
         P_off=0;
    LU_off < LU.length;
    LU_off += N*N,
     P_off += N
  )
  {
    for( let i=0; i < N; i++ ) P[P_off+i] = i;
    for( let i=0; i < N; i++ )
    {
      const row_i = LU_off + i*N;
      // ROW PIVOTING
      {
        let p=i;
        for( let j=i+1; j < N; j++ )
          if( Math.abs( LU[LU_off + N*j+i] )
            > Math.abs( LU[LU_off + N*p+i] ) )
            p=j;

        if( i != p )
        {
          const   P_p = P[P_off+i]; P[P_off+i] = P[P_off+p]; P[P_off+p] = P_p; // KEEP TRACK OF ROW SWAPS
          const row_p = LU_off + p*N;
          // SWAP ROWS
          for( let j=0; j < N; j++ ) {
            const tmp = LU[row_i+j]; LU[row_i+j] = LU[row_p+j]; LU[row_p+j] = tmp;
          }
        }
      }
      // ELIMINATE ELEMENTS BELOW PIVOT
      for( let j=i+1; j < N; j++ )
      {
        const
          row_j = LU_off + j*N,
          scale = LU[row_j+i] / LU[row_i+i];
        LU[row_j+i] = scale;
        for( let k=i+1; k < N; k++ )
          LU[row_j+k] -= scale * LU[row_i+k];
      }
    }
  }

  return [
    new NDArray(A.shape,           LU),
    new NDArray(A.shape.slice(0,-1),P)
  ];
}


export function lu_solve(LU,P,y)
{
  if( undefined == y ) { y=P; [LU,P] = LU; }
  LU = asarray(LU); if( LU.ndim < 2 ) throw new Error('LU must be at least 2D.');
  P  = asarray(P ); if( P .ndim < 1 ) throw new Error( 'P must be at least 1D.');
   y = asarray( y); if(  y.ndim < 2 ) throw new Error( 'y must be at least 2D.');

  const
    [N,M] = LU.shape.slice(-2),
    [I,J] =  y.shape.slice(-2);
  if( M != N ) throw new Error('Last two dimensions of LU must be quadratic.')
  if( M != I ) throw new Error("LU and y don't match.");
  if( M != P.shape.slice(-1) )
    throw new Error("LU and P don't match.")

  const
    ndim = Math.max(LU.ndim, P.ndim+1, y.ndim),
    shape = Int32Array.from({ length: ndim }, () => 1 );
  shape[ndim-2] = I;
  shape[ndim-1] = J;

  // FIND COMMON (BROADCASTED) SHAPE
  for( let arr of [LU,y] )
    for( let i=ndim-2, j=arr.ndim-2; i-- > 0 && j-- > 0; )
      if( 1 === shape[i] )
        shape[i] = arr.shape[j];
      else if( shape[i] != arr.shape[j] && arr.shape[j] != 1 )
        throw new Error('LU and y are not broadcast-compatible.');

  for( let i=ndim-2, j=P.ndim-1; i-- > 0 && j-- > 0; )
    if( 1 === shape[i] )
      shape[i] = P.shape[j]
    else if( shape[i] != P.shape[j] && P.shape[j] != 1 )
      throw new Error('P is not broadcast-compatible.');

  // GENERATE RESULT DATA
  const x_dat = new ARRAY_TYPES[[LU,P,y].every(a => a.dtype==='float32') ? 'float32' : 'float64']( shape.reduce((a,b) => a*b, 1) );

  const
    LU_dat = LU.data,
     P_dat =  P.data,
     y_dat =  y.data;
  let
    LU_off = 0, LU_stride = 1,
     P_off = 0,  P_stride = 1,
     y_off = 0,  y_stride = 1,
     x_off = 0;

  function solv(d) {
    if( d === ndim-2 ) {
      LU_stride = N*N;
       P_stride = N;
       y_stride = N*J;

      // COPYING PERMUTED y
      for( let i=0; i < I; i++ ) {
        const row_x = x_off + J * i,
              row_y = y_off + J * P_dat[P_off+i];
        for( let j=0; j < J; j++ )
          x_dat[row_x+j] = y_dat[row_y+j];
      }

      // FORWARD SUBSTITUTION
      for( let i=0; i < I; i++ )
      for( let j=0; j < J; j++ )
      for( let k=0; k < i; k++ )
        x_dat[x_off+i*J+j] -= LU_dat[LU_off+N*i+k] * x_dat[x_off+k*J+j]

      // BACKWARD SUBSTITUTION
      _triu_solve(I,I,J, LU_dat,LU_off, x_dat,x_off);

      LU_off += LU_stride;
       P_off +=  P_stride;
       y_off +=  y_stride;
       x_off +=  y_stride;

      return;
    }
    for( let l=shape[d]; ; l-- ) {
      solv(d+1);
      if( l == 1 ) break
      // RESET ALONG BROADCAST AXES
      if( ! (LU.shape[ d - ndim + LU.ndim   ] > 1) ) LU_off -= LU_stride
      if( ! ( P.shape[ d - ndim +  P.ndim+1 ] > 1) )  P_off -=  P_stride
      if( ! ( y.shape[ d - ndim +  y.ndim   ] > 1) )  y_off -=  y_stride
    }
    LU_stride *= LU.shape[ d - ndim + LU.ndim   ] || 1
     P_stride *=  P.shape[ d - ndim +  P.ndim+1 ] || 1
     y_stride *=  y.shape[ d - ndim +  y.ndim   ] || 1
  }
  solv(0);

  return new NDArray(shape,x_dat);
}
