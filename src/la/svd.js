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

import math from '../math'
import {asarray, NDArray} from '../nd_array'
import {ARRAY_TYPES, eps} from '../dt'
import {SingularMatrixSolveError} from './singular_matrix_solve_error'


export {svd_jac_2sided as svd_decomp} from './svd_jac_2sided'


export function svd_rank( sv )
{
  sv = asarray(sv)

  const   N = sv.shape[sv.ndim-1],
    r_shape = sv.shape.slice(0,-1),
        EPS = math.sqrt(eps(sv.dtype))

  sv = sv.data
  const r = new ARRAY_TYPES['int32'](sv.length/N)

  for( let off=0; off < r.length; off++ )
  {
    const T = EPS * math.abs(sv[N*off])
    if( r[off] !== 0 )
      throw new Error('Assertion failed.')
    for( ; r[off] < N; r[off]++ )
    {
      const sv_r = math.abs(sv[N*off + r[off]])
      if( ! isFinite(sv_r) )
        throw new Error('svd_rank(): NaN or Infinity encountered.')
      if( sv_r <= T )
        break
    }
  }

  return new NDArray(r_shape,r)
}


export function svd_solve(U,sv,V, y)
{
  U = asarray(U)
  sv= asarray(sv)
  V = asarray(V)

  if( U.shape[U.ndim-2] !== V.shape[V.ndim-1] )
    throw new Error('rrqr_solve(Q,R,P, y): System not square.')

  const x = svd_lstsq(U,sv,V, y),
      EPS = math.sqrt( eps(sv.dtype) ),
        N = sv.shape[sv.ndim-1]

  sv = sv.data

  for( let sv_off = 0;
           sv_off < sv.length;
           sv_off += N )
  {
    const T = EPS * math.abs(sv[sv_off])

    for( let r; r < N; r++ )
    {
      const sv_r = math.abs(sv[sv_off + r])
      if( ! isFinite(sv_r) )
        throw new Error('svd_solve(): NaN or Infinity encountered.')
      if( sv_r <= T )
        throw new SingularMatrixSolveError(x)
    }
  }

  return x
}


export function svd_lstsq(U,sv,V, y)
{
  if( y == undefined )
  {
    if( V != undefined )
      throw new Error('svd_lstsq(Q,R,P, y): Either 2 ([Q,R,P], y) or 4 arguments (Q,R,P, y) expected.')
    y = R
    ([U,sv,V] = Q)
  }

  U  = asarray(U ); if( U .ndim < 2 ) throw new Error('svd_lstsq(U,sv,V, y): U.ndim must be at least 2.' )
  sv = asarray(sv); if( sv.ndim < 1 ) throw new Error('svd_lstsq(U,sv,V, y): sv.ndim must be at least 1.')
  V  = asarray(V ); if( V .ndim < 2 ) throw new Error('svd_lstsq(U,sv,V, y): V.ndim must be at least 2.' )
  y  = asarray(y ); if( y .ndim < 2 ) throw new Error('svd_lstsq(U,sv,V, y): y.ndim must be at least 2.' )

  const
    [N,M] = U.shape.slice(-2),
    [I]   = V.shape.slice(-1),
    [J]   = y.shape.slice(-1)
  if( N !=  y.shape[ y.ndim-2] ) throw new Error("svd_lstsq(U,sv,V, y): U and y don't match.")
  if( M != sv.shape[sv.ndim-1] ) throw new Error("svd_lstsq(U,sv,V, y): U and sv don't match.")
  if( M !=  V.shape[ V.ndim-2] ) throw new Error("svd_lstsq(U,sv,V, y): V and sv don't match.")

  const ndim = Math.max(U.ndim, sv.ndim+1, V.ndim, y.ndim),
       shape = Int32Array.from({length: ndim}, () => 1);
  shape[ndim-2] = I;
  shape[ndim-1] = J;

  // FIND COMMON (BROADCASTED) SHAPE
  for( let arr of [U,V,y] )
    for( let i=ndim-2, j=arr.ndim-2; i-- > 0 && j-- > 0; )
      if( 1 === shape[i] )
        shape[i] = arr.shape[j];
      else if( shape[i] != arr.shape[j] && arr.shape[j] != 1 )
        throw new Error('svd_lstsq(U,sv,V, y): U,sv,V,y not broadcast-compatible.');

  for( let i=ndim-2, j=sv.ndim-1; i-- > 0 && j-- > 0; )
    if( 1 === shape[i] )
      shape[i] = sv.shape[j];
    else if( shape[i] != sv.shape[j] && sv.shape[j] != 1 )
      throw new Error('svd_lstsq(U,sv,V, y): U,sv,V,y not broadcast-compatible.');

  // GENERATE RESULT DATA
  const    EPS = math.sqrt(eps(sv.dtype)),
    DTypeArray = ARRAY_TYPES[ [U,sv,V,y].every( a => a.dtype==='float32' ) ? 'float32' : 'float64' ],
       tmp = new DTypeArray(M*J),
     x_dat = new DTypeArray( shape.reduce((a,b) => a*b, 1) ),
     U_dat =  U.data,
    sv_dat = sv.data,
     V_dat =  V.data,
     y_dat =  y.data;
  let
     U_off = 0,  U_stride = 1,
    sv_off = 0, sv_stride = 1,
     V_off = 0,  V_stride = 1,
     y_off = 0,  y_stride = 1,
     x_off = 0;

  function solv(d) {
    if( d === ndim-2 ) {
       U_stride = N*M;
      sv_stride =   M;
       V_stride = M*I;
       y_stride = N*J;

      const R = function(){
        const T = EPS * math.abs(sv_dat[sv_off])

        for( let r=0; r < M; r++ )
        {
          const sv_r = math.abs(sv_dat[sv_off + r])
          if( ! isFinite(sv_r) )
            throw new Error('svd_solve(): NaN or Infinity encountered.')
          if( sv_r <= T )
            return r
        }
        return M
      }()
      // SEE: Gene H. Golub, Charles F. Van Golub
      //     "Matrix Computations", 4th edition
      //      Page 260f, Chap. 5.3.1 (Implications of Full Rank)

      // tmp = U.T @ y
      tmp.fill(0)
      for( let k=0; k < N; k++ )
      for( let i=0; i < R; i++ )
      for( let j=0; j < J; j++ )
        tmp[J*i+j] += U_dat[U_off + M*k+i] * y_dat[y_off + J*k+j]

      // tmp \= diag(sv)
      for( let i=0; i < R; i++ )
      for( let j=0; j < J; j++ )
        tmp[J*i+j] /= sv_dat[sv_off + i]

      // x = V.T @ tmp
      for( let k=0; k < R; k++ )
      for( let i=0; i < I; i++ )
      for( let j=0; j < J; j++ )
        x_dat[x_off + J*i+j] += V_dat[V_off + I*k+i] * tmp[J*k+j]

       U_off +=  U_stride;
      sv_off += sv_stride;
       V_off +=  V_stride;
       y_off +=  y_stride;
       x_off += I*J;

      return;
    }
    for( let l=shape[d]; ; l-- ) {
      solv(d+1);
      if( l == 1 ) break;
      if( ! ( U.shape[ d - ndim   +  U.ndim ] > 1) )  U_off -=  U_stride;
      if( ! (sv.shape[ d - ndim+1 + sv.ndim ] > 1) ) sv_off -= sv_stride;
      if( ! ( V.shape[ d - ndim   +  V.ndim ] > 1) )  V_off -=  V_stride;
      if( ! ( y.shape[ d - ndim   +  y.ndim ] > 1) )  y_off -=  y_stride;
    }
     U_stride *=  U.shape[ d - ndim   +  U.ndim ] || 1;
    sv_stride *= sv.shape[ d - ndim+1 + sv.ndim ] || 1;
     V_stride *=  V.shape[ d - ndim   +  V.ndim ] || 1;
     y_stride *=  y.shape[ d - ndim   +  y.ndim ] || 1;
  }
  solv(0);

  return new NDArray(shape,x_dat);

  throw new Error('Not yet implemented!')
}
