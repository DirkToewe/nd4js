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

import {schur_decomp,
        schur_eigen,
        schur_eigenvals} from './schur'
import {asarray, NDArray} from '../nd_array'
import {ARRAY_TYPES} from '../dt'
import {zip_elems} from '../zip_elems'
import math from '../math'


export function eigen(A)
{
  const [D,B]= eigen_balance_pre(A,2),
        [Q,T]= schur_decomp(B),
        [Λ,V]= schur_eigen(Q,T)

  if( ! V.dtype.startsWith('complex') ) throw new Error('Assertion failed.')
  if( ! D.dtype.startsWith('float'  ) ) throw new Error('Assertion failed.')

  const N = A.shape[A.ndim-1],
        norm_sum = new Float64Array(N),
        norm_max = new Float64Array(N),
        V_dat = V.data._array,
        D_dat = D.data

  // UNDO BALANCING
  
  for( let D_off=0,  V_off=0;
                     V_off < V_dat.length;
           D_off+=N, V_off += 2*N*N )
  {
    norm_sum.fill(0)
    norm_max.fill(0)
    // SCALE ROWS & COMPUTE SCALED COLUMN NORMS
    for( let i=0; i <   N; i++ ) { const D_i = D_dat[D_off + i]
    for( let j=0; j < 2*N; j++ ) {
      // scale row
      const V_ij = Math.abs(V_dat[V_off + 2*N*i+j] *= D_i)
      // update norm
      if(   V_ij > 0 ) {
        const k = j>>1
        if( V_ij > norm_max[k] ) {
          norm_sum[k] *= (norm_max[k]/V_ij)**2; norm_max[k]  = V_ij
        } norm_sum[k] += (V_ij/norm_max[k])**2
      }
    }}
    for( let j=0; j < N; j++ ) {
      let max = norm_max[j];
      norm_sum[j] = isFinite(max) ? Math.sqrt(norm_sum[j])*max : max;
    }

    // NORMALIZE COLUMNS
    for( let i=0; i <   N; i++ )
    for( let j=0; j < 2*N; j++ )
      V_dat[V_off + 2*N*i+j] /= norm_sum[j>>1]
  }

  return [Λ,V]
}


export function eigenvals(A)
{
  const [D,B]= eigen_balance_pre(A,2),
        [Q,T]= schur_decomp(B);
  return schur_eigenvals(T);
}


export function eigen_balance_pre(A,p)
{
  // SEE: RODNEY JAMES, JULIEN LANGOU, AND BRADLEY R. LOWERY,
  //      "ON MATRIX BALANCING AND EIGENVECTOR COMPUTATION"
  if( p == null ) p = 2;
  if( p > Number.MAX_VALUE ) return eigen_balance_pre_inf(A);

  const N = A.shape[A.ndim-1],
        DTypeArray = ARRAY_TYPES[A.dtype==='float32' ? 'float32' : 'float64'],
        A_shape = A.shape,
        D_shape = A_shape.slice(0,-1);
  if( ! (p >= 1) ) throw new Error(`Invalid norm p=${p};`);
  if( A.shape[A.ndim-2] != N ) throw new Error('A is not square');
  A = DTypeArray.from(A.data); // <- protection copy

  const D = new DTypeArray(A.length/N),
      TOL = 0.95 ** (1/p);
  D.fill(1.0);

  for( let A_off=A.length,
           D_off=D.length; A_off > 0; ) { A_off -= N*N;
                                          D_off -= N;
    for( let i,   done=false; !done; )
    for(     i=0, done=true; i < N; i++ )
    { let c=0.0, c_max=0.0,
          r=0.0, r_max=0.0;
      // COMPUTE ROW AND COLUMN NORM
      for( let j=0; j < N; j++ )
//        if(true)
        if( i !== j )
        {{const A_ij = Math.abs(A[A_off + N*i+j]);
          if(   A_ij > 0 ) {
            if( A_ij > r_max ) {
              const scale = r_max / A_ij ; r_max = A_ij; r *= scale**p;
            } const ratio =         A_ij / r_max;        r += ratio**p;
          }
        }{const A_ji = Math.abs(A[A_off + N*j+i]);
          if(   A_ji > 0 ) {
            if( A_ji > c_max ) {
              const scale = c_max / A_ji ; c_max = A_ji; c *= scale**p;
            } const ratio =         A_ji / c_max;        c += ratio**p;
          }
        }}
      r = ! isFinite(r) ? r : r**(1/p) * r_max;
      c = ! isFinite(c) ? c : c**(1/p) * c_max;

      if( r*c == 0.0 ) continue;
      if( ! isFinite(r*c) ) throw new Error('NaN encountered.');
      let old_norm = ( c >= r
        ? ( 1 + (r/c)**p )**(1/p) * c
        : ( 1 + (c/r)**p )**(1/p) * r
      );

      // FIND A POWER OF TWO THAT MORE OR LESS BALANCES r AND c
      let scale = 1.0
      while( r >= c*2 ) { c*=2; r/=2; scale*=2; }
      while( c >= r*2 ) { c/=2; r*=2; scale/=2; }

      let new_norm = ( c >= r
        ? ( 1 + (r/c)**p )**(1/p) * c
        : ( 1 + (c/r)**p )**(1/p) * r
      );

      if( new_norm >= TOL*old_norm ) continue

      done = false;
      D[D_off + i] *= scale;
      for( let j=0; j < N; j++ ) {
        A[A_off + N*i+j] /= scale;
        A[A_off + N*j+i] *= scale;
      }
    }
  }

  return [
    new NDArray(D_shape, D),
    new NDArray(A_shape, A)
  ];
}


function eigen_balance_pre_inf(A)
{
  // SEE: RODNEY JAMES, JULIEN LANGOU, AND BRADLEY R. LOWERY,
  //      "ON MATRIX BALANCING AND EIGENVECTOR COMPUTATION"
  const N = A.shape[A.ndim-1],
        DTypeArray = ARRAY_TYPES[A.dtype==='float32' ? 'float32' : 'float64'],
        A_shape = A.shape,
        D_shape = A_shape.slice(0,-1);
  if( A.shape[A.ndim-2] != N ) throw new Error('A is not square');
  A = DTypeArray.from(A.data);

  const D = new DTypeArray(A.length/N);
  D.fill(1.0);

  for( let A_off=A.length,
           D_off=D.length; A_off > 0; ) { A_off -= N*N;
                                          D_off -= N;
    for( let i,   done=false; !done; )
    for(     i=0, done=true; i < N; i++ )
    { let c=0.0,
          r=0.0;
      // COMPUTE ROW AND COLUMN NORM
      for( let j=0; j < N; j++ )
        if( i !== j ) {
          const A_ij = Math.abs(A[A_off + N*i+j]); r = Math.max(r,A_ij);
          const A_ji = Math.abs(A[A_off + N*j+i]); c = Math.max(c,A_ji);
        }

      if( r*c === 0.0 ) continue;
      let old_norm = Math.max(c,r);
      if( old_norm == 0.0 ) continue;
      if( ! isFinite(old_norm) ) throw new Error('NaN encountered.');

      // FIND A POWER OF TWO THAT MORE OR LESS BALANCES r AND c
      let scale = 1.0
      while( r >= c*2 ) { c*=2; r/=2; scale*=2; }
      while( c >= r*2 ) { c/=2; r*=2; scale/=2; }

      let new_norm = Math.max(c,r);
      if( new_norm >= old_norm ) continue;

      done = false;
      D[D_off + i] *= scale;
      for( let j=0; j < N; j++ ) {
        A[A_off + N*i+j] /= scale;
        A[A_off + N*j+i] *= scale;
      }
    }
  }

  return [
    new NDArray(D_shape,D),
    new NDArray(A_shape,A)
  ];
}


export function eigen_balance_post(D,V)
{
  if( V.ndim < 2 ) throw new Error('eigen_balance_post(D,V): V.ndim must be at least 2.')

  const [M,N]= V.shape.slice(-2)
  if( M !== N ) throw new Error('eigen_balance_post(D,V): V must be square.')

  const DTypeArray = ARRAY_TYPES['float64'],
      V_arr = zip_elems([V,D.reshape(...D.shape,1)], 'complex128', math.mul)
  V = V_arr.data

  const norm_sum = new DTypeArray(N),
        norm_max = new DTypeArray(N)

  // NORMALIZE COLUMNS
  for( let off=0; off < V.length; off += M*N )
  {
    norm_sum.fill(0)
    norm_max.fill(0)
    // COMPUTE COLUMN NORMS
    for( let i=0; i < N; i++ )
    for( let j=0; j < N; j++ ) {
      const V_ij = math.abs(V[off + N*i+j])
      if(   V_ij > 0 ) {
        if( V_ij > norm_max[j] ) {
          norm_sum[j] *= (norm_max[j]/V_ij)**2; norm_max[j]  = V_ij
        } norm_sum[j] += (V_ij/norm_max[j])**2
      }
    }
    for( let j=0; j < N; j++ ) {
      let max = norm_max[j];
      norm_sum[j] = isFinite(max) ? Math.sqrt(norm_sum[j])*max : max;
    }

    // NORMALIZE COLUMNS
    for( let i=0; i < N; i++ )
    for( let j=0; j < N; j++ )
      V[off + N*i+j] = math.div(V[off + N*i+j], norm_sum[j])
  }

  return V_arr
}
