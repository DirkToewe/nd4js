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
import {matmul2} from './matmul'
import {qr_decomp} from './qr'
import {ARRAY_TYPES, eps} from '../dt'
import {transpose_inplace} from './transpose_inplace'

import math from '../math';


export function svd_jac_2sided(A)
{
  A = asarray(A);
  if( A.dtype.startsWith('complex') )
    throw new Error('svd_jac_1sided(A): A.dtype must be float.');
  const shape = A.shape,
    N = shape[shape.length-2];

  // QR DECOMPOSE RECTANGULAR MATRICES AND WORK ON (QUADRATIC) R
  {
    const M = shape[shape.length-1];
    // if A is not square use QR Decomposition
    if( N > M ) {
      const [Q,R] = qr_decomp(A),
         [U,sv,V] = svd_jac_2sided(R)
      return [matmul2(Q,U), sv, V]
    }
    if( N < M ) {
      const [Q,R] = qr_decomp(A.T),
         [U,sv,V] = svd_jac_2sided(R)
      transpose_inplace(V)
      return [V, sv, matmul2(Q,U).T]
    }
  }
  // ALLOCATE RESULT DATA
  const DType = A.dtype==='float32' ? 'float32' : 'float64',
        DTypeArray = ARRAY_TYPES[DType],
        TOL = eps(DType),
        U = DTypeArray.from(A.data); A = undefined; // <- potentially allow GC
  const D = new DTypeArray(N*N), // <- tempory storage for decomposition
        V = new DTypeArray(U.length),
        sv= new DTypeArray(U.length/N);

  if( 1 >  N ) throw new Error('Assertion failed.');
  if( 1 == N ) {
    for( let i=U.length; i-- > 0; )
      if( U[i] < +0.0 ) {
          U[i] *= -1.0;
         sv[i]  = -1.0;
      }
      else sv[i] = +1.0;
    return [
      new NDArray(shape,            sv ),
      new NDArray(shape.slice(0,-1), U ),
      new NDArray(shape,      V.fill(1))
    ];
  }


  for( let UV_off=0,
           sv_off=0; sv_off < sv.length; UV_off += N*N,
                                         sv_off += N )
  {
    // MOVE FROM U TO D
    for( let i=0; i < N; i++ )
    for( let j=0; j < N; j++ ) { D[N*i+j] = U[UV_off + N*i+j]; U[UV_off + N*i+j] = i != j ? 0 : 1 };
    // INIT V TO IDENTITY
    for( let i=0; i < N; i++ )
    for( let j=0; j < N; j++ ) V[UV_off + N*i+j] = i != j ? 0 : 1;

     //
    // JACOBI SVD ITERATIONS
   //
    let finished = false;
    while( ! finished )
    {
      finished = true;

      for( let k=0; k < N; k++ )
      for( let l=0; l < k; l++ )
      {
        const D_kk = D[N*k+k], D_kl = D[N*k+l],
              D_lk = D[N*l+k], D_ll = D[N*l+l];
        if( ! ( Math.hypot(D_kl,D_lk) / Math.max( Math.abs(D_kk), Math.abs(D_ll) ) > TOL ) )
          continue; // <- TODO check if really a good stopping criterion (there may be smaller off-diag. entries larger relative to their respective diag. entries)

        finished = false;

        // determine rotation angles
        // ┌                ┐ ┌            ┐ ┌                ┐   ┌        ┐
        // │ cos(β) -sin(β) │ │ D_kk  D_kl │ │ cos(α) -sin(α) │ ! │ d1   0 │
        // │                │ │            │ │                │ = │        │
        // │ sin(β)  cos(β) │ │ D_lk  D_ll │ │ sin(α)  cos(α) │   │  0  d2 │
        // └                ┘ └            ┘ └                ┘   └        ┘
        //
        // such that: |d1| >= |d2|
        //
        // => 0 = cos(α)*{D_kl*cos(β) - D_ll*sin(β)} + sin(α)*{D_kk*cos(β) - D_lk*sin(β)} = ½⋅{(D_ll - D_kk)⋅sin(α-β) - (D_kk + D_ll)⋅sin(α+β) + (D_kl + D_lk)⋅cos(α-β) + (D_kl - D_lk)⋅cos(α+β)}
        //    0 = cos(α)*{D_kk*sin(β) + D_lk*cos(β)} - sin(α)*{D_kl*sin(β) + D_ll*cos(β)} = ½⋅{(D_ll - D_kk)⋅sin(α-β) + (D_kk + D_ll)⋅sin(α+β) + (D_kl + D_lk)⋅cos(α-β) - (D_kl - D_lk)⋅cos(α+β)}
        //   d1 = cos(α)*{D_kk*cos(β) - D_lk*sin(β)} - sin(α)*{D_kl*cos(β) - D_ll*sin(β)}
        //   d2 = cos(α)*{D_kl*sin(β) + D_ll*cos(β)} + sin(α)*{D_kk*sin(β) + D_lk*cos(β)}
        //
        // => 0 = {D_kl⋅cos(α) - D_kk⋅sin(α)}⋅cos(β) + {D_lk⋅sin(α) - D_ll⋅cos(α)}⋅sin(β)
        //    0 = {D_kl⋅sin(α) + D_kk⋅cos(α)}⋅sin(β) + {D_lk⋅cos(α) + D_ll⋅sin(α)}⋅cos(β)
        //   d1 = {D_kl⋅sin(α) + D_kk⋅cos(α)}⋅cos(β) + {D_lk⋅cos(α) - D_ll⋅sin(α)}⋅sin(β)
        //   d2 = {D_kl⋅cos(α) - D_kk⋅sin(α)}⋅sin(β) + {D_lk⋅sin(α) + D_ll⋅cos(α)}⋅cos(β)
        //
        // => 0 = (D_kk+D_ll)⋅sin(α-β) + (D_kl-D_lk)⋅cos(α-β)  
        //    0 = (D_kk-D_ll)⋅sin(α+β) + (D_kl+D_lk)⋅cos(α+β)  
        const [cα,sα,cβ,sβ] = function(){
          let Cα,Sα,Cβ,Sβ;
  
          // ANGLE COMPUTATION VARIANT A (SEE ASCII ART ABOVE)
          {
            const m = Math.atan2(D_lk - D_kl, D_ll + D_kk),// = α - β
                  p = Math.atan2(D_lk + D_kl, D_ll - D_kk),// = α + β
                  α = (p+m)/2,
                  β = (p-m)/2;
                 Cα = Math.cos(α); Sα = Math.sin(α);
                 Cβ = Math.cos(β); Sβ = Math.sin(β);
          }
  
          // tan is 180°-periodical so lets try all possible solutions
          for( const [cα,sα,cβ,sβ] of [
            [ Cα, Sα,   Cβ, Sβ], // [p,m]+[  0°,  0°]
            [-Sα, Cα,   Sβ,-Cβ], // [p,m]+[  0°,180°]
            [-Sα, Cα,  -Sβ, Cβ], // [p,m]+[180°,  0°]
            [-Cα,-Sα,   Cβ, Sβ]  // [p,m]+[180°,180°]
          ]) {
            const d00 = (D_kl*sα + D_kk*cα)*cβ + (D_lk*cα - D_ll*sα)*sβ,
                  d11 = (D_kl*cα - D_kk*sα)*sβ + (D_lk*sα + D_ll*cα)*cβ;
            // ROTATE IN A WAY THAT ENSURES DESCENDING ORDER
            if( d11 >= Math.abs(d00) ) {
              if( d00 >= -0.0 )
                return [cα,sα,cβ,sβ];
              Cα = cα; Sα = sα;
              Cβ = cβ; Sβ = sβ;
            }
          }
          return [Cα,Sα,Cβ,Sβ];
        }();
  
        // ROTATE COLUMNS IN D
        for( let i=N; i-- > 0; ) {
          const D_ik = D[N*i+k],
                D_il = D[N*i+l];
          D[N*i+k] = D_ik*cα - D_il*sα;
          D[N*i+l] = D_il*cα + D_ik*sα;
        }
  
        // ROTATE ROWS IN D
        for( let i=N; i-- > 0; ) {
          const D_ki = D[N*k+i],
                D_li = D[N*l+i];
          D[N*k+i] = D_ki*cβ - D_li*sβ;
          D[N*l+i] = D_li*cβ + D_ki*sβ;
        }
  
        if( ! math.is_close(0, D[N*k+l]) ) throw new Error(`Assertion failed: 0 =/= ${D[N*k+l]}.`)
        if( ! math.is_close(0, D[N*l+k]) ) throw new Error(`Assertion failed: 0 =/= ${D[N*l+k]}.`)
        // ENTRIES (k,l) AND (l,k) ARE REMAINDERS (CANCELLATION ERROR) FROM ELIMINATION => SHOULD BE SAFELY ZEROABLE
        D[N*k+l] = 0.0;
        D[N*l+k] = 0.0;
  
        // ROTATE ROWS IN U (TRANSPOSED FOR CACHE LOCALITY REASONS)
        for( let i=0; i < N; i++ ) {
          const ki = UV_off + N*k+i, U_ki = U[ki],
                li = UV_off + N*l+i, U_li = U[li];
          U[ki] = U_ki*cβ - U_li*sβ;
          U[li] = U_li*cβ + U_ki*sβ;
        }
  
        // ROTATE ROWS IN V
        for( let i=0; i < N; i++ ) {
          const ki = UV_off + N*k+i, V_ki = V[ki],
                li = UV_off + N*l+i, V_li = V[li];
          V[ki] = V_ki*cα - V_li*sα;
          V[li] = V_li*cα + V_ki*sα;
        }
      }
    }

    // MOVE D TO SV
    for( let i=0; i < N; i++ ) sv[sv_off + i] = D[N*i+i];

    // USE INSERTION SORT SV (should take as most O(N) because it is heavily pre-sorted by choosing the rotation angles accordingly)
    for( let i=0; i < N; i++ )
    {
      let sv_i = sv[sv_off + i];
      // flip sign
      if( sv_i < 0.0 ) {
        sv[sv_off + i] = (sv_i *= -1);
        for( let k=0; k < N; k++ ) U[UV_off + N*i+k] *= -1;
      }

      // insertion sort
      for( let j=i; j-- > 0; ) {
        const sv_j = sv[sv_off + j];
        if( sv_j >= sv_i ) break;
        // swap sv
        sv[sv_off + j+1] = sv[sv_off + j];
                           sv[sv_off + j] = sv_i;
        // swap U and V rows
        const rowJ = UV_off + N*j,
              rowI =  rowJ  + N;
        for( let k=0; k < N; k++ ) { const U_ik = U[rowI+k]; U[rowI+k] = U[rowJ+k]; U[rowJ+k] = U_ik; }
        for( let k=0; k < N; k++ ) { const V_ik = V[rowI+k]; V[rowI+k] = V[rowJ+k]; V[rowJ+k] = V_ik; }
      }
    }

    // TRANSPOSE U IN-PLACE (TRANSPOSED FOR CACHE LOCALITY REASONS)
    for( let i=0  ; i < N-1; i++ )
    for( let j=1+i; j < N  ; j++ ) {
      const
        ij = UV_off + N*i+j,
        ji = UV_off + N*j+i, U_ij = U[ij];
                                    U[ij] = U[ji];
                                            U[ji] = U_ij;
    }
  }

  return [
    new NDArray(shape,             U),
    new NDArray(shape.slice(0,-1),sv),
    new NDArray(shape,             V)
  ]
}
