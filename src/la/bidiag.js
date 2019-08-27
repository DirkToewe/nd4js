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
import {super_dtype, ARRAY_TYPES} from '../dt'


// TODO consider Householder reflections for eliminating right side of a row
export function bidiag_decomp(A)
{
  A = asarray(A);
  if( A.ndim < 2 ) throw new Error('A must be at least 2D.');
  const
    DTypeArray = ARRAY_TYPES[ A.dtype === 'float32' ? 'float32' : 'float64' ], // <- ensure at least double precision
    U_shape = Int32Array.from(A.shape),
    B_shape = Int32Array.from(U_shape),
    V_shape = Int32Array.from(U_shape),
    [N,M] = A.shape.slice(  -2),
    len   = A.shape.slice(0,-2).reduce( (a,b) => a*b, 1 ),
    I     = Math.min(N,M),   // #rows    of bidiag. matrix B
    J     = N >= M ? I : I+1;// #columns of bidiag. matrix B
  U_shape[U_shape.length-1] = I;
  B_shape[B_shape.length-2] = I;
  B_shape[B_shape.length-1] = J;
  V_shape[V_shape.length-2] = J;

  let U,V;
  // copy A to U or V, whichever is large enough
  if( N >= M ) {
    U  =     DTypeArray.from(A.data); A = U;
    V  = new DTypeArray( len*M*M );
  } else {
    V  = new DTypeArray(len*(N+1)*M);
    A = A.data;
    for( let i=A.length,
             j=V.length; i-- > 0; ) V[--j] = A[i];
    A = V;
    U  = new DTypeArray(len*N*N);
  }

  const
    B  = new DTypeArray( len*I*J ),                                //<- resulting bidiagonal matrix
    csU= new DTypeArray(         I * (N-1 + N-I  ) ),              //<- caches ALL row    rotations later used to compute U
    csV= new DTypeArray( N < M ? I * (M-2 + M-I-1) : (I-1)*(M-2) ),//<- caches ALL column rotations later used to compute V
    uv = new DTypeArray( Math.max(N-I,M) );                        //<- cached remainder of a row of U or a column of V to compute U/V

  for(
    let A_off= N >= M ? 0 : len*M,
        U_off=0,
        B_off=0,
        V_off=0; B_off < B.length;
        A_off += N*M,
        U_off += N*I,
        B_off += I*J,
        V_off += J*M
  )
  {
    let
      csiU = 0,
      csiV = 0;
    for( let i=0; i < I; i++ )
    {
      // ELIMINATE THE ELEMENTS IN (i+1)-th COLUMN BELOW DIAG (IN B)
      for( let j=i+1; j < N; j++ )
      {
        const     A_ji = A[A_off + M*j+i]; if( 0.0 == A_ji ) { csU[csiU++]=1; csU[csiU++]=0; continue; }
        const     A_ii = A[A_off + M*i+i],
                         norm = Math.hypot(A_ii,A_ji),
              c = A_ii / norm,
              s = A_ji / norm;
        A[A_off + M*i+i] = norm;
        A[A_off + M*j+i] = 0;
        // ROTATE ROWS IN A
        for( let k=i; ++k < M; ) {
          const
            ik = A_off + M*i+k, A_ik = A[ik],
            jk = A_off + M*j+k, A_jk = A[jk];
          A[ik] = A_jk*s + A_ik*c;
          A[jk] = A_jk*c - A_ik*s;
        }
        csU[csiU++] = c;
        csU[csiU++] = s;
      }
      if( i+3 > M ) continue;
      // ELIMINATE (i+1)-th ROW RIGHT OF BIDIAG
      for( let j=i+2; j < M; j++ ) {
        const ii = A_off + M*i+(i+1), A_ii = A[ii],
              ij = A_off + M*i+ j   , A_ij = A[ij], norm = Math.hypot(A_ii,A_ij);
        A[ii] = norm;
        A[ij] = 0;
        csV[csiV++] = norm == 0 ? 1 : A_ii / norm;
        csV[csiV++] = norm == 0 ? 0 : A_ij / norm;
      }
      // APPLY ABOVE ELIMINATION TO REMAINING ROWS RIGHT OF BIDIAG (ROW-BY-ROW IMPROVES CACHE LOCALITY)
      for( let k=i+1; k < N; k++ ) { csiV -= 2*(M-i-2);
      for( let j=i+2; j < M; j++ ) {
          const
            c = csV[csiV++],
            s = csV[csiV++],
            ki = A_off + M*k+(i+1), A_ki = A[ki],
            kj = A_off + M*k+ j   , A_kj = A[kj];
          A[ki] = A_kj*s + A_ki*c;
          A[kj] = A_kj*c - A_ki*s;
      }}
    }

    if( csiU != csU.length ) throw new Error('Assertion failed: (csiU='+csiU+") != (csU.length="+csU.length+")" );
    if( csiV != csV.length ) throw new Error('Assertion failed: (csiV='+csiV+") != (csV.length="+csV.length+")" );

    // MOVE A -> B
    for( let i=0; i < I; i++ )
    for( let j=0; j < J; j++ )
      B[B_off + J*i+j] = A[A_off + M*i+j];

     //
    // COMPUTE U
   //
    // ROTATE U
    for( let k = 0; k < N; k++ )
    { csiU = 0
      for( let i=0; i < N-I; i++ )          uv[i] = k != (i+I) ? 0.0 : 1.0;
      for( let i=0; i < I  ; i++ ) U[U_off+I*k+i] = k !=  i    ? 0.0 : 1.0;
      for( let i=0; i < I  ; i++ ) {
        const skip = Math.max(1,k-i);
        csiU += 2*(skip-1)
        for( let j = i+skip; j < I; j++ ) {
          const
            c = csU[csiU++],
            s = csU[csiU++],
            ki = U_off + I*k+i, U_ki = U[ki],
            kj = U_off + I*k+j, U_kj = U[kj];
          U[ki] = U_kj*s + U_ki*c;
          U[kj] = U_kj*c - U_ki*s;              
        }
        for( let j = Math.max(0,i+skip-I); j < (N-I); j++ ) {
          const
            c = csU[csiU++],
            s = csU[csiU++], ki = U_off + I *k+i, U_ki = U[ki];
          U[ki] = uv[j]*s + U_ki*c;
          uv[j] = uv[j]*c - U_ki*s;
        }
      }
    }
     //
    // COMPUTE V
   //
    // ROTATE V
    for( let k = 0; k < M; k++ )
    { csiV = 0
      // uv CACHES (k+1)-th COLUMN OF V
      uv.fill(0.0, 0,M)
      uv[k] = 1
      for( let i=0; i < I; i++ )
      {
        const skip = Math.max(2,k-i);
        csiV += 2*(skip-2);
        for( let j=i+skip; j < M; j++ ) {
          const
            c = csV[csiV++],
            s = csV[csiV++],
            v_i = uv[i+1],
            v_j = uv[j  ];
          uv[i+1] = v_j*s + v_i*c;
          uv[j  ] = v_j*c - v_i*s;
        }
      }
      // WRITE (k+1)-th COLUMN
      for( let i=0; i < J; i++ ) V[V_off + M*i+k] = uv[i];
    }
  }

  return [
    new NDArray(U_shape, U),
    new NDArray(B_shape, B),
    new NDArray(V_shape, V)
  ];
}
