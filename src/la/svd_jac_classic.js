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

import {ARRAY_TYPES, eps} from '../dt'
import {asarray, NDArray} from '../nd_array'
import {matmul2} from './matmul'
import {qr_decomp} from './qr'
import {transpose_inplace} from './transpose_inplace'
import {_svd_jac_angles,
        _svd_jac_post } from './_svd_jac_utils'
import {_giv_rot_rows,
        _giv_rot_cols} from './_giv_rot'


export function svd_jac_classic(A)
{
  A = asarray(A);
  if( A.dtype.startsWith('complex') )
    throw new Error('svd_jac_1sided(A): A.dtype must be float.');
  const
    shape = A.shape,
    N = shape[shape.length-2];

  // QR DECOMPOSE RECTANGULAR MATRICES AND WORK ON (QUADRATIC) R
  {
    const M = shape[shape.length-1];
    // if A is not square use QR Decomposition
    if( N > M ) {
      const [Q,R] = qr_decomp(A),
         [U,sv,V] = svd_jac_classic(R)
      return [matmul2(Q,U), sv, V]
    }
    if( N < M ) {
      const [Q,R] = qr_decomp(A.T),
         [U,sv,V] = svd_jac_classic(R)
      transpose_inplace(V)
      return [V, sv, matmul2(Q,U).T]
    }
  }
  // ALLOCATE RESULT DATA
  const DType = A.dtype==='float32' ? 'float32' : 'float64',
        DTypeArray = ARRAY_TYPES[DType],
        TOL = eps(DType) * N,
        U = DTypeArray.from(A.data); A = undefined; // <- potentially allow GC
  const S = new DTypeArray(N*N), // <- tempory storage for decomposition
        V = new DTypeArray(U.length),
       sv = new DTypeArray(U.length/N),
     diag = new DTypeArray(N),
      ord = Int32Array.from({length: N}, (_,i) => i);

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

   //
  // BUILD TRIANGLE TREE
 //
  // the size of each level of the tree
  const treeSizes = Int32Array.from( function*(){
    for( let n=N; n > 1; ) {
      n = n+1 >>> 1;
      yield n;
    }
  }() );
  const treeData = new DTypeArray( treeSizes.reduce( (len,N) => len + (N*N+N >>> 1), 0 ) );

//  const piv = (i,j) => {
//    const S_ij = S[N*i+j],
//          S_ji = S[N*j+i];
//    return S_ij*S_ij + S_ji*S_ji;
//  };

  const find_pivot = () => {
    let k=0,
        l=0,
        i=0,
        j=0;

    for( let off=treeData .length,
              h =treeSizes.length; h-- > 0; )
    {
      const val = (k,l) => {
        k += i;
        l += j; return treeData[off + (k*k+k >>> 1) + l];
      };
      const n = treeSizes[h];
      off -= n*n+n >>> 1;       
      let max = val(0,0); k=i;
                          l=j;
      if( i+1 < n ) { const v10 = val(1,0); if( v10 > max ) { max = v10; k=i+1; l=j  ; }
                      const v11 = val(1,1); if( v11 > max ) { max = v11; k=i+1; l=j+1; } }
      if( i != j )  { const v01 = val(0,1); if( v01 > max ) { max = v01; k=i;   l=j+1; } }
      i = 2*k;
      j = 2*l;
    }
    let max = -Infinity;
    const hyp = (s,t) => {
      s += i;
      t += j; const S_st = S[N*s+t],
                    S_ts = S[N*t+s]; return S_st*S_st + S_ts*S_ts;
    };
    if( i+1 < N ) { const h10 = hyp(1,0); if( h10 > max ) { max=h10; k=i+1; l=j  ; } }
    if( j+1 < N ) { const h01 = hyp(0,1); if( h01 > max ) { max=h01; k=i  ; l=j+1; } }
    if( i != j  ) { const h00 = hyp(0,0); if( h00 > max ) { max=h00; k=i  ; l=j  ; }
    if( i+1 < N ) { const h11 = hyp(1,1); if( h11 > max ) {          k=i+1; l=j+1; } } }

    return [l,k];
  };

  /** updates the specified row in the triangle tree.
   */
  const update_row = row =>
  {
    row = row >>> 1 << 1; // <- round down to pow2
    // build bottom tree level
    for( let i=row; i < row+2 && i < N; i++ ) { const r = i>>>1, off = r*r+r >>> 1;
    for( let j=0;   j < row+1         ; j++ ) {
      const x = S[N*i+j],
            y = S[N*j+i], S_ij = i===j ? -Infinity : x*x + y*y,
            k = off + (j >>> 1);
      treeData[k] = i%2 || j%2 ? Math.max(treeData[k],S_ij) : S_ij;
    }}
    // build remaining tree levels
    for( let off=0, h=1; h < treeSizes.length; h++ )
    {
      const N = treeSizes[h-1],
          OFF = off; off += N*N+N >>> 1;
      row = row >>> 2 << 1;
      for( let R=row; R < row+2 && R < N; R++ ) { const ROFF = OFF + (R*R+R >>> 1), r = R >>> 1,
                                                        roff = off + (r*r+r >>> 1);
      for( let C= 0 ; C <= R            ; C++ ) {
        const k =          roff +(C >>> 1),
              T = treeData[ROFF + C];
        treeData[k] = R%2 || C%2 ? Math.max(treeData[k],T) : T;
      }}
    }
  }

  /** updates the specified column in the triangle tree.
   */
  const update_col = col =>
  {
    col = col >>> 1 << 1; // <- round down to pow2
    // build bottom tree level
    const J = Math.min(col+2,N);
    for( let i=col; i < N; i++ ) { const r = i>>>1, off = r*r+r >>> 1;
    for( let j=col; j < J; j++ ) {
      const x = S[N*i+j],
            y = S[N*j+i], S_ij = i===j ? -Infinity : x*x + y*y,
            k = off + (j >>> 1);
      treeData[k] = i%2 || j%2 ? Math.max(treeData[k],S_ij) : S_ij;
    }}
    // build remaining tree levels
    for( let off=0, h=1; h < treeSizes.length; h++ )
    {
      col = col >>> 2 << 1;
      const N = treeSizes[h-1],
            J = Math.min(col+2,N),
          OFF = off; off += N*N+N >>> 1;
      for( let R=col; R < N          ; R++ ) { const ROFF = OFF + (R*R+R >>> 1), r = R >>> 1,
                                                     roff = off + (r*r+r >>> 1);
      for( let C=col; C < J && C <= R; C++ ) {
        const k =          roff +(C >>> 1),
              T = treeData[ROFF + C];
        treeData[k] = R%2 || C%2 ? Math.max(treeData[k],T) : T;
      }}
    }
  }

  for( let UV_off=0,
           sv_off=0; sv_off < sv.length; UV_off += N*N,
                                         sv_off += N )
  {
    // MOVE FROM U TO S
    for( let i=0; i < N; i++ )
    for( let j=0; j < N; j++ ) { S[N*i+j] = U[UV_off + N*i+j];
                                            U[UV_off + N*i+j] = +(i===j) };
    // INIT V TO IDENTITY
    for( let i=0; i < N; i++ )
    for( let j=0; j < N; j++ ) V[UV_off + N*i+j] = +(i===j);
    // INIT TRIANGLE TREE
    for( let i=0; i < N; i += 2 ) update_row(i);

     //
    // (CLASSICAL) JACOBI SVD ITERATIONS
   //
    for( let s=-1,
             t=-1;; )
    {
      // FIND THE OFF DIAGONAL PAIR WITH THE LARGEST HYPOTHENUSE
      const [l,k] = find_pivot();
      if( l >= k ) throw new Error('Assertion failed.')
      if( s == k &&
          t == l ) break; // <- DRY PRINCIPLE
      s=k;
      t=l;

//      // CHECK THAT THIS IS TRULY THE MAXIMUM (TODO: COMMENT OUT)
//      for( let i=1; i < N; i++ )
//      for( let j=0; j < i; j++ )
//        if( ! (piv(i,j) <= piv(k,l)) ) throw new Error(`Assertion failed: ${i}, ${j}.`);

      const
        S_kk = S[N*k+k],
        S_kl = S[N*k+l],
        S_lk = S[N*l+k],
        S_ll = S[N*l+l];
      // stopping criterion inspiredy by:
      //  "Jacobi's Method is More Accurate than QR"
      //   by James Demmel
      //   SIAM J. Matrix Anal. Appl, vol. 13, pp. 1204-1245, 1992
      if( ! ( Math.max( Math.abs(S_kl), Math.abs(S_lk) ) > Math.sqrt(Math.abs(S_kk*S_ll)) * TOL ) )
        break;
 
      const [cα,sα,cβ,sβ] = _svd_jac_angles(S_ll, S_lk,
                                            S_kl, S_kk);

      // ROTATE S
      _giv_rot_rows(S, N, N*l,N*k, cα,sα);
      _giv_rot_cols(S, N,   l,  k, cβ,sβ);

      // ENTRIES (k,l) AND (l,k) ARE REMAINDERS (CANCELLATION ERROR) FROM ELIMINATION => SHOULD BE SAFELY ZEROABLE
      S[N*k+l] = 0.0;
      S[N*l+k] = 0.0;

      // UPDATE TRIANGLE TREE ROWS AND COLUMNS
      update_row(k); update_row(l);
      update_col(k); update_col(l);

      // ROTATE U & V
      _giv_rot_rows(U, N, UV_off + N*l,
                          UV_off + N*k, cα, sα);
      _giv_rot_rows(V, N, UV_off + N*l,
                          UV_off + N*k, cβ,-sβ);
    }

    _svd_jac_post( N, U,S,V, UV_off, sv, sv_off, ord );
  }

  return [
    new NDArray(shape,             U),
    new NDArray(shape.slice(0,-1),sv),
    new NDArray(shape,             V)
  ]
}
