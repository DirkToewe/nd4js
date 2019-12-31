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

import {eps, ARRAY_TYPES} from '../dt'
import math from '../math'
import {asarray, NDArray} from '../nd_array'

import {_giv_rot_rows} from './_giv_rot'
import {FrobeniusNorm} from './norm'
import {_norm,
        _norm_update} from './rrqr'
import {_transpose_inplace} from './transpose_inplace'


export function srrqr_decomp_full( A, opt={} )
{
  // SEE: Ming Gu, Stanley C. Eisenstat,
  //     "EFFICIENT ALGORITHMS FOR COMPUTING A STRONG RANK-REVEALING QR FACTORIZATION"
  //      https://math.berkeley.edu/~mgu/MA273/Strong_RRQR.pdf

  A = asarray(A)
  if( A.ndim < 2 ) throw new Error('srrqr_decomp_full(A,opt): A must be at least 2D.')
  const
    dtype = A.dtype==='float32' ? 'float32' : 'float64',
    DTypeArray = ARRAY_TYPES[dtype],
    R_shape =                 A.shape,
    Q_shape = Int32Array.from(R_shape),
    P_shape =                 Q_shape.slice(0,-1),
    [M,N]   =                 R_shape.slice(  -2),
       K    = Math.min(M,N); // <- M not M-1, because the last row still needs pivotization
  Q_shape[Q_shape.length-1] = M;
  P_shape[P_shape.length-1] = N;

  const R = DTypeArray.from(A.data); // <- we could encourage GC by setting `A = undefined` after this line
//                            A = undefined
  const P = new Int32Array(R.length/M), // <- tracks column permutations
        Q = new DTypeArray(R.length/N*M),
       AB = new DTypeArray(M*N), // 
     norm = new DTypeArray(N<<1), // <─ underflow-safe representation of the column norm
     NORM = new FrobeniusNorm();

  let {
    dtol = 1.1 * K**0.1,
    ztol = eps(dtype) * R.reduce((x,y) => Math.max(Math.abs(y),x), 0)
  } = opt;
                       if(!(dtol >=1)) throw new Error(`srrqr_decomp_full(A,opt): invalid opt.dtol: ${dtol}. Must be 1 or a greater number.`);
  dtol = Number(dtol); if(!(ztol >=0)) throw new Error(`srrqr_decomp_full(A,opt): invalid opt.ztol: ${ztol}. Must be a non-negative number.`);

  for(
    let Q_off=0,
        R_off=0,
        P_off=0; Q_off < Q.length; Q_off += M*M,
                                   R_off += M*N,
                                   P_off +=   N
  )
  {
    // INIT P
    for( let i=0; i < N; i++ ) P[P_off + i] = i;

    // INIT Q (TO IDENTITY)
    for( let i=0; i < M; i++ ) Q[Q_off + M*i+i] = 1;

    // INIT COLUMN NORM
    norm.fill(0);
    for( let i=0; i < M; i++ )
      _norm_update(norm, R, R_off + N*i, 0);

    // INIT AB
    AB.fill(0);

    let swapped = false;
//    let predict_detA = undefined;

    // ELIMINATE COLUMN BY COLUMN OF R
    outer_loop:for( let k=0; k < K; )
    { // FIND PIVOT COLUMN THAT (HOPEFULLY) ENSURES RANK REVEAL (RRQR is inherently not guaranteed to do that)
      let    p = -1,
        norm_p = -Infinity;
      for( let j=k; j < N; j++ ) {
        const  norm_j =_norm(norm,j);
        if(    norm_p < norm_j ) {
          p=j; norm_p = norm_j
        }
      }
      // swap pivot to column k
      if( p !== k ) {
        // swap columns in R
        for( let j=0; j < M; j++ ) {
          const R_jk = R[R_off + N*j+k];
                       R[R_off + N*j+k] = R[R_off + N*j+p];
                                          R[R_off + N*j+p] = R_jk;
        }
        // swap columns in A\B
        for( let j=0; j < k; j++ ) {
          const AB_jk = AB[N*j+k];
                        AB[N*j+k] = AB[N*j+p];
                                    AB[N*j+p] = AB_jk;
        }
        // swap P
        const P_i = P[P_off+k];
                    P[P_off+k] = P[P_off+p];
                                 P[P_off+p] = P_i;
      }

      if( Math.sqrt(N-k)*norm_p <= ztol ) // <- rank-deficient case TODO: consider finishing the triangulation
        break outer_loop;

      inner_loop:while(true)
      {
        // RESET COLUMN NORM (INDEX k IS SET TO ZERO FOR THE NEXT RRQR)
        norm.fill(0.0, (k+1)<<1);

        // ELIMINATE COLUMN k BELOW DIAGONAL (USING GIVEN ROTATIONS)
        let count = 0;
        for( let j=k; ++j < M; )
        { const       kk = R_off + N*k+k,
                      jk = R_off + N*j+k,
                    R_jk = R[jk];
          if( 0 !== R_jk )
          {   const R_kk = R[kk],
                           norm = Math.hypot(R_kk,R_jk),
                c = R_kk / norm,
                s = R_jk / norm; R[jk] = 0;
            if( s !== 0 ) {      R[kk] = norm;
              _giv_rot_rows(R, N-1-k, kk+1,
                                      jk+1, c,s);
              _giv_rot_rows(
                Q, swapped ? M : j+1,
                Q_off + M*k,
                Q_off + M*j, c,s
              );
            }
          }
          _norm_update(norm, R, R_off + N*j, k+1); // <- keep track of column norm
        }

        if( k >= K-1 ) break outer_loop;

//        const detA = function(){
//          let detA = 0;
//          for( let i=0; i <= k; i++ )
//            detA += Math.log2(Math.abs(R[R_off + N*i+i]));
//          return detA;
//        }();

        // UPDATE inv(A)
        { const                                   R_kk = - R[R_off + N*k+k];
                                   AB[N*k+k]=-1 / R_kk;
          for( let i=k; i-- > 0; ) AB[N*i+k]    /=R_kk;
        }

        // UPDATE A\B
        for( let i=0;   i <= k; i++ )
        for( let j=k; ++j <  N;     )
          AB[N*i+j] += AB[N*i+k] * R[R_off + N*k+j];

        k += 1;

//        // ASSERTIONS
//        // check triangularity of R
//        for( let i=1; i < M;          i++ )
//        for( let j=0; j < k && j < i; j++ )
//          if( !(Math.abs(R[R_off + N*i+j]) <= 1e-8) )
//            throw new Error('Assertion failed.');
//        // check triangularity of inv(A)
//        for( let i=1; i < M;          i++ )
//        for( let j=0; j < k && j < i; j++ )
//          if( !(Math.abs(AB[N*i+j]) <= 1e-6) )
//            throw new Error('Assertion failed.');
//        // check orthogonality of Q
//        for( let i=0; i < M; i++ )
//        for( let j=0; j < M; j++ ) {
//          let sum = 0;
//          for( let h=0; h < M; h++ )
//            sum += Q[Q_off + M*i+h] * Q[Q_off + M*j+h];
//          if( !(Math.abs(sum - (i===j)) <= 1e-8) )
//            throw new Error(`Assertion failed: ${sum} != ${(i===j)*1}.`);
//        }
//        // check Q @ R = Input (keep in mind that Q ist still column major)
//        for( let i=0; i < M; i++ )
//        for( let j=0; j < N; j++ )
//        {
//          const J = P[P_off + j];
//
//          let sum = 0;
//          for( let h=0; h < M; h++ )
//            sum += Q[Q_off + M*h+i] * R[R_off + N*h+j];
//
//          if( !(Math.abs(sum - A.data[R_off + N*i+J]) <= 1e-8) )
//            throw new Error(`Assertion failed: ${sum} =/= ${A.data[R_off + N*i+J]}.`);
//        }
//        // check inv(A)
//        for( let i=0; i < k; i++ )
//        for( let j=0; j < k; j++ ) {
//          let sum = 0;
//          for( let h=0; h < k; h++ )
//            sum += R[R_off + N*i+h] * AB[N*h+j];
//          if( !(Math.abs(sum - (i===j)) <= 1e-8) ) throw new Error(`Assertion failed: ${sum} != ${(i===j)*1}.`);
//        }
//        // check A\B
//        for( let i=0; i < k; i++ )
//        for( let j=k; j < N; j++ ) {
//          let sum = 0;
//          for( let h=0; h < k; h++ )
//            sum += AB[N*i+h] * R[R_off + N*h+j];
//
//          if( !(Math.abs(sum - AB[N*i+j]) <= 1e-8) )
//            throw new Error(`Assertion failed: ${sum} =/= ${AB[N*i+j]}.`);
//        }

//        // check determinant prediction
//        if( undefined !== predict_detA )
//          if( !(Math.abs(detA - predict_detA) <= 1e-8) )
//            throw new Error(`Assertion failed detA=${detA} != predict_detA=${predict_detA}.`);

        // SEARCH BEST COLUMN SWAPS
        let p = -1,
            q = -1,
            F = -Infinity;

        for( let i=0; i < k; i++ )
        {
          NORM.reset();
          for( let j=i; j < k; j++ )
            NORM.include(AB[N*i+j]);
                                       const row_norm =  NORM.result;
          for( let j=k; j < N; j++ ) { const col_norm = _norm(norm,j);
            const   f = Math.hypot( AB[N*i+j], row_norm*col_norm );
            if( F < f )
               [F,p,q] = [f,i,j];
          }
        }

        // IF NO GOOD COLUMN SWAP FOUND, START NEXT COLUMN
        if( !(F > dtol) ) { // <- handles NaN... I hope...
//          predict_detA = undefined;
          break inner_loop;
        }

        swapped = true;
        k -= 1; // <- go back one step since the swappend column needs retriangulation

        predict_detA = detA + Math.log2(F);

        // MOVE COLUMN p TO k (VIA CYCLIC PERMUTATION) TODO use triangulary property to reduces Ops
        // CYCLE COLUMNS of R
        for( let i=0; i <= k; i++ ) {
          const                             R_ip = R[R_off + N*i+p];
          for( let j=Math.max(p,i-1); j < k; j++ ) R[R_off + N*i+j] = R[R_off + N*i+(j+1)];
                                                   R[R_off + N*i+k] = R_ip;
        }
        // CYCLE ROWS OF inv(A) TODO: adjust to memory alignment!
        for( let j=p; j < N; j++ ) norm[j]=AB[N*p+j]; // move (row p) -> norm

        for( let i=p; i < k; i++ )
        for( let j=i; j < N; j++ ) AB[N*i+j] = AB[N*(i+1)+j];

        for( let j=p; j < N; j++ ) AB[N*k+j] = norm[j]; // move norm -> (row k)
        // CYCLE P P
        { const                P_p = P[P_off + p];
          for( let j=p; j < k; j++ ) P[P_off + j] = P[P_off + j+1];
                                     P[P_off + k] = P_p;
        }

        // RETRIANGULATE USING GIVENS ROTATIONS
        // (since cyclic permutation turned R from triangular to Hessenberg)
        for( let i=p; i < k; i++ )
        { const       ii = R_off + N* i   +i,
                      ji = R_off + N*(i+1)+i,
                    R_ji = R[ji];
          if( 0 !== R_ji )
          {   const R_ii = R[ii],
                           norm = Math.hypot(R_ii,R_ji),
                c = R_ii / norm,
                s = R_ji / norm; R[ji] = 0;
            if( s !== 0 ) {      R[ii] = norm;
              _giv_rot_rows(R, N-1-i, ii+1,
                                      ji+1,        c,s);
              _giv_rot_rows(Q, M, Q_off + M* i,
                                  Q_off + M*(i+1), c,s);
              // Givens rotate columns of inv(A)
              for( let h=0; h <= i; h++ ) { // <- TODO try an tighten loop bounds
                const AB_hi = AB[N*h+ i   ],
                      AB_hj = AB[N*h+(i+1)];
                AB[N*h+ i   ] =  c*AB_hi + s*AB_hj; // <- we know inv(A) is going to be triangular again so remove cancellation errors
                AB[N*h+(i+1)] = -s*AB_hi + c*AB_hj;
              }
              // rotate last row of inv(A)
              AB[N*k+(i+1)] = -s*AB[N*k+i] + c*AB[N*k+(i+1)];
              AB[N*k+ i   ] = 0; // <- aside from rounding, inv(A) is going to be triangular again
            }
          }
        }

        // DOWNDATE A\B
        for( let i=0;   i <= k; i++ )
        for( let j=k; ++j <  N;     )
          AB[N*i+j] -= AB[N*i+k] * R[R_off + N*k+j];

        // DOWNDATE inv(A)
        { const                                 R_kk = - R[R_off + N*k+k];
          for( let i=k; i-- > 0; ) AB[N*i+k] *= R_kk;
        }
        AB.fill(0, N*k+k, N*(k+1));

        // SWAP COLUMN k AND q
        // swap columns of R
        for( let i=0; i < M; i++ ) {
          const R_ik = R[R_off + N*i+k];
                       R[R_off + N*i+k] = R[R_off + N*i+q];
                                          R[R_off + N*i+q] = R_ik;
        }
        // swap columns of A\B
        for( let i=0; i < k; i++ ) {
          const AB_ik = AB[N*i+k];
                        AB[N*i+k] = AB[N*i+q];
                                    AB[N*i+q] = AB_ik;
        }
        // swap P
        const P_k = P[P_off + k];
                    P[P_off + k] = P[P_off + q];
                                   P[P_off + q] = P_k;
      }
    }
    _transpose_inplace(M, Q,Q_off);
  }

  return [
    new NDArray(Q_shape, Q),
    new NDArray(R_shape, R),
    new NDArray(P_shape, P)
  ];
}
