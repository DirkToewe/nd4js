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

import {SingularMatrixSolveError} from './singular_matrix_solve_error'
import {asarray, NDArray} from '../nd_array'
import {eps, ARRAY_TYPES} from '../dt'
import {unpermute_rows} from './permute'
import {FrobeniusNorm} from './norm'
import {_giv_rot_rows} from './_giv_rot'
import {_transpose_inplace} from './transpose_inplace'
import {_triu_solve} from './tri'


export function _norm_update(norm, Q, i,j)
{
  // https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Math/hypot#Polyfill
  i |= 0;
  j <<= 1;

  for( ; j < norm.length; j = j+2 | 0 ) {
    let s = Math.abs(Q[i+(j>>>1)]);
    if( s !== 0 ) {
      if(         norm[j] < s ) {
        const r = norm[j] / s; norm[j+1] *= r*r;
                  norm[j] = s;
      }       s/= norm[j]
                  norm[j+1] += s*s;
    }
  }
}


export function _norm(norm, i)
{
  // https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Math/hypot#Polyfill
  i <<= 1;
  const max = norm[i];
  return isFinite(max) ? Math.sqrt(norm[i+1])*max : max;
}


export function _rrqr_rank(M,N, R,R_off, tmp)
{
  // SEE: Gene H. Golub, Charles F. Van Golub
  //      "Matrix Computations", 4th edition
  //      Page 276f, Chap. 5.4.2 (QR with Column Pivoting) & 5.4.3 (Numerical Rank and AΠ=QR)
  const L = Math.min(M,N)
  if( ! (tmp.length >= L) ) throw new Error('Assertion failed.')

  const norm = new FrobeniusNorm();

  for( let i=L; i-- > 0; )
  {
    for( let j=N; j-- > i; )
      norm.include( R[R_off + N*i+j] );

    tmp[i] = norm.result;

    if( ! isFinite(tmp[i]) )
      throw new Error('Infinity or NaN encountered during rank estimation.')
  }

  const dtype = tmp instanceof Float32Array ? 'float32' : 'float64',
            T = eps(dtype)*2 * Math.max(M,N) * tmp[0] // <- threshold

  let    r = L
  while( r > 0 && tmp[r-1] <= T ) // <- TODO use binary search here for slightly improved performance
       --r
  return r
}


export function rrqr_decomp_full(A)
{
  A = asarray(A)
  if( A.ndim < 2 ) throw new Error('A must be at least 2D.')
  const
    DTypeArray = ARRAY_TYPES[A.dtype==='float32' ? 'float32' : 'float64'], // <- ensure at least double precision
    R_shape =                 A.shape,
    Q_shape = Int32Array.from(R_shape),
    P_shape =                 Q_shape.slice(0,-1),
    [M,N]   =                 R_shape.slice(  -2),
       K    = Math.min(M,N); // <- M not M-1, because the last row still needs pivotization
  Q_shape[Q_shape.length-1] = M;
  P_shape[P_shape.length-1] = N;

  const R = DTypeArray.from(A.data); // <- we could encourage GC by setting `A = undefined` after this line
  A = undefined
  const P = new Int32Array(R.length/M), // <- tracks column permutations
        Q = new DTypeArray(R.length/N*M),
     norm = new DTypeArray(N<<1); // <─ underflow-safe representation of the column norm (https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Math/hypot#Polyfill)

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
//*DEBUG*/    if( norm.some(x => x!==0) )
//*DEBUG*/      throw new Error('Assertion failed.');
    for( let j=0; j < M; j++ )
      _norm_update(norm, R, R_off + N*j, 0);

    // ELIMINATE COLUMN BY COLUMN OF R
    for( let i=0; i < K; i++ )
    { // FIND PIVOT COLUMN THAT (HOPEFULLY) ENSURES RANK REVEAL (RRQR is inherently not guaranteed to do that)
      let    p = -1,
        norm_p = -Infinity;
      for( let j=i; j < N; j++ ) {
        const  norm_j =_norm(norm,j);
        if(    norm_p < norm_j ) {
          p=j; norm_p = norm_j
        }
      }
      // swap pivot to column i
      if( p !== i ) {
        for( let j=0; j < M; j++ ) {
          const ji = R_off + N*j+i,
                jp = R_off + N*j+p, R_ji = R[ji];
                                           R[ji] = R[jp];
                                                   R[jp] = R_ji
        }
        const P_i = P[P_off+i];
                    P[P_off+i] = P[P_off+p];
                                 P[P_off+p] = P_i
      }

      // RESET COLUMN NORM (INDEX i IS SET TO ZERO FOR THE NEXT RRQR)
      norm.fill(0.0, i<<1);

      if( 0 === norm_p ) break; // <- clearly rank-deficient case

      // ELIMINATE COLUMN BELOW DIAGONAL
                                 const ii = R_off + N*i+i;
      for( let j=i; ++j < M; ) { const ji = R_off + N*j+i;
        const     R_ji = R[ji];
        if( 0 !== R_ji )
        {   const R_ii = R[ii],
                         norm = Math.hypot(R_ii,R_ji),
              c = R_ii / norm,
              s = R_ji / norm; R[ji] = 0;
          if( s !== 0 ) {      R[ii] = norm;
            _giv_rot_rows(R, N-1-i, ii+1,
                                    ji+1,        c,s);
            _giv_rot_rows(Q,   1+j, Q_off + M*i,
                                    Q_off + M*j, c,s);
          }
        }
        _norm_update(norm, R, R_off + N*j, i+1);
      }
      R[ii] = (R[ii] < 0 || Object.is(-0,R[ii]) ? -1 : +1) * norm_p;
    }
    _transpose_inplace(M, Q,Q_off);
  }

  return [
    new NDArray(Q_shape, Q),
    new NDArray(R_shape, R),
    new NDArray(P_shape, P)
  ];
}


export function _rrqr_decomp_inplace( M,N,L, A,A_off, Y,Y_off, P,P_off, norm )
{
  // using this this method could be used to implement rrqr_decomp_full
  // BUT it would be inefficient because, due to the special structure of Q,
  // Givens rotations of Q can be made more efficient in rrqr_decomp_full

  if( 0 !== M%1 ) throw new Error('Assertion failed.');
  if( 0 !== N%1 ) throw new Error('Assertion failed.');
  if( 0 !== L%1 ) throw new Error('Assertion failed.');
  if( 0 !== A_off%1 ) throw new Error('Assertion failed.');
  if( 0 !== Y_off%1 ) throw new Error('Assertion failed.');
  if( 0 !== P_off%1 ) throw new Error('Assertion failed.');

  if( ! (0 < M) ) throw new Error('Assertion failed.');
  if( ! (0 < N) ) throw new Error('Assertion failed.');
  if( ! (0 < L) ) throw new Error('Assertion failed.');
  if( ! (0 <= A_off) ) throw new Error('Assertion failed.');
  if( ! (0 <= Y_off) ) throw new Error('Assertion failed.');
  if( ! (0 <= P_off) ) throw new Error('Assertion failed.');

  if( ! (M*N <= A.length-A_off) ) throw new Error('Assertion failed.');
  if( ! (M*L <= Y.length-Y_off) ) throw new Error('Assertion failed.');
  if( ! (  N <= P.length-P_off) ) throw new Error('Assertion failed.');

  M |= 0;
  N |= 0;
  L |= 0;
  A_off |= 0;
  Y_off |= 0;
  P_off |= 0;

  if( norm.length !== N<<1 ) throw new Error('Assertion failed.');

  // INIT COLUMN NORM
  norm.fill(0.0, 0,N<<1)
  for( let j=0; j < M; j++ )
    _norm_update(norm, A, A_off + N*j, 0);

  const K = Math.min(M,N); // <- M not M-1, because the last row still needs pivotization

  // ELIMINATE COLUMN BY COLUMN OF R
  for( let i=0; i < K; i++ )
  { // FIND PIVOT COLUMN THAT (HOPEFULLY) ENSURES RANK REVEAL (RRQR is inherently not guaranteed to do that)
    let    p = -1,
      norm_p = -Infinity;
    for( let j=i; j < N; j++ ) {
      const  norm_j =_norm(norm,j);
      if(    norm_p < norm_j ) {
        p=j; norm_p = norm_j;
      }
    }
    // swap pivot to column i
    if( p !== i ) {
      for( let j=0; j < M; j++ ) {
        const ji = A_off + N*j+i,
              jp = A_off + N*j+p, A_ji = A[ji];
                                         A[ji] = A[jp];
                                                 A[jp] = A_ji;
      }
      const P_i = P[P_off+i];
                  P[P_off+i] = P[P_off+p];
                               P[P_off+p] = P_i;
    }

    // RESET COLUMN NORM (INDEX i IS SET TO ZERO FOR THE NEXT RRQR)
    norm.fill(0.0, i<<1);

    if( 0 === norm_p ) break; // <- clearly rank-deficient case

    // ELIMINATE COLUMN BELOW DIAGONAL
                               const ii = A_off + N*i+i;
    for( let j=i; ++j < M; ) { const ji = A_off + N*j+i;
      const     A_ji = A[ji];
      if( 0 !== A_ji )
      {   const A_ii = A[ii],
                       norm = Math.hypot(A_ii,A_ji),
            c = A_ii / norm,
            s = A_ji / norm; A[ji] = 0;
        if( s !== 0 ) {      A[ii] = norm;
          _giv_rot_rows(A, N-1-i, ii+1,
                                  ji+1,    c,s);
          _giv_rot_rows(Y, L, Y_off + L*i,
                              Y_off + L*j, c,s);
        }
      }
      _norm_update(norm, A, A_off + N*j, i+1);
    }
    A[ii] = (A[ii] < 0 || Object.is(-0,A[ii]) ? -1 : +1) * norm_p;
  }
}


export function rrqr_decomp(A)
{
  A = asarray(A)
  if( A.ndim < 2 ) throw new Error('A must be at least 2D.')
  const
    DTypeArray = ARRAY_TYPES[A.dtype==='float32' ? 'float32' : 'float64'], // <- ensure at least double precision
    Q_shape =                 A.shape,
    R_shape = Int32Array.from(Q_shape),
    P_shape =                 Q_shape.slice(0,-1),
    [N,M]   =                 Q_shape.slice(  -2);
  R_shape[R_shape.length-2] = M;
  P_shape[P_shape.length-1] = M;

  if( N <= M ) return rrqr_decomp_full(A)

  const Q = DTypeArray.from(A.data); // <- we could encourage GC by setting `A = undefined` after this line
  A = undefined
  const P = new Int32Array(Q.length/N), // <- tracks column permutations
        R = new DTypeArray(Q.length/N*M), // <- cache cos() and sin() values to apply M column rotations to Q at once
     norm = new DTypeArray(M<<1); // <─ underflow-safe representation of the column norm (https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Math/hypot#Polyfill)

  for(
    let Q_off=0,
        R_off=0,
        P_off=0; Q_off < Q.length; Q_off += N*M,
                                   R_off += M*M,
                                   P_off +=   M
  )
  {
    // INIT P
    for( let i=0; i < M; i++ ) P[P_off + i] = i;

    // INIT COLUMN NORM
//*DEBUG*/    if( norm.some(x => x!==0) )
//*DEBUG*/      throw new Error('Assertion failed.');
    for( let j=0; j < N; j++ )
      _norm_update(norm, Q, Q_off + M*j, 0);

    // ELIMINATE COLUMN BY COLUMN OF R (WHICH IS CURRENTLY STORED IN Q)
    for( let i=0; i < M; i++ )
    { // FIND PIVOT COLUMN THAT (HOPEFULLY) ENSURES RANK REVEAL (RRQR is inherently not guaranteed to do that)
      let    p = -1,
        norm_p = -Infinity;
      for( let j=i; j < M; j++ ) {
        const  norm_j =_norm(norm,j);
        if(    norm_p < norm_j ) {
          p=j; norm_p = norm_j
        }
      }
      // swap pivot to column i
      if( p !== i ) {
        for( let j=0; j < N; j++ ) {
          const ji = Q_off + M*j+i,
                jp = Q_off + M*j+p, A_ji = Q[ji];
                                           Q[ji] = Q[jp];
                                                   Q[jp] = A_ji
        }
        const P_i = P[P_off+i];
                    P[P_off+i] = P[P_off+p];
                                 P[P_off+p] = P_i
      }

      // RESET COLUMN NORM (INDEX i IS SET TO ZERO FOR THE NEXT RRQR)
      norm.fill(0.0, i<<1);

      if( 0 === norm_p ) break; // <- clearly rank-deficient case

      // ELIMINATE COLUMN BELOW DIAGONAL
                                 const ii = Q_off + M*i+i;
      for( let j=i; ++j < N; ) { const ji = Q_off + M*j+i;
        const     A_ji = Q[ji];
        if( 0 !== A_ji )
        {   const A_ii = Q[ii];
          let            norm = Math.hypot(A_ii,A_ji),
              c = A_ii / norm,
              s = A_ji / norm;
          if( s !== 0 ) {
            if( c < 0 ) {
                c *= -1;
                s *= -1;
             norm *= -1;
            }
            // rotate i and j
            _giv_rot_rows(Q, M-1-i, ii+1,
                                    ji+1, c,s);
            Q[ii] = norm;
          } Q[ji] = s;
        }
        _norm_update(norm, Q, Q_off + M*j, i+1);
      }
      Q[ii] = (Q[ii] < 0 || Object.is(Q[ii],-0) ? -1 : +1) * norm_p;
    }

    // MOVE R FROM Q -> R
    for( let i=0; i < M; i++ ) {
      const R_ii = R[R_off + M*i+i],
        s = R_ii < 0 || Object.is(R_ii,-0) ? -1 : +1;
      for( let j=i; j < M; j++ ) {
        R[R_off + M*i+j] = Q[Q_off + M*i+j] * s;
                           Q[Q_off + M*i+j] = s*(i===j);
      }
    }

    // COMPUTE Q
    for( let i=M; i-- > 0; )
    for( let j=N; --j > i; ) {
      const s = Q[Q_off + M*j+i]; if(0 === s) continue;
                Q[Q_off + M*j+i]  =  0; 
      const c = Math.sqrt( (1-s)*(1+s) ); 
      _giv_rot_rows(Q, M-i, Q_off + M*j+i,
                            Q_off + M*i+i, c,s);
    }
  }

  return [
    new NDArray(Q_shape, Q),
    new NDArray(R_shape, R),
    new NDArray(P_shape, P)
  ];
}


export function rrqr_rank(R)
{
  R = asarray(R)
  const [M,N] = R.shape.slice(  -2),
      r_shape = R.shape.slice(0,-2)
  R = R.data
  const r = new ARRAY_TYPES['int32'](R.length/M/N),
      tmp = new ARRAY_TYPES[R.dtype==='float32' ? 'float32' : 'float64']( Math.min(M,N) )

  for( let r_off=0,
           R_off=0; r_off < r.length;
           r_off ++,
           R_off += M*N )
    r[r_off] = _rrqr_rank(M,N, R,R_off, tmp)

  return new NDArray(r_shape, r)
}


export function rrqr_solve(Q,R,P, y)
{
  Q = asarray(Q)
  R = asarray(R)
  const N = Q.shape[Q.ndim-2]
  if( N !== R.shape[R.ndim-1] )
    throw new Error('rrqr_solve(Q,R,P, y): Q @ R not square.')

  const x = rrqr_lstsq(Q,R,P, y),
      tmp = new ARRAY_TYPES[R.dtype==='float32' ? 'float32' : 'float64'](N)

  for( let R_off = 0;
           R_off < R.data.length;
           R_off += N*N )
  {
    const rank = _rrqr_rank(N,N, R.data,R_off, tmp)
    if( rank < N ) // FIXME: what if (Qᵀy)[i >= rank] ≈ 0 then the system would still be solvable
      throw new SingularMatrixSolveError(x) 
  }

  return x
}


export function rrqr_lstsq(Q,R,P, y)
{
  if( y == undefined )
  {
    if( P != undefined )
      throw new Error('rrqr_lstsq(Q,R,P, y): Either 2 ([Q,R,P], y) or 4 arguments (Q,R,P, y) expected.')
    y = R
    ([Q,R,P] = Q)
  }

  Q = asarray(Q); if( Q.ndim < 2 ) throw new Error('rrqr_lstsq(Q,R,P, y): Q.ndim must be at least 2.')
  R = asarray(R); if( R.ndim < 2 ) throw new Error('rrqr_lstsq(Q,R,P, y): R.ndim must be at least 2.')
  P = asarray(P); if( P.ndim < 1 ) throw new Error('rrqr_lstsq(Q,R,P, y): P.ndim must be at least 1.')
  y = asarray(y); if( y.ndim < 2 ) throw new Error('rrqr_lstsq(Q,R,P, y): y.ndim must be at least 2.')
  if( P.dtype !== 'int32' ) throw new Error('rrqr_lstsq(Q,R,P, y): P.dtype must be "int32".')

  //  ________________   ______                   ___________
  // |                | |\(MxI)|                 |           |
  // |                | | \ R  |  ___________    |           |
  // |     (NxM)      | |  \   | |   (IxJ)   |   |   (NxJ)   |
  // |       Q        | |   \  | |     X     | = |     Y     |
  // |                | |    \ | |           |   |           |
  // |                | |     \| |___________|   |           |
  // |                | |  0   |                 |           |
  // |________________| |______|                 |___________|
  const
    [N,M] = Q.shape.slice(-2),
    [I]   = R.shape.slice(-1),
    [J]   = y.shape.slice(-1)
  if( N != y.shape[y.ndim-2] ) throw new Error("rrqr_lstsq(Q,R,P,y): Q and y don't match.")
  if( M != R.shape[R.ndim-2] ) throw new Error("rrqr_lstsq(Q,R,P,y): Q and R don't match.")
  if( I != P.shape[P.ndim-1] ) throw new Error("rrqr_lstsq(Q,R,P,y): R and P don't match.")

  const ndim = Math.max(Q.ndim, R.ndim, P.ndim+1, y.ndim),
       shape = Int32Array.from({length: ndim}, () => 1);
  shape[ndim-2] = I;
  shape[ndim-1] = J;

  // FIND COMMON (BROADCASTED) SHAPE
  for( let arr of [Q,R,y] )
    for( let i=ndim-2, j=arr.ndim-2; i-- > 0 && j-- > 0; )
      if( 1 === shape[i] )
        shape[i] = arr.shape[j];
      else if( shape[i] != arr.shape[j] && arr.shape[j] != 1 )
        throw new Error('rrqr_lstsq(Q,R,P,y): Q,R,P,y not broadcast-compatible.');

  for( let i=ndim-2, j=P.ndim-1; i-- > 0 && j-- > 0; )
    if( 1 === shape[i] )
      shape[i] = P.shape[j];
    else if( shape[i] != P.shape[j] && P.shape[j] != 1 )
      throw new Error('rrqr_lstsq(Q,R,P,y): Q,R,P,y not broadcast-compatible.');

  // GENERATE RESULT DATA
  const
    DTypeArray = ARRAY_TYPES[ [Q,R,y].every( a => a.dtype==='float32' ) ? 'float32' : 'float64' ],
    tmp_rank = new DTypeArray( Math.min(M,I) ),
    tmp_perm = new Int32Array(I),
    x_dat = new DTypeArray( shape.reduce((a,b) => a*b, 1) ),
    Q_dat = Q.data,
    R_dat = R.data,
    P_dat = P.data,
    y_dat = y.data;
  let
    Q_off = 0, Q_stride = 1,
    R_off = 0, R_stride = 1,
    P_off = 0, P_stride = 1,
    y_off = 0, y_stride = 1,
    x_off = 0;

  function solv(d) {
    if( d === ndim-2 ) {
      Q_stride = N*M;
      R_stride = M*I;
      P_stride =   I;
      y_stride = N*J;

      const R = _rrqr_rank(M,I, R_dat,R_off, tmp_rank)

      // Q.T @ y
      for( let k=0; k < N; k++ )
      for( let i=0; i < R; i++ )
      for( let j=0; j < J; j++ )
        x_dat[x_off+i*J+j] += Q_dat[Q_off+k*M+i] * y_dat[y_off+k*J+j]

      _triu_solve(R,I,J, R_dat,R_off, x_dat,x_off);

      // APPLY P TO X (PERMUTE ROWS)
      // https://www.geeksforgeeks.org/reorder-a-array-according-to-given-indexes/
      for( let i=I; i-- > 0; )
        tmp_perm[i] = P_dat[P_off + i]

      for( let i=I; i-- > 0; )
      {
        let k
        while( (k = tmp_perm[i]) !== i )
        {
          if( tmp_perm[k] === k )
            throw new Error("rrqr_lstsq(Q,R,P,y): Invalid indices in P.")
          tmp_perm[i] = tmp_perm[k]
          tmp_perm[k] = k

          for( let j=J; j-- > 0; )
          {
            const x_ij = x_dat[x_off + J*i+j]
                         x_dat[x_off + J*i+j] = x_dat[x_off + J*k+j]
                                                x_dat[x_off + J*k+j] = x_ij
          }
        }
      }

      Q_off += Q_stride;
      R_off += R_stride;
      P_off += P_stride;
      y_off += y_stride;
      x_off += I*J;

      return;
    }
    for( let l=shape[d]; ; l-- ) {
      solv(d+1);
      if( l == 1 ) break;
      if( ! (Q.shape[ d - ndim   + Q.ndim ] > 1) ) Q_off -= Q_stride;
      if( ! (R.shape[ d - ndim   + R.ndim ] > 1) ) R_off -= R_stride;
      if( ! (P.shape[ d - ndim+1 + P.ndim ] > 1) ) P_off -= P_stride;
      if( ! (y.shape[ d - ndim   + y.ndim ] > 1) ) y_off -= y_stride;
    }
    Q_stride *= Q.shape[ d - ndim   + Q.ndim ] || 1;
    R_stride *= R.shape[ d - ndim   + R.ndim ] || 1;
    P_stride *= P.shape[ d - ndim+1 + P.ndim ] || 1;
    y_stride *= y.shape[ d - ndim   + y.ndim ] || 1;
  }
  solv(0);

  return new NDArray(shape,x_dat);
}
