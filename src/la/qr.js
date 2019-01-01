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


export function qr_decomp_full(A)
{
  A = asarray(A);
  if( A.ndim < 2 ) throw new Error('A must be at least 2D.');
  const
    DTypeArray = ARRAY_TYPES[A.dtype === 'float32' ? 'float32' : 'float64'], // <- ensure at least double precision
    R_shape = A.shape,
    Q_shape = Int32Array.from(R_shape),
    [N,M]   = R_shape.slice(-2),
     R = DTypeArray.from(A.data);
  A = undefined
  Q_shape[Q_shape.length-1] = N;
  const Q = new DTypeArray(R.length/M*N);
  Q.fill(0); // <- in case of an object array

  for(
    let Q_off=0,
        R_off=0;
    Q_off < Q.length;
    Q_off += N*N,
    R_off += N*M
  )
  {
    // INIT Q TO IDENTITY MATRIX
    for( let i=0; i < N; i++ ) Q[Q_off+N*i+i] = 1.0;

    for( let i=1; i < N; i++ ) { const I = Math.min(i,M);
    for( let j=0; j < I; j++ )
    {
      // USE GIVENS ROTATION TO ELIMINATE ELEMENT R_ji
      const R_ij = R[R_off+M*i+j]; if( R_ij == 0.0 ) continue;
      const R_jj = R[R_off+M*j+j],
                   norm = Math.hypot(R_jj,R_ij),
        c = R_jj / norm,
        s = R_ij / norm;
      R[R_off+M*i+j] = 0;
      R[R_off+M*j+j] = norm;
      // ROTATE ROW i AND j IN R
      for( let k=j; ++k < M; ) {
        const ik = R_off+M*i+k, R_ik = R[ik],
              jk = R_off+M*j+k, R_jk = R[jk];
        R[ik] = c*R_ik - s*R_jk;
        R[jk] = s*R_ik + c*R_jk;
      }
      // ROTATE COL i AND j IN Q (Q TRANSPOSED FOR CACHE LOCALITY REASONS) 
      for( let k=0; k <= i; k++ ) {
        const ik = Q_off+N*i+k, Q_ik = Q[ik],
              jk = Q_off+N*j+k, Q_jk = Q[jk];
        Q[ik] = c*Q_ik - s*Q_jk;
        Q[jk] = s*Q_ik + c*Q_jk;
      }
    }}
    // TRANSPOSE Q (Q TRANSPOSED FOR CACHE LOCALITY REASONS)
    for( let i=0;   i < N; i++ )
    for( let j=i+1; j < N; j++ ) {
      const
        ij = Q_off+N*i+j,
        ji = Q_off+N*j+i,
        Q_ij = Q[ij]; Q[ij] = Q[ji]; Q[ji] = Q_ij;
    }
  }

  return [
    new NDArray(Q_shape, Q),
    new NDArray(R_shape, R)
  ];
}


export function qr_decomp(A)
{
  A = asarray(A);
  if( A.ndim < 2 ) throw new Error('qr_decomp(A): A.ndim must be at least 2.');
  const
    DTypeArray = ARRAY_TYPES[A.dtype === 'float32' ? 'float32' : 'float64'], // <- ensure at least double precision
    Q_shape = A.shape,
    R_shape = Int32Array.from(Q_shape),
    [N,M] = Q_shape.slice(-2)
  R_shape[R_shape.length-2] = M;

  if( N <= M ) return qr_decomp_full(A);

  const Q = DTypeArray.from(A.data); // <- we could encourage GC by setting `A = undefined` after this line
  A = undefined
  const R = new DTypeArray(Q.length/N*M),
       cs = new DTypeArray( N*M - (M*(M+1) >>> 1) ),
        r = function(){
          try      { return    cs.subarray(M); }
          catch(e) { return new DTypeArray(M); }
        }();  // <- additional space to temp. store rows of R not contained in the result

  for(
    let R_off=0,
        Q_off=0; Q_off < Q.length; Q_off += N*M,
                                   R_off += M*M
  )
  {
    let csi=0;

    // COMPUTE R (inside of Q)
    for( let i=1; i < N; i++ ) { const I = Math.min(i,M);
    for( let j=0; j < I; j++ )
    { // USE GIVENS ROTATION TO ELIMINATE ELEMENT R_ji
      const R_ij = Q[Q_off + M*i+j]; if( R_ij === 0.0 ) { cs[csi++] = 0.0; continue; }
      const R_jj = Q[Q_off + M*j+j];
      let          norm = Math.hypot(R_jj,R_ij),
        c = R_jj / norm,
        s = R_ij / norm;
      if( c < 0 ) {
           c *= -1;
           s *= -1;
        norm *= -1;
      }
      cs[csi++] = s;
      Q[Q_off + M*i+j] = 0;
      Q[Q_off + M*j+j] = norm;
      // ROTATE ROW i AND j IN R
      for( let k=j; ++k < M; ) {
        const ik = Q_off + M*i+k, R_ik = Q[ik],
              jk = Q_off + M*j+k, R_jk = Q[jk];
        Q[ik] = c*R_ik - s*R_jk;
        Q[jk] = s*R_ik + c*R_jk;
      }
    }}

    if( csi != cs.length ) throw new Error('Assertion failed!');

    // MOVE R FROM Q -> R AND INIT Q TO I
    for( let i=0; i < M; i++ )
    for( let j=i; j < M; j++ ) {
      R[R_off + M*i+j] = Q[Q_off + M*i+j];
                         Q[Q_off + M*i+j] = i !== j ? 0.0 : 1.0;
    }

    // COMPUTE Q
    for( let i=N; --i > 0; ) { const I = Math.min(i,M);
    for( let j=I; j-- > 0; )
    {
      const s = cs[--csi]; if( 0.0 === s ) continue;
      const c = Math.sqrt( (1-s)*(1+s) );
      // ROTATE ROW i AND j IN Q
      for( let k=j; k < M; k++ ) {
        const ik = Q_off + M*i+k, R_ik = Q[ik],
              jk = Q_off + M*j+k, R_jk = Q[jk];
        Q[ik] = s*R_jk + c*R_ik;
        Q[jk] = c*R_jk - s*R_ik;
      }
    }}

    if( csi != 0 ) throw new Error('Assertion failed!');
  }

  return [
    new NDArray(Q_shape, Q),
    new NDArray(R_shape, R)
  ];
}


export function qr_lstsq(Q,R, y)
{
  if( undefined == y ) { y=R; [Q,R] = Q; }
  Q = asarray(Q); if( Q.ndim < 2 ) throw new Error('qr_lstsq(Q,R,y): Q.ndim must be at least 2.')
  R = asarray(R); if( R.ndim < 2 ) throw new Error('qr_lstsq(Q,R,y): R.ndim must be at least 2.')
  y = asarray(y); if( y.ndim < 2 ) throw new Error('qr_lstsq(Q,R,y): y.ndim must be at least 2.')

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
    [J]   = y.shape.slice(-1),
     L    = Math.min(M,I)
  if( N != y.shape[y.ndim-2] ) throw new Error("qr_lstsq(Q,R,y): Q and y don't match.")
  if( M != R.shape[R.ndim-2] ) throw new Error("qr_lstsq(Q,R,y): Q and R don't match.")

  if( I > N ) throw new Error("qr_lstsq(Q,R,y): Under-determined systems not supported. Use rrqr instead.")

  const ndim = Math.max(Q.ndim, R.ndim, y.ndim),
       shape = Int32Array.from({length: ndim}, () => 1);
  shape[ndim-2] = I;
  shape[ndim-1] = J;

  // FIND COMMON (BROADCASTED) SHAPE
  for( let arr of [Q,R,y] )
    for( let i=ndim-2, j=arr.ndim-2; i-- > 0 && j-- > 0; )
      if( 1 === shape[i] )
        shape[i] = arr.shape[j];
      else if( shape[i] != arr.shape[j] && arr.shape[j] != 1 )
        throw new Error('Q, R, y are not broadcast-compatible.');

  // GENERATE RESULT DATA
  const
    DTypeArray = ARRAY_TYPES[ [Q,R,y].every( a => a.dtype==='float32' ) ? 'float32' : 'float64' ],
    x_dat = new DTypeArray( shape.reduce((a,b) => a*b, 1) ),
    Q_dat = Q.data,
    R_dat = R.data,
    y_dat = y.data;
  let
    Q_off = 0, Q_stride = 1,
    R_off = 0, R_stride = 1,
    y_off = 0, y_stride = 1,
    x_off = 0;

  function solv(d) {
    if( d === ndim-2 ) {
      Q_stride = N*M;
      R_stride = M*I;
      y_stride = N*J;

      // Q.T @ y
      for( let i=0; i < L; i++ )
      for( let j=0; j < J; j++ )
      for( let k=0; k < N; k++ )
        x_dat[x_off+i*J+j] += Q_dat[Q_off+k*M+i] * y_dat[y_off+k*J+j]

      // BACKWARD SUBSTITUTION
      for( let i=L; i-- > 0; )
      for( let j=J; j-- > 0; ) {
        for( let k=L; --k > i; )
          x_dat[x_off+i*J+j] -= R_dat[R_off+I*i+k] * x_dat[x_off+k*J+j]
        x_dat[x_off+i*J+j] /= R_dat[R_off+I*i+i]
      }

      Q_off += Q_stride;
      R_off += R_stride;
      y_off += y_stride;
      x_off += I*J;

      return;
    }
    for( let l=shape[d]; ; l-- ) {
      solv(d+1);
      if( l == 1 ) break;
      if( ! (Q.shape[ d - ndim + Q.ndim ] > 1) ) Q_off -= Q_stride;
      if( ! (R.shape[ d - ndim + R.ndim ] > 1) ) R_off -= R_stride;
      if( ! (y.shape[ d - ndim + y.ndim ] > 1) ) y_off -= y_stride;
    }
    Q_stride *= Q.shape[ d - ndim + Q.ndim ] || 1;
    R_stride *= R.shape[ d - ndim + R.ndim ] || 1;
    y_stride *= y.shape[ d - ndim + y.ndim ] || 1;
  }
  solv(0);

  return new NDArray(shape,x_dat);
}
