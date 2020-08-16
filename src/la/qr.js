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
import {_giv_rot_qr,
        _giv_rot_rows} from './_giv_rot'
import {_transpose_inplace} from './transpose_inplace'
import {_triu_solve} from './tri'


export function qr_decomp_full(A)
{
  A = asarray(A);
  if( A.ndim < 2 ) throw new Error('A must be at least 2D.');
  const DType = A.dtype === 'float32' ? 'float32' : 'float64',
        DTypeArray = ARRAY_TYPES[DType], // <- ensure at least double precision
    B = DType === 'float32' ? 64/4 : 64/8,
    R_shape = A.shape,
    Q_shape = Int32Array.from(R_shape),
    [M,N]   = R_shape.slice(-2),
     R = DTypeArray.from(A.data);
  A = undefined
  Q_shape[Q_shape.length-1] = M;
  const Q = new DTypeArray(R.length/N*M);
//  Q.fill(0); // <- in case of an object array

  for(
    let Q_off=0,
        R_off=0;
        Q_off < Q.length;
        Q_off += M*M,
        R_off += M*N
  )
  {
    // INIT Q TO IDENTITY MATRIX
    for( let i=0; i < M; i++ ) Q[Q_off+M*i+i] = 1.0;
                                                   // The idea of the blocked loop is that the top B rows
    for( let J=0; J < N; J += B )                  // are cached and used to elimate B columns in the remaining
    for( let I=J; I < M; I += B )                  // rows. This should reduce the number of cache misses.
    for( let i=I; i < I+B && i < M         ; i++ ) // Some quick and dirty benchmarks indicate a ~15% perfomance
    for( let j=J; j < J+B && j < N && j < i; j++ ) // improvement for matrix sizes of [300,300] and above.
    {
      // USE GIVENS ROTATION TO ELIMINATE ELEMENT R_ji
      const ij = R_off+N*i+j, R_ij = R[ij]; if(0 === R_ij) continue;
      const jj = R_off+N*j+j, R_jj = R[jj],
        [c,s,norm] =_giv_rot_qr(R_jj,R_ij);
          R[ij]= 0; if(0 === s) continue;
          R[jj]= norm;
      _giv_rot_rows( R, N-1-j, jj+1,
                               ij+1,      c,s);
      _giv_rot_rows( Q,   1+i, Q_off+M*j,
                               Q_off+M*i, c,s);
    }
    _transpose_inplace(M, Q,Q_off);
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
  const R = new DTypeArray(Q.length/N*M);  // <- additional space to temp. store rows of R not contained in the result

  for(
    let R_off=0,
        Q_off=0; Q_off < Q.length; Q_off += N*M,
                                   R_off += M*M
  )
  {
    // COMPUTE R (inside of Q)
    for( let i=1; i < N; i++ ) { const I = Math.min(i,M);
    for( let j=0; j < I; j++ )
    { // USE GIVENS ROTATION TO ELIMINATE ELEMENT R_ji
      const ij = Q_off + M*i+j, R_ij = Q[ij]; if(0 === R_ij) continue;
      const jj = Q_off + M*j+j, R_jj = Q[jj];
      let [c,s,norm] =_giv_rot_qr(R_jj,R_ij);
      if( s !== 0 ) {
        if( c < 0 ) {
            c *= -1;
            s *= -1;
         norm *= -1;
        }
        _giv_rot_rows(Q, M-1-j, jj+1,
                                ij+1, c,s);
        Q[jj] = norm;
      } Q[ij] = s;
    }}

    // MOVE R FROM Q -> R AND INIT Q TO I
    for( let i=0; i < M; i++ )
    for( let j=i; j < M; j++ ) {
      R[R_off + M*i+j] = Q[Q_off + M*i+j];
                         Q[Q_off + M*i+j] = +(i===j);
    }

    // COMPUTE Q
    for( let i=N; --i > 0; ) { const I = Math.min(i,M);
    for( let j=I; j-- > 0; )
    {
      const s = Q[Q_off + M*i+j]; if(0 === s) continue;
                Q[Q_off + M*i+j]  =  0;
      const c = Math.sqrt( (1-s)*(1+s) );
      _giv_rot_rows(Q, M-j, Q_off + M*i+j,
                            Q_off + M*j+j, c,s);
    }}
  }

  return [
    new NDArray(Q_shape, Q),
    new NDArray(R_shape, R)
  ];
}


export function _qr_decomp_inplace(M,N,L, A,A_off, Y,Y_off)
{
  if( M%1 !== 0 ) throw new Error('Assertion failed.');
  if( N%1 !== 0 ) throw new Error('Assertion failed.');
  if( L%1 !== 0 ) throw new Error('Assertion failed.');
  if( A_off%1 !== 0 ) throw new Error('Assertion failed.');
  if( Y_off%1 !== 0 ) throw new Error('Assertion failed.');
  if( A.length%1 !== 0 ) throw new Error('Assertion failed.');
  if( Y.length%1 !== 0 ) throw new Error('Assertion failed.');

  if( !(0 <= M) ) throw new Error('Assertion failed.');
  if( !(0 <= N) ) throw new Error('Assertion failed.');
  if( !(0 <= L) ) throw new Error('Assertion failed.');
  if( !(0 <= A_off) ) throw new Error('Assertion failed.');
  if( !(0 <= Y_off) ) throw new Error('Assertion failed.');
  if( !(0 <= A.length) ) throw new Error('Assertion failed.');
  if( !(0 <= Y.length) ) throw new Error('Assertion failed.');

  if( !(M*N <= A.length - A_off) ) throw new Error('Assertion failed.');
  if( !(M*L <= Y.length - Y_off) ) throw new Error('Assertion failed.');

  for( let i=1; i < M         ; i++ )
  for( let j=0; j < N && j < i; j++ )
  {
    // USE GIVENS ROTATION TO ELIMINATE ELEMENT R_ji
    const ij = A_off+N*i+j, A_ij = A[ij]; if(0 === A_ij) continue;
    const jj = A_off+N*j+j, A_jj = A[jj],
      [c,s,norm] =_giv_rot_qr(A_jj,A_ij);
        A[ij]= 0; if(0 === s) continue;
        A[jj]= norm;
    _giv_rot_rows(A,N-1-j,  jj+1,
                            ij+1, c,s);
    _giv_rot_rows(Y,L, Y_off+L*j,
                       Y_off+L*i, c,s);
  }
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

      _triu_solve(L,I,J, R_dat,R_off, x_dat,x_off);

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
