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

import {ARRAY_TYPES} from '../dt'
import {asarray, NDArray} from '../nd_array'


export function tril(m, k=0)
{
  m = asarray(m);
  if( m.ndim < 2 ) throw new Error('Input must be at least 2D.');
  return m.mapElems( m.dtype, (x,...indices) => {
    const [i,j] = indices.slice(-2);
    return i < j-k ? 0 : x
  });
}


export function triu(m, k=0)
{
  m = asarray(m);
  if( m.ndim < 2 ) throw new Error('Input must be at least 2D.');
  return m.mapElems( m.dtype, (x,...indices) => {
    const [i,j] = indices.slice(-2);
    return i > j-k ? 0 : x
  });
}


export function _tril_solve(M,N,O, L,L_off, X,X_off)
{
  if( !(M <= N) ) throw new Error('Assertion failed.');
  M |= 0
  N |= 0
  O |= 0
  L_off |= 0
  X_off |= 0

  // FORWARD SUBSTITUTION
  for( let i=0; i < M; i++ )
  {
    for( let k=0; k < i; k++ )
    for( let j=0; j < O; j++ )
      X[X_off+O*i+j] -= L[L_off+N*i+k] * X[X_off+O*k+j];

    for( let j=0; j < O; j++ )
      X[X_off+O*i+j] /= L[L_off+N*i+i];
  }
}


export function tril_solve(L,Y)
{
  L = asarray(L); if( L.ndim < 2 ) throw new Error('tril_solve(L,Y): L.ndim must be at least 2.');
  Y = asarray(Y); if( Y.ndim < 2 ) throw new Error('tril_solve(L,Y): Y.ndim must be at least 2.');

  const [M,N] = Y.shape.slice(-2);
  
  if( L.shape[L.ndim-2] !== M ) throw new Error("tril_solve(L,Y): L and Y don't match.");
  if( L.shape[L.ndim-1] !== M ) throw new Error('tril_solve(L,Y): Last two dimensions of L must be quadratic.');

  const
    ndim = Math.max(L.ndim, Y.ndim),
    L_shape = L.shape,
    Y_shape = Y.shape,
    X_shape = new Int32Array(ndim);
  X_shape.fill(1, 0,-2);
  X_shape[ndim-2] = M;
  X_shape[ndim-1] = N;

  // FIND COMMON (BROADCASTED) SHAPE
  for( let shape of [L_shape,Y_shape] )
    for( let i=ndim-2, j=shape.length-2; i-- > 0 && j-- > 0; )
      if( 1 === X_shape[i] )
        X_shape[i] = shape[j];
      else if( X_shape[i] != shape[j] && shape[j] != 1 )
        throw new Error('tril_solve(L,Y): L and Y not broadcast-compatible.');

  // GENERATE RESULT DATA
  const dtype = [L,Y].every(A => A.dtype==='float32') ? 'float32' : 'float64';
        L = L.data;
        Y = Y.data;
  const X = new ARRAY_TYPES[dtype]( X_shape.reduce((a,b) => a*b) );
    
  let
    L_off = 0, L_stride = 1,
    Y_off = 0, Y_stride = 1,
    X_off = 0;

  function solv(d) {
    if( d === ndim-2 ) {
      L_stride = M*M;
      Y_stride = M*N;

      // COPYING y
      for( let i=Y_stride; i-- > 0; )
        X[X_off+i] = Y[Y_off+i];

      _tril_solve(M,M,N, L,L_off, X,X_off);

      L_off += L_stride;
      Y_off += Y_stride;
      X_off += Y_stride;
    }
    else {
      for( let l=X_shape[d]; ; l-- ) {
        solv(d+1);
        if( l == 1 ) break;
        if( ! (L_shape[ d - ndim + L_shape.length ] > 1) ) L_off -= L_stride;
        if( ! (Y_shape[ d - ndim + Y_shape.length ] > 1) ) Y_off -= Y_stride;
      }
      L_stride *= L_shape[ d - ndim + L_shape.length ] || 1;
      Y_stride *= Y_shape[ d - ndim + Y_shape.length ] || 1;
    }
  }
  solv(0);

  return new NDArray(X_shape, X);
}


export function _triu_solve(M,N,O, U,U_off, X,X_off)
{
  if( !(M <= N) ) throw new Error('Assertion failed.');
  M |= 0
  N |= 0
  O |= 0
  U_off |= 0
  X_off |= 0

  if( !(0 <= U_off) ) throw new Error('Assertion failed.');
  if( !(0 <= X_off) ) throw new Error('Assertion failed.');

  if( !(M*N <= U.length - U_off) ) throw new Error('Assertion failed.');
  if( !(M*O <= X.length - X_off) ) throw new Error('Assertion failed.');

  // BACKWARD SUBSTITUTION
  for( let i=M; i-- > 0; )
  for( let j=O; j-- > 0; )
  {
    for( let k=M; --k > i; )
      X[X_off + O*i+j] -= U[U_off + N*i+k] * X[X_off + O*k+j]

    X[X_off + O*i+j] /= U[U_off + N*i+i];
  }
}


export function triu_solve (U,Y)
{
  U = asarray(U); if( U.ndim < 2 ) throw new Error('triu_solve(U,Y): U.ndim must be at least 2.');
  Y = asarray(Y); if( Y.ndim < 2 ) throw new Error('triu_solve(U,Y): Y.ndim must be at least 2.');

  const [M,N] = Y.shape.slice(-2);
  
  if( U.shape[U.ndim-2] !== M ) throw new Error("triu_solve(U,Y): U and Y don't match.");
  if( U.shape[U.ndim-1] !== M ) throw new Error('triu_solve(U,Y): Last two dimensions of U must be quadratic.');

  const
    ndim = Math.max(U.ndim, Y.ndim),
    U_shape = U.shape,
    Y_shape = Y.shape,
    X_shape = new Int32Array(ndim);
  X_shape.fill(1, 0,-2);
  X_shape[ndim-2] = M;
  X_shape[ndim-1] = N;

  // FIND COMMON (BROADCASTED) SHAPE
  for( let shape of [U_shape,Y_shape] )
    for( let i=ndim-2, j=shape.length-2; i-- > 0 && j-- > 0; )
      if( 1 === X_shape[i] )
        X_shape[i] = shape[j];
      else if( X_shape[i] != shape[j] && shape[j] != 1 )
        throw new Error('triu_solve(U,Y): U and Y not broadcast-compatible.');

  // GENERATE RESULT DATA
  const dtype = [U,Y].every(A => A.dtype==='float32') ? 'float32' : 'float64';
        U = U.data;
        Y = Y.data;
  const X = new ARRAY_TYPES[dtype]( X_shape.reduce((a,b) => a*b) );
    
  let
    U_off = U.length, U_stride = 1,
    Y_off = Y.length, Y_stride = 1,
    X_off = X.length;

  function solv(d) {
    if( d === ndim-2 ) {
      U_stride = M*M;
      Y_stride = M*N;

      U_off -= U_stride;
      Y_off -= Y_stride;
      X_off -= Y_stride;

      // COPYING y
      for( let i=Y_stride; i-- > 0; )
        X[X_off+i] = Y[Y_off+i];

      _triu_solve(M,M,N, U,U_off, X,X_off);
    }
    else {
      for( let l=X_shape[d]; ; l-- ) {
        solv(d+1);
        if( l == 1 ) break;
        if( ! (U_shape[ d - ndim + U_shape.length ] > 1) ) U_off += U_stride;
        if( ! (Y_shape[ d - ndim + Y_shape.length ] > 1) ) Y_off += Y_stride;
      }
      U_stride *= U_shape[ d - ndim + U_shape.length ] || 1;
      Y_stride *= Y_shape[ d - ndim + Y_shape.length ] || 1;
    }
  }
  solv(0);

  return new NDArray(X_shape, X);
}
