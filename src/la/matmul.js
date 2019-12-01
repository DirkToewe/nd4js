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

import {NDArray, asarray} from '../nd_array'
import {ARRAY_TYPES, super_dtype} from '../dt'
import math from '../math'


const [
  matmul2_RR,
  matmul2_RC,
  matmul2_CR,
  matmul2_CC,
  matmul2_ELSE
] = function*(){
  const mk_matmul2 = addProduct => new Function(
    'A_shape', 'A',
    'B_shape', 'B',
    'C_shape', 'C',
    `
      const I = A_shape[A_shape.length-2],
            K = A_shape[A_shape.length-1],
            J = B_shape[B_shape.length-1]; 
      let
        a = 0, aStride = 1,
        b = 0, bStride = 1,
        c = 0;

      const loops = new Int32Array(C_shape.length-2);
      for( let d=0; d >= 0; )
        if( d === C_shape.length-2 ) {
          aStride = I*K;
          bStride = K*J;
          for( const cEnd = c + I*J; c < cEnd; c += J, b -= bStride ) {
          for( const aEnd = a + K  ; a < aEnd; c -= J, a++ ) {
          for( const bEnd = b +   J; b < bEnd; c += 1, b++ ) {
            ${addProduct}
          }}}
          b += bStride;
          d -= 1;
        }
        else
        {
          if( loops[d]++ > 0 ) {
            if( loops[d] > C_shape[d] ) {
              aStride *= A_shape[ d - C_shape.length + A_shape.length ] || 1;
              bStride *= B_shape[ d - C_shape.length + B_shape.length ] || 1;
              loops[d--] = 0;
              continue;
            }
            if( ! (A_shape[ d - C_shape.length + A_shape.length ] > 1) ) a -= aStride;
            if( ! (B_shape[ d - C_shape.length + B_shape.length ] > 1) ) b -= bStride;
          }
          ++d;
        }
    `
  )

  yield mk_matmul2('C[c] += A[a]*B[b];')
  yield mk_matmul2(`
    C[2*c+0] += A[a] * B[2*b+0];
    C[2*c+1] += A[a] * B[2*b+1];
  `)
  yield mk_matmul2(`
    C[2*c+0] += B[b] * A[2*a+0];
    C[2*c+1] += B[b] * A[2*a+1];
  `)
  yield mk_matmul2(`
    C[2*c+0] += B[2*b+0] * A[2*a+0]  -  B[2*b+1] * A[2*a+1];
    C[2*c+1] += B[2*b+0] * A[2*a+1]  +  B[2*b+1] * A[2*a+0];
  `)
  yield mk_matmul2('C[c] = this.add(C[c], this.mul(A[a],B[b]));').bind(math)
}()


export function matmul2(a,b)
{
  a = asarray(a)
  b = asarray(b)
  if( a.ndim < 2 ) throw new Error('A must be at least 2D.');
  if( b.ndim < 2 ) throw new Error('B must be at least 2D.');

  const
    [I,K] = a.shape.slice(-2),
     J    = b.shape[b.ndim-1];
  if( b.shape[b.ndim-2] != K )
    throw new Error('The last dimension of A and the 2nd to last dimension of B do not match.');

  const
    ndim = Math.max(a.ndim, b.ndim),
    shape = Int32Array.from({ length: ndim }, () => 1 );
  shape[ndim-2] = I;
  shape[ndim-1] = J;

  // FIND COMMON (BROADCASTED) SHAPE
  for( let arr of [a,b] )
    for( let i=ndim-2, j=arr.ndim-2; i-- > 0 && j-- > 0; )
      if( 1 === shape[i] )
        shape[i] = arr.shape[j];
      else if( shape[i] != arr.shape[j] && arr.shape[j] != 1 )
        throw new Error('Shapes are not broadcast-compatible.');

  // GENERATE RESULT DATA
  const DTypeArray = ARRAY_TYPES[super_dtype(a.dtype,b.dtype)];

  // INTEGER TYPE ID (0: scalar, 1: complex, 2: object)
  const c = new DTypeArray( shape.reduce((m,n) => m*n, 1) );
  if( c instanceof Array )
    c.fill(0)

  let pairing = 0
  for( const ab of [a,b] ) {
    pairing *= 3
    switch(ab.dtype) {
      default          : pairing += 1
      case 'complex128': pairing += 1
      case      'int32':
      case    'float32':
      case    'float64':
    }
  }
  switch(pairing)
  {
    case 0 : matmul2_RR  (a.shape, a.data       , b.shape, b.data       , shape, c       ); break
    case 1 : matmul2_RC  (a.shape, a.data       , b.shape, b.data._array, shape, c._array); break
    case 3 : matmul2_CR  (a.shape, a.data._array, b.shape, b.data       , shape, c._array); break
    case 4 : matmul2_CC  (a.shape, a.data._array, b.shape, b.data._array, shape, c._array); break
    default: matmul2_ELSE(a.shape, a.data       , b.shape, b.data       , shape, c       ); break
  }

  return new NDArray(shape,c);
}


export function matmul(...matrices)
{
  matrices = matrices.map(a => asarray(a))
  if( matrices.length == 1 ) return matrices[0];
  if( matrices.length == 2 ) return matmul2(...matrices);

  /** Returns the number of floating point operations necessary to
   *  matrix multiply two arrays of the given shapes.
   */
  function nOps( shapeA, shapeB ) {
    const [I,K] = shapeA.slice(-2),
           J    = shapeB[shapeB.length-1];
    if( shapeB[shapeB.length-2] != K )
      throw new Error('Shape mismatch.');

    const
      ndim = Math.max(shapeA.length, shapeB.length),
      shape = Int32Array.from({ length: ndim }, () => 1 );
    shape[ndim-2] = I;
    shape[ndim-1] = J;

    // FIND COMMON (BROADCASTED) SHAPE
    for( let shp of [shapeA,shapeB] )
      for( let i=ndim-2, j=shp.length-2; i-- > 0 && j-- > 0; )
        if( 1 === shape[i] )
          shape[i] = shp[j];
        else if( shape[i] != shp[j] && shp[j] != 1 )
          throw new Error('Shapes are not broadcast-compatible.');

    return [ shape.reduce( (a,b) => a*b, 1 )*K, shape ]
  };
  // https://en.wikipedia.org/wiki/Matrix_chain_multiplication
  // https://www.geeksforgeeks.org/matrix-chain-multiplication-dp-8/
  
  // opMap[i][j] (j >= i) caches the optimal number of FLOPs required to multiply matrices[i] up to (including) matrices[j]
  const opMap = Array.from({ length: matrices.length }, () => [] );

  // initialized opMap
  for( let i=0; i < matrices.length; i++ ) opMap[i][i] = [ 0, matrices[i].shape ]
  // compute remaining opMap
  for( let len=2; len <= matrices.length;   len++ )
  for( let  i =0;  i  <= matrices.length-len; i++ )
  {
    let
      minFlops = Infinity,
      minShape;
    for( let j=1; j < len; j++ )
    {
      const [lFlops,lShape] = opMap[i  ][i+ j -1]
      const [rFlops,rShape] = opMap[i+j][i+len-1]
      let [flops,shape] = nOps(lShape,rShape)
      flops += lFlops + rFlops;
      if( flops < minFlops ) {
        minFlops = flops;
        minShape = shape; // <- the shape should always be the same so this is not strictly necessary
      }
    }
    if( minShape === undefined ) throw new Error('Integer overflow (too many FLOPs).');
    opMap[i][i+len-1] = [ minFlops, minShape ]
  }

  // compute the result using the minimal number of FLOPs using opMap
  function product(from, to)
  {
    if( from == to ) return matrices[from];
    let
      minFlops = Infinity,
      minIdx;
    for( let i=from; i < to; i++ )
    {
      const [lFlops,lShape] = opMap[from][i]
      const [rFlops,rShape] = opMap[i+1][to]
      let [flops,_] = nOps(lShape,rShape)
      flops += lFlops + rFlops;
      if( flops < minFlops ) {
        minFlops = flops;
        minIdx   = i;
      }
    }
    return matmul2(
      product(from, minIdx),
      product(minIdx+1, to)
    );
  }

  return product(0, matrices.length-1)
}
