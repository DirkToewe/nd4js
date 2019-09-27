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
import {transpose_inplace} from './transpose_inplace'
import {_svd_jac_angles} from './_svd_jac_utils'
import {_giv_rot_rows} from './_giv_rot'
import {FrobeniusNorm} from './norm'
import {root1d_bisect} from '../opt/root1d_bisect'


// TODO:
//   * [_svd_dc_neves]      Use matrix multiplication to update U and V
//   * [_svd_dc_neves] Optimize matrix multiplication to update U and V (exploit sparsity)

/* Computes the SVD of the 1x2 matrix B
 */
export function _svd_dc_1x2( N, U,U_off, B,B_off, V,V_off )
{
  N     |= 0;
  U_off |= 0;
  B_off |= 0;
  V_off |= 0;

  let c = B[B_off  ],
      s = B[B_off+1];
  const  norm = Math.hypot(c,s);
  if(0!==norm) {
    c /= norm;
    s /= norm;

    B[B_off  ] = norm;
    B[B_off+1] = NaN;

    V[V_off      ] =  c;
    V[V_off+    1] =  s;
    V[V_off+N*1  ] = -s;
    V[V_off+N*1+1] =  c;
  }
  else {
    V[V_off      ] = 1;
    V[V_off+N*1+1] = 1;
  }
  U[U_off] = 1;
}


/* Computes the SVD of the 2x3 matrix B
 *     ┌            ┐
 *     │ b1  b2   0 │
 * B = │            │
 *     │ 0   b3  b4 │
 *     └            ┘
 */
export function _svd_dc_2x3( N, U,U_off, B,B_off, V,V_off )
{
  U_off  |= 0;
  B_off  |= 0;
  V_off  |= 0;
  N      |= 0; const
  M = N-1|  0;

  let b1 = B[B_off  ],
      b2 = B[B_off+1],
      b3 = B[B_off+2],
      b4 = B[B_off+3];
  // STEP 1:
  //   RQ-decompose down to 2x2 problem
  if( 0 !== b4 )
  { // eliminate B[1,2] (stored in b4)
    const     norm = Math.hypot(b3,b4),
      ca = b3/norm,
      sa = b4/norm;
    V[V_off + N*1+1] =  ca;
    V[V_off + N*1+2] =  sa;
    b3 = norm;
    b4 = -sa*b2;
    b2 =  ca*b2;
    if( 0 !== b4 )
    { // eliminate B[0,2] (stored in b4)
      const     norm = Math.hypot(b1,b4),
        cb = b1/norm,
        sb = b4/norm;
           b1 = norm;
      V[V_off      ] =  cb;
      V[V_off+    1] =  sb*-sa;
      V[V_off+    2] =  sb* ca;
      V[V_off+N*2  ] = -sb;
      V[V_off+N*2+1] =  cb*-sa;
      V[V_off+N*2+2] =  cb* ca;
    }
    else {
      V[V_off      ] =   1;
      V[V_off+N*2+1] = -sa;
      V[V_off+N*2+2] =  ca;
    }
  }
  else {
    V[V_off      ] = 1;
    V[V_off+N*1+1] = 1;
    V[V_off+N*2+2] = 1;
  }

  // STEP 2;
  //   USE SINGLE JACOBI ROTATION TO COMPUTE SVD
  const [ca,sa, cb,sb] = _svd_jac_angles(b1,b2,
                                         0 ,b3);
  // WRITE SINGULAR VALUES TO B  
  B[B_off+1] = NaN; // <- TODO remove after testing
  B[B_off+3] = NaN; // <- TODO remove after testing

  B[B_off] = ca*b1*cb - ( ca*b2 + sa*b3)*sb;
                   b3 = (-sa*b2 + ca*b3)*cb - sa*b1*sb;
  const        s = b3 < 0 ? -1 : +1;
  B[B_off+2] = s * b3;

  // WRITE U
  U[U_off      ] =  ca;
  U[U_off+    1] = -sa*s;
  U[U_off+M*1  ] =  sa;
  U[U_off+M*1+1] =  ca*s;

  // ROTATE V
  V[V_off+N*1] = sb*V[V_off];
  V[V_off    ]*= cb;
  const       V_01 = V[V_off+1]*cb - V[V_off+N*1+1]*sb;
    V[V_off+N*1+1] = V[V_off+1]*sb + V[V_off+N*1+1]*cb;
    V[V_off+    1] = V_01;
  const       V_02 = V[V_off+2]*cb - V[V_off+N*1+2]*sb;
    V[V_off+N*1+2] = V[V_off+2]*sb + V[V_off+N*1+2]*cb;
    V[V_off+    2] = V_02;
}


/* Computes the singular values of a matrix shaped like a backwards seven.
 *     ┌                        ┐
 *     │ d1                     │
 *     │                        │
 *     │     d2                 │
 *     │       .                │
 * B = │         .              │
 *     │           .            │
 *     │            d[n-1]      │
 *     │                        │
 *     │ z1  . . .  z[n-1] z[n] │
 *     └                        ┘
 * See:
 *   - "INTRODUCTION OF DOUBLE DIVIDE AND CONQUER AND THE RECENT PROGRESS", TARO KONDA & YOSHIMASA NAKAMURA
 *
 * Keep in mind that B is stored in a sparse memory layout.
 * 
 * B = Array(d1, z1, d2, z2, ..., d[n]=0, z[n])
 * 
 * where d[i] >= d[i+1]
 */
export function _svd_dc_neves( N, n, U,U_off, B,B_off, V,V_off, outerOrder )
{
  U_off |= 0;
  B_off |= 0;
  V_off |= 0;
  N     |= 0; const M = N-1 | 0;
  n     |= 0; const m = n-1 | 0;

  if( n <= 2 )
    throw new Error('Assertion failed.');

  //     DIAGONAL ELEMENTS: d[i] = B[B_off + 2*i  ]
  // OFF-DIAGONAL ELEMENTS: z[i] = B[B_off + 2*i+1]

  // TODO remove checks after testing
  // ensure d[-1] === 1
  if( B[B_off + 2*m-2] !== 0 )
    throw new Error('Assertion failed.');
  // ensure d[i] >= d[i+1]
  for( let i=1; i < m; i++ )
    if( B[B_off + 2*(i-1)] < B[B_off + 2*i] )
      throw new Error('Assertion failed.');

  const TOL = eps('float64'), // <- FIXME make general purpose (independent of dtype)
       NORM = new FrobeniusNorm();

  for( let i=0; i < m; i++ )
    NORM.include(B[B_off + 2*i+1]);
  const zNorm = NORM.result;


  // STEP 1: DEFLATION
  //   TEMPORARILY MOVE VALUES WITH z[j] ≈ 0 TO THE LEFT SUCH THAT
  //
  //       ┌                                   ┐
  //       │ d1                                │
  //       │    .                              │
  //       │      .                            │
  //       │        .                          │
  //       │        d[i]                       │
  //       │                                   │
  //   B = │            d[i+1]                 │
  //       │                 .                 │
  //       │                   .               │
  //       │                     .             │
  //       │                       d[n-1]      │
  //       │                                   │
  //       │ 0 . . . 0  z[i+1] ... z[n-1] z[n] │
  //       └                                   ┘
  //
  //   This way only the lower right quadrant needs to be solved
  //   iteratively (using the secular equations). The actual
  //   implementation does not more the deflated values to μ
  //   and not to the left.
  const      μ = new Float64Array(m),
    innerOrder = new   Int32Array(m);

  const n0 = function(){
    let n0 = 0;
    for( let j=m-1,
             i=m-1; i-- > 0; )
    { const di = B[B_off + 2*i  ],
            zi = B[B_off + 2*i+1],
            oi = outerOrder[i];
      if( Math.abs(zi) <= di*TOL ) {
                 μ[n0] = di;
        innerOrder[n0] = oi; // <- used as temp. for outerOrder
                 ++n0;
      }
      else {      --j;
        B[B_off + 2*j  ] = di;
        B[B_off + 2*j+1] = B[B_off + 2*i+1];
           innerOrder[j] = oi; // <- used as temp. for outerOrder
      }
    }
    return n0;
  }();

  // innerOrder just used as temp. memory, move to outerOrder
  for( let i=0; i < n0; i++ )
    outerOrder[i] = innerOrder[i];


  // STEP 2:
  //   HANDLE DUPLICATE VALUES ON THE DIAGONAL
  //
  // Lets say d[i] == d[i+1] and z[i] != 0 with 0 < i < n-1
  //     ┌                                    ┐
  //     │ d1                                 │
  //     │   .                                │
  //     │     .                              │
  //     │       .                            │
  //     │        d[i]                        │
  //     │                                    │
  // B = │             d[i+1]                 │
  //     │                  .                 │
  //     │                    .               │
  //     │                      .             │
  //     │                        d[n-1]      │
  //     │                                    │
  //     │ z1 ... z[i] z[i+1] ... z[n-1] z[n] │
  //     └                                    ┘
  //
  // We can now use a left and right Givens rotation
  // to eliminate z[i+1].
  // n := √(z²ᵢ + z²ᵢ₊₁)
  // c := zᵢ   / n
  // s := zᵢ₊₁ / n
  //
  //     ┌                            ┐
  //     │ 1                          │
  //     │   .                        │
  //     │     .                      │
  //     │       .                    │
  //     │         1                  │
  //     │                            │
  //     │            c -s            │
  // G = │                            │
  //     │            s  c            │
  //     │                            │
  //     │                  1         │
  //     │                    .       │
  //     │                      .     │
  //     │                        .   │
  //     │                          1 │
  //     └                            ┘
  //
  //                           ┌                                   ┐
  //                           │ d1                                │
  //                           │   .                               │
  //                           │     .                             │
  //                           │       .                           │
  //                           │        d[i]                       │
  //                           │                                   │
  // B = GᵀG⋅B⋅GᵀG = Gᵀ⋅B'⋅G = │            d[i+1]                 │
  //                           │                 .                 │
  //                           │                   .               │
  //                           │                     .             │
  //                           │                       d[n-1]      │
  //                           │                                   │
  //                           │ z1 ...  0    n  . . . z[n-1] z[n] │
  //                           └                                   ┘
  //
  // Thus z'[i] is now 0 and can be moved to the left side (similar to Step 1)
  const W = new Float64Array(m*m), // <- temp. storage to update U and V
     rotJ = new Int32Array(m);

  const n1 = function(){
    let n1 = n0;
    for( let j=m-1,
             i=m-1; i-- > n0; )
    { const di = B[B_off + 2*i  ],
            dj = B[B_off + 2*j  ],  d = (di+dj) / 2,
            zi = B[B_off + 2*i+1],
            zj = B[B_off + 2*j+1], oi = innerOrder[i];
      if( d === di ||
          d === dj ) {
        const      z = Math.hypot(zi,zj),
          c = zj / z,
          s = zi / z;
        W[2*n1  ] = c;
        W[2*n1+1] = s;
        B[B_off + 2*j+1] = z;
        B[B_off + 2*j  ] = μ[n1] = di;
                  outerOrder[n1] = oi; // <- used as temp. for outerOrder
                        rotJ[n1] =  j;
                           ++n1;
      }
      else {      --j;
        B[B_off + 2*j  ] = di;
        B[B_off + 2*j+1] = B[B_off + 2*i+1];
           outerOrder[j] = oi;
      }
    }
    return n1;
  }();
  for( let i = 2*n0;
           i < 2*n1; i++ ) {
    B[B_off + i] = W[i];
                   W[i] = 0;
  }


  // STEP 3:
  //   SOLVE THE SECULAR EQUATIONS
  //
  //      !            zₗ²
  //   0  =  1 + Σₗ ────────      k = 0, ..., m-1
  //                dₗ² - σₖ²
  //
  for( let i=n1; i < m; i++ )
  {
    let sLo =          B[B_off + 2*i  ],
        sHi = n1 < i ? B[B_off + 2*i-2] : (sLo + zNorm);

    if( sLo > sHi ) throw new Error('Assertion failed.');

    const shift = function(){
      if( i > n1 )
      { const mid = (sLo + sHi) / 2;
        let   sum = 1;
        for( let i=n1; i < m; i++ ) {
          const  di = B[B_off + 2*i  ],
                 zi = B[B_off + 2*i+1];
          sum += zi / (di-mid)
              * (zi / (di+mid));
        }
        if( ! isFinite(sum) ) throw new Error('Assertion failed.');
        if( sum > 0 ) { const s=sHi; sLo = sLo-sHi; sHi = -Number.MIN_VALUE; return s; }
      }                 const s=sLo; sHi = sHi-sLo; sLo = +Number.MIN_VALUE; return s;
    }();

    for(;;)
    { // bisect
      const s  = (sLo + sHi) / 2;
      if(   s === sLo
         || s === sHi ) {
        μ[i] = s; // FIXME at this point sLo and sHi still have to be compared
        break;
      }
      // evalue the secular equation
      let sum = 1;
      for( let i=n1; i < m; i++ ) {
        const  di = B[B_off + 2*i  ],
               zi = B[B_off + 2*i+1];
        sum += zi / (di-shift-s)
            * (zi / (di+shift+s));
      }
      if( ! isFinite(sum) ) throw new Error('Assertion failed.');
      // adjust bisection bounds
      if( sum <= 0 ) sLo = s;
      if( sum >= 0 ) sHi = s;
    }
  }
  if( Math.abs(B[B_off + 2*m-1]) === 0 ) {
    B[B_off + 2*m-2] = 0;
    B[B_off + 2*m-1] = 0; μ[m-1] = 0;
  }


  // STEP 4:
  //   RECOMPUTE z TO ORTHOGONALIZE U & V
  //   (as originally suggested by Gu and Eisenstat)
  {
    const σn = μ[m-1],
          sn = B[B_off + 2*(m-1) - 2*(σn < 0)]; // <- shift
    for( let i=n1; i < m; i++ )
    {
      const di = B[B_off + 2*i];
      let   zi = (sn-di+σn)
               * (sn+di+σn);

      for( let j=n1; j < i; j++ )
      { const σj = μ[j],
              sj = B[B_off + 2*j - 2*(σj < 0)], // <- shift
              dj = B[B_off + 2*j];
        zi *= ( (sj-di+σj) / (dj-di) )
           *  ( (sj+di+σj) / (dj+di) );
      }

      for( let j=i; j < m-1; j++ )
      { const σj = μ[j],
              sj = B[B_off + 2*j - 2*(σj < 0)], // <- shift
              dj = B[B_off + 2*j+2];
        zi *= ( (sj-di+σj) / (dj-di) )
           *  ( (sj+di+σj) / (dj+di) );
      }

      B[B_off + 2*i+1] = Math.sign(B[B_off + 2*i+1]) * Math.sqrt(zi);
    }
  }

  // triple-merge-sort singular values from deflation
  // and the ones from the secular equation solutions
  for( let h=n0-1,
           i=n1-1,
           j=n1,
           k= 0; k < m; k++ )
  { let val = -Infinity,
       best = 4;
    if( j <  m )                    { best=2; val = μ[j] + B[B_off + 2*j - 2*(μ[j] < 0)]; }
    if( i >= n0 && ! (μ[i] < val) ) { best=1; val = μ[i]; }
    if( h >=  0 && ! (μ[h] < val) ) { best=0; val = μ[h]; }
    switch(best){
      case 0: innerOrder[h--] = k; continue;
      case 1: innerOrder[i--] = k; continue;
      case 2: innerOrder[j++] = k; continue;
      default: throw new Error('Assertion failed.');
    }
  }
  Object.freeze(innerOrder.buffer);

  // STEP 5:
  //   UPDATE U
  for( let i=n1; i < m; i++ )
  {
    const σi = μ[i],
          si = B[B_off + 2*i - 2*(σi < 0)]; // <- shift
    NORM.reset();
    for( let j=n1; j < m-1; j++ ) {
      const dj = B[B_off + 2*j  ],
            zj = B[B_off + 2*j+1],
          W_ij = ( zj / (dj-si-σi) )
               * ( dj / (dj+si+σi) );
      NORM.include(W[m*i+j] = W_ij);
    }
    NORM.include(W[m*i+m-1] = -1);
    const norm = NORM.result;

    if( ! (0 < norm) ) throw new Error('Assertion failed.');

    for( let j=n1; j < m; j++ )
      W[m*i+j] /= norm;
  }

  // rotate W
  for( let i=n0; i < n1; i++ )
    W[m*i+i] = 1;

  for( let i=n1; i-- > n0; )
  { const j = rotJ[i];
    if(   j < m-1 )
    { const c = B[B_off + 2*i  ],
            s = B[B_off + 2*i+1];
      for( let k=i; k < m; k++ ) // <- k should be started at i
      { const      W_ki = W[m*k+i],
                   W_kj = W[m*k+j];
        W[m*k+i] = W_ki* c  +  W_kj*s;
        W[m*k+j] = W_ki*-s  +  W_kj*c;
      }
    }
  }


  // TEST ONLY: write W -> U
  for( let i=0; i < n0; i++ )
    U[U_off + M*outerOrder[i]
            +   innerOrder[i]] = 1;

  for( let i=n0; i < m; i++ )
  for( let j=n0; j < m; j++ )
    U[U_off + M*outerOrder[i]
            +   innerOrder[j]] = W[m*j+i];


  // STEP 6:
  //   UPDATE V
  W.fill(0.0);

  for( let i=n1; i < m; i++ )
  {
    const σi = μ[i],
          si = B[B_off + 2*i - 2*(σi < 0)];
    NORM.reset();
    for( let j=n1; j < m; j++ ) {
      const dj = B[B_off + 2*j  ],
            zj = B[B_off + 2*j+1],
          W_ij = zj / (dj-si-σi)
                    / (dj+si+σi);
      NORM.include(W[m*i+j] = W_ij);
    }
    const norm = NORM.result;

    if( ! (0 < norm || i === m-1) ) throw new Error('Assertion failed.');

    for( let j=n1; j < m; j++ )
      W[m*i+j] /= norm;
  }

  if( 0 === B[B_off + 2*m-1] ) {
    for( let i=n1; i < m-1; i++ )
      W[m*(m-1)+i] = 0;
    W[m*(m-1)+(m-1)] = 1;
  }


  // rotate W
  for( let i=n0; i < n1; i++ )
    W[m*i+i] = 1;

  for( let i=n1; i-- > n0; )
  {
    const j = rotJ[i],
          c = B[B_off + 2*i  ],
          s = B[B_off + 2*i+1];
    for( let k=i; k < m; k++ ) // <- k should be started at i
    { const      W_ki = W[m*k+i],
                 W_kj = W[m*k+j];
      W[m*k+i] = W_ki* c  +  W_kj*s;
      W[m*k+j] = W_ki*-s  +  W_kj*c;
    }
  }


  // TEST ONLY: write W -> V
  for( let i=0; i < n0; i++ )
    V[V_off + N*innerOrder[i]
            +   outerOrder[i]] = 1;

  for( let i=n0; i < m; i++ )
  for( let j=n0; j < m; j++ )
    V[V_off + N*innerOrder[i]
            +   outerOrder[j]] = W[m*i+j];

  V[V_off + N*m+m] = 1;


  // STEP 7:
  //   GENERAL POSTPROCESSING
  for( let i=n1; i < m; i++ )
    μ[i] += B[B_off + 2*i - 2*(μ[i] < 0)];

  for( let i=0; i < m; i++ ) {
    const j = innerOrder[i];
    B[B_off + 2*j  ] = μ[i];
    B[B_off + 2*j+1] = NaN; // <- TEST ONLY
  }
}


/* Computes the svd of an upper bidiagonal matrix inplace.
 * The bidiagonal matrix is stored in a memory efficient
 * banded memory format where B[B_off+2*i] is the diagonal
 * entry B(i,i) and B[B_off+2*i+1] is the off-diagonal
 * element B(i,i+1).
 */
function _svd_dc( N, n, U,U_off, B,B_off, V,V_off )
{
  if(     n > N) throw new Error('Assertion failed.');
  if(0 >= n    ) throw new Error('Assertion failed.');

  if(2===n) return _svd_1x2(N, U,U_off, B,B_off, V,V_off);
  if(3===n) return _svd_2x3(N, U,U_off, B,B_off, V,V_off);

  // V1ᵀ[-1] ≙ Last  row of V1ᵀ
  // V2ᵀ[ 0] ≙ First row of V2ᵀ
  // S1'     ≙ The square part of S1
  // S2'     ≙ The square part of S2
  //
  // The first step is to divide the bidiagonal matrix 
  //     ┏                     ┓     ┏                     ┓
  //     ┃          ┆          ┃     ┃          ┆          ┃
  //     ┃    B1    ┆          ┃     ┃ U1⋅S1⋅V1 ┆          ┃
  //     ┃          ┆          ┃     ┃          ┆          ┃
  //     ┃┄┄┄┄┄┄┄┬┄┄┼┄┄┐       ┃(SVD)┃┄┄┄┄┄┄┄┬┄┄┼┄┄┐       ┃
  // B = ┃       ┆b1┆b2┆       ┃  =  ┃       ┆b1┆b2┆       ┃
  //     ┃       └┄┄┼┄┄┴┄┄┄┄┄┄┄┃     ┃       └┄┄┼┄┄┴┄┄┄┄┄┄┄┃
  //     ┃          ┆          ┃     ┃          ┆          ┃
  //     ┃          ┆    B2    ┃     ┃          ┆ U2⋅S2⋅V2 ┃
  //     ┃          ┆          ┃     ┃          ┆          ┃
  //     ┗                     ┛     ┗                     ┛
  //                                              ┏                     ┓                                                       
  //     ┏               ┓ ┏                    ┓ ┃          ┆          ┃
  //     ┃      ┆        ┃ ┃        ┆ ┆         ┃ ┃          ┆          ┃
  //     ┃  U1  ┆        ┃ ┃   S1'  ┆0┆         ┃ ┃    V1    ┆          ┃
  //     ┃      ┆        ┃ ┃        ┆ ┆         ┃ ┃          ┆          ┃
  //     ┃┄┄┄┄┄┄┼┄┐      ┃ ┃┄┄┄┄┄┄┄┄┴┄┼┄┄┄┄┄┄┄┄┄┃ ┃          ┆          ┃
  //   = ┃      ┆1┆      ┃ ┃b1⋅V1ᵀ[-1]┆b2⋅V2ᵀ[0]┃ ┃┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┄┄┄┄┄┃
  //     ┃      └┄┼┄┄┄┄┄┄┃ ┃┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┄┄┬┄┃ ┃          ┆          ┃
  //     ┃        ┆      ┃ ┃          ┆       ┆ ┃ ┃          ┆          ┃
  //     ┃        ┆  U2  ┃ ┃          ┆   S2' ┆0┃ ┃          ┆    V2    ┃
  //     ┃        ┆      ┃ ┃          ┆       ┆ ┃ ┃          ┆          ┃
  //     ┗               ┛ ┗                    ┛ ┃          ┆          ┃
  //                                              ┗                     ┛
  // Let's now get the middle matrix into neves (backwards seven) shape,
  // stating with permuting rows and columns.
  //
  // W1 := V1[:-1]
  //       ┏             ┓
  //       ┃   V2[ -2]   ┃
  //       ┃┄┄┄┄┄┄┄┄┄┄┄┄┄┃
  // W2 := ┃             ┃
  //       ┃   V2[:-2]   ┃
  //       ┃             ┃
  //       ┗             ┛
  // w1 := V1[-1]
  // w2 := V2[-1]
  //
  // C2 =   b1⋅V1ᵀ[-1,:-1]
  //         ┏          ╷                  ┓
  // C2 = b2⋅┃ V2[-2,0] ┆    V2ᵀ[0,:-2]    ┃
  //         ┗          ╵                  ┛
  // c1 = b1⋅V1[-1,0]
  // c2 = b2⋅V2[-1,0]                            ┏                     ┓                                                       
  //     ┏               ┓ ┏                   ┓ ┃          ┆          ┃
  //     ┃ ┆      ┆      ┃ ┃  C1  ┆  C2  ┆c1┆c2┃ ┃    W1    ┆          ┃
  //     ┃ ┆  U1  ┆      ┃ ┃┄┄┄┄┄┄┼┄┄┄┄┄┄┴┄┄┴┄┄┃ ┃          ┆          ┃
  //     ┃ ┆      ┆      ┃ ┃      ┆            ┃ ┃┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┄┄┄┄┄┃
  //     ┃┄┼┄┄┄┄┄┄┤      ┃ ┃  S1' ┆            ┃ ┃          ┆          ┃
  // B = ┃1┆      ┆      ┃ ┃      ┆            ┃ ┃          ┆    W2    ┃
  //     ┃┄┘      ├┄┄┄┄┄┄┃ ┃┄┄┄┄┄┄┼┄┄┄┄┄┄┐     ┃ ┃          ┆          ┃
  //     ┃        ┆      ┃ ┃      ┆      ┆     ┃ ┃┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┄┄┄┄┄┃
  //     ┃        ┆  U2  ┃ ┃      ┆  S2' ┆     ┃ ┃    w1    ┆          ┃
  //     ┃        ┆      ┃ ┃      ┆      ┆     ┃ ┃┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┄┄┄┄┄┃
  //     ┗               ┛ ┗                   ┛ ┃          ┆    w2    ┃
  //                                             ┗                     ┛
  // We can no use a single Given's rotation to turn the matrix into Neves shape.
  //           h := √( c1² + c2²)
  // c := c1 / h
  // s := c2 / h                                ┏                     ┓                                                       
  //     ┏               ┓ ┏                  ┓ ┃          ┆          ┃
  //     ┃ ┆      ┆      ┃ ┃  C1  ┆  C2  ┆h┆  ┃ ┃    W1    ┆          ┃
  //     ┃ ┆  U1  ┆      ┃ ┃┄┄┄┄┄┄┼┄┄┄┄┄┄┴┄┘  ┃ ┃          ┆          ┃
  //     ┃ ┆      ┆      ┃ ┃      ┆           ┃ ┃┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┄┄┄┄┄┃
  //     ┃┄┼┄┄┄┄┄┄┤      ┃ ┃  S1' ┆           ┃ ┃          ┆          ┃
  // B = ┃1┆      ┆      ┃ ┃      ┆           ┃ ┃          ┆    W2    ┃
  //     ┃┄┘      ├┄┄┄┄┄┄┃ ┃┄┄┄┄┄┄┼┄┄┄┄┄┄┐    ┃ ┃          ┆          ┃
  //     ┃        ┆      ┃ ┃      ┆      ┆    ┃ ┃┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┄┄┄┄┄┃
  //     ┃        ┆  U2  ┃ ┃      ┆  S2' ┆    ┃ ┃   c*w1   ┆  s*w2    ┃
  //     ┃        ┆      ┃ ┃      ┆      ┆    ┃ ┃┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┄┄┄┄┄┃
  //     ┗               ┛ ┗                  ┛ ┃  -s*w1   ┆  c*w2    ┃
  //                                            ┗                     ┛
  throw new Error('Not yet implemented!');
}


export function svd_dc(A)
{
  // SEE:
  //  - "INTRODUCTION OF DOUBLE DIVIDE AND CONQUER AND THE RECENT PROGRESS", TARO KONDA & YOSHIMASA NAKAMURA
  //  - "A Divide-and-Conquer Approach for Solving Singular Value Decomposition on a Heterogeneous System", Ding Liu & Ruixuan Li & David J. Lilja & Weijun Xiao
  A = asarray(A);
  if( A.dtype.startsWith('complex') )
    throw new Error('svd_dc(A): A.dtype must be float.');
  const [M,N] = A.shape.slice(-2);

  if( M > N ) {
    const [U,sv,V] = svd_dc(A)
    transpose_inplace(V);
    return [U.T,sv,V];
  }

  
}
