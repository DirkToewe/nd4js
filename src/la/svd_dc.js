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
import {transpose_inplace} from './transpose_inplace'
import {_giv_rot_rows} from './_giv_rot'
import {_svd_jac_angles} from './_svd_jac_utils'
import {FrobeniusNorm} from './norm'
import {_bidiag_decomp_horiz} from './bidiag'

// WARNING:
//   The following methods expect V to be stored in transposed/colum-major fashion.
//     * _svd_dc_1x2
//     * _svd_dc_2x3
//     * _svd_dc_neves
//     * _svd_dc_bidiag

/* Computes the SVD of the 1x2 matrix B
 */
export const _svd_dc_1x2 = ( N, U,U_off, F,B_off,V_off ) =>
{
  N     |= 0;
  U_off |= 0;
  B_off |= 0;
  V_off |= 0;

  let c = F[B_off  ],
      s = F[B_off+1];
  const  norm = Math.hypot(c,s);
  if(0!==norm) {
    c /= norm;
    s /= norm;

    F[B_off  ] = norm;
    F[B_off+1] = NaN;

    F[V_off      ] =  c;
    F[V_off+    1] = -s;
    F[V_off+N*1  ] =  s;
    F[V_off+N*1+1] =  c;
  }
  else {
    F[V_off      ] = 1;
    F[V_off+N*1+1] = 1;
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
export const _svd_dc_2x3 = ( N, U,U_off, F,B_off,V_off ) =>
{
  U_off  |= 0;
  B_off  |= 0;
  V_off  |= 0;
  N      |= 0; const
  M = N-1|  0;

  let b1 = F[B_off  ],
      b2 = F[B_off+1],
      b3 = F[B_off+2],
      b4 = F[B_off+3];
  // STEP 1:
  //   RQ-decompose down to 2x2 problem
  if( 0 !== b4 )
  { // eliminate B[1,2] (stored in b4)
    const     norm = Math.hypot(b3,b4),
      ca = b3/norm,
      sa = b4/norm;
    F[V_off + N*1+1] =  ca;
    F[V_off + N*2+1] =  sa;
    b3 = norm;
    b4 = -sa*b2;
    b2 =  ca*b2;
    if( 0 !== b4 )
    { // eliminate B[0,2] (stored in b4)
      const     norm = Math.hypot(b1,b4),
        cb = b1/norm,
        sb = b4/norm;
           b1 = norm;
      F[V_off      ] =  cb;
      F[V_off+N*1  ] =  sb*-sa;
      F[V_off+N*2  ] =  sb* ca;
      F[V_off+    2] = -sb;
      F[V_off+N*1+2] =  cb*-sa;
      F[V_off+N*2+2] =  cb* ca;
    }
    else {
      F[V_off      ] =   1;
      F[V_off+N*1+2] = -sa;
      F[V_off+N*2+2] =  ca;
    }
  }
  else {
    F[V_off      ] = 1;
    F[V_off+N*1+1] = 1;
    F[V_off+N*2+2] = 1;
  }

  // STEP 2;
  //   USE SINGLE JACOBI ROTATION TO COMPUTE SVD
  const [ca,sa, cb,sb] = _svd_jac_angles(b1,b2,
                                         0 ,b3);
  // WRITE SINGULAR VALUES TO B  
  F[B_off+1] = NaN; // <- TODO remove after testing
  F[B_off+3] = NaN; // <- TODO remove after testing

  F[B_off] = ca*b1*cb - ( ca*b2 + sa*b3)*sb;
                   b3 = (-sa*b2 + ca*b3)*cb - sa*b1*sb;
  const        s = b3 < 0 ? -1 : +1;
  F[B_off+2] = s * b3;

  // WRITE U
  U[U_off      ] =  ca;
  U[U_off+    1] = -sa*s;
  U[U_off+M*1  ] =  sa;
  U[U_off+M*1+1] =  ca*s;

  // ROTATE V
  F[V_off+1] = sb*F[V_off];
  F[V_off  ]*= cb;
  const         V1 = F[V_off+N*1]*cb - F[V_off+N*1+1]*sb;
    F[V_off+N*1+1] = F[V_off+N*1]*sb + F[V_off+N*1+1]*cb;
    F[V_off+N*1  ] = V1;
  const         V2 = F[V_off+N*2]*cb - F[V_off+N*2+1]*sb;
    F[V_off+N*2+1] = F[V_off+N*2]*sb + F[V_off+N*2+1]*cb;
    F[V_off+N*2  ] = V2;
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
 * 
 * Keep in mind that V is stored in a column-major order in this method.
 */
export const _svd_dc_neves = (N,n, U,U_off, F,B_off,V_off, I) =>
{
  U_off |= 0;
  B_off |= 0;
  V_off |= 0;
  N     |= 0; const M = N-1 | 0;
  n     |= 0; const m = n-1 | 0;

  if( I.length < M*3                                   ) throw new Error('Assertion failed: Integer work matrix I too small.');
  if( F.length < M*(M+2) + M*2/*B*/ + (M+1)*(M+1)/*V*/ ) throw new Error('Assertion failed: Scalar work matrix F too small.');

  // Amount of temp. float memory:
  //   - m*m entries to store the matrix to update U and V
  //   -   m entries to store the (shifted) singular values
  //   -   m entries to compute the matrix multiplication (row by row)
  const σ_off = M*(M+2) - m,
       mm_off = M*(M+2) - m*2,
        W_off = M*(M+2) - m*(m+2); // <- mm as in "matrix multiplication"

  // Amount of temp. int memory:
  //   - m entries for the outer (row   ) order of U & V update matrices
  //   - m entries for the inner (column) order of U & V update matrices
  //   - m entries for the rotation pairings from step 2
  const rot_off = m*2,
        inn_off = m, // <- inner order
        out_off = 0; // <- outer order

  if( n < 2 ) throw new Error('Assertion failed.');

  //     DIAGONAL ELEMENTS: d[i] = F[B_off + 2*i  ]
  // OFF-DIAGONAL ELEMENTS: z[i] = F[B_off + 2*i+1]

  // TODO remove checks after testing
  // ensure d[-1] === 1
  if( F[B_off + 2*m-2] !== 0 )
    throw new Error('Assertion failed.');
  // ensure d[i] >= d[i+1]
  for( let i=1; i < m; i++ )
    if( F[B_off + 2*(i-1)] < F[B_off + 2*i] )
      throw new Error('Assertion failed.');

  const NORM = new FrobeniusNorm();

  const [zNorm,scale] = function(){
    for( let i=0; i < m; i++ )
      NORM.include(F[B_off + 2*i+1]);

    let zNorm = NORM.result;

    for( let i=0; i < m; i++ )
      NORM.include(F[B_off + 2*i]);

    let scale = NORM.result;
    if( scale===0 )
        scale = 1;

    return [zNorm/scale, scale];
  }();

  // normalize
  for( let i=0; i < 2*m; i++ )
    F[B_off + i] /= scale;

  const TOL = eps(
    U instanceof Float32Array ? 'float32'
                              : 'float64'
  );


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
  //   implementation moves the deflated values to σ and not to
  //   the left.
  const n0 = function(){
    let n0 = 0;
    for( let j=m-1,
             i=m-1; i-- > 0; )
    { const di = F[  B_off + 2*i  ],
            zi = F[  B_off + 2*i+1],
            oi = I[out_off +   i  ];
      if( Math.abs(zi)/TOL <= di ) { // <- FIXME this could be estimated more accurately
        F[  σ_off + n0] = di;
        I[inn_off + n0] = oi; // <- used as temp. for outerOrder
                  ++n0;
      }
      else {        --j;
        F[  B_off + 2*j  ] = di;
        F[  B_off + 2*j+1] = F[B_off + 2*i+1];
        I[inn_off +   j  ] = oi; // <- used as temp. for outerOrder
      }
    }
    return n0;
  }();

  // innerOrder just used as temp. memory, move to outerOrder
  for( let i=0; i < n0; i++ )
    I[out_off + i] = I[inn_off + i];


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
  const n1 = function(){
    let n1 = n0;
    for( let j=m-1,
             i=m-1; i-- > n0; )
    { const di = F[B_off + 2*i  ],
            dj = F[B_off + 2*j  ], oi = I[inn_off + i],
            zi = F[B_off + 2*i+1],
            zj = F[B_off + 2*j+1], z = Math.hypot(zi,zj);
      if( (di-dj)/TOL <= di || ! isFinite( Math.sqrt(m) / (di-dj) ) ) // <- TODO find better criteria
      { const c = zj / z,
              s = zi / z;
        F[W_off + 2*n1  ] = c;
        F[W_off + 2*n1+1] = s;
        F[B_off + 2*j+1] = z;
        F[B_off + 2*j  ] = F[  σ_off + n1] = di;
                           I[out_off + n1] = oi; // <- used as temp. for outerOrder
                           I[rot_off + n1] =  j;
                                     ++n1;
      }
      else {      --j;
        F[B_off + 2*j  ] = di;
        F[B_off + 2*j+1] = F[B_off + 2*i+1];
        I[out_off + j  ] = oi;
      }
    }
    return n1;
  }();
  for( let i = 2*n0;
           i < 2*n1; i++ ) {
    F[B_off + i] = F[W_off + i];
                   F[W_off + i] = 0;
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
    let sLo =          F[B_off + 2*i  ],
        sHi = n1 < i ? F[B_off + 2*i-2] : (sLo + zNorm);

    if( sLo > sHi ) throw new Error('Assertion failed.');

    const shift = function(){
      if( i > n1 )
      { const mid = (sLo + sHi) / 2;
        let   sum = 1;
        for( let i=n1; i < m; i++ )
        { const  di = F[B_off + 2*i  ],
                 zi = F[B_off + 2*i+1];
          sum += zi / (di-mid)
              * (zi / (di+mid));
        }
        if( ! isFinite(sum) ) throw new Error('Assertion failed.');
        if( sum < 0 ) { const s=sHi; sLo = sLo-sHi; sHi = -Number.MIN_VALUE; return s; }
      }                 const s=sLo; sHi = sHi-sLo; sLo = +Number.MIN_VALUE; return s;
    }();

    // TODO: Make bisection more accurate (see nd.opt.root1d_bisect)
    for(;;)
    { // bisect
      const s  = (sLo + sHi) / 2;
      if(   s === sLo
         || s === sHi )
      { F[σ_off + i] = s; // FIXME at this point sLo and sHi still have to be compared
        break;
      }
      // evalue the secular equation
      let sum = 1;
      for( let i=n1; i < m; i++ )
      { const  di = F[B_off + 2*i  ],
               zi = F[B_off + 2*i+1];
        sum += zi / (di-shift-s)
            * (zi / (di+shift+s));
      }
      if( ! isFinite(sum) ) throw new Error('Assertion failed.');
      // adjust bisection bounds
      if( sum <= 0 ) sLo = s;
      if( sum >= 0 ) sHi = s;
    }
  }
  if( Math.abs(F[B_off + 2*m-1]) === 0 ) {
    F[B_off + 2*m-2] = 0;
    F[B_off + 2*m-1] = 0; F[σ_off + m-1] = 0;
  }


  // STEP 4:
  //   RECOMPUTE z TO ORTHOGONALIZE U & V
  //   (as originally suggested by Gu and Eisenstat)
  {
    const σn = F[σ_off + m-1],
          sn = F[B_off + 2*(m-1 - (σn < 0))]; // <- shift
    for( let i=n1; i < m; i++ )
    { const di = F[B_off + 2*i];
      let   zi = (sn-di+σn)
               * (sn+di+σn);

      for( let j=n1; j < i; j++ )
      { const σj = F[σ_off + j],
              sj = F[B_off + 2*(j - (σj < 0))], // <- shift
              dj = F[B_off + 2*j];
        zi *= ( (sj-di+σj) / (dj-di) )
           *  ( (sj+di+σj) / (dj+di) );
      }

      for( let j=i; j < m-1; j++ )
      { const σj = F[σ_off + j],
              sj = F[B_off + 2*(j - (σj < 0))], // <- shift
              dj = F[B_off + 2*j+2];
        zi *= ( (sj-di+σj) / (dj-di) )
           *  ( (sj+di+σj) / (dj+di) );
      }

      F[B_off + 2*i+1] = Math.sign(F[B_off + 2*i+1]) * Math.sqrt(zi);
    }
  }

  // triple-merge-sort singular values from deflation
  // and the ones from the secular equation solutions
  for( let h=n0-1,
           i=n1-1,
           j=n1,
           k= 0; k < m; k++ )
  { let val = -Infinity,
       best = 3;
    if( j <  m  ) { const σj = F[σ_off + j];                      best=2; val = σj + F[B_off + 2*(j - (σj < 0))]; }
    if( i >= n0 ) { const σi = F[σ_off + i]; if( ! (σi < val) ) { best=1; val = σi; } }
    if( h >=  0 ) { const σh = F[σ_off + h]; if( ! (σh < val) ) { best=0; val = σh; } }
    switch(best){
      case 0: I[inn_off + h--] = k; continue;
      case 1: I[inn_off + i--] = k; continue;
      case 2: I[inn_off + j++] = k; continue;
      default: throw new Error('Assertion failed.');
    }
  }


  // STEP 5:
  //   UPDATE U
  for( let i=n1; i < m; i++ )
  { const σi = F[σ_off + i],
          si = F[B_off + 2*(i - (σi < 0))]; // <- shift
    NORM.reset();
    for( let j=n1; j < m-1; j++ )
    { const dj = F[B_off + 2*j  ],
            zj = F[B_off + 2*j+1],
          W_ij = ( zj / (dj-si-σi) )
               * ( dj / (dj+si+σi) );
      NORM.include(F[W_off + m*i+j] = W_ij);
    }
    NORM.include(F[W_off + m*i+m-1] = -1);
    const norm = NORM.result;

    if( ! (0 < norm) ) throw new Error('Assertion failed.');

    for( let j=n1; j < m; j++ )
      F[W_off + m*i+j] /= norm;
  }

  // transpose dense part of W
  for( let i=n1;   i < m; i++ )
  for( let j=i ; ++j < m;     )
  { const W_ij = F[W_off + m*i+j];
                 F[W_off + m*i+j] = F[W_off + m*j+i];
                                    F[W_off + m*j+i] = W_ij;
  }

  if( n0 < n1 )
  { // init deflated region in W
    for( let i=n0; i < n1; i++ ) {
      F.fill(0.0, W_off + m*i+n0,
                  W_off + m*i+m);
      F[W_off + m*i+i] = 1;
    }
    for( let i=n1; i < m; i++ )
      F.fill(0.0, W_off + m*i+n0,
                  W_off + m*i+n1);
    // rotate W
    for(  let i = n1; i-- > n0; )
    {   const j = I[rot_off + i];
      if(     j < m-1 )
      { const c = F[B_off + 2*i  ],
              s = F[B_off + 2*i+1];
        _giv_rot_rows(F, m-i, W_off + m*i+i,
                              W_off + m*j+i, c,s);
      }
    }
  }

  // UPDATE U: U = U⋅Wᵀ
  for( let r=0; r < m; r++ )
  { // compute dense part of row (matrix multiplication)
    // U is fairly sparse so the matrix multiplication
    // is designed to benefit from that fact
    F.fill(0.0, mm_off + n0,
                mm_off + m);
                     for( let i=n0; i < m; i++ ) { const            U_ri = U[U_off + M*r + I[out_off + i]];
    if( 0 !== U_ri ) for( let j=n0; j < m; j++ ) { F[mm_off + j] += U_ri * F[W_off + m*i+j]; }}

    // compute deflated part of row
    for( let i=0; i < n0; i++ )
    { const                         c = I[out_off + i];
      F[mm_off + i] = U[U_off + M*r+c];
    }

    // write back row
    for( let i=0; i < m; i++ )
    { const         c  = I[inn_off + i];
      U[U_off + M*r+c] = F[ mm_off + i];
    }
  }


  // STEP 6:
  //   UPDATE V
  for( let i=n1; i < m; i++ )
  { const σi = F[σ_off + i],
          si = F[B_off + 2*(i - (σi < 0))];
    NORM.reset();
    for( let j=n1; j < m; j++ ) {
      const dj = F[B_off + 2*j  ],
            zj = F[B_off + 2*j+1],
          W_ij = zj / (dj-si-σi)
                    / (dj+si+σi);
      NORM.include(F[W_off + m*i+j] = W_ij);
    }
    const norm = NORM.result;

    if( ! (0 < norm || i === m-1) ) throw new Error('Assertion failed.');

    for( let j=n1; j < m; j++ )
      F[W_off + m*i+j] /= norm;
  }

  if( 0 === F[B_off + 2*m-1] )
  { for( let i=n1; i < m-1; i++ )
      F[W_off + m*(m-1)+i] = 0;
    F[W_off + m*(m-1)+(m-1)] = 1;
  }

  // transpose dense part of W
  for( let i=n1;   i < m; i++ )
  for( let j=i ; ++j < m;     ) {
    const W_ij = F[W_off + m*i+j];
                 F[W_off + m*i+j] = F[W_off + m*j+i];
                                    F[W_off + m*j+i] = W_ij;
  }

  if( n0 < n1 ) {
    // init deflated region in W
    for( let i=n0; i < n1; i++ ) {
      F.fill(0.0, W_off + m*i+n0,
                  W_off + m*i+m);
      F[W_off + m*i+i] = 1;
    }
    for( let i=n1; i < m; i++ )
      F.fill(0.0, W_off + m*i+n0,
                  W_off + m*i+n1);
    // rotate W
    for( let i=n1; i-- > n0; )
    { const  j = I[rot_off + i],
             c = F[B_off + 2*i  ],
             s = F[B_off + 2*i+1];
      _giv_rot_rows(F, m-i, W_off + m*i+i,
                            W_off + m*j+i, c,s);
    }
  }


  // UPDATE V: V = V⋅Wᵀ
  for( let r=0; r < n; r++ )
  { // compute dense part of row (matrix multiplication)
    // U is fairly sparse so the matrix multiplication
    // is designed to benefit from that fact
    F.fill(0.0, mm_off + n0,
                mm_off + m);
                     for( let i=n0; i < m; i++ ) { const            V_ri = F[V_off + N*r + I[out_off + i]];
    if( 0 !== V_ri ) for( let j=n0; j < m; j++ ) { F[mm_off + j] += V_ri * F[W_off + m*i+j]; }}

    // compute deflated part of row
    for( let i=0; i < n0; i++ ) {
      const                         c = I[out_off + i];
      F[mm_off + i] = F[V_off + N*r+c];
    }

    // write back row
    for( let i=0; i < m; i++ ) {
      const         c  = I[inn_off + i];
      F[V_off + N*r+c] = F[ mm_off + i];
    }
  }


  // STEP 7:
  //   GENERAL POSTPROCESSING
  for( let i=n1; i < m; i++ )
    F[σ_off + i] += F[B_off + 2*(i - (F[σ_off + i] < 0))];

  for( let i=0; i < m; i++ ) {
    const       j =  I[inn_off + i];
    F[B_off + 2*j  ] = F[σ_off + i] * scale;
    F[B_off + 2*j+1] = NaN; // <- TEST ONLY
  }
}


/* Computes the svd of an upper bidiagonal matrix inplace.
 * The bidiagonal matrix is stored in a memory efficient
 * banded memory format where B[B_off+2*i] is the diagonal
 * entry B(i,i) and B[B_off+2*i+1] is the off-diagonal
 * element B(i,i+1).
 */
export function _svd_dc_bidiag( N, n, U,U_off, F,B_off,V_off, I )
{
  // TODO:
  //   Replacing this recursion with a loop would likely
  //   increase performance by a significant amount
  // TODO:
  //   For small enough bidiagonal matrices (around m <= 64),
  //   the QR method could be used to compute the SVD not just
  //   fast but also more accurately, see:
  //     * "Accurate Singular Values of Bidiagonal Matrices"
  //       by James Demmel & W. Kahan
  //       Lapack Working Note (LAWN) No. 3
  //       http://www.netlib.org/lapack/lawnspdf/lawn03.pdf
  //     * "Computing Small Singular Values of Bidiagonal Matrices with Guaranteed High Relative Accuracy"
  //       by James Demmel & W. Kahan
  if(     n > N) throw new Error('Assertion failed.');
  if(1 >= n    ) throw new Error('Assertion failed.');

  if(2===n) return _svd_dc_1x2(N, U,U_off, F,B_off,V_off);
  if(3===n) return _svd_dc_2x3(N, U,U_off, F,B_off,V_off);

  const M  = N-1,
        m  = n-1,
        n0 = n >>> 1,
        m0 = n0-1;

  if( I.length < M  ) throw new Error('Assertion failed: Integer work matrix I too small.');
  if( F.length < M*2) throw new Error('Assertion failed: Scalar work matrix F too small.');


  // STEP 1: DIVIDE
  // --------------
  // The idea is to divide the bidiagonal matrix into an upper left block B1 and
  // a lower right block B2 with one row of the bidiagonal with the non-zero entries
  // b1 and b2 separating the two blocks.
  //       ┏                     ┓     ┏                     ┓
  //       ┃          ┆          ┃     ┃          ┆          ┃
  //       ┃    B1    ┆          ┃     ┃ U1⋅S1⋅V1 ┆          ┃
  //       ┃          ┆          ┃     ┃          ┆          ┃
  //       ┃┄┄┄┄┄┄┄┬┄┄┼┄┄┐       ┃(SVD)┃┄┄┄┄┄┄┄┬┄┄┼┄┄┐       ┃
  //   B = ┃       ┆b1┆b2┆       ┃  =  ┃       ┆b1┆b2┆       ┃
  //       ┃       └┄┄┼┄┄┴┄┄┄┄┄┄┄┃     ┃       └┄┄┼┄┄┴┄┄┄┄┄┄┄┃
  //       ┃          ┆          ┃     ┃          ┆          ┃
  //       ┃          ┆    B2    ┃     ┃          ┆ U2⋅S2⋅V2 ┃
  //       ┃          ┆          ┃     ┃          ┆          ┃
  //       ┗                     ┛     ┗                     ┛
  _svd_dc_bidiag(N,  n0, U,U_off          , F,B_off,      V_off          , I);
  _svd_dc_bidiag(N,n-n0, U,U_off + M*n0+n0, F,B_off+2*n0, V_off + N*n0+n0, I);
                         U[U_off + M*m0+m0] = 1;

  const b1 = F[B_off + 2*m0  ],
        b2 = F[B_off + 2*m0+1];

//  if( ! isFinite(b1) ) throw new Error('Assertion failed.');
//  if( ! isFinite(b2) ) throw new Error('Assertion failed.');

  // The orthogonal matrix U1,U2,V1,V2 can moved to separate left and right orthogonal
  // matrices. V1ᵀ[-1] is the last row of the transpose of V1. V2ᵀ[0] is the first row of
  // the transpose of V2.                                                   ┏                     ┓
  //     ┏                     ┓   ┏               ┓ ┏                    ┓ ┃          ┆          ┃
  //     ┃          ┆          ┃   ┃      ┆        ┃ ┃        ┆ ┆         ┃ ┃          ┆          ┃
  //     ┃ U1⋅S1⋅V1 ┆          ┃   ┃  U1  ┆        ┃ ┃   S1'  ┆0┆         ┃ ┃    V1    ┆          ┃
  //     ┃          ┆          ┃   ┃      ┆        ┃ ┃        ┆ ┆         ┃ ┃          ┆          ┃
  //     ┃┄┄┄┄┄┄┄┬┄┄┼┄┄┐       ┃   ┃┄┄┄┄┄┄┼┄┐      ┃ ┃┄┄┄┄┄┄┄┄┴┄┼┄┄┄┄┄┄┄┄┄┃ ┃          ┆          ┃
  // B = ┃       ┆b1┆b2┆       ┃ = ┃      ┆1┆      ┃ ┃b1⋅V1ᵀ[-1]┆b2⋅V2ᵀ[0]┃ ┃┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┄┄┄┄┄┃
  //     ┃       └┄┄┼┄┄┴┄┄┄┄┄┄┄┃   ┃      └┄┼┄┄┄┄┄┄┃ ┃┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┄┄┬┄┃ ┃          ┆          ┃
  //     ┃          ┆          ┃   ┃        ┆      ┃ ┃          ┆       ┆ ┃ ┃          ┆          ┃
  //     ┃          ┆ U2⋅S2⋅V2 ┃   ┃        ┆  U2  ┃ ┃          ┆   S2' ┆0┃ ┃          ┆    V2    ┃
  //     ┃          ┆          ┃   ┃        ┆      ┃ ┃          ┆       ┆ ┃ ┃          ┆          ┃
  //     ┗                     ┛   ┗               ┛ ┗                    ┛ ┃          ┆          ┃
  //                                                                        ┗                     ┛
  for( let i= 0; i < m0; i++ ) F[B_off + 2*i +1] = b1 * F[V_off + N*m0 + i];
  for( let i=n0; i < m ; i++ ) F[B_off + 2*i +1] = b2 * F[V_off + N*n0 + i];

  // With the following definitions:
  //                 ┏          ╷    ┓
  //   b1⋅V1ᵀ[-1] =: ┃    R1    ┆ r1 ┃
  //                 ┗          ╵    ┛
  //                 ┏          ╷    ┓
  //   b2⋅V2ᵀ[ 0] =: ┃    R2    ┆ r2 ┃
  //                 ┗          ╵    ┛
  //   h := √(r1² + r2²)
  //   c = r1 / h
  //   s = r2 / h
  //         ┏             ┓
  //         ┃             ┃
  //         ┃      W1     ┃
  //   V1 =: ┃             ┃
  //         ┃┄┄┄┄┄┄┄┄┄┄┄┄┄┃
  //         ┃      w1     ┃
  //         ┗             ┛
  //         ┏             ┓
  //         ┃             ┃
  //         ┃      W2     ┃
  //   V2 =: ┃             ┃
  //         ┃┄┄┄┄┄┄┄┄┄┄┄┄┄┃
  //         ┃      w2     ┃
  //         ┗             ┛
  // And a single Givens rotation, we can make the right-most column of the middle matrix zero.
  //                                                  ┏                     ┓                                            ┏                     ┓
  //       ┏               ┓ ┏                      ┓ ┃          ┆          ┃   ┏               ┓ ┏                    ┓ ┃          ┆          ┃
  //       ┃      ┆        ┃ ┃        ┆  ┆          ┃ ┃          ┆          ┃   ┃      ┆        ┃ ┃        ┆ ┆         ┃ ┃     W1   ┆          ┃
  //       ┃  U1  ┆        ┃ ┃   S1'  ┆ 0┆          ┃ ┃    V1    ┆          ┃   ┃  U1  ┆        ┃ ┃   S1'  ┆0┆         ┃ ┃          ┆          ┃
  //       ┃      ┆        ┃ ┃        ┆  ┆          ┃ ┃          ┆          ┃   ┃      ┆        ┃ ┃        ┆ ┆         ┃ ┃┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┄┄┄┄┄┃
  //       ┃┄┄┄┄┄┄┼┄┐      ┃ ┃┄┄┄┄┄┄┄┄┼┄┄┼┄┄┄┄┄┄┄┬┄┄┃ ┃          ┆          ┃   ┃┄┄┄┄┄┄┼┄┐      ┃ ┃┄┄┄┄┄┄┄┄┼┄┼┄┄┄┄┄┄┄┬┄┃ ┃   c*w1   ┆  s*w2    ┃
  //   B = ┃      ┆1┆      ┃ ┃   R1   ┆r1┆  R2   ┆r2┃ ┃┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┄┄┄┄┄┃ = ┃      ┆1┆      ┃ ┃   R1   ┆h┆  R2   ┆ ┃ ┃┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┄┄┄┄┄┃
  //       ┃      └┄┼┄┄┄┄┄┄┃ ┃┄┄┄┄┄┄┄┄┴┄┄┼┄┄┄┄┄┄┄┼┄┄┃ ┃          ┆          ┃   ┃      └┄┼┄┄┄┄┄┄┃ ┃┄┄┄┄┄┄┄┄┴┄┼┄┄┄┄┄┄┄┤ ┃ ┃          ┆          ┃
  //       ┃        ┆      ┃ ┃           ┆       ┆  ┃ ┃          ┆          ┃   ┃        ┆      ┃ ┃          ┆       ┆0┃ ┃          ┆    W2    ┃
  //       ┃        ┆  U2  ┃ ┃           ┆   S2' ┆ 0┃ ┃          ┆    V2    ┃   ┃        ┆  U2  ┃ ┃          ┆   S2' ┆ ┃ ┃          ┆          ┃
  //       ┃        ┆      ┃ ┃           ┆       ┆  ┃ ┃          ┆          ┃   ┃        ┆      ┃ ┃          ┆       ┆ ┃ ┃┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┄┄┄┄┄┃
  //       ┗               ┛ ┗                      ┛ ┃          ┆          ┃   ┗               ┛ ┗                    ┛ ┃  -s*w1   ┆  c*w2    ┃
  //                                                  ┗                     ┛                                            ┗                     ┛
  const [c,s,h] = function(){
    let   c = b1 * F[V_off + N*m0 + m0],
          s = b2 * F[V_off + N*n0 + m ];
    const h = Math.hypot(c,s);
    return [c/h, s/h, h];
  }();

//  if( ! isFinite(h) ) throw new Error('Assertion failed.'); 

  F[B_off + 2*m0  ] = 0;
  F[B_off + 2*m0+1] = h;

  if( 0 !== h ) {
    for( let i= 0; i < n0; i++ ) { F[V_off + N*i + m ] = F[V_off + N*i + m0] * -s;
                                                         F[V_off + N*i + m0] *= c; }
    for( let i=n0; i < n ; i++ ) { F[V_off + N*i + m0] = F[V_off + N*i + m ] *  s;
                                                         F[V_off + N*i + m ] *= c; }
  }


  // STEP 2: CONQUER
  // ---------------
  // After step 1, row and column swaps can be used to turn the middle matrix into a Neves matrix.
  // See _svd_dc_neves() for more information about the shape of the Neves matrix. An integer array
  // called the "outer order" is used to keep track of the row and column swaps. The outer order
  // indicates where the rows and column of the middle matrix have originaally been.
  I[m-1] = m0;
  // merge sort diagonal entries
  for( let i = 0,
           j =n0,
           k = 0; k < m-1; k++ )
    I[k] = j >= m || (i < m0 && F[B_off + 2*i]
                             >= F[B_off + 2*j]) ? i++ : j++;

//  if( I.slice(0,m).sort().some((i,j) => i !== j) ) throw new Error('Assertion failed.');

  for( let i=0; i < m; i++ ) {
    const                     j = I[i];
    F[2*i  ] = F[B_off + 2*j  ];
    F[2*i+1] = F[B_off + 2*j+1];
  }

  for( let i=0; i < m; i++ ) {
    F[B_off + 2*i  ] = F[2*i  ];
    F[B_off + 2*i+1] = F[2*i+1];
  }

  _svd_dc_neves(N,n, U,U_off, F,B_off,V_off, I);
}


export function _svd_dc( M,N, U,U_off, sv,sv_off, V,V_off, I, F )
{
  if( M > N ) throw new Error('Assertion failed.');

  if( I.length < M*3                                   ) throw new Error('Assertion failed: Integer work matrix I too small.');
  if( F.length < M*(M+2) + M*2 + (M+1)*(M+1) + (M+1)*N ) throw new Error( 'Assertion failed: Scalar work matrix F too small.');

  const B_off =          M   *(M+2),
       V1_off =  B_off + M*2,
       V2_off = V1_off +(M+1)*(M+1);

  F.fill(0.0, V1_off,
              V2_off + N);

  // COPY A FROM V TO V2
  for( let i=0; i < M; i++ )
  for( let j=0; j < N; j++ )
    F[V2_off+N + N*i+j] = V[V_off + N*i+j];

  V.fill(0.0,  V_off,
               V_off + M*N);

  _bidiag_decomp_horiz(M,N, U,U_off, F,0, F,V2_off);

  for( let i=0; i < M; i++ ) {
    F[B_off + 2*i  ] = F[(M+1)*i + i  ];
    F[B_off + 2*i+1] = F[(M+1)*i + i+1];
  }

  _svd_dc_bidiag(
    M+1, M+1,
    V, V_off,
    F, B_off,
      V1_off,
    I
  );

  for( let i=0; i < M; i++ )
    sv[sv_off + i] = F[B_off + 2*i];

  // update U = U @ U2 (U2 is stored in V)
  for( let i=0; i < M; i++ )
  { F.fill(0.0, 0,M);
    for( let k=0; k < M; k++ )
    for( let j=0; j < M; j++ )                    F[j] += U[U_off + M*i+k] * V[V_off + M*k+j];
    for( let j=0; j < M; j++ ) U[U_off + M*i+j] = F[j];
  }

  // update V = V1 @ V2 (keep in mind that VV was computed in a transposed/column-major fashion)
  V.fill(0.0, V_off,V_off + M*M);
  for( let k=0; k < M+1; k++ )
  for( let i=0; i < M;   i++ )
  for( let j=0; j < N;   j++ ) V[V_off + N*i+j] += F[V1_off + (M+1)*k+i] * F[V2_off + N*k+j];
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
    const [U,sv,V] = svd_dc(A.T)
    transpose_inplace(U);
    return [V.T,sv,U];
  }

  const V_shape = A.shape,
        U_shape = V_shape.slice(),
       sv_shape = V_shape.slice(0, -1);
  U_shape[U_shape.length - 1] = M;

  const DType = A.dtype === 'float32' ? 'float32' : 'float64',
        DTypeArray = ARRAY_TYPES[DType];
  A = A.data;

  const len = A.length / (M*N),
          V = DTypeArray.from(A);

  A = undefined;

  const U = new DTypeArray(len * M*M),
       sv = new DTypeArray(len *   M),
        F = new DTypeArray(M*(M+2)/*F*/ + M*2/*B*/ + (M+1)*(M+1)/*V1*/ + (M+1)*Math.max(1+M,N)/*V2*/),
        I = new Int32Array(M*3);

  for(
    let U_off=0,
       sv_off=0,
        V_off=0; sv_off < sv.length;
        U_off += M*M,
       sv_off +=   M,
        V_off +=   M*N
  ) _svd_dc(M,N, U,U_off, sv,sv_off, V,V_off, I, F);

  return [
    new NDArray( U_shape, U),
    new NDArray(sv_shape,sv),
    new NDArray( V_shape, V)
  ];
}
