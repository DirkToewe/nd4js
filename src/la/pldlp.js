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
import {KahanSum} from '../kahan_sum'
import {asarray, NDArray} from '../nd_array'


// REFERENCES
// ----------
// .. [1] "DSYTF2.f"
//         Reference-LAPACK v3.9.0
//         https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/SRC/dsytf2.f
// .. [2] "Matrix Computations" 4th Edition
//         Chapter 4   "Special Linear Systems"
//         Section 4.4 "Symmetric Indefinite Systems"
//         pp. 186ff
//         Hindustan Book Agency, 2015


const _pldlp_decomp_1x1 = (M,N, LD,LD_off, k) =>
{
  if( 0 !== M%1 ) throw new Error('Assertion failed.');
  if( 0 !== N%1 ) throw new Error('Assertion failed.');
  if( 0 !== k%1 ) throw new Error('Assertion failed.');
  if( 0 !== LD_off   %1 ) throw new Error('Assertion failed.');
  if( 0 !== LD.length%1 ) throw new Error('Assertion failed.');
  if( ! (0 <= M) ) throw new Error('Assertion failed.');
  if( ! (M <= N) ) throw new Error('Assertion failed.');
  if( ! (k >= 0) ) throw new Error('Assertion failed.');
  if( ! (k <  M) ) throw new Error('Assertion failed.');

  const                                                                       D_kk = LD[LD_off + N*k+k];
  for( let i=k; ++i < M; ) {               const LD_ik = LD[LD_off + N*i+k] / D_kk;
  for( let j=k; j++ < i; ) LD[LD_off + N*i+j] -= LD_ik * LD[LD_off + N*j+k]; }

  for( let i=k; ++i < M; ) LD[LD_off + N*i+k] /= D_kk;
}


const _pldlp_decomp_2x2 = (M,N, LD,LD_off, k) =>
{
  if( 0 !== M%1 ) throw new Error('Assertion failed.');
  if( 0 !== N%1 ) throw new Error('Assertion failed.');
  if( 0 !== k%1 ) throw new Error('Assertion failed.');
  if( 0 !== LD_off   %1 ) throw new Error('Assertion failed.');
  if( 0 !== LD.length%1 ) throw new Error('Assertion failed.');
  if( ! ( 0 <= M ) ) throw new Error('Assertion failed.');
  if( ! ( M <= N ) ) throw new Error('Assertion failed.');
  if( ! ( k >= 0 ) ) throw new Error('Assertion failed.');
  if( ! ( k < M-1) ) throw new Error('Assertion failed.');

  // INVERSE OF DIAGONAL 2X2 BLOCK
  const tmp = LD[LD_off + N*(k+1)+ k   ],
        D00 = LD[LD_off + N*(k+1)+(k+1)] / tmp,
        D11 = LD[LD_off + N* k   + k   ] / tmp,
        D10 =          1 / (D00*D11 - 1) / tmp; // <- TODO: understand this underflow-safe Vodoo magic

  for( let j=k+2; j < M; j++ )
  {
    const s0 = D10 * ( D00 * LD[LD_off + N*j + k  ] - LD[LD_off + N*j + k+1] ),
          s1 = D10 * ( D11 * LD[LD_off + N*j + k+1] - LD[LD_off + N*j + k  ] );

    for( let i=j; i < M; i++ )
      LD[LD_off + N*i+j] -= LD[LD_off + N*i+k  ]*s0 +
                            LD[LD_off + N*i+k+1]*s1;

    LD[LD_off + N*j+k  ] = s0;
    LD[LD_off + N*j+k+1] = s1;
  }
}


const _pldlp_swap = (M,N, LD,LD_off, P,P_off, k,l) =>
{
  if( ! (k <= l) ) throw new Error('Assertion failed: ' + JSON.stringify({k,l}));
  if( ! (M <= N) ) throw new Error('Assertion failed.');

  const P_k = P[P_off + k];
              P[P_off + k] = P[P_off + l];
                             P[P_off + l] = P_k;

  for( let i=0; i < k; i++ ) {
    const LD_ki = LD[LD_off + N*k+i];
                  LD[LD_off + N*k+i] = LD[LD_off + N*l+i];
                                       LD[LD_off + N*l+i] = LD_ki;
  }

  for( let i=k; ++i < l; ) {
    const LD_ik = LD[LD_off + N*i+k];
                  LD[LD_off + N*i+k] = LD[LD_off + N*l+i];
                                       LD[LD_off + N*l+i] = LD_ik;
  }

  const LD_kk = LD[LD_off + N*k+k];
                LD[LD_off + N*k+k] = LD[LD_off + N*l+l];
                                     LD[LD_off + N*l+l] = LD_kk;

  for( let i=l; ++i < M; ) {
    const LD_ik = LD[LD_off + N*i+k];
                  LD[LD_off + N*i+k] = LD[LD_off + N*i+l];
                                       LD[LD_off + N*i+l] = LD_ik;
  }
}


export function _pldlp_decomp(M,N, LD,LD_off, P,P_off)
{
  if( ! (M <= N) ) throw new Error('Assertion failed.');

  for( let i=0; i < M; i++ )
    P[P_off + i] = i;

  const α = (Math.sqrt(17) + 1) / 8;

  for( let k=0; k < M; k++ )
  {
    let is1x1 = true;

    if( k < M-1 )
    {
      let r,
        A_rk = -Infinity,
        A_kk = Math.abs(LD[LD_off + N*k+k]);
  
      for( let i=k; ++i < M; )
      { const         A_ik = Math.abs(LD[LD_off + N*i+k]);
        if( !(A_rk >= A_ik) ) { // <- handles NaN
              A_rk  = A_ik;
                r   =   i;
        }
      }

      if( !(0 < Math.max(A_rk,A_kk)) )
        throw new Error('_pldlp_decomp(M,N, LD,LD_off, P,P_off): Zero column or NaN encountered.');

      if( A_kk < α*A_rk )
      { let s,
          A_sr = -Infinity;
        for( let i=k; i < r; i++ )
        { const         A_ri = Math.abs(LD[LD_off + N*r+i]);
          if( !(A_sr >= A_ri) ) { // <- handles NaN
                A_sr  = A_ri;
                  s   =    i;
          }
        }
        for( let i=r; ++i < M; )
        { const         A_ir = Math.abs(LD[LD_off + N*i+r]);
          if( !(A_sr >= A_ir) ) { // <- handles NaN
                A_sr  = A_ir;
                  s   =   i;
          }
        }
        if( r===s )
          throw new Error('Assertion failed.');

        if( A_kk < α * A_rk * (A_rk / A_sr) )
        {
          const A_rr = Math.abs(LD[LD_off + N*r+r])
          if(   A_rr < α*A_sr ) {
            is1x1 = false;
            P[P_off + r] ^= -1;
            ++k;
          }
          if( k !== r )
            _pldlp_swap(M,N, LD,LD_off, P,P_off, k,r);
        }
      }
    }

    if(is1x1) _pldlp_decomp_1x1(M,N, LD,LD_off, k  );
    else      _pldlp_decomp_2x2(M,N, LD,LD_off, k-1);
  }
}


export function pldlp_decomp(S)
{
  S = asarray(S);
  const
    DTypeArray = ARRAY_TYPES[S.dtype === 'float32' ? 'float32' : 'float64'],
    LD_shape =  S.shape,
     P_shape = LD_shape.slice(0,-1),
       [N,M] = LD_shape.slice(  -2);
  S = S.data;

  if( N !== M )
    throw new Error('Last two dimensions must be quadratic.')

  const LD = new DTypeArray(S.length  ),
         P = new Int32Array(S.length/N);

  for( let LD_off=0; LD_off < LD.length; LD_off += N*N )
  {
    for( let i=0; i < N; i++ )
    for( let j=0; j <=i; j++ )
      LD[LD_off + N*i+j] = S[LD_off + N*i+j];

    const P_off = LD_off / N;

    _pldlp_decomp(N,N, LD,LD_off, P,P_off);
  }

  return [
    new NDArray(LD_shape,LD),
    new NDArray( P_shape, P)
  ];
}


export function pldlp_l( LD, P )
{
  if( null == P ) [LD,P] = LD;

  LD = asarray(LD);
   P = asarray( P);
  if(LD.ndim < 2 ) throw new Error('pldlp_solve(LD,P,y): LD must be at least 2D.');
  if( P.ndim < 1 ) throw new Error('pldlp_solve(LD,P,y): P must be at least 1D.');

  const N = LD.shape[LD.ndim-2];
  if(   N!==LD.shape[LD.ndim-1] ) throw new Error('pldlp_solve(LD,P,y): LD must be square.');
  if(   N!== P.shape[ P.ndim-1] ) throw new Error('pldlp_solve(LD,P,y): LD and P shape mismatch.');

  const
    ndim = Math.max(LD.ndim, P.ndim+1),
    shape = Int32Array.from({ length: ndim }, ()=>1);
  shape[ndim-2] = N;
  shape[ndim-1] = N;

  // FIND COMMON (BROADCASTED) SHAPE
  for( let shp of [
   LD.shape.slice(0,-2),
    P.shape.slice(0,-1)
  ])
    for( let i=ndim-2, j=shp.length; i-- > 0 && j-- > 0; )
      if( 1 === shape[i] )
        shape[i] = shp[j];
      else if( shape[i] !== shp[j] && shp[j] !== 1 )
        throw new Error('Shapes are not broadcast-compatible.');

  // GENERATE RESULT DATA
  const         DType = [LD,P].every(x => x.dtype === 'float32') ? 'float32' : 'float64',
                DTypeArray = ARRAY_TYPES[DType],
    L_dat = new DTypeArray( shape.reduce((a,b) => a*b, 1) ),
   LD_dat =LD.data,
    P_dat = P.data;
  let
   LD_off = 0,LD_stride = 1,
    P_off = 0, P_stride = 1,
    L_off = 0;

  function solv(d) {
    if( d === ndim-2 ) {
     LD_stride = N*N;
      P_stride = N;

      for( let i=0; i < N; i++ )
      {
        const             I = i - (P_dat[P_off + i] < 0);
        for( let j=0; j < I; j++ )
          L_dat[L_off + N*i+j] = LD_dat[LD_off + N*i+j];

        L_dat[L_off + N*i+i] = 1;
        
//        L_dat[L_off + N*i+i] = LD_dat[LD_off + N*i+i];
//        if( P[P_off + i] < 0 )
//          L_dat[ L_off + N* i   + i-1] =
//          L_dat[ L_off + N*(i-1)+ i  ] =
//         LD_dat[LD_off + N* i   + i-1];
      }

      LD_off +=LD_stride;
       L_off +=LD_stride;
       P_off += P_stride;

      return;
    }
    for( let l=shape[d]; ; l-- ) {
      solv(d+1);
      if( l===1 ) break;
      if( ! (LD.shape[ d - ndim +LD.ndim  ] > 1) ) LD_off -=LD_stride;
      if( ! ( P.shape[ d - ndim + P.ndim+1] > 1) )  P_off -= P_stride;
    }
    LD_stride *=LD.shape[ d - ndim +LD.ndim  ] || 1;
     P_stride *= P.shape[ d - ndim + P.ndim+1] || 1;
  }
  solv(0);

  return new NDArray(shape,L_dat);
}


export function pldlp_d( LD, P )
{
  if( null == P ) [LD,P] = LD;

  LD = asarray(LD);
   P = asarray( P);
  if(LD.ndim < 2 ) throw new Error('pldlp_solve(LD,P,y): LD must be at least 2D.');
  if( P.ndim < 1 ) throw new Error('pldlp_solve(LD,P,y): P must be at least 1D.');

  const N = LD.shape[LD.ndim-2];
  if(   N!==LD.shape[LD.ndim-1] ) throw new Error('pldlp_solve(LD,P,y): LD must be square.');
  if(   N!== P.shape[ P.ndim-1] ) throw new Error('pldlp_solve(LD,P,y): LD and P shape mismatch.');

  const
    ndim = Math.max(LD.ndim, P.ndim+1),
    shape = Int32Array.from({ length: ndim }, ()=>1);
  shape[ndim-2] = N;
  shape[ndim-1] = N;

  // FIND COMMON (BROADCASTED) SHAPE
  for( let shp of [
   LD.shape.slice(0,-2),
    P.shape.slice(0,-1)
  ])
    for( let i=ndim-2, j=shp.length; i-- > 0 && j-- > 0; )
      if( 1 === shape[i] )
        shape[i] = shp[j];
      else if( shape[i] !== shp[j] && shp[j] !== 1 )
        throw new Error('Shapes are not broadcast-compatible.');

  // GENERATE RESULT DATA
  const         DType = [LD,P].every(x => x.dtype === 'float32') ? 'float32' : 'float64',
                DTypeArray = ARRAY_TYPES[DType],
    D_dat = new DTypeArray( shape.reduce((a,b) => a*b, 1) ),
   LD_dat =LD.data,
    P_dat = P.data;
  let
   LD_off = 0,LD_stride = 1,
    P_off = 0, P_stride = 1,
    D_off = 0;

  function solv(d) {
    if( d === ndim-2 ) {
     LD_stride = N*N;
      P_stride = N;

      for( let i=0; i < N; i++ )
      {
        D_dat[D_off + N*i+i] = LD_dat[LD_off + N*i+i];

        if( P_dat[ P_off + i] < 0 )
            D_dat[ D_off + N* i   + i-1] =
            D_dat[ D_off + N*(i-1)+ i  ] = LD_dat[LD_off + N* i   + i-1];
      }

      LD_off +=LD_stride;
       D_off +=LD_stride;
       P_off += P_stride;

      return;
    }
    for( let l=shape[d]; ; l-- ) {
      solv(d+1);
      if( l===1 ) break;
      if( ! (LD.shape[ d - ndim +LD.ndim  ] > 1) ) LD_off -=LD_stride;
      if( ! ( P.shape[ d - ndim + P.ndim+1] > 1) )  P_off -= P_stride;
    }
    LD_stride *=LD.shape[ d - ndim +LD.ndim  ] || 1;
     P_stride *= P.shape[ d - ndim + P.ndim+1] || 1;
  }
  solv(0);

  return new NDArray(shape,D_dat);
}


export function pldlp_p( LD, P )
{
  if( null == P ) [LD,P] = LD;

  LD = asarray(LD);
   P = asarray( P);
  if(LD.ndim < 2 ) throw new Error('pldlp_solve(LD,P,y): LD must be at least 2D.');
  if( P.ndim < 1 ) throw new Error('pldlp_solve(LD,P,y): P must be at least 1D.');

  const N = LD.shape[LD.ndim-2];
  if(   N!==LD.shape[LD.ndim-1] ) throw new Error('pldlp_solve(LD,P,y): LD must be square.');
  if(   N!== P.shape[ P.ndim-1] ) throw new Error('pldlp_solve(LD,P,y): LD and P shape mismatch.');

  for( let i=LD.ndim-2,
           j= P.ndim-1; i-- > 0 &&
                        j-- > 0; )
    if( P.shape[j] !== 1 &&
       LD.shape[i] !== 1 &&
       LD.shape[i] !== P.shape[j] )
      throw new Error('Shapes are not broadcast-compatible.');

  const P_shape =  P.shape;
  P = Int32Array.from(P.data);
  LD= undefined;

  const test = new Uint8Array(N);

  for( let P_off=P.length; (P_off-=N) >= 0; )
  {
    test.fill(false);

    for( let i=N; i-- > 0; )
    {
      let p = P[P_off + i];
      if( p < 0 ) {
        if( !(0 < i) || P[P_off + i-1] < 0 )
          throw new Error('pldlp_solve(LD,P,y): P is invalid.');
        p ^= -1;
      }

      if( !(p >=0) ) throw new Error('pldlp_solve(LD,P,y): P contains out-of-bounds indices.');
      if( !(p < N) ) throw new Error('pldlp_solve(LD,P,y): P contains out-of-bounds indices.');

      if( test[p] )
        throw new Error('pldlp_solve(LD,P,y): P contains duplicate indices.');

      test[p] = true;
      P[P_off + i] = p;
    }

    if( test.some(x => !x) )
      throw new Error('pldlp_solve(LD,P,y): P is invalid.');
  }

  return new NDArray(P_shape, P);
}


export function _pldlp_solve(M,N,O, LD,LD_off, P,P_off, X,X_off, tmp)
{
  if( 0 !== LD.length%1 ) throw new Error('Assertion failed.');
  if( 0 !==  P.length%1 ) throw new Error('Assertion failed.');
  if( 0 !==  X.length%1 ) throw new Error('Assertion failed.');
  if( 0 !==tmp.length%1 ) throw new Error('Assertion failed.');
  if( 0 !==LD_off%1 ) throw new Error('Assertion failed.');
  if( 0 !== P_off%1 ) throw new Error('Assertion failed.');
  if( 0 !== X_off%1 ) throw new Error('Assertion failed.');
  if( 0 !== M%1 ) throw new Error('Assertion failed.');
  if( 0 !== N%1 ) throw new Error('Assertion failed.');
  if( 0 !== O%1 ) throw new Error('Assertion failed.');

  if( LD.length -LD_off < M*N ) throw new Error('Assertion failed.');
  if(  P.length - P_off < M   ) throw new Error('Assertion failed.');
  if(  X.length - X_off < M*O ) throw new Error('Assertion failed.');
  if(tmp.length         < M*O ) throw new Error('Assertion failed.');

  if( !(0 < M) ) throw new Error('Assertion failed.');
  if( !(0 < N) ) throw new Error('Assertion failed.');
  if( !(0 < O) ) throw new Error('Assertion failed.');

  if( !(M <= N) ) throw new Error('Assertion failed.');

  // PERMUTE INPUT
  for( let i=0; i < M; i++ ) { let I = P[P_off + i]; I ^= -(I<0);
  for( let j=0; j < O; j++ )
    tmp[O*i+j] = X[X_off + O*I+j];
  }

  // FORWARD SUBSTITUTION
  for( let i=1; i < M; i++ ) { const K = i - (P[P_off + i] < 0);
  for( let k=0; k < K; k++ )
  for( let j=0; j < O; j++ )
    tmp[O*i+j] -= LD[LD_off + N*i+k] * tmp[O*k+j];
  }

  // SCALING
  for( let i=0; i < M; i++ )
  {
    if( i < M-1 && P[P_off + i+1] < 0 )
    { // 2x2 BLOCK
      const s = LD[LD_off + N*(i+1) + i  ],
          D00 = LD[LD_off + N*(i+1) + i+1] / s,
          D11 = LD[LD_off + N* i    + i  ] / s,
          D10 =          1 / (D00*D11 - 1) / s; // <- TODO: understand this underflow-safe Vodoo magic (see [1])

      for( let j=0; j < O; j++ ) {
        const y0 = tmp[O* i   +j],
              y1 = tmp[O*(i+1)+j];
        tmp[O* i   +j] = D10 * (D00*y0 - y1);
        tmp[O*(i+1)+j] = D10 * (D11*y1 - y0);
      }

      ++i;
    }
    else
    { // 1x1 BLOCK
      for( let j=0; j < O; j++ )
        tmp[O*i+j] /= LD[LD_off + N*i+i];
    }
  }

  // BACKWARD SUBSTITUTION
  for( let k=M; k-- > 1; ) { const I = k - (P[P_off + k] < 0);
  for( let i=I; i-- > 0; )
  for( let j=O; j-- > 0; )
    tmp[O*i+j] -= LD[LD_off + N*k+i] * tmp[O*k+j]
  }

  // (UN)PERMUTE OUTPUT
  for( let i=0; i < M; i++ ) { let I = P[P_off + i]; I ^= -(I<0);
  for( let j=0; j < O; j++ )
    X[X_off + O*I+j] = tmp[O*i+j];
  }
}


export function pldlp_solve(LD,P,y)
{
  if( null == y ) {
    if( null == P ) throw new Error('Assertion failed.');
    y = P;
    [LD,P] = LD;
  }
  LD = asarray(LD);
   P = asarray( P);
   y = asarray( y);
  if(LD.ndim < 2 ) throw new Error('pldlp_solve(LD,P,y): LD must be at least 2D.');
  if( P.ndim < 1 ) throw new Error('pldlp_solve(LD,P,y): P must be at least 1D.');
  if( y.ndim < 2 ) throw new Error('pldlp_solve(LD,P,y): y must be at least 2D.');

  const N =LD.shape[LD.ndim-2],
        O = y.shape[ y.ndim-1];
  if( N !==LD.shape[LD.ndim-1] ) throw new Error('pldlp_solve(LD,P,y): LD must be square.');
  if( N !== P.shape[ P.ndim-1] ) throw new Error('pldlp_solve(LD,P,y): LD and P shape mismatch.');
  if( N !== y.shape[ y.ndim-2] ) throw new Error('pldlp_solve(LD,P,y): LD and y shape mismatch.');

  const
    ndim = Math.max(LD.ndim, P.ndim+1, y.ndim),
    shape = Int32Array.from({ length: ndim }, ()=>1);
  shape[ndim-2] = N;
  shape[ndim-1] = O;

  // FIND COMMON (BROADCASTED) SHAPE
  for( let shp of [
   LD.shape.slice(0,-2),
    P.shape.slice(0,-1),
    y.shape.slice(0,-2)
  ])
    for( let i=ndim-2, j=shp.length; i-- > 0 && j-- > 0; )
      if( 1 === shape[i] )
        shape[i] = shp[j];
      else if( shape[i] !== shp[j] && shp[j] !== 1 )
        throw new Error('Shapes are not broadcast-compatible.');

  // GENERATE RESULT DATA
  const         DType = [LD,P,y].every(x => x.dtype === 'float32') ? 'float32' : 'float64',
                DTypeArray = ARRAY_TYPES[DType],
      tmp = new DTypeArray(N*O),
    x_dat = new DTypeArray( shape.reduce((a,b) => a*b, 1) ),
   LD_dat =LD.data,
    P_dat = P.data,
    y_dat = y.data;
  let
   LD_off = 0,LD_stride = 1,
    P_off = 0, P_stride = 1,
    y_off = 0, y_stride = 1,
    x_off = 0;

  function solv(d) {
    if( d === ndim-2 ) {
     LD_stride = N*N;
      P_stride = N;
      y_stride = N*O;

      // COPY y
      for( let i=0; i < y_stride; i++ )
        x_dat[x_off+i] = y_dat[y_off+i];

      _pldlp_solve(N,N,O, LD_dat,LD_off, P_dat,P_off, x_dat,x_off, tmp);

      LD_off +=LD_stride;
       P_off += P_stride;
       y_off += y_stride;
       x_off += y_stride;

      return;
    }
    for( let l=shape[d]; ; l-- ) {
      solv(d+1);
      if( l===1 ) break;
      if( ! (LD.shape[ d - ndim +LD.ndim  ] > 1) ) LD_off -=LD_stride;
      if( ! ( P.shape[ d - ndim + P.ndim+1] > 1) )  P_off -= P_stride;
      if( ! ( y.shape[ d - ndim + y.ndim  ] > 1) )  y_off -= y_stride;
    }
    LD_stride *=LD.shape[ d - ndim +LD.ndim  ] || 1;
     P_stride *= P.shape[ d - ndim + P.ndim+1] || 1;
     y_stride *= y.shape[ d - ndim + y.ndim  ] || 1;
  }
  solv(0);

  return new NDArray(shape,x_dat);
}
