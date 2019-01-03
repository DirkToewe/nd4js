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
import {rrqr_decomp,
       _rrqr_rank} from './rrqr'
import {root1d_bisect} from '../opt/root1d_bisect'


export function svd_jac_1sided(A)
{
  // SEE:
  //   http://www.netlib.org/lapack/lawnspdf/lawn169.pdf
  //   http://www.netlib.org/lapack/lawnspdf/lawn170.pdf
  A = asarray(A);
  const
     U_shape =                 A.shape,
     V_shape = Int32Array.from(U_shape),
     sv_shape=                 U_shape.slice(0,-1),
    [N,M]    =                 U_shape.slice(  -2);

  if( N < M ) {
    const [U,sv,V] = svd_jac_1sided(A.T);
    // TODO: transpose V in-place
    return [V.T, sv, U.T];
  }

  const
    DTypeArray = ARRAY_TYPES[A.dtype==='float32' ? 'float32' : 'float64'],
    TOL = Number.EPSILON * N,// <- FIXME what about float32 ?
    sqrt2 = Math.sqrt(2);
  sv_shape[sv_shape.length-1] = M;
   V_shape[ V_shape.length-2] = M;

  // TODO: check if QR preconditioning is sensible for N==M
  let [Q,R,P] = rrqr_decomp( N >= M ? A : A.T );
  A = undefined;

  Q = Q.data;
  R = R.data;
  P = P.data;

  const
    U  = new DTypeArray(M*M), // <- temporarily stores U
    sv = new DTypeArray(Q.length/N),
    tmp= new DTypeArray(M) // <- temporary storage for matrix multiplication

  function* indices_241(M) {
    // SEE: W. F. Mascarenhas "On the Convergence of the Jacobi Method for Arbitrary Orderings."
    // ┏    ╷    ╷         ┓
    // ┃ B1 ┊ B3 ┊         ┃
    // ┃┄┄┄┄┼┄┄┄┄┤   B7    ┃
    // ┃    ┊ B2 ┊         ┃
    // ┃    └┄┄┄┄┼┄┄┄┄┬┄┄┄┄┃ = RRᵀ
    // ┃         ┊ B4 ┊ B6 ┃
    // ┃         └┄┄┄┄┼┄┄┄┄┃
    // ┃              ┊ B5 ┃
    // ┗              ╵    ┛
    // ORDER: {1,2,3,4,5,6,1,2,4,5,7}
    const
      m = M >>> 2,
      IJ = new Int32Array(2);
    for( let k=0; k < 2; k++ )
    {
      // B1
      for( IJ[0]=0      ; IJ[0] <   m; IJ[0]++ )
      for( IJ[1]=1+IJ[0]; IJ[1] <   m; IJ[1]++ ) yield IJ;
      // B2
      for( IJ[0]=  m    ; IJ[0] < 2*m; IJ[0]++ )
      for( IJ[1]=1+IJ[0]; IJ[1] < 2*m; IJ[1]++ ) yield IJ;
      // B3
      if( 0==k )
      for( IJ[0]=0      ; IJ[0] <   m; IJ[0]++ )
      for( IJ[1]=  m    ; IJ[1] < 2*m; IJ[1]++ ) yield IJ;
      // B4
      for( IJ[0]=2*m    ; IJ[0] < 3*m; IJ[0]++ )
      for( IJ[1]=1+IJ[0]; IJ[1] < 3*m; IJ[1]++ ) yield IJ;
      // B5
      for( IJ[0]=3*m    ; IJ[0] <   M; IJ[0]++ )
      for( IJ[1]=1+IJ[0]; IJ[1] <   M; IJ[1]++ ) yield IJ;
      // B6
      if( 0==k )
      for( IJ[0]=2*m    ; IJ[0] < 3*m; IJ[0]++ )
      for( IJ[1]=3*m    ; IJ[1] <   M; IJ[1]++ ) yield IJ;
    }
    // B7
    for( IJ[0]=  0; IJ[0] < 2*m; IJ[0]++ )
    for( IJ[1]=2*m; IJ[1] <   M; IJ[1]++ ) yield IJ;
  }

  for( let Q_off    = 0,
           R_off    = 0,
           P_sv_off = 0; R_off < R.length;
           Q_off    += N*M,
           R_off    += M*M,
           P_sv_off += M
  )
  {
    // [STEP 1] determine numerical rank of R
    const RANK = _rrqr_rank(N,M, R,R_off, tmp)

    // [STEP 2] eliminate the right-most, rank-deficient columns of R, then only perform the Jacobi SVD on the upper left square part
    for( let i=M   ; i-- > RANK; ) { // <- columns that are to be eliminated
    for( let j=RANK; j-- > 0   ; ) {
    for( let k=RANK; --k > j   ; ) {
        const
          s    = U[        M*i+k],
          c    = R[R_off + M*i+k],
          R_ji = R[R_off + M*j+i],
          R_jk = R[R_off + M*j+k];
        R[R_off + M*j+i] = c*R_ji - s*R_jk;
        R[R_off + M*j+k] = s*R_ji + c*R_jk;
    }
      const
        R_ji = R[R_off + M*j+i],
        R_jj = R[R_off + M*j+j],
        norm = Math.hypot(R_ji,R_jj),
        c = R_jj / norm,
        s = R_ji / norm;
      R[R_off + M*j+i] = c*R_ji - s*R_jj;
      R[R_off + M*j+j] = s*R_ji + c*R_jj;
      R[R_off + M*i+j] = c;
      U[        M*i+j] = s;
      if( ! (Math.abs(R[R_off + M*j+i]) < 1e-8) ) throw new Error('Assertion failed.');
      R[R_off + M*j+i] = 0;
    }}

    // [STEP 3]COPY R -> U
    for( let i=0   ; i < RANK; i++ )
    for( let j=0   ; j < RANK; j++ ) U[M*i+j] = R[R_off + M*i+j];

    // [STEP 4] PERFORM JACOBI SWEEPS
    let sweeps = 0;
    for( let maxErr=Infinity; maxErr > TOL*TOL; )
    {
      if( ++sweeps > 128 ) throw new Error('Too many iterations.')
      maxErr=0.0;
      for( const [i,j] of indices_241(RANK) ) // "2.41"-order with locally super-quadratic convergence (in theory)
//      for( let i=0  ; i < RANK; i++ ) // ◀─┬─ alternative, simpler (De Rijk) order
//      for( let j=1+i; j < RANK; j++ ) // ◀─╯
      {
        if( j <= i ) throw new Error('Assertion failed');
        let R_ii = 0.0, // ◀─╮
            R_ij = 0.0, // ◀─┼─ TODO: make this underflow-saf(ish)er
            R_jj = 0.0; // ◀─╯ 
        { // compute dot products
          const    IK = R_off + M*i + RANK;
          for( let ik = R_off + M*i,
                   jk = R_off + M*j; ik < IK; )
          {
            const R_ik = R[ik++]; R_ii += R_ik*R_ik;
            const R_jk = R[jk++]; R_ij += R_ik*R_jk;
                                  R_jj += R_jk*R_jk;
          }
        }
        if( R_ii <= 0 ) throw new Error('Unexpected underflow.');
        if( R_jj <= 0 ) throw new Error('Unexpected underflow.');
        maxErr = Math.max( maxErr, (R_ij / R_ii) * (R_ij / R_jj) );// <- FIXME: check whether underflow makes a problem
        if( maxErr/TOL <= TOL ) continue;

        // determine rotation angle using binary search. TODO: use sth. faster than binary search
        let
          s = root1d_bisect(
            s => {
              const cs = Math.sqrt( (1+s)*(1-s) ) * s;
              return cs*R_ii + R_ij*(1 + sqrt2*s)*(1 - sqrt2*s) - cs*R_jj; // <- theoretically better approximation (underflow-safe?)
            },
            /*s_min=*/-sqrt2/2,
            /*s_max=*/+sqrt2/2
          ),
          c = Math.sqrt( (1+s)*(1-s) );
        

        // CHOOSE (ANGLE) SOLUTION THAT KEEPS MAIN DIAGONAL IN ORDER, i.e.:
        // c*(c*R_ii - s*R_ij) - s*(c*R_ij - s*R_jj) = r_ii < r_jj = s*(s*R_ii + c*R_ij) + c*(s*R_ij + c*R_jj)
        if( R_ii*(1 + sqrt2*s)*(1 - sqrt2*s) < R_jj + 2*c*s*R_ij ) { const C = c; c = -s; s = C; }

        { // perform rotation
          const    ik0= R_off + M*i;
          for( let ik = R_off + M*i+RANK,
                   jk = R_off + M*j+RANK; ik > ik0; )
          { const R_ik = R[--ik],
                  R_jk = R[--jk];
            R[ik] = c*R_ik - s*R_jk;
            R[jk] = c*R_jk + s*R_ik;
          }
        }
      }
    }

    // [STEP 5] COMPUTE SINGULAR VALUES AND NORMALIZE ROWS OF V
    // TODO: V might also have to be reorthogonalized
    for( let i=0; i < RANK; i++ )
    {
      let sum = 0.0,
          max = 0.0;
      for( let j=0; j < RANK; j++ )
      {
        const R_ij = Math.abs(R[R_off + M*i+j]);
        if(   R_ij != 0 ) { // <- handles NaN (by making the result NaN)
          if( R_ij > max ) {
            const scale = max / R_ij; max = R_ij;
            sum *= scale*scale;
          }
          const ratio = R_ij / max;
          sum += ratio*ratio;
        }
      }
      const norm = isFinite(max) ? Math.sqrt(sum)*max : max;
      if( norm == 0 ) throw new Error('Unhandled underflow. Ask the developer nicely and bribe him with a couple of beers, then he may address this issue.');
      for( let j=RANK; j-- > 0; ) R[R_off + M*i+j] /= norm;
      sv[P_sv_off + i] = norm;
    }

    // [STEP 6] COMPUTE U ( Using U = R @ inv(V) @ inv(SV) = R @ V.T @ diag(1/sv) )
    //   ( should be numerically stable since V is orthogonal, i.e. well-conditioned ... right? )
    for( let i=0; i < RANK; i++ ) {
      for( let j=0; j < RANK; j++ ) { tmp[j] = U[M*i+j]; U[M*i+j]=0; }

      for( let j=0; j < RANK; j++ ) {
        for( let k=0; k < RANK; k++ )
          U[M*i+j] += R[R_off + M*j+k]*tmp[k];
        U[M*i+j] /= sv[P_sv_off + j];
      }
    }

    // [STEP 7] UNDO [STEP 2] 
    for( let i=RANK; i < M   ; i++ ) {
    for( let j=0   ; j < i   ; j++ ) {
    for( let k=0   ; k < RANK; k++ ) {
        const
          s    = U[        M*i+k],
          c    = R[R_off + M*i+k],
          R_ji = R[R_off + M*j+i],
          R_jk = R[R_off + M*j+k];
        R[R_off + M*j+i] =  c*R_ji + s*R_jk;
        R[R_off + M*j+k] = -s*R_ji + c*R_jk;
    }}
      // LAST REMAINING ROW
      let R_ii = 1.0;
      for( let k=0; k < RANK; k++ ) {
        const
          s = U[        M*i+k],
          c = R[R_off + M*i+k];
              R[R_off + M*i+k] = -s*R_ii; R_ii *= c;
      }
      R[R_off + M*i+i] = R_ii;
    }

    // [STEP 8] MULTIPLY Q = Q @ U
    for( let i=0; i < N; i++ )
    {
      for( let j=0; j < RANK; j++ ) { tmp[j] = Q[Q_off + M*i+j]; Q[Q_off + M*i+j]=0; }

      for( let k=0; k < RANK; k++ )
      for( let j=0; j < RANK; j++ ) Q[Q_off + M*i+j] += tmp[k] * U[M*k+j];
    }

    // [STEP 9] APPLY P TO R (COLUMN PERMUTATIONS)
    for( let i=0; i < M; i++ )
    {
      for( let j=0; j   < M; j++ ) tmp[P[P_sv_off + j]] = R[R_off + M*i+j];
      for( let j=M; j-- > 0;     )                        R[R_off + M*i+j] = tmp[j];
    }

    // TODO SORT SV???
  }

  return [
    new NDArray( U_shape, Q ),
    new NDArray(sv_shape, sv),
    new NDArray( V_shape, R )
  ];
};


//la._svd_decomp_jac_1sided = A => {
//  // SEE:
//  //   http://www.netlib.org/lapack/lawnspdf/lawn169.pdf
//  //   http://www.netlib.org/lapack/lawnspdf/lawn170.pdf
//  A = nd.asarray(A);
//  const
//     U_shape =                 A.shape,
//     V_shape = Int32Array.from(U_shape),
//     sv_shape=                 U_shape.slice(0,-1),
//    [N,M]    =                 U_shape.slice(  -2);
//
//  if( N < M ) {
//    const [U,sv,V] = la._svd_decomp_jac_1sided(A.T);
//    // TODO: transpose V in-place
//    return [V.T,sv,U.T];
//  }
//
//  const
//    DTypeArray = nd.dtypes[nd.super_dtype(A.dtype,'float64')],
//    TOL = Number.EPSILON * N; // <- FIXME what about float32 ?
//  sv_shape[sv_shape.length-1] = M;
//   V_shape[ V_shape.length-2] = M;
//
//  // TODO: check if QR preconditioning is sensible for N==M
//  let [Q,R,P] = la.rrqr_decomp( N >= M ? A : A.T ); // <- TODO: determine rank and perform second QR decomp.
//  A = undefined;
//
//  Q = Q.data;
//  R = R.data;
//  P = P.data;
//
//  const
//    U  = new DTypeArray(M*M), // <- temporarily stores U
//    sv = new DTypeArray(Q.length/N),
//    tmp= new DTypeArray(M); // <- temporary storage for matrix multiplication
//
//  function* indices_241() {
//    // SEE: W. F. Mascarenhas "On the Convergence of the Jacobi Method for Arbitrary Orderings."
//    // ┏    ╷    ╷         ┓
//    // ┃ B1 ┊ B3 ┊         ┃
//    // ┃┄┄┄┄┼┄┄┄┄┤   B7    ┃
//    // ┃    ┊ B2 ┊         ┃
//    // ┃    └┄┄┄┄┼┄┄┄┄┬┄┄┄┄┃ = RRᵀ
//    // ┃         ┊ B4 ┊ B6 ┃
//    // ┃         └┄┄┄┄┼┄┄┄┄┃
//    // ┃              ┊ B5 ┃
//    // ┗              ╵    ┛
//    // ORDER: {1,2,3,4,5,6,1,2,4,5,7}
//    const
//      m = M >>> 2,
//      IJ = new Int32Array(2);
//    for( let k=0; k < 2; k++ )
//    {
//      // B1
//      for( IJ[0]=0      ; IJ[0] <   m; IJ[0]++ )
//      for( IJ[1]=1+IJ[0]; IJ[1] <   m; IJ[1]++ ) yield IJ;
//      // B2
//      for( IJ[0]=  m    ; IJ[0] < 2*m; IJ[0]++ )
//      for( IJ[1]=1+IJ[0]; IJ[1] < 2*m; IJ[1]++ ) yield IJ;
//      // B3
//      if( 0==k )
//      for( IJ[0]=0      ; IJ[0] <   m; IJ[0]++ )
//      for( IJ[1]=  m    ; IJ[1] < 2*m; IJ[1]++ ) yield IJ;
//      // B5
//      for( IJ[0]=3*m    ; IJ[0] <   M; IJ[0]++ )
//      for( IJ[1]=1+IJ[0]; IJ[1] <   M; IJ[1]++ ) yield IJ;
//      // B4
//      for( IJ[0]=2*m    ; IJ[0] < 3*m; IJ[0]++ )
//      for( IJ[1]=1+IJ[0]; IJ[1] < 3*m; IJ[1]++ ) yield IJ;
//      // B6
//      if( 0==k )
//      for( IJ[0]=2*m    ; IJ[0] < 3*m; IJ[0]++ )
//      for( IJ[1]=3*m    ; IJ[1] <   M; IJ[1]++ ) yield IJ;
//    }
//    // B7
//    for( IJ[0]=  0; IJ[0] < 2*m; IJ[0]++ )
//    for( IJ[1]=2*m; IJ[1] <   M; IJ[1]++ ) yield IJ;
//  }
//
//  for( let Q_off    = 0,
//           R_off    = 0,
//           P_sv_off = 0; R_off < R.length;
//           Q_off    += N*M,
//           R_off    += M*M,
//           P_sv_off += M
//  )
//  {
//    // determine the rank
//    let RANK = M;
//    while( R[R_off + (M+1)*(--RANK)] == 0 );
//    ++RANK;
//    console.log(`RANK: ${RANK}/${M}`);
//
//    // TODO eliminate the right-most, rank-deficient columns of R, then only perform the Jacobi SVD on the upper left square part
//
//    // COPY R -> U
//    for( let i=M*M; i-- > 0; ) U[i] = R[R_off + i];
//
//    // PERFORM JACOBI SWEEPS
//    let sweeps = 0;
//    for( let maxErr=Infinity; maxErr/TOL > TOL; )
//    {
//      if( ++sweeps > 128 ) throw new Error('Too many iterations.')
//      maxErr=0.0;
//      for( const [i,j] of indices_241() ) // "2.41"-order with locally super-quadratic convergence (in theory)
////      for( let i=0  ; i < M; i++ ) // ◀─┬─ alternative, simpler (De Rijk) order
////      for( let j=1+i; j < M; j++ ) // ◀─╯
//      {
//        if( j <= i ) throw new Error('Assertion failed');
//        let R_ii = 0.0,
//            R_ij = 0.0,
//            R_jj = 0.0;
//        { // compute dot products
//          const    IK = R_off + M*(i+1);
//          for( let ik = R_off + M* i,
//                   jk = R_off + M* j; ik < IK; )
//          {
//            const R_ik = R[ik++]; R_ii += R_ik*R_ik;
//            const R_jk = R[jk++]; R_ij += R_ik*R_jk;
//                                  R_jj += R_jk*R_jk;
//          }
//        }
//        if( R_ii <= 0 ) throw new Error('Assertion failed!');
//        if( R_jj <= 0 ) throw new Error('Assertion failed!');
//        maxErr = Math.max( maxErr, (R_ij / R_ii) * (R_ij / R_jj) );// <- FIXME: check whether underflow makes a problem
//
//        // determine rotation angle using binary search. TODO: use faster approximation than binary search
//        let c,s; {
//          let s_min = -Math.sqrt(2)/2,
//              s_max = +Math.sqrt(2)/2;
//          const
//            sign = R_ii > R_jj ? +1.0 : -1.0, // <- re-orient the function to be montonously increasing
//            F = s => {
//              const c = Math.sqrt( (1+s)*(1-s) );
////              return sign * ( s*c*(R_ii - R_jj) + 2*R_ij*(0.5 + s)*(0.5 - s) ); // <- theoretically better approximation (underflow-safe)
//              return sign * ( // <- approximates actually performed numerical operations
//                  s*(c*R_ii - s*R_ij)
//                + c*(c*R_ij - s*R_jj)
//              )
//            };
//          // begin binary search
//          for(;;) {
//            const s = (s_min + s_max)/2;
//            if( s_min == s ||
//                s_max == s ) break;
//            const f = F(s);
//            if( isNaN(f) ) throw new Error('NaN encountered.');
//            if( f <= 0 ) s_min = s;
//            if( f >= 0 ) s_max = s;
//          }
//          // choose between α_min and α_max
//          s = Math.abs( F(s_min) ) < Math.abs( F(s_max) ) ? s_min : s_max,
//          c = Math.sqrt( (1+s)*(1-s) );
//        }
//
//        { // always choose the rotation that keeps the main diagonal in order
//          const r_ii = c*(c*R_ii - s*R_ij) - s*(c*R_ij - s*R_jj),
//                r_jj = s*(s*R_ii + c*R_ij) + c*(s*R_ij + c*R_jj);
//          if( r_ii < r_jj ) { const C = c; c = -s; s = C; }
//        }
//
//        { // perform rotation
//          const    ik0= R_off + M* i;
//          for( let ik = R_off + M*(i+1),
//                   jk = R_off + M*(j+1); ik > ik0; )
//          {
//            const R_ik = R[--ik],
//                  R_jk = R[--jk];
//            R[ik] = c*R_ik - s*R_jk;
//            R[jk] = c*R_jk + s*R_ik;
//          }
//        }
//      }
//    }
////    console.log(`${sweeps} sweeps.`)
//
//    // COMPUTE SINGULAR VALUES AND NORMALIZE ROWS OF V
//    for( let i=0; i < M; i++ )
//    {
//      let sum = 0.0,
//          max = 0.0;
//      for( let j=0; j < M; j++ )
//      {
//        const R_ij = Math.abs(R[R_off + M*i+j]);
//        if(   R_ij != 0 ) { // <- handles NaN (by making the result NaN)
//          if( R_ij > max ) {
//            const scale = max / R_ij; max = R_ij;
//            sum *= scale*scale;
//          }
//          const ratio = R_ij / max;
//          sum += ratio*ratio;
//        }
//      }
//      const norm = isFinite(max) ? Math.sqrt(sum)*max : max;
//      for( let j=M; j-- > 0; ) R[R_off + M*i+j] /= norm;
//      sv[P_sv_off + i] = norm;
//    }
//
//    // COMPUTE U ( Using U = R @ inv(V) @ inv(SV) = R @ V.T @ diag(1/sv) )
//    //   ( should be numerically stable since V is orthogonal, i.e. well-conditioned ... right? )
//    for( let i=0; i < M; i++ ) {
//      for( let j=0; j < M; j++ ) { tmp[j] = U[M*i+j]; U[M*i+j]=0; }
//
//      for( let j=0; j < M; j++ ) {
//        for( let k=0; k < M; k++ )
//          U[M*i+j] += R[R_off + M*j+k]*tmp[k];
//        U[M*i+j] /= sv[P_sv_off + j];
//      }
//    }
//
//    // MULTIPLY Q = Q @ U
//    for( let i=0; i < N; i++ )
//    {
//      for( let j=0; j < M; j++ ) { tmp[j] = Q[Q_off + M*i+j]; Q[Q_off + M*i+j]=0; }
//
//      for( let k=0; k < M; k++ )
//      for( let j=0; j < M; j++ ) Q[Q_off + M*i+j] += tmp[k] * U[M*k+j];
//    }
//
//    // APPLY P TO R (COLUMN PERMUTATIONS)
//    for( let i=0; i < M; i++ )
//    {
//      for( let j=0; j   < M; j++ ) tmp[P[P_sv_off + j]] = R[R_off + M*i+j];
//      for( let j=M; j-- > 0;     )                        R[R_off + M*i+j] = tmp[j];
//    }
//  }
//
//  return [
//    new nd.Array( U_shape, Q ),
//    new nd.Array(sv_shape, sv),
//    new nd.Array( V_shape, R )
//  ];
//};
