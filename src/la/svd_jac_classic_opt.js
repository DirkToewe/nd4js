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
import {qr_decomp} from './qr'
import {transpose_inplace} from './transpose_inplace'
import {_svd_jac_rot_rows,
        _svd_jac_rot_cols,
        _svd_jac_angles,
        _svd_jac_post } from './_svd_jac_utils'


function _find_branch( T, T_off, n )
{
  let i=0,
      j=0;
  if( n > 1 )
    [i,j] = _find_branch(T, T_off+2*n*(n+1), n+1 >>> 1 );

  let ij=0;
  T_off += 2*i*(i+1) + 4*j;

  for( let k=1; k < 4; k++ )
    if( T[T_off+k] > T[T_off+ij] )
      ij = k;

  return [
    (i<<1) + (ij>>>1),
    (j<<1) + (ij % 2)
  ];
}


export function svd_jac_classic_opt(A)
{
  A = asarray(A);
  if( A.dtype.startsWith('complex') )
    throw new Error('svd_jac_1sided(A): A.dtype must be float.');
  const
    shape = A.shape,
    N = shape[shape.length-2];

  // QR DECOMPOSE RECTANGULAR MATRICES AND WORK ON (QUADRATIC) R
  {
    const M = shape[shape.length-1];
    // if A is not square use QR Decomposition
    if( N > M ) {
      const [Q,R] = qr_decomp(A),
         [U,sv,V] = svd_jac_classic_opt(R)
      return [matmul2(Q,U), sv, V]
    }
    if( N < M ) {
      const [Q,R] = qr_decomp(A.T),
         [U,sv,V] = svd_jac_classic_opt(R)
      transpose_inplace(V)
      return [V, sv, matmul2(Q,U).T]
    }
  }
  // ALLOCATE RESULT DATA
  const DType = A.dtype==='float32' ? 'float32' : 'float64',
        DTypeArray = ARRAY_TYPES[DType],
        TOL = (N * eps(DType))**2,
        U =     DTypeArray.from(A.data); A = undefined; // <- potentially allow GC
  const S = new DTypeArray(N*N), // <- tempory storage for decomposition
        V = new DTypeArray(U.length),
       sv = new DTypeArray(U.length/N),
      ord = Int32Array.from({length: N}, (_,i) => i);

  if( 1 >  N ) throw new Error('Assertion failed.');
  if( 1 == N ) {
    for( let i=U.length; i-- > 0; )
      if( U[i] < +0.0 ) {
          U[i] *= -1.0;
         sv[i]  = -1.0;
      }
      else sv[i] = +1.0;
    return [
      new NDArray(shape,            sv ),
      new NDArray(shape.slice(0,-1), U ),
      new NDArray(shape,      V.fill(1))
    ];
  }

  // stopping criterion inspiredy by:
  //  "Jacobi's Method is More Accurate than QR"
  //   by James Demmel
  //   SIAM J. Matrix Anal. Appl, vol. 13, pp. 1204-1245, 1992
  const piv = (i,j) => {
    const  S_ij = S[N*i+j],
           S_ji = S[N*j+i];
    return i < N && j < N ? S_ij*S_ij + S_ji*S_ji : 0;
  };

  // TRIANGULAR TREE
  const T = new DTypeArray(
    function(){
      for( let len = 0,
                 n = N+1 >>> 1;; )
      {          n = n+1 >>> 1;
        len += 2*n*(n+1);
        if( 1 >= n ) return len;
      }
    }()
  );

  const update = (x,y) => {
    const idx = Int32Array.of(x,y);

    // UPDATE LEAVES
    for( let s=idx.length; s-- > 0; )
    {
      const k = (idx[s] >>>= 1) << 1; // <- round down to multiple of 2
      for( let l=0; l < N; l += 2 )
      {
        let max = Math.max(
          piv(k,1+l),
          piv(k+1,l)
        );
        if( k !== l ) {
          max = Math.max(max,  piv(k  ,l  ));
          max = Math.max(max,  piv(k+1,l+1));
        }
        let i = Math.max(k,l) >>> 1,
            j = Math.min(k,l) >>> 1,
            O = 2*(i%2) + (j%2);
        i >>>= 1;
        j >>>= 1;
        T[O + 2*i*(i+1) + 4*j] = max;
      }
    }

    let OFF = 0;

    for( let  n = N+1 >>> 1;; )
    {         n = n+1 >>> 1;
      if(1 >= n) break;
      const off = OFF + 2*n*(n+1);

      for( let s=idx.length; s-- > 0; )
      {
        const p = (idx[s] >>>= 1);

        for( let q=0; q < n; q++ )
        {
          let   i = Math.max(p,q),
                j = Math.min(p,q);
          const O = OFF + 2*i*(i+1) + 4*j;
          let   o = off + 2*(i%2) + (j%2);
          i >>>= 1;
          j >>>= 1;
          o += 2*i*(i+1) + 4*j;

          T[o] = T[O];
              
          for( let k=1; k < 4; k++ )
            T[o] = Math.max(T[o], T[O+k]);
        }
      }

      OFF = off;
    }
  };

  const find_piv = () => {
    let [i,j] = _find_branch(T,0, N+3 >>> 2);

    i <<= 1;
    j <<= 1;

    let p = -1,
        q = -1,
       mx = -Infinity;

    for( let k=i; k < i+2 && k < N; k++ )
    for( let l=j; l < j+2 && l < k; l++ )
    {
      let piv_kl = piv(k,l);
      if( piv_kl > mx ) {
        q = k,
        p = l, mx = piv_kl;
      }
    }

    return [p,q];
  };

  const debug_print = () => {
    console.log('\n')
    console.log('DEBUG')
    console.log('-----')
    console.log([...T].map(f => f.toFixed(3)))

    console.log('')
    console.log('PIVOTS')
    console.log('------')
    for( let i=0; i < N; i++ )
    {
      let row = [];
      for( let j=0; j < i; j++ )
        row.push( `${ piv(i,j).toFixed(3) }`.padEnd(8) );
      console.log( row.join('') );
    }

    console.log('')
    console.log('TREE')
    console.log('----')
    for( let off = 0,
               n = N+1 >>> 1;; off += 2*n*(n+1) )
    {          n = n+1 >>> 1;
      for( let r=0; r < n; r++ )
      for( const dr of [0,2] )
      {
        const row = [];
        for( let c=0; c <= r; c++ )
        for( const dc of [0,1] )
        {
          const t = T[off + 2*r*(r+1) + 4*c + dr + dc];
          row.push(`${t.toFixed(3)}`.replace('Infinity','∞').padEnd(8));
        }
        console.log(row.join(''))
      }

      console.log('-----');

      if( n <= 1 ) break;
    }
  };

  for( let UV_off=0,
           sv_off=0; sv_off < sv.length; UV_off += N*N,
                                         sv_off += N )
  {
    // MOVE FROM U TO S
    for( let i=0; i < N; i++ )
    for( let j=0; j < N; j++ ) { S[N*i+j] = U[UV_off + N*i+j];
                                            U[UV_off + N*i+j] = i != j ? 0 : 1 };
    // INIT V TO IDENTITY
    for( let i=0; i < N; i++ )
    for( let j=0; j < N; j++ ) V[UV_off + N*i+j] = i != j ? 0 : 1;

    // INIT TRIANGLE TREE
    for( let i=0; i < N; i += 4 ) update(i, Math.min(i+2,N-1)); // <- FIXME this can be done in half as many steps

     //
    // (CLASSICAL) JACOBI SVD ITERATIONS
   //
    for(;;)
    {
      // FIND THE OFF DIAGONAL PAIR WITH THE LARGEST HYPOTHENUSE
      const [p,q] = find_piv();

      // CHECK THAT THIS IS TRULY THE MAXIMUM (TODO: COMMENT OUT)
      if( Math.random()*N < 16 )
      for( let i=1; i < N; i++ )
      for( let j=0; j < i; j++ )
        if( ! (piv(i,j) <= piv(q,p)) ) {
          debug_print();
          throw new Error(`Assertion failed: ${p},${q} > ${i},${j}.`);
        }

      const
        S_pp = S[N*p+p],
        S_pq = S[N*p+q],
        S_qp = S[N*q+p],
        S_qq = S[N*q+q];

      if( ! ( S_pq*S_pq + S_qp*S_qp > Math.abs(S_pp*S_qq) * TOL ) )
        break;
 
      const [cα,sα,cβ,sβ] = _svd_jac_angles(S_pp, S_pq,
                                            S_qp, S_qq);

      // ROTATE S
      _svd_jac_rot_rows(S, N, N*p,N*q, cα,sα);
      _svd_jac_rot_cols(S, N,   p,  q, cβ,sβ);

      // ENTRIES (k,l) AND (l,k) ARE REMAINDERS (CANCELLATION ERROR) FROM ELIMINATION => SHOULD BE SAFELY ZEROABLE
      S[N*p+q] = 0.0;
      S[N*q+p] = 0.0;

      // UPDATE TRIANGLE TREE ROWS AND COLUMNS
      update(p,q);

      // ROTATE U & V
      _svd_jac_rot_rows(U, N, UV_off + N*p,
                              UV_off + N*q, cα, sα);
      _svd_jac_rot_rows(V, N, UV_off + N*p,
                              UV_off + N*q, cβ,-sβ);
    }

    _svd_jac_post( N, U,S,V, UV_off, sv, sv_off, ord );
  }

  return [
    new NDArray(shape,             U),
    new NDArray(shape.slice(0,-1),sv),
    new NDArray(shape,             V)
  ]
}
