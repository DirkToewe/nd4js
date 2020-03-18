'use strict';

/* This file is part of ND4JS.
 *
 * ND4JS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ND4JS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ND4JS. If not, see <http://www.gnu.org/licenses/>.
 */

import {_cholesky_decomp} from '../la/cholesky'
import {_tril_solve} from '../la/tri'

// References
// ----------
// .. [1] "REPRESENTATIONS OF QUASI-NEWTON MATRICES AND THEIR USE IN LIMITED MEMORY METHODS"
//         by Richard H. Byrd, Jorge Nocedal and Robert B. Schnabel
//         Technical Report NAM-03, June 1992; revised January 21, 1996
//         http://users.iems.northwestern.edu/~nocedal/PDFfiles/representations.pdf
//
// .. [2] "A LIMITED MEMORY ALGORITHM FOR BOUND CONSTRAINED OPTIMIZATION"
//         by Richard H. Byrd, Peihuang Lu, Jorge Nocedal and Ciyou Zhu
//         Technical Report NAM-08; Revised: May 1994
//         http://users.iems.northwestern.edu/~nocedal/PDFfiles/limited.pdf


export class LBFGSSolver
{
  constructor( M, N )
  {
    if( M%1 !== 0 ) throw new Error('Assertion failed.');
    if( N%1 !== 0 ) throw new Error('Assertion failed.');
    if( ! (0 < M) ) throw new Error('Assertion failed.');
    if( ! (0 < N) ) throw new Error('Assertion failed.');

    this.m =   0;
    this.M = M|0;
    this.N = N|0;

    this._scale = 1;

    this.dX  = new Float64Array(M*N);
    this.dG  = new Float64Array(M*N);

    this.dXdX= new Float64Array(M*M);
    this.dXdG= new Float64Array(M*M);
    this.dGdG= new Float64Array(M*M);

    this.D   = new Float64Array(M);

    this.J   = new Float64Array(M*M);
  }


  _computeJ()
  {
    const {m,M,N, dXdX,
                  dXdG, J, scale} = this;

    for( let i=0; i <  m; i++ )
    for( let j=0; j <= i; j++ )
      J[M*i+j] = scale * dXdX[M*i+j];

    for( let i=0; i <  m; i++ )
    for( let j=0; j <= i; j++ )
    for( let k=0; k <  j; k++ )
      J[M*i+j] += dXdG[M*i+k] / dXdG[M*k+k] * dXdG[M*j+k];

    _cholesky_decomp(m,M, J,0);
  }


  get scale() {
    return this._scale;
  }


  set scale( s )
  {
    s *= 1;
    if( isNaN(s) ) throw new Error('Assertion failed.');
    if( !(0 < s) ) throw new Error('Assertion failed.');

    if( this._scale!==s ) {
        this._scale = s;
        this._computeJ();
    }
  }


  update( dx, dg )
  {
    const {M,N, dX,
                dG, dXdX,
                    dXdG,
                    dGdG, D} = this;

    if( N !== dx.length ) throw new Error('Assertion failed.');
    if( N !== dg.length ) throw new Error('Assertion failed.');

    if( this.m >= M )
    {
      if( this.m > M ) throw new Error('Assertion failed.');

      // MAKE SPACE IN LAST ROW BY OVERWRITING 1ST ROW
      for( const A of [dX,dG] )
        for( let i=1; i < M; i++ )
        for( let j=0; j < N; j++ )
          A[N*(i-1)+j] = A[N*i+j];

      for( const A of [dXdX,dXdG,dGdG] )
        for( let i=0; ++i < M; )
        for( let j=0; ++j < M; )
          A[M*(i-1)+(j-1)] = A[M*i+j];

      for( let i=0; ++i < M; )
        D[i-1] = D[i];
    }
    else
      this.m++;

    const {m} = this,
           i = m-1;

    for( let j=0; j < N; j++ ) dX[N*i+j] = dx[j];
    for( let j=0; j < N; j++ ) dG[N*i+j] = dg[j];

    // UPDATE dX @ dX.T
    for( let j=0; j < m; j++ )
    { let                        dxdx = 0;
      for( let k=0; k < N; k++ ) dxdx += dX[N*i+k]*dX[N*j+k];
      dXdX[M*i+j] =
      dXdX[M*j+i] = dxdx;
    }

    // UPDATE dX @ dG.T
    for( let j=0; j < m; j++ )
    {
      dXdG[M*i+j] = 0; for( let k=0; k < N; k++ ) dXdG[M*i+j] += dX[N*i+k]*dG[N*j+k];
      dXdG[M*j+i] = 0; for( let k=0; k < N; k++ ) dXdG[M*j+i] += dX[N*j+k]*dG[N*i+k];
    }

    // UPDATE sqrt( diag(dX @ dG.T) )
    D[i] = Math.sqrt( dXdG[M*i+i] );

    // UPDATE dG @ dG.T
    for( let j=0; j < m; j++ )
    { let                        dgdg = 0;
      for( let k=0; k < N; k++ ) dgdg += dG[N*i+k]*dG[N*j+k];
      dGdG[M*i+j] =
      dGdG[M*j+i] = dgdg;
    }

    // UPDATE J
    this._computeJ();
  }


  /** This method "sort of" computes half an inner product over B.
   *  To compute an inner product u.T @ B @ v, You can use:
   *
   *    compute_ub(u, ub);
   *    compute_ub(v, vb);
   *    compute_ubbv(ub, dot(u,v), bv)
   */
  compute_bv( v, bv )
  { //             ┏                 ┓
    //             ┃        ┊        ┃
    //             ┃        ┊        ┃ ┏                 ┓ -1  ┏                 ┓   ┏                 ┓ -1  ┏                         ┓
    //             ┃        ┊        ┃ ┃        ┊        ┃     ┃        ┊        ┃   ┃        ┊        ┃     ┃                         ┃
    //             ┃        ┊        ┃ ┃   √D   ┊-√D \ Lᵀ┃     ┃   -I   ┊    0   ┃   ┃   √D   ┊   0    ┃     ┃           dG            ┃
    //             ┃        ┊        ┃ ┃        ┊        ┃     ┃        ┊        ┃   ┃        ┊        ┃     ┃                         ┃
    // uᵀ⋅B⋅v = uᵀ⋅┃   dGᵀ  ┊ B0⋅dXᵀ ┃⋅┃┈┈┈┈┈┈┈┈┼┈┈┈┈┈┈┈┈┃  ⋅  ┃┈┈┈┈┈┈┈┈┼┈┈┈┈┈┈┈┈┃ ⋅ ┃┈┈┈┈┈┈┈┈┼┈┈┈┈┈┈┈┈┃  ⋅  ┃┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┃⋅v
    //             ┃        ┊        ┃ ┃        ┊        ┃     ┃        ┊        ┃   ┃        ┊        ┃     ┃                         ┃
    //             ┃        ┊        ┃ ┃    0   ┊   Jᵀ   ┃     ┃    0   ┊    I   ┃   ┃-L / √D ┊   J    ┃     ┃          dX⋅B0          ┃
    //             ┃        ┊        ┃ ┃        ┊        ┃     ┃        ┊        ┃   ┃        ┊        ┃     ┃                         ┃
    //             ┃        ┊        ┃ ┗                 ┛     ┗                 ┛   ┗                 ┛     ┗                         ┛
    //             ┃        ┊        ┃
    //             ┗                 ┛
    //                 ┏                 ┓                  ┏                 ┓ -1  ┏                      ┓
    //                 ┃        ┊        ┃                  ┃        ┊        ┃     ┃                      ┃
    //                 ┃   -I   ┊    0   ┃                  ┃   √D   ┊   0    ┃     ┃         dG           ┃
    //                 ┃        ┊        ┃                  ┃        ┊        ┃     ┃                      ┃
    //        =: uᵀ⋅bᵀ⋅┃┈┈┈┈┈┈┈┈┼┈┈┈┈┈┈┈┈┃⋅b⋅v    =>    b = ┃┈┈┈┈┈┈┈┈┼┈┈┈┈┈┈┈┈┃  ⋅  ┃┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┃
    //                 ┃        ┊        ┃                  ┃        ┊        ┃     ┃                      ┃
    //                 ┃    0   ┊    I   ┃                  ┃-L / √D ┊   J    ┃     ┃        dX⋅B0         ┃
    //                 ┃        ┊        ┃                  ┃        ┊        ┃     ┃                      ┃
    //                 ┗                 ┛                  ┗                 ┛     ┗                      ┛
    // Where:
    //   L = tril(dX⋅dGᵀ, -1)
    //   D = diag(dX⋅dGᵀ)
    //   I = eye(m)
    //   J = cholesky( dX⋅B0⋅dXᵀ + (L/D)⋅Lᵀ )
    //  B0 = I⋅scale
    const {m,M,N, dX,
                  dG, dXdG, J, D, scale} = this;

    if( v.length !==  N ) throw new Error('Assertion failed.');
    if(bv.length !==2*M ) throw new Error('Assertion failed.');

    bv.fill(0.0, 0,2*m);

    for( let i=0; i < m; i++ )
    for( let j=0; j < N; j++ )
      bv[i] += dG[N*i+j] * v[j];

    for( let i=1; i < m; i++ )
    for( let j=0; j < i; j++ )
      bv[m+i] += dXdG[M*i+j] * bv[j] / dXdG[M*j+j];

    for( let i=0; i < m; i++ )
      bv[i] /= D[i];

    for( let i=0; i < m; i++ )
    for( let j=0; j < N; j++ )
      bv[m+i] += dX[N*i+j] * v[j] * scale;

    _tril_solve(m,M,1, J,0, bv,m);
  }


  compute_be( j, hej )
  {
    const {m,M,N, dX,
                  dG, dXdG, J, D, scale} = this;

    if( j%1 !== 0 ) throw new Error('Assertion failed.');

    j |= 0;

    if( !(j >= 0) ) throw new Error('Assertion failed.');
    if( !(j <  N) ) throw new Error('Assertion failed.');
    if( hej.length !== 2*M ) throw new Error('Assertion failed.');

    hej.fill(0.0, 0,2*m);

    for( let i=0; i < m; i++ )
      hej[i] += dG[N*i+j];

    for( let i=1; i < m; i++ )
    for( let j=0; j < i; j++ )
      hej[m+i] += dXdG[M*i+j] * hej[j] / dXdG[M*j+j];

    for( let i=0; i < m; i++ )
      hej[i] /= D[i];

    for( let i=0; i < m; i++ )
      hej[m+i] += dX[N*i+j] * scale;

    _tril_solve(m,M,1, J,0, hej,m);
  }


  compute_ubbv( hu, uv, hv )
  {
    const {m,M, scale} = this;

    uv *= scale;
    if( isNaN(uv) ) throw new Error('Assertion failed.');

    if( hu.length !== 2*M ) throw new Error('Assertion failed.');
    if( hv.length !== 2*M ) throw new Error('Assertion failed.');

    let                        result  = 0;
    for( let i=2*m; i-- > m; ) result -= hu[i] * hv[i];
    for( let i=  m; i-- > 0; ) result += hu[i] * hv[i];
    return                     result +  uv;
  }
}
