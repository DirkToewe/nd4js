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
import {_ldl_decomp,
        _ldl_solve} from '../la/ldl'
import {_tril_solve} from '../la/tri'

import {_heap_sort} from './_opt_utils'

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

// Compares to [2], the following symbols have been renamed in the source:
//
//   Y -> dG
//   S -> dX
//
// dG is the delta in gradients between iterations. dX is the delta in X.
// The names dG and dX therefore seemed to be more intuitive.

/** Computes the dot product of two arrays.
 */
const dot = (u,v) => {
  if( u.length%1 !== 0       ) throw new Error('Assertion failed.');
  if( v.length%1 !== 0       ) throw new Error('Assertion failed.');
  if( u.length   !== v.length) throw new Error('Assertion failed.');

  let uv = 0;
  for( let i=u.length; i-- > 0; )
    uv += u[i]*v[i];

  return uv;
};


export class L_BFGS_B_Solver
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

    this.travels = new Float64Array(N);
    this.Bdx = new Float64Array(2*M);
    this.Bg  = new Float64Array(2*M);
    this.Bei = new Float64Array(2*M);

    this.M_WAAW = new Float64Array(4*M*M);
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
    const {M,N, dX, dXdX,
                dG, dXdG,
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


  compute_cauchyGeneralized( G, X, bounds, indices )
  {
    const {m,N, Bg,Bdx,Bei, travels} = this;

    if(      G.length !==  N) throw new Error('Assertion failed.');
    if(      X.length !==  N) throw new Error('Assertion failed.');
//    if(     dX.length !==  N) throw new Error('Assertion failed.');
    if( bounds.length !==2*N) throw new Error('Assertion failed.');
    if(indices.length !==  N) throw new Error('Assertion failed.');

/*DEBUG*/    if( ! indices.slice().sort((i,j) => i-j).every((i,I) => i===I) )
/*DEBUG*/      throw new Error('Assertion failed.');

    for( let i=N; i-- > 0; ) {
      if( !(X[i] >= bounds[2*i+0]) ) throw new Error('Assertion failed.');
      if( !(X[i] <= bounds[2*i+1]) ) throw new Error('Assertion failed.');
    }

    for( let i=N; i-- > 0; )
      if( 0===G[i] )
        travels[i] = Infinity;
      else {
        const                           end = bounds[2*i + (G[i] < 0)];
                   travels[i] = (X[i] - end) / G[i];
        if( !(0 <= travels[i]) )
          throw new Error('Assertion failed: ' + travels[i]);
      }

    let   t = 0,
      n_fix = 0;

    this.compute_bv(G, Bg);
    Bdx.fill(0.0, 0,2*m);

    let fp  = dot(G,G),
        fpp = this.compute_ubbv(Bg, fp, Bg);

    // Make sure zero gradient entries are the very last
    // TODO: entries with travels[i]===0 could be skipped before _heap_sort()
    const order = _heap_sort(indices, (i,j) => 0 !== G[j] && travels[i] < travels[j]);

    loop:for( const i of order )
    {
      if( 0 === G[i] ) {
        for( let j=n_fix; j < N; j++ ) {
          const                            k = indices[j];
          if( !(Number.MAX_VALUE < travels[k]) ) throw new Error('Assertion failed.');
          if(              0 !== G[travels[k]] ) throw new Error('Assertion failed.');
        }
        break loop;
      }

      const end = bounds[2*i + (G[i] < 0)],
             dt = travels[i] - t;

      if( ! (0 <= dt) )
        throw new Error('Assertion failed.');

      if( 0===dt ) { // <- TODO: entries with travels[i]===0 could be skipped before this loop and _heap_sort()
        X[i] = end;
        G[i] = 0;
        t = travels[i];
        ++n_fix;
        continue loop;
      }

      if( !(0 <= fpp) ) throw new Error('Assertion failed.'); // <- B is supposed to be positive definite

      const cp = fp / fpp;

      if( cp < dt ) {
        t += Math.max(0,cp);
        break loop;
      }

      t = travels[i];

      for( let j=2*m; j-- > 0; )
        Bdx[j] -= dt * Bg[j];

      this.compute_be(i, Bei);

      const eBx = this.compute_ubbv(Bei,-t*G[i], Bdx),
            eBe = this.compute_ubbv(Bei,  1    , Bei),
            eBg = this.compute_ubbv(Bei,   G[i], Bg );

      fp  -=  dt * fpp  +  G[i] * (G[i] + eBx);
      fpp -=  G[i] * (2*eBg - eBe*G[i]);

      for( let j=2*m; j-- > 0; )
        Bg[j] -= G[i] * Bei[j];

      X[i] = end;
      G[i] = 0;
      ++n_fix;
    }

    for( let  j=n_fix; j < N; j++ )
    { const   i = indices[j];
            X[i] -= t*G[i];
      if( !(X[i] >= bounds[2*i+0]) ) { if( !(X[i] > bounds[2*i+0]-1e-8) ) throw new Error('Assertion failed.'); X[i] = bounds[2*i+0]; }
      if( !(X[i] <= bounds[2*i+1]) ) { if( !(X[i] < bounds[2*i+1]+1e-8) ) throw new Error('Assertion failed.'); X[i] = bounds[2*i+1]; }
    }

    return n_fix;
  }


  compute_subspace_Hv( v, Hv, indices,n_fix )
  {
    // For more information, see [2].
    // The goal is to solve the following subspace_problem:
    //        _
    // Hv = Z⋅B⋅Zᵀ⋅v
    //       
    // Where v and Hv are the input and output parameters of
    // of compute_subspace_Hv (see above).
    // _
    // B is the subspace Hessian (approximation):
    // _
    // B = Zᵀ⋅B⋅Z
    //                _
    // The inverse of B can be written as:
    // _                     _           _       _
    // B⁻¹ = I/scale + Zᵀ⋅W⋅(M⁻¹ + scale⋅Wᵀ⋅A⋅Aᵀ⋅W)⋅Wᵀ⋅Z
    //
    //       ┏                                  ┓
    //       ┃                ┊                 ┃
    //       ┃                ┊                 ┃
    //       ┃-(D+dGᵀdG/scale)┊       -Rᵀ       ┃
    //       ┃                ┊                 ┃
    //       ┃                ┊                 ┃
    // _     ┃                ┊                 ┃
    // M⁻¹ = ┃┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┼┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┃
    //       ┃                ┊                 ┃
    //       ┃                ┊                 ┃
    //       ┃                ┊                 ┃
    //       ┃       -R       ┊        0        ┃
    //       ┃                ┊                 ┃
    //       ┃                ┊                 ┃
    //       ┗                                  ┛
    //
    // R = triu(dXᵀdG)
    // D = diag(dXᵀdG)
    //
    const {m,M,N, dX, dXdG,
                  dG, dGdG, M_WAAW, Bei: u, scale} = this;

    if(      v.length !== N ) throw new Error('Assertion failed.');
    if(     Hv.length !== N ) throw new Error('Assertion failed.');
    if(indices.length !== N ) throw new Error('Assertion failed.');
    if( n_fix%1 !== 0 ) throw new Error('Assertion failed.');
    if( !(n_fix >= 0) ) throw new Error('Assertion failed.');
    if( !(n_fix <= N) ) throw new Error('Assertion failed.');

/*DEBUG*/    if( ! indices.slice().sort((i,j) => i-j).every((i,I) => i===I) )
/*DEBUG*/      throw new Error('Assertion failed.');

/*DEBUG*/    for( let i=0; i < n_fix; i++ )
/*DEBUG*/      if( v[indices[i]] !== 0 ) throw new Error('Assertion failed.');

    //
    // RIGHT PART
    //
    for( let i=0; i < m; i++ )
    { let sum = 0;
      for( let k=n_fix; k < N; k++ )
      { const                j = indices[k];
        sum += dG[N*i+j] * v[j];
      }
      u[i] = sum / scale;
    }

    for( let i=0; i < m; i++ )
    { let sum = 0;
      for( let k=n_fix; k < N; k++ )
      { const                j = indices[k];
        sum += dX[N*i+j] * v[j];
      }
      u[m+i] = sum;
    }

    //              _
    // MIDDLE PART (M PART)
    //
    if( 0 < n_fix )
    {
      // UPPER LEFT BLOCK
      for( let i=0; i < m; i++ )
      for( let j=0; j <=i; j++ )
      { let sum = 0;
        for( let l=0; l < n_fix; l++ )
        { const         k = indices[l];
          sum += dG[N*i+k] * dG[N*j+k];
        }
        M_WAAW[2*m*i+j] = sum;
      }
      
      // LOWER LEFT BLOCK
      for( let i=0; i < m; i++ )
      for( let j=0; j < m; j++ )
      { let sum = 0;
        for( let l=0; l < n_fix; l++ )
        { const         k = indices[l];
          sum += dX[N*i+k] * dG[N*j+k];
        }
        M_WAAW[2*m*(m+i) + j] = sum;
      }
      
      // LOWER RIGHT BLOCK
      for( let i=0; i < m; i++ )
      for( let j=0; j <=i; j++ )
      { let sum = 0;
        for( let l=0; l < n_fix; l++ )
        { const         k = indices[l];
          sum += dX[N*i+k] * dX[N*j+k];
        }
        M_WAAW[2*m*(m+i)+(m+j)] = sum * scale;
      }
    }
    else {
      for( let i=0; i < m*2; i++ )
      for( let j=0; j <=i  ; j++ )
        M_WAAW[2*m*i+j] = 0;
    }

    // UPPER LEFT BLOCK
    for( let i=0; i < m; i++ )
    for( let j=0; j <=i; j++ ) {
      M_WAAW[2*m*i+j] -= dGdG[M*i+j];
      M_WAAW[2*m*i+j] /= scale
    }

    for( let i=0; i < m; i++ )
      M_WAAW[2*m*i+i] -= dXdG[M*i+i];

    // LOWER LEFT BLOCK
    for( let i=0; i < m; i++ )
    for( let j=i; j < m; j++ )
      M_WAAW[2*m*(m+i) + j] -= dXdG[M*i+j];

    // FIXME: It remains to be seen whether or not LDL is accurate enough.
    //        If it isn't, we have to switch to PLDLP (Bunch-Kaufman).
    _ldl_decomp(2*m,2*m,   M_WAAW,0);
    _ldl_solve (2*m,2*m,1, M_WAAW,0, u,0);

    //
    // LEFT PART
    //
    Hv.fill(0.0, 0,N);

    for( let j=0    ; j < m; j++ )
    for( let k=n_fix; k < N; k++ )
    { const           i = indices[k];
      Hv[i] += dG[N*j+i] * u[j];
    }

    for( let k=n_fix; k < N; k++ )
    { const i = indices[k];
         Hv[i] /= scale;
    }

    for( let j=0    ; j < m; j++ )
    for( let k=n_fix; k < N; k++ )
    { const           i = indices[k];
      Hv[i] += dX[N*j+i] * u[m+j];
    }

    for( let k=n_fix; k < N; k++ )
    { const      i = indices[k];
      Hv[i] += v[i] / scale;
    }
  }
}