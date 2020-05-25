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

import { LineSearchNoProgressError } from "./line_search/line_search_error";


export class LBFGS_Solver
{
  constructor( M, N )
  {
    this.dXdG= new Float64Array(M);
    this.dX  = new Float64Array(M*N);
    this.dG  = new Float64Array(M*N);
    this.tmp = new Float64Array(M);

    Object.assign(this, {M,N});
    this.m     = 0;
    this.off   = 0;

    Object.seal(this);
  }


  update( dx, dg, )
  {
    const {M,N, dG,
              dXdG,
              dX} = this;

    if( dx.length !== N ) throw new Error('Assertion failed.');
    if( dg.length !== N ) throw new Error('Assertion failed.');

    let                      dxdg = 0;
    for( let j=N; j-- > 0; ) dxdg += dx[j]*dg[j];

    if( !(0 < dxdg) )
      throw new LineSearchNoProgressError();

    if( this.m === M ) this.off = (this.off+1) % M;
    else               this.m++;

    const i = (this.off + this.m-1) % M;

    dXdG[i] = dxdg;
    for( let j=N; j-- > 0; ) dX[N*i+j] = dx[j];
    for( let j=N; j-- > 0; ) dG[N*i+j] = dg[j];
  }


  forget( k )
  {
    if( 0 !== k%1      ) throw new Error('Assertion failed.');
    if( !(k >= 0     ) ) throw new Error('Assertion failed.');
    if( !(k <= this.m) ) throw new Error('Assertion failed.');

    this.off = (this.off+k) % this.M;
    this.m -= k;
  }


  compute_Hv_phase1( v )
  {
    const {M,N, dG,
           m, dXdG,
              dX, off, tmp: α} = this;

    if( v.length !== N ) throw new Error(`Assertion failed.`);

    for( let I=m; I-- > 0; ) { let i = off+I;
                                   i -= (i>=M)*M;
      let                      αi = 0;
      for( let j=N; j-- > 0; ) αi += dX[N*i+j] * v[j];
                       α[i] = (αi /= dXdG[i]);
      for( let j=N; j-- > 0; )
        v[j] -= αi * dG[N*i+j];
    }
  }


  compute_Hv_phase2( v )
  {
    const {M,N, dG,
           m, dXdG,
              dX, off, tmp: α} = this;

    if( v.length !== N ) throw new Error(`Assertion failed.`);

    for( let I=0; I < m; I++ ) { let i = off+I;
                                     i-= (i>=M)*M;
      let                      βi = 0;
      for( let j=N; j-- > 0; ) βi += dG[N*i+j] * v[j];
                               βi /= dXdG[i];
      for( let j=N; j-- > 0; )
        v[j] += (α[i] - βi) * dX[N*i+j];
    }
  }
}
