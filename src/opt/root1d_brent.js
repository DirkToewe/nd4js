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

import {nextDown,
        nextUp} from "../dt/float64_utils";


// REFERENCES
// ----------
// .. [1] https://en.wikipedia.org/wiki/Brent's_method
// .. [2] http://www.netlib.no/netlib/go/zeroin.f
// .. [3] https://blogs.mathworks.com/cleve/2015/10/26/zeroin-part-2-brents-version/?s_tid=blogs_rc_1


const min = (x,y) => Math. min(x,y),
      max = (x,y) => Math. max(x,y),
      abs =  x    => Math. abs(x),
     sign =  x    => Math.sign(x);
      

export function root1d_brent(F, a,b)
{
  a *= 1;
  b *= 1;
  if( isNaN(a) ) throw new Error('root1d_bisect(F,a,b): a must be a number.');
  if( isNaN(b) ) throw new Error('root1d_bisect(F,a,b): b must be a number.');

  if( ! (F instanceof Function) ) throw new Error('root1d_bisect(F,a,b): F must be a function.');

  let fa = F(a),
      fb = F(b);

  if( ! (sign(fa)*fb <= 0) )
    throw new Error('root1d_falsi(F,a,b): F(a) and F(b) must have opposite signs.');

  loop:for(;;)
  {
    // interval [b,c] contains the root, where a is a previous iterate or c
    let c     = a,
       fc     =  fa,
       Δb     = b-a,
       Δb_old = Δb;

    do
    {
      // fb should be the better than fa (b as in better)
      if( abs(fc) < abs(fb) ) {
         a =  b;  b =  c;  c =  a;
        fa = fb; fb = fc; fc = fa;
      }

      // [l,r] = (b,c) is the interval in which we are looking for a new iterate
      const l = nextUp  ( min(b,c) ),
            r = nextDown( max(b,c) );
      if( r < l || fb === 0 ) break loop;

      const bc = c-b; // <- b+xm would be the bisection point, i.e. mid point of [b,c]

      // see if a bisection is forced (either if no progress or interpolation bad)
      if( abs(Δb_old) < Number.EPSILON * abs(b) || abs(fa) <= abs(fb) )
        Δb_old = Δb = bc/2;
      else
      {
        let p,q,s = fb/fa;

        if( a === c )
        {
          // LINEAR INTERPOLATION
          p = bc*s;
          q =  1-s;
        }
        else
        {
          // INVERSE QUADRATIC INTERPOLATION
          const r = fb/fc;
                q = fa/fc;
          p = s * ( bc*q*(q-r) - (b-a)*(r-1) );
          q = (q-1)*(r-1)*(s-1);
        }

        if( p <= 0 ) p = -p;
        else         q = -q;

        s = Δb_old;
            Δb_old = Δb;

        // same as: if( |Δb_new| >= 3/4*|c-b|  OR  |Δb_new| >= |0.5*s| )
        if( 4*p >= 3*bc*q || p >= abs(0.5*s*q) )
          Δb_old = Δb = bc/2;
        else
          Δb = p/q;
      }

       a = b;
      fa =fb;

      b += Δb;
      if( !(b>=l) ) b = l;
      if( !(b<=r) ) b = r;

      fb = F(b);
    }
    while( sign(fb)*fc <= 0 )
  }

  return b;
}
