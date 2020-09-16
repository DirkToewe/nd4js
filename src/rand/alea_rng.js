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

import {ARRAY_TYPES} from '../dt'
import {NDArray} from '../nd_array';

import {_giv_rot_qr,
        _giv_rot_rows} from "../la/_giv_rot";
import { RNG } from '../_test_data_generators';


const mul32 = 2**+32,
      div32 = 2**-32,
      div53 = 2**-53;


// Slightly modified version of the Mash hashing algorithm as implemented by Johanne Baagøe.
// Sadly Baagoe's original blog posts (http://baagoe.com/en/RandomMusings/hash/avalanche.xhtml)
// are not available anymore. This hashing algorithm has the desireable "Avalanche Effect" property
// which mean each change in an input bit on average changes 50% of the output (hash value) bits.
// See: https://en.wikipedia.org/wiki/Avalanche_effect.
function mash( str, seed )
{
  str = str.toString();

  for( let i=0; i < str.length; i++ )
  {
    seed += str.charCodeAt(i);
    let     temp = 0.02519603282416938 * seed;
    seed  = temp >>> 0;
    temp -= seed;
    temp *= seed;
    seed  = temp >>> 0;
    temp -= seed;
    seed += temp * mul32;
  }

  return seed;
}

// Slightly modernized version of the Alea pseudo-random number generator by Johanne Baagøe.
// Sadly Baagoe's original blog posts (http://baagoe.com/en/RandomMusings/javascript/) are not
// available anymore. There is however a partial mirror of his site on GitHub:
//
//   https://github.com/nquinlan/better-random-numbers-for-javascript-mirror
//
export class AleaRNG
{
  constructor( seed )
  {
    if( seed == null )
      throw new Error('Assertion failed.');

    seed = seed.toString();
    let s0 = mash(' ', 0xefc8249d),
        s1 = mash(' ', s0),
        s2 = mash(' ', s1),
        t0 = mash(seed,s2),
        t1 = mash(seed,t0),
        t2 = mash(seed,t1);
    s0 >>>= 0; s0 -= t0 >>> 0; s0 *= div32; s0 += s0 < 0;
    s1 >>>= 0; s1 -= t1 >>> 0; s1 *= div32; s1 += s1 < 0;
    s2 >>>= 0; s2 -= t2 >>> 0; s2 *= div32; s2 += s2 < 0;
    // this#nextNormal = NaN;
    Object.assign(this, {s0,s1,s2, c: 1, __nextNormal: NaN});
    Object.seal(this);
  }


  __next() {
    const t = 2091639 * this.s0 + this.c * div32;
    this.s0 = this.s1;
              this.s1 = this.s2;
    return              this.s2 = t - (this.c = t|0);
  }


  bool() {
    return this.uniform() < 0.0;
  }


  int( from, until )
  {
    if( until === undefined ) {
        until = from;
                from = 0;
    }
    if( from%1 !== 0 ) throw new Error('AlreaRNG::int(from,until): from must be a valid int.');
    if(until%1 !== 0 ) throw new Error('AlreaRNG::int(from,until): until must be a valid int.');
    if( !(from < until) ) throw new Error('AlreaRNG::int(from,until): from must be less than until.');
    return Math.floor( this.uniform(from,until) );
  }


  shuffle( array, from, until )
  {
    if(  array.length%1 !== 0 ) throw new Error('AlreaRNG::normal(array,from=0,until=array.length): array must be Array-like.');
    if(!(array.length   >=  0)) throw new Error('AlreaRNG::normal(array,from=0,until=array.length): array must be Array-like.');

    if( null == from ) from = 0;
    if( null == until) until= array.length;

    if(0!==from %1        ) throw new Error('AlreaRNG::normal(array,from=0,until=array.length): from must be int.');
    if(0!==        until%1) throw new Error('AlreaRNG::normal(array,from=0,until=array.length): until must be int.');
    if( ! (from         >= 0    ) ) throw new Error('AlreaRNG::normal(array,from=0,until=array.length): from must be non-negative.');
    if( ! (from         <= until) ) throw new Error('AlreaRNG::normal(array,from=0,until=array.length): from must not be greater than until.');
    if( ! (array.length >= until) ) throw new Error('AlreaRNG::normal(array,from=0,until=array.length): until must not be greater than array.length.');

    // https://en.wikipedia.org/wiki/Fisher-Yates_shuffle
    for( let i=from; i < until-1; i++ )
    { const j = this.int(i,until),
           aj = array[j];
                array[j] = array[i];
                           array[i] = aj;
    }
  }


  uniform( lo=-1, hi=+1 )
  {
    lo*=1; if( ! isFinite(lo) ) throw new Error('Assertion failed.');
    hi*=1; if( ! isFinite(hi) ) throw new Error('Assertion failed.');
    const s = this.__next() + 
             (this.__next() * 0x200000 | 0) * div53;
    return lo*(1-s) + s*hi;
  }


  normal( mean=0, sigma=1 )
  {
    mean *= 1; if( ! isFinite(mean ) ) throw new Error('AlreaRNG::normal(mean,sigma): mean must be valid number.');
    sigma*= 1; if( ! isFinite(sigma) ) throw new Error('AlreaRNG::normal(mean,sigma): sigma must be valid number.');
    const       nxt  =  this.__nextNormal;
    if( ! isNaN(nxt) ){ this.__nextNormal = NaN;
      return    nxt*sigma + mean;
    }
    // https://en.wikipedia.org/wiki/Marsaglia_polar_method
    let x,y,r;
    do {
      x = this.uniform();
      y = this.uniform();
      r = x*x + y*y;
    }
    while(r > 1 || r == 0);

    const               z = Math.sqrt(-2 * Math.log(r) / r);
    this.__nextNormal = z*x;
    return       mean + z*y*sigma;
  }


  ortho( dtype, ...shape )
  {
    // REFERENCES:
    //   - https://blogs.sas.com/content/iml/2012/03/28/generating-a-random-orthogonal-matrix.html
    if( ! (dtype in ARRAY_TYPES) ) {
      shape.unshift(dtype);
      dtype = 'float64';
    }

    if( shape.length  <  1 ) throw new Error(`AlreaRNG::orthogonal([dtype,] ...shape): shape.length must be at least 1.`);
    if( shape.length === 1 )
        shape.push(shape[0]);

    for( let i=shape.length; i-- > 0; )
      if( shape[i] % 1 !== 0 || shape[i] < 1 )
        throw new Error(`AlreaRNG::orthogonal([dtype,] ...shape): shape[${i}] = ${shape[i]} not a positive integer.`);
    shape = Int32Array.from(shape);

    const [M,N] = shape.subarray(-2),
           K    = Math.max(M,N),
             L  = Math.min(M,N);

    const     DTypeArray = ARRAY_TYPES[dtype],
      U = new DTypeArray( shape.reduce((x,y) => x*y) ),
      Q = new Float64Array(K*L);

    for( let U_off=U.length; (U_off-=M*N) >= 0; )
    {
      // INIT Q TO IDENTITY
      for( let i=K; i-- > 0; )
      for( let j=L; j-- > 0; )
        Q[L*i+j] = i!==j ? 0 : (this.bool() ? -1 : +1); // <- randomly flips sign of det(Q)

      // PSEUDO QR DECOMPOSITION
      // -----------------------
      // A vector with random entries drawn from a normal distribution multiplied by
      // an orthogonal matrix results in another vector with random normal entries.
      // This fact allows us to draw the entries of the QR decomposed matrix "on-the-fly",
      // avoiding the memory and computation overhead of the decomposed matrix.
      for( let j=0; j < K; j++ )
      {
        let A_jj = this.normal();

        for( let i=j; ++i < K; )
        {
          const   A_ij = this.normal(),
            [c,s,norm] =_giv_rot_qr(A_jj,A_ij);
          if( 0 === s ) continue;        A_jj = norm;
          _giv_rot_rows(Q, Math.min(i+1,L), L*j,L*i, c,s);
        }
      }

      // COPY Q -> U (potentially transposing the entries)
      for( let i=K; i-- > 0; )
      for( let j=L; j-- > 0; )
        U[U_off + (M < N ? N*j+i : N*i+j)] = Q[L*i+j];
    }

    return new NDArray(shape,U);
  }
}
