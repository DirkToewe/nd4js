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

export * from './min_max'


export function* linspace( start, end, num )
{
  if( ! (num%1 === 0)   ) throw new Error('linspace(start,end,num): num must be an integer greater than 1.');
  if( ! (num    >  1)   ) throw new Error('linspace(start,end,num): num must be an integer.greater than 1.');
  if( ! isFinite(start) ) throw new Error('linspace(start,end,num): start must be a finite number.');
  if( ! isFinite(end  ) ) throw new Error('linspace(start,end,num): end must be a finite number.');

  for( let i=0; i < num; i++ )
  {
    const s = i / (num-1);
    yield start * (1-s) + s * end;
  }
}


export function* range(start, stop, step=1)
{
  if( start%1 !== 0 ) throw new Error('range(start, stop, step=1): start not a valid integer.');
  if( stop %1 !== 0 ) throw new Error('range(start, stop, step=1): stop not a valid integer.');
  if( step %1 !== 0 ) throw new Error('range(start, stop, step=1): step not a valid integer.');
  if( step    === 0 ) throw new Error('range(start, stop, step=1): step must not be 0.');

  if( 0 < step ) for( let i=start; i < stop; i+=step ) yield i;
  else           for( let i=start; i > stop; i+=step ) yield i;
}


export function* cartesian_prod( ...seqs )
{
  seqs = seqs.map(
    seq => (
      seq instanceof        Array ||
      seq instanceof Float32Array ||
      seq instanceof Float64Array ||
      seq instanceof    Int8Array ||
      seq instanceof   Int16Array ||
      seq instanceof   Int32Array ||
      seq instanceof   Uint8Array ||
      seq instanceof  Uint16Array ||
      seq instanceof  Uint32Array
    ) ? seq.slice() : [...seq]
  );

  const result = new Array(seqs.length);

  function* iter( i ) {
    if( i < seqs.length )
      for( const x of seqs[i] ) {
        result[i] = x;
        yield* iter(i+1);
      }
    else
      yield result.slice();
  }

  yield* iter(0);
}


export function* enumerate( seq )
{
  let i=0;
  for( const x of seq )
    yield [i++, x];
}


export function* zip( ...seqs )
{
  if( ! (0 < seqs.length) )
    throw new Error('zip(...seqs): seqs.length must be at least 1.');

  const N = seqs.length;

  for( let i=N; i-- > 0; )
    seqs[i] = seqs[i][Symbol.iterator]();

  for(;;) {
    const next = [];

    for( let i=0; i < N; i++ ) {
      const     item = seqs[i].next();
             if(item.done ) return;
      next.push(item.value);
    }

    yield next;
  }
}


export function* repeat( n, seq )
{
  if( undefined === seq ) {
    seq = n;
          n = Infinity;
  }

  if( n !== Infinity ) {
    if( n%1 !== 0 ) throw new Error('repeat(n,seq): n must be a non-negative int or Infinity.');
    if( ! (n >= 0)) throw new Error('repeat(n,seq): n must be a non-negative int or Infinity.');
  }

  seq = (
    seq instanceof        Array ||
    seq instanceof Float32Array ||
    seq instanceof Float64Array ||
    seq instanceof    Int8Array ||
    seq instanceof   Int16Array ||
    seq instanceof   Int32Array ||
    seq instanceof   Uint8Array ||
    seq instanceof  Uint16Array ||
    seq instanceof  Uint32Array
  ) ? seq.slice() : [...seq];

  for( let i=n; i-- > 0; )
    yield* seq;
}
