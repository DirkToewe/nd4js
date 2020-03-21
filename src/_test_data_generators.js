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


import {tabulate} from './tabulate'
import {matmul} from './la/matmul'
import {rand_ortho} from './la/rand_ortho'


export function _rand_int(from,until)
{
  if( 0  !==  from %1 ) throw new Error('Assertion failed.');
  if( 0  !==  until%1 ) throw new Error('Assertion failed.');
  if( !(from < until) ) throw new Error('Assertion failed.');

  return Math.floor( Math.random() * (until-from) ) + from;
}


///** Generates a series of shapes that are broadcast compatible except for their tails.
// */
//export function _rand_shapes( ...tails ) // <- TODO
//{
//  throw new Error('Not yet implemented.');
//}


export function _shuffle( array, from, until )
{
  if( array.length%1 !== 0 )
    throw new Error('Assertion failed.');

  if( null == from ) from = 0;
  if( null == until) until= array.length;

  if(0!==from %1        ) throw new Error('Assertion failed.');
  if(0!==        until%1) throw new Error('Assertion failed.');
  if( ! (from <= until) ) throw new Error('Assertion failed.');

  // https://en.wikipedia.org/wiki/Fisher-Yates_shuffle
  for( let i=from; i < until-1; i++ )
  { const j = _rand_int(i,until),
         aj = array[j];
              array[j] = array[i];
                         array[i] = aj;
  }
}


export function _shuffled( array, from, until )
{
   array = array.slice();
  _shuffle(array,from,until);
  return   array;
}


export function _rand_rows0(...shape)
{
  if( !(2 <= shape.length) ) throw new Error('Assertion failed.');

  const AA = tabulate(shape, 'float64', () => Math.random()*8-4);

  const N = shape.pop(),
        M = shape.pop();

  const A = AA.data;

  const rows = new Int32Array(M);

  for( let A_off=A.length; (A_off -= M*N) >= 0; )
  {
    for( let i=M; i-- > 0; ) rows[i] = i;
    const rank = _rand_int(0,M);

    for( let m=M; m > rank; )
    {
      const      k = _rand_int(0,m--),
        i = rows[k];
            rows[k] = rows[m];

      for( let j=N; j-- > 0; )
        A[A_off + N*i+j] = 0;
    }
  }

  Object.freeze(A.buffer);
  Object.freeze(AA);
  return AA;
}


export function _rand_cols0(...shape)
{
  if( !(2 <= shape.length) ) throw new Error('Assertion failed.');

  const AA = tabulate(shape, 'float64', () => Math.random()*8-4);

  const N = shape.pop(),
        M = shape.pop();

  const A = AA.data;

  const cols = new Int32Array(N);

  for( let A_off=A.length; (A_off -= M*N) >= 0; )
  {
    for( let i=N; i-- > 0; ) cols[i] = i;
    const rank = _rand_int(0,N);

    for( let n=N; n > rank; )
    {
      const      k = _rand_int(0,n--),
        j = cols[k];
            cols[k] = cols[n];

      for( let i=M; i-- > 0; )
        A[A_off + N*i+j] = 0;
    }
  }

  Object.freeze(A.buffer);
  Object.freeze(AA);
  return AA;
}


export function _rand_rankdef(...shape)
{
  if( !(2 <= shape.length) ) throw new Error('Assertion failed.');

  const N = shape.pop(),
        M = shape.pop(),
        L = Math.min(M,N);

  // use random, mostly rank-deficient SVD to generate test matrix A
  const ranks = tabulate(shape, 'int32', () => _rand_int(0,L+1) ), // <- ranks
            U = rand_ortho('float64', ...shape,M,L),
            V = rand_ortho('float64', ...shape,L,N),
            S = tabulate([...shape,L,L], 'float64', (...idx) => {
              const j = idx.pop(),
                    i = idx.pop(), rank = ranks(...idx);

              if( i !== j || rank <= i ) return 0; 

              return Math.random()*8 + 0.1
            }),
            A = matmul(U,S,V);

  // add random scaling
  for( let S_off=S.data.length; (S_off -= L*L) >= 0; )
  {
    const scale = Math.random()*1e9 + 1;

    for( let i=L; i-- > 0; )
      S.data[S_off + L*i+i] *= scale;
  }

  Object.freeze(ranks.data.buffer);
  Object.freeze(A.data.buffer);
  return [A, ranks];
}
