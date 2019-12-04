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

import {array, asarray, NDArray} from '../nd_array'
import {lstsq} from '../la/lstsq'


export function* root_newton_gen( fJ, x0 )
{
  let x = x0 instanceof NDArray
    ? x0.mapElems('float64') // <- protection copy
    : array('float64', x0);
  x0 = undefined;

  if( x.ndim !== 1 )
    throw new Error('root_newton_gen(fJ, x0): x0.ndim must be 1.');

  const [N] = x.shape;

  for(;;)
  {
    const [f,J] = fJ(x).map( a => asarray(a) );

    yield [x,f,J];

    if( J.ndim !== 2 || J.shape[0] !== N
                     || J.shape[1] !== N
     || f.ndim !== 1 || f.shape[0] !== N )
      throw new Error('root_newton_gen(fJ, x0: float[N]): fJ return type must be [float[N],float[N,N]].');

    const X =  x.data,
         dx = lstsq(J, f.reshape(N,1)).reshape(N),
         DX = dx.data;

    for( let i=DX.length; i-- > 0; )
      DX[i] = X[i] - DX[i];

    x = dx;
  }
}
