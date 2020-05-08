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


export class Comparator extends Function
{
  constructor( compare_fn )
  {
    const self = (x,y) =>
    {
      const len = x.length;
      if(   len!==y.length )
        throw new Error('arrays.Comparator(x,y): x and y must have same length.');

      for( let i=0; i < len; i++ ) {
        const    c = compare_fn(x[i],y[i]);
        if(0 !== c)
          return c;
      }

      return 0;
    }
    self.prototype = Comparator.prototype;
    return Object.freeze(self);
  }
}


export const compare = new Comparator( (x,y) => (x > y) - (x < y) )
