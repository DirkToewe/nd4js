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

export class LineSearchError extends Error
{
}

export class LineSearchNoProgressError extends LineSearchError
{
}

export class LineSearchBisectionError extends LineSearchError
{
  constructor(x,f,g)
  {
    super();
    Object.assign(this, {x,f,g});
  }
}

export class LineSearchBoundReachedError extends LineSearchError
{
  constructor(x,f,g)
  {
    super();
    Object.assign(this, {x,f,g});
  }
}
