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

import {albaali_fletcher} from './albaali_fletcher'


export const strong_wolfe = (opt={}) => {
  opt = {...opt};

  for( const [old_name, new_name] of [
    'c1', 'fRed',
    'c2', 'gRed',
    'c3', 'grow'
  ])
    if( !(new_name in opt) && old_name in opt ) {
      console.warn(`strong_wolfe(opt): opt.${old_name} is deprecated, use opt.${new_name} instead.`);
      opt[new_name] = opt[old_name];
      delete          opt[old_name];
    }

  return albaali_fletcher(opt);
}
