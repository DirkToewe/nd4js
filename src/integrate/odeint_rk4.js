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

import {asarray} from '../nd_array'
import {zip_elems} from '../zip_elems'


export function odeint_rk4( dy, y0, t0, dt )
{
  y0 = asarray(y0);

  // https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
  const k1 = dy(           y0,                          t0       ).mapElems(dy => dy*dt),
        k2 = dy(zip_elems([y0,k1], (y0,k1) => y0+k1/2), t0 + dt/2).mapElems(dy => dy*dt),
        k3 = dy(zip_elems([y0,k2], (y0,k2) => y0+k2/2), t0 + dt/2).mapElems(dy => dy*dt),
        k4 = dy(zip_elems([y0,k3], (y0,k3) => y0+k3  ), t0 + dt  ).mapElems(dy => dy*dt);

  return zip_elems(
    [y0,k1,k2,k3,k4],
    (y0,k1,k2,k3,k4) => y0 + (k1 + 2*(k2+k3) + k4)/6
  );
}
