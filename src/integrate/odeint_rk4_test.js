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

import {odeint_rk4} from './odeint_rk4'

import {array} from '../nd_array'
import {CUSTOM_MATCHERS, forEachItemIn} from '../jasmine_utils'
import {solve} from '../la/solve'
import {tabulate} from '../tabulate'


describe('odeint_rk4', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  forEachItemIn([0.125, 0.25]).it("integrates Jansen's linkage correctly", dt => {
    const y0 = array([
        0.0,             0.0,
      -35.0353991829,   19.5071987762,
       17.7855809088,   37.4956412365,
       50.7538954493,  - 0.0954512769988,
        2.77183874468, -39.2021288959,
      -17.0054756067,  -84.0335669394,
      -28.0419102143,  -19.2671627536
    ]);
    Object.freeze(y0);
    Object.freeze(y0.data.buffer);

    const edges = [
      [0, 1],
      [0, 2],
      [0, 4],
      [1, 2],
      [1, 6],
      [2, 3],
      [3, 4],
      [4, 5],
      [4, 6],
      [5, 6]
    ];
    Object.freeze(edges);

    const angVel = Math.PI/180;

    const dy = (y,t) =>
    {
      const A = tabulate([14,14], () => 0),
            v = tabulate([14, 1], () => 0);

      let row=0;

      // linkage constraints
      for( const [i,j] of edges )
      {
        let dx = y(2*i  ) - y(2*j  ),
            dy = y(2*i+1) - y(2*j+1);
        const hyp = Math.hypot(dx,dy);
        dx /= hyp;
        dy /= hyp;
        A.set([row, 2*i  ],  dx);
        A.set([row, 2*i+1],  dy);
        A.set([row, 2*j  ], -dx);
        A.set([row, 2*j+1], -dy);
        v.set([row,0], 0);
        row++;
      }

      // fixed bearing BC
      A.set([row, 2*0  ], 1); v.set([row,0], 0); row++;
      A.set([row, 2*0+1], 1); v.set([row,0], 0); row++;

      // rotation BC
      let dx = y(2*3  ) -38.0,
          dy = y(2*3+1) - 7.8;
      
      A.set([row, 2*3  ], 1); v.set([row,0],  dy*angVel); row++;
      A.set([row, 2*3+1], 1); v.set([row,0], -dx*angVel); row++;

      const  result = solve(A,v).reshape(-1);
      return result;
    }

    const dist = (y, [i,j]) => {
      const xi = y(2*i), yi = y(2*i+1),
            xj = y(2*j), yj = y(2*j+1);
      return Math.hypot(
        xi-xj,
        yi-yj
      );
    };

    let y = y0;
    for( let t=0; t < 360; t += dt )
    {
      y = odeint_rk4(dy, y,t, dt);

      // check rod linkage distances
      for( const ij of edges )
        expect( dist(y, ij) ).toBeAllCloseTo( dist(y0, ij) );
    }

    expect(y).toBeAllCloseTo(y0);
  })
})
