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

import {forEachItemIn, CUSTOM_MATCHERS} from '../jasmine_utils'
import {tabulate} from '../tabulate'
import {_shuffle} from "../_test_data_generators";

import {norm} from '../la/norm'

import {fit_lin} from './fit_lin'
import {num_grad} from './num_grad'


describe('fit_lin', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  forEachItemIn(
    function*(){

      for( let run=512; run-- > 0; )
      {
        let funcs = [
          x => 1,
          x => x,
          x => x*x,
          x => x*x*x,
          x => 1/(1 + x*x),
          x => Math.exp(x),
          x => Math.sign(x)*Math.sqrt(Math.abs(x)),
          x => Math.sin(x)
        ];
        _shuffle(funcs);

        const N = Math.random()*funcs.length + 1 | 0;

        funcs = funcs.slice(0,N);
        const coeffs = Float64Array.from(funcs, () => Math.random()*2 - 1);
        yield [funcs, coeffs];
      }
    }()
  ).it(`fits random univariate functions correctly.`, ([funcs, coeffs]) => {
    const N = Math.random()*1024 + 128 | 0;

    const x = tabulate([N], () => Math.random()*8 - 4),
          y = x.mapElems(
            x => funcs.reduce((y,f,i) => y + coeffs[i]*f(x), 0)
          );

    for( const args of [
      [funcs],
      [0, funcs],
      [1e-8, funcs]
    ])
    {
      const f = fit_lin(x,y, ...args),
            Y = x.mapElems( x => f(x) );

      expect(f.coeffs).toBeAllCloseTo(coeffs);
      expect(Y).toBeAllCloseTo(y);
    }
  })


  forEachItemIn(
    function*(){

      for( let run=137; run-- > 0; )
      {
        let funcs = [
          ([x,y]) => 1,
          ([x,y]) => x,
          ([x,y]) => x*y,
          ([x,y]) =>   y,
          ([x,y]) => Math.hypot(x,y),
          ([x,y]) => Math.atan2(x,y),
          ([x,y]) => Math.cos(x),
          ([x,y]) => Math.sin(y),
        ];
        _shuffle(funcs);

        const N = Math.random()*funcs.length + 1 | 0;

        funcs = funcs.slice(0,N);
        const coeffs = Float64Array.from(funcs, () => Math.random()*2 - 1);
        yield [funcs, coeffs];
      }
    }()
  ).it(`fits random bivariate functions correctly.`, ([funcs, coeffs]) => {
    const N = Math.random()*512 + 128 | 0;

    const x = tabulate([N,2], () => Math.random()*8 - 4),
          y = [];

    for( const {data: xy} of x )
      y.push(
        funcs.reduce((y,f,i) => y + coeffs[i]*f(xy), 0)
      );

    for( const args of [
      [funcs],
      [0, funcs],
      [1e-8, funcs]
    ])
    {
      const f = fit_lin(x,y, ...args),
            Y = [];

      for( const {data: xy} of x )
        Y.push( f(xy) );

      expect(f.coeffs).toBeAllCloseTo(coeffs);
      expect(Y).toBeAllCloseTo(y);
    }
  })
});

