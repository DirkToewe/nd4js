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

import {forEachItemIn, CUSTOM_MATCHERS} from '../jasmine_utils'
import {FrobeniusNorm, norm} from './norm'
import {tabulate} from '../tabulate'


describe('norm', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from

      for( let run=16*1024; run-- > 0; )
      {
        const values = Float64Array.from({length: randInt(0,128)}, () => Math.random()*2 - 1);
        Object.freeze(values.buffer);
        yield values;
      }
    }()
  ).it('FrobeniusNorm works on random examples', values => {
    const ref = Math.hypot(...values);

    const norm = new FrobeniusNorm();
    for( const x of values )
      norm.include(x);

    expect(norm.result).toBeCloseTo(ref, {atol: 0, rtol: 1e-15});

    norm.reset();
    for( const x of values )
      norm.include(x);

    expect(norm.result).toBeCloseTo(ref, {atol: 0, rtol: 1e-15});
  });

  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from;

      for( let run=16*1024; run-- > 0; )
      {
        const values = new Float64Array( randInt(0,128) );

        for( let i=randInt(0,values.length); i >= 0; --i )
          values[i] = Math.random()*1e-162;

        // SHUFFLE
        for( let i=values.length; i > 0; )
        {
          const j = randInt(0,i--);
          const tmp = values[i];
                      values[i] = values[j];
                                  values[j] = tmp;
        }

        Object.freeze(values.buffer);
        yield values;
      }
    }()
  ).it('FrobeniusNorm is underflow-safe', values => {
    const ref = Math.hypot(...values);

    const norm = new FrobeniusNorm();
    for( const x of values )
      norm.include(x);

    expect(norm.result).toBeCloseTo(ref, {atol: 0, rtol: 1e-15});

    norm.reset();
    for( const x of values )
      norm.include(x);

    expect(norm.result).toBeCloseTo(ref, {atol: 0, rtol: 1e-15});
  });

  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from;

      for( let run=16*1024; run-- > 0; )
      {
        const shape = Int32Array.from({length: randInt(0,5)}, () => randInt(1,9));
        Object.freeze(shape.buffer);
        const A = tabulate(shape, () => Math.random()*2-1);
        Object.freeze(A.data.buffer);
        yield A;
      }
    }()
  ).it("norm(_,'fro') works for random examples", A => {
    const ref = Math.hypot(...A.data);

    expect( norm(A                ) ).toBeCloseTo(ref, {atol: 0, rtol: 1e-14});
    expect( norm(A,'fro'          ) ).toBeCloseTo(ref, {atol: 0, rtol: 1e-14});
    expect( norm(A,'fro',null     ) ).toBeCloseTo(ref, {atol: 0, rtol: 1e-14});
    expect( norm(A,'fro',undefined) ).toBeCloseTo(ref, {atol: 0, rtol: 1e-14});
  });

  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from;

      for( let run=16*1024; run-- > 0; )
      {
        const shape = Int32Array.from({length: randInt(0,5)}, () => randInt(1,9));
        Object.freeze(shape.buffer);

        const  A = tabulate(shape, () => 0),
          values = A.data;

        for( let i=randInt(0,values.length); i >= 0; --i )
          values[i] = Math.random()*1e-162;

        // SHUFFLE
        for( let i=values.length; i > 0; )
        {
          const j = randInt(0,i--);
          const tmp = values[i];
                      values[i] = values[j];
                                  values[j] = tmp;
        }

        Object.freeze(values.buffer);
        yield A;
      }
    }()
  ).it("norm(_,'fro') is underflow-safe", A => {
    const ref = Math.hypot(...A.data);

    expect( norm(A                ) ).toBeCloseTo(ref, {atol: 0, rtol: 1e-14});
    expect( norm(A,'fro'          ) ).toBeCloseTo(ref, {atol: 0, rtol: 1e-14});
    expect( norm(A,'fro',null     ) ).toBeCloseTo(ref, {atol: 0, rtol: 1e-14});
    expect( norm(A,'fro',undefined) ).toBeCloseTo(ref, {atol: 0, rtol: 1e-14});
  });
})

