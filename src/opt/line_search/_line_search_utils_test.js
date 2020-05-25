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

import {CUSTOM_MATCHERS, forEachItemIn} from '../../jasmine_utils'

import {linspace} from "../../iter";

import {num_grad} from "../num_grad";

import {_min1d_interp_ffg,
        _min1d_interp_ffgg,
        _min1d_interp_gg} from './_line_search_utils'


describe('_min1d_interp_ffgg', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  forEachItemIn(
    function*(){
      const N = 13,
            M = Math.PI;

        for( const a of linspace(-M,+M, N) ) if( 0 !== a )
        for( const b of linspace(-M,+M, N) )
        for( const c of linspace(-M,+M, N) )
        {    const d = Math.random()*M*2 - M;

          const f = x => ((a/3*x - a/2*(b+c))*x + a*b*c)*x + d,
                g = x => a*(x-b)*(x-c);

          for( const x1 of linspace(-M,+M, N) )
          for( const x2 of linspace(-M,+M, N) )
            if( x1 !== x2 )
            {
              expect( g(x1) ).toBeAllCloseTo( num_grad(f)(x1) )
              expect( g(x2) ).toBeAllCloseTo( num_grad(f)(x2) )

              const l = Math.min(b,c),
                    r = Math.max(b,c);

              if( a < 0 ) expect( f(l) ).toBeLessThanOrEqual( f(r) );
              else        expect( f(r) ).toBeLessThanOrEqual( f(l) );

              yield Object.freeze([
                [x1,x2, f(x1), f(x2), g(x1), g(x2)],
                a < 0 ? l : r
              ]);
            }
        }
    }()
  ).it('works on generated examples', ([args, xMin]) => {
    expect( _min1d_interp_ffgg(...args) ).toBeAllCloseTo(xMin, {atol: 1e-5});
  });

  forEachItemIn(
    function*(){
        const N = 17,
              M = Math.PI;

        for( let run=0; run++ < 773; )
        {
          const a = (Math.random()*M + Number.EPSILON) * (Math.random() < 0.5 ? +1 : -1),
                b =  Math.random()*M*2 - M,
                c =  Math.random()*M*2 - M,
                d =  Math.random()*M*2 - M;

          const f = x => ((a/3*x - a/2*(b+c))*x + a*b*c)*x + d,
                g = x => a*(x-b)*(x-c);

          for( const x1 of linspace(-M,+M, N) )
          for( const x2 of linspace(-M,+M, N) )
            if( x1 !== x2 )
            {
              expect( g(x1) ).toBeAllCloseTo( num_grad(f)(x1) )
              expect( g(x2) ).toBeAllCloseTo( num_grad(f)(x2) )

              const l = Math.min(b,c),
                    r = Math.max(b,c);

              if( a < 0 ) expect( f(l) ).toBeLessThanOrEqual( f(r) );
              else        expect( f(r) ).toBeLessThanOrEqual( f(l) );

              yield Object.freeze([
                [x1,x2, f(x1), f(x2), g(x1), g(x2)],
                a < 0 ? l : r
              ]);
            }
        }
    }()
  ).it('works given random examples', ([args, xMin]) => {
    expect( _min1d_interp_ffgg(...args) ).toBeAllCloseTo(xMin, {atol: 1e-5});
  });
})


describe('_min1d_interp_gg', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  forEachItemIn(
    function*(){
      const N = 37,
            M = Math.PI;

      for( const b of linspace(Number.EPSILON,M, N) )
      for( const z of linspace(-M,+M, N) )
      {
        const g = x => 2*b*(x-z);

        for( const x1 of linspace(-M,+M, N) )
        for( const x2 of linspace(-M,+M, N) )
          if( x1 !== x2 )
            yield Object.freeze([ [x1,x2, g(x1), g(x2)], z ]);
      }
    }()
  ).it('works on generated examples.', ([args, xMin]) => {
    expect( _min1d_interp_gg(...args) ).toBeAllCloseTo(xMin);
  });
})


describe('_min1d_interp_ffg', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })

  forEachItemIn(
    function*(){
      const N = 4,
            M = 2;

      for( let a = -M; a <= +M; a += 1/N )
      for( let b =1/N; b <= +M; b += 1/N )
      for( let z = -M; z <= +M; z += 1/N )
      {
        const f = x => a + b*(x-z)*(x-z),
              g = x =>   2*b*(x-z);

        for( let x1 = -M; x1 <= +M; x1 += 1/N )
        for( let x2 = -M; x2 <= +M; x2 += 1/N )
          if( x1 !== x2 )
            yield Object.freeze([ [x1,x2, f(x1),f(x2), g(x1)], z ]);
      }
    }()
  ).it('works on generated examples.', ([args, xMin]) => {
    expect( _min1d_interp_ffg(...args) ).toBeAllCloseTo(xMin);
  });

  forEachItemIn(
    function*(){
      for( let i=1024; i-- > 0; )
      {
        const a = Math.random()*8 - 4,
              b = Math.random()*8 + 1/4,
              z = Math.random()*8 - 4;

        const f = x => a + b*(x-z)*(x-z),
              g = x =>   2*b*(x-z);

        for( let j=1024; j-- > 0; )
        {
          const x1 =  Math.random()*8 - 4,
           x2 = x1 + (Math.random()*4 + 1e-4) * (Math.random() < 0.5 ? +1 : -1);
          yield Object.freeze([ [x1,x2, f(x1),f(x2), g(x1)], z ]);
        }
      }
    }()
  ).it('works on random examples.', ([args, xMin]) => {
    expect( _min1d_interp_ffg(...args) ).toBeAllCloseTo(xMin, {atol: 1e-6});
  });
})
