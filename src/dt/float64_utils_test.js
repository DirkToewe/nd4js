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

import {nextDown,
        nextUp} from './float64_utils'


describe('nextUp', () => {
  it('nextUp(NaN) is NaN', () => {
    expect( nextUp(NaN) ).toEqual(NaN);
  });

  it('nextUp(-∞) === -MAX_VALUE', () => {
    expect( nextUp(-Infinity) ).toBe( -Number.MAX_VALUE );
  });

  it('nextUp( -MIN_VALUE * (2**54-1) ) === -MIN_VALUE * (2**54-2)', () => {
    expect(
      nextUp( -Number.MIN_VALUE * (2**54-1) )
    ).toBe(
              -Number.MIN_VALUE * (2**54-2)
    );
  });

  it('nextUp(-MIN_VALUE) === 0', () => {
    expect( nextUp(-Number.MIN_VALUE) ).toBe( 0 );
  });

  it('nextUp(-0) === MIN_VALUE', () => {
    expect( nextUp(-0.0) ).toBe( Number.MIN_VALUE );
  })

  it('nextUp(+0) === MIN_VALUE', () => {
    expect( nextUp(+0.0) ).toBe( Number.MIN_VALUE );
  });

  it('nextUp(+MIN_VALUE) === MIN_VALUE*2', () => {
    expect( nextUp(Number.MIN_VALUE) ).toBe( Number.MIN_VALUE*2 );
  });

  it('nextUp( MIN_VALUE*(2**51  ) ) === MIN_VALUE*(2**51+1)', () => {
    expect(
      nextUp( Number.MIN_VALUE*(2**51  ) )
    ).toBe(
              Number.MIN_VALUE*(2**51+1)
    );
  });

  it('nextUp(-MIN_VALUE*(2**51+1) ) ===-MIN_VALUE*(2**51  )', () => {
    expect(
      nextUp(-Number.MIN_VALUE*(2**51+1) )
    ).toBe(
             -Number.MIN_VALUE*(2**51  )
    );
  });

  it('nextUp( MIN_VALUE*(2**53-1) ) === MIN_VALUE*(2**53  )', () => {
    expect(
      nextUp( Number.MIN_VALUE*(2**53-1) )
    ).toBe(
              Number.MIN_VALUE*(2**53  )
    );
  });

  it('nextUp( MIN_VALUE* 2**53    ) === MIN_VALUE* 2**53 * (1+ϵ)', () => {
    expect(
      nextUp( Number.MIN_VALUE*2**53 )
    ).toBe(
              Number.MIN_VALUE*2**53 * (1+Number.EPSILON)
    );
  });

  it('nextUp(+MAX_VALUE) === ∞', () => {
    expect( nextUp(Number.MAX_VALUE) ).toBe( Infinity );
  });

  for( const sgn of [-1,+1] )
  {
    const sgn_label = sgn < 0 ? 'negative': 'positive';

    forEachItemIn(
      function*(){
        const TINY = Number.MIN_VALUE * 2**54;

        for( let run=0; run++ < 4*1024*1024; )
          yield sgn * Math.random() * TINY;
      }()
    ).it(`nextUp() works given random, ${sgn_label}, tiny inputs`, x => {
      const X = nextUp(x);
  
      expect(X).toBeGreaterThan(x);
  
      expect([x,X]).toContain(     (x+X)/2 );
      expect([x,X]).toContain( x + (X-x)/2 );
    });

    forEachItemIn(
      function*(){
        for( let run=0; run++ < 4*1024*1024; )
          yield sgn * Math.random() * 1337;
      }()
    ).it(`nextUp() works given random, ${sgn_label} inputs`, x => {
      const X = nextUp(x);
  
      expect(X).toBeGreaterThan(x);
  
      expect([x,X]).toContain(     (x+X)/2 );
      expect([x,X]).toContain( x + (X-x)/2 );
    });
  
    forEachItemIn(
      function*(){
        for( let run=0; run++ < 4*1024*1024; )
          yield sgn * Math.random() * (Number.MAX_VALUE * (1-Number.EPSILON/2));
      }()
    ).it(`nextUp() works given random, ${sgn_label}, huge inputs`, x => {
      const X = nextUp(x);
  
      expect(X).toBeGreaterThan(x);
  
      expect([x,X,sgn*Infinity]).toContain( (x+X)/2 );
      expect([x,X]).toContain( x + (X-x)/2 );
    });
  }
})


describe('nextDown', () => {
  it('nextDown(NaN) is NaN', () => {
    expect( nextDown(NaN) ).toEqual(NaN);
  });

  it('nextDown(+∞) === MAX_VALUE', () => {
    expect( nextDown(+Infinity) ).toBe( +Number.MAX_VALUE );
  });

  it('nextDown( MIN_VALUE * (2**54-1) ) === MIN_VALUE * (2**54-2)', () => {
    expect(
      nextDown( Number.MIN_VALUE * (2**54-1) )
    ).toBe(
                Number.MIN_VALUE * (2**54-2)
    );
  });

  it('nextDown(+MIN_VALUE) === 0', () => {
    expect( nextDown(+Number.MIN_VALUE) ).toBe( 0 );
  });

  it('nextDown(-0) === -MIN_VALUE', () => {
    expect( nextDown(-0.0) ).toBe( -Number.MIN_VALUE );
  })

  it('nextDown(+0) === -MIN_VALUE', () => {
    expect( nextDown(+0.0) ).toBe( -Number.MIN_VALUE );
  });

  it('nextDown(-MIN_VALUE) === -MIN_VALUE*2', () => {
    expect( nextDown(-Number.MIN_VALUE) ).toBe( -Number.MIN_VALUE*2 );
  });

  it('nextDown(-MIN_VALUE*(2**51  ) ) ===-MIN_VALUE*(2**51+1)', () => {
    expect(
      nextDown( -Number.MIN_VALUE*(2**51  ) )
    ).toBe(
                -Number.MIN_VALUE*(2**51+1)
    );
  });

  it('nextDown( MIN_VALUE*(2**51+1) ) === MIN_VALUE*(2**51  )', () => {
    expect(
      nextDown( Number.MIN_VALUE*(2**51+1) )
    ).toBe(
                Number.MIN_VALUE*(2**51  )
    );
  });

  it('nextDown(-MIN_VALUE*(2**53-1) ) ===-MIN_VALUE*(2**53  )', () => {
    expect(
      nextDown(-Number.MIN_VALUE*(2**53-1) )
    ).toBe(
              -Number.MIN_VALUE*(2**53  )
    );
  });

  it('nextDown(-MIN_VALUE*2**53     ) ===-MIN_VALUE* 2**53 * (1+ϵ) )', () => {
    expect(
      nextDown(-Number.MIN_VALUE*2**53 )
    ).toBe(
               -Number.MIN_VALUE*2**53 * (1+Number.EPSILON)
    );
  });

  it('nextDown(-MAX_VALUE) === -∞', () => {
    expect( nextDown(-Number.MAX_VALUE) ).toBe(-Infinity );
  });

  for( const sgn of [-1,+1] )
  {
    const sgn_label = sgn < 0 ? 'negative': 'positive';

    forEachItemIn(
      function*(){
        const TINY = Number.MIN_VALUE * 2**54;

        for( let run=0; run++ < 4*1024*1024; )
          yield sgn * Math.random() * TINY;
      }()
    ).it(`nextDown() works given random, ${sgn_label}, tiny inputs`, x => {
      const X = nextDown(x);
  
      expect(X).toBeLessThan(x);
  
      expect([x,X]).toContain(     (x+X)/2 );
      expect([x,X]).toContain( x + (X-x)/2 );
    });

    forEachItemIn(
      function*(){
        for( let run=0; run++ < 4*1024*1024; )
          yield sgn * Math.random() * 1337;
      }()
    ).it(`nextDown() works given random, ${sgn_label} inputs`, x => {
      const X = nextDown(x);
  
      expect(X).toBeLessThan(x);
  
      expect([x,X]).toContain(     (x+X)/2 );
      expect([x,X]).toContain( x + (X-x)/2 );
    });
  
    forEachItemIn(
      function*(){
        for( let run=0; run++ < 4*1024*1024; )
          yield sgn * Math.random() * (Number.MAX_VALUE * (1-Number.EPSILON/2));
      }()
    ).it(`nextDown() works given random, ${sgn_label}, huge inputs`, x => {
      const X = nextDown(x);
  
      expect(X).toBeLessThan(x);
  
      expect([x,X,sgn*Infinity]).toContain( (x+X)/2 );
      expect([x,X]).toContain( x + (X-x)/2 );
    });
  }
})
