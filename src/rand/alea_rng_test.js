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

import {forEachItemIn, CUSTOM_MATCHERS} from "../jasmine_utils";

import {binary_rangesearch} from '../arrays/binary_search'

import {eye} from "../la/eye";
import {matmul2} from "../la/matmul";

import {AleaRNG} from './alea_rng';
import { SEEDS } from "./alea_rng_test_data";
import {alea as AleaRef, alea} from 'seedrandom';
import {cartesian_prod,
        linspace,
        range,
        repeat,
        zip} from '../iter'


// const PROGRESS = [
//   " ",
//   "▏",
//   "▎",
//   "▍",
//   "▌",
//   "▋",
//   "▊",
//   "▉",
//   "█"
// ];


// function histogram( values, { nBuckets=32, width=56 }={} )
// {
//   const xMin = values.reduce( (x,y) => Math.min(x,y) ),
//         xMax = values.reduce( (x,y) => Math.max(x,y) ),
//         size = (xMax - xMin) / nBuckets;

//   const buckets = new Uint32Array(nBuckets);

//   for( const v of values )
//   {
//     const   i = Math.min( nBuckets-1, Math.floor( (v - xMin) / size ) );
//     buckets[i]++;
//   }

//   const yMax = buckets.reduce( (x,y) => Math.max(x,y) );

//   let hist = '';

//   hist += '╻'.padStart(11) + '\n';

//   for( let i=0; i < nBuckets; i++ )
//   {
//     if( 0 < i )
//       hist += '\n';
//     if( i < nBuckets-1 )
//       hist += ' <' + ( (xMin+size*(i+1)).toFixed(3) ).padStart(7);
//     else
//       hist += ''.padStart(9);
//     hist += '╶╂ '

//     const h = Math.ceil( buckets[i] * 8*width / yMax );

//     hist += ''.padStart(h>>>3, '█');
//     hist += PROGRESS[h%8];
//     hist += buckets[i];
//   }

//   hist += '\n' + '╹'.padStart(11);

//   return hist;
// }


// {
//   const rng = new AleaRNG( Math.random() ),
//         nrm = Float64Array.from({length: 64*1024}, () => rng.normal () ),
//         uni = Float64Array.from({length: 64*1024}, () => rng.uniform() ),
//         int =   Int32Array.from({length: 64*1024}, () => rng.int(-8,+9) );

//   console.log( '\n\nAleaRNG::normal\n---------------\n'       + histogram(nrm) + '\n' );
//   console.log( '\n\nAleaRNG::uniform\n----------------\n'     + histogram(uni) + '\n' );
//   console.log( '\n\nAleaRNG::int(-8,9)\n------------------\n' + histogram(int,{nBuckets: 17}) + '\n' );
// }


describe('AleaRNG', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  });


  forEachItemIn(
    function*(){
      const lens = function*(){
        for(;;)
        for( const len of linspace(18*1024,36*1024, 8) )
          yield len|0;
      }();

      yield* zip(SEEDS, lens);
    }()
  ).it('normal() passes jarque bera test', ([seed,length]) => {
    // https://de.wikipedia.org/wiki/Jarque-Bera-Test
    const rng = new AleaRNG(seed);

    const vals = Float64Array.from({length}, () => rng.normal());

    const E = vals.reduce((E,x) => E +  x      , 0) / length,
          s = vals.reduce((s,x) => s + (x-E)**2, 0) / length,
          S = vals.reduce((S,x) => S + (x-E)**3, 0) / length / s**1.5,
          C = vals.reduce((C,x) => C + (x-E)**4, 0) / length / s**2;

    // (I am not a statistician so it is highly probability that I am talking out of my ass is statistically significant)
    // If we run `n` tests, the probability that we falsly reject our normality hypothesis at least once is:
    //
    //   p(reject|normal) = 1 - (1-α)**n
    //
    // For `n=281` and `α≈0.0001` that should imply a significance of `p≈0.0278`... right?
    const  JB = length/6 * (S*S + (C-3)*(C-3)/4);
    // expect(JB).toBeLessThan( 9.210); // α ≈ 0.01
    // expect(JB).toBeLessThan(10.597); // α ≈ 0.005
    // expect(JB).toBeLessThan(13.816); // α ≈ 0.001
    // expect(JB).toBeLessThan(15.202); // α ≈ 0.0005
    expect(JB).toBeLessThan(18.421); // α ≈ 0.0001
  });


  forEachItemIn(
    function*(){
      const args = function*(){
        for(;;)
        for( const from  of range(-3,+3) )
        for( const until of range(from+2,from+5) )
          yield [from,until];
      }();

      const lens = function*(){
        for(;;)
        for( const len of linspace(128*1024,384*1024, 13) )
          yield len|0;
      }();

      yield *zip(SEEDS, lens, args);
    }()
  ).it('int(from,until) passes χ²-test', ([seed, len, [from,until]]) => {
    // χ²-table; α = 0.001;
    const chiq = [// DF
      10.828,     //  1
      13.816,     //  2
      16.266,     //  3
      18.467,     //  4
      20.515,     //  5
      22.458,     //  6
      24.322,     //  7
      26.124,     //  8
      27.877,     //  9
      29.588,     // 10
      31.264,     // 11
      32.909,     // 12
      34.528,     // 13
      36.123,     // 14
      37.697,     // 15
      39.252      // 16
    ];

    const rng = new AleaRNG(seed);

    const counts = new Uint32Array(until-from);
    for( let i=len; i-- > 0; )
    {
      const nxt = rng.int(from,until) - from;
      if(   nxt%1 !== 0         ) throw new Error('Assertion failed.');
      if( !(nxt   >=  0)        ) throw new Error('Assertion failed.');
      if( !(nxt < counts.length)) throw new Error('Assertion failed.');
      counts[nxt]++;
    }

    const E = len / (until-from);

    let chi2 = 0.0;

    for( let i=counts.length; i-- > 0; )
    { const   d = counts[i] - E; 
      chi2 += d*d / E;
    }

    expect(chi2).toBeLessThan( chiq[counts.length-2] );
  });


  forEachItemIn(
    SEEDS
  ).it("uniform(0,1) matches reference implementation given lorem ipsum sentence seeds", seed => {
    const tst = new AleaRNG(seed),
          ref = new AleaRef(seed);

    for( let i=0; i++ < 8*1024; )
    {
      const t = tst.uniform(0,1),
            r = ref.double();

      expect(t).toBe(r);
    }
  });


  forEachItemIn(
    function*(){
      const lens = function*(){
        for(;;)
        for( const len of linspace(256*1024,512*1024, 7) )
          yield len|0;
      }();

      yield* zip(SEEDS.slice(0,73), lens);
    }()
  ).it('normal() Has σ=1 and E=0', ([seed,len]) => {
    const rng = new AleaRNG(seed);

    let E = 0,
        s = 0;

    for( let i=0; i < len; i++ )
    {
      const x = rng.normal();
      E +=  x;
      s += x*x;
    }

    E = E/len;
    s = s/len - E*E; 

    expect(E).toBeAllCloseTo(0, {atol: 1e-2, rtol: 0});
    expect(s).toBeAllCloseTo(1, {atol: 1e-2, rtol: 0});
  });


  forEachItemIn(
    function*(){
      const args = function*(){
        for(;;)
        for( const mean of linspace(-2.0,+2.0, 7) )
        for( const sigm of linspace(+0.1,+2.0, 3) )
          yield [mean,sigm];
      }();

      const lens = function*(){
        for(;;)
        for( const len of linspace(128*1024,256*1024, 7) )
          yield len|0;
      }();

      yield *zip(SEEDS, lens, args);
    }()
  ).it('normal(mean,sigma) Has σ=sigma and E=mean', ([seed, len, [mean,sigma]]) => {
    const rng = new AleaRNG(seed);

    let E = 0,
        s = 0;

    for( let i=0; i < len; i++ )
    {
      const x = rng.normal(mean,sigma);
      E +=  x;
      s += x*x;
    }

    E = E/len;
    s = s/len - E*E; 

    expect(E).toBeAllCloseTo(mean,    {atol: 1e-2, rtol: 1e-2});
    expect(s).toBeAllCloseTo(sigma**2,{atol: 1e-2, rtol: 1e-2});
  });
});


describe(`AleaRNG::shuffle`, () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  });

  // χ²-table; α = 0.001;
  const chiq = [
    10.828,
    13.816,
    16.266,
    18.467,
    20.515,
    22.458,
    24.322,
    26.124,
    27.877,
    29.588,
    31.264,
    32.909,
    34.528,
    36.123,
    37.697,
    39.252,
    40.790,
    42.312,
    43.820,
    45.315,
    46.797,
    48.268,
    49.728,
    51.179,
    52.620,
    54.052,
    55.476,
    56.892,
    58.301,
    59.703,
    61.098,
    62.487
  ];


  for( let N=1; N++ < 7; )
  {
    const               description = `shuffle(a: Float64Array) passes χ²-test [a.length =${N.toString().padStart(2)}]`,
      rng = new AleaRNG(description);

    it(description, () => {
      expect(N).toBeLessThan(31);
      const  M = 1e6;

      const array = new Float64Array(N),
            freqs = new Float64Array(N*N);

      for( let repeat=0; repeat++ < M; )
      {
        const off = rng.int(-(1<<30), +(1<<30)),
            scale = rng.int(1,1<<30) * (rng.bool() ? -1 : +1);

        // init array
        for( let i=N; i-- > 0; )
          array[i] = off + scale*i;

        rng.shuffle(array);

        let bits = 0;
        for( let i=N; i-- > 0; )
        {
          const j = (array[i] - off) / scale;
          if(   j%1 !== 0  ) throw new Error();
          if( !(j   >=  0) ) throw new Error();
          if( !(j   <   N) ) throw new Error();
          if( bits>>>j & 1 ) throw new Error();
          bits = bits | 1<<j;

          freqs[N*i+j] += 1;
        }
        if( bits !== (1<<N) - 1 )
          throw new Error();
      }

      const E = M/N;

      for( let i=N; i-- > 0; )
      {
        let chi = 0;
        for( let j=N; j-- > 0; ) {
          const   d = freqs[N*i+j] - E;
          chi += d*d / E;
        }
        expect(chi).toBeLessThan(chiq[N-2]);
      }
    });
  }
});


for( const dtype of [[], ['float64'], ['object']] )
describe(`AleaRNG::ortho(${dtype.map(s => `'${s}', `)}...)`, () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  });


  for( const [rndo_name, rndo_method] of [
    [ 'N,N', (rng,N) => rng.ortho(...dtype,N,N) ],
    [ 'N  ', (rng,N) => rng.ortho(...dtype,N  ) ]
  ])
    forEachItemIn(
      function*(){
        const steps_per_binade = 3;

        const sizes = function*(){
          for(;;)
          for( let run=0*steps_per_binade; run <= 6*steps_per_binade; run++ )
            yield Math.round(2**(run/steps_per_binade));
        }();

        yield* zip(repeat(4,SEEDS), sizes);
      }()
    ).it(`ortho(${dtype.map(s => `'${s}',`)}  ${rndo_name}) generates orthogonal square matrices`, ([seed,N]) => {
      const               rng = new AleaRNG(seed);
      let Q = rndo_method(rng,N);

      expect( Q.dtype ).toBe( [...dtype,'float64'][0] );
      expect( Q.shape ).toEqual( Int32Array.of(N,N) );

      if( ! Q.dtype.startsWith('float') )
            Q   = Q.mapElems('float64');
      const Q_T = Q.T,                        I = eye(N);
      expect( matmul2(Q,Q_T) ).toBeAllCloseTo(I);
      expect( matmul2(Q_T,Q) ).toBeAllCloseTo(I);
    });


    forEachItemIn(
      cartesian_prod(
        SEEDS,
        range(1,8),
        range(1,8)
      )
    ).it(`ortho(${dtype.map(s => `'${s}',`)}  M,N) generates orthogonal rectangular matrices`, ([seed,M,N]) => {
      const   rng = new AleaRNG(seed);
      let Q = rng.ortho(...dtype,M,N);

      expect( Q.dtype ).toBe( [...dtype,'float64'][0] );
      expect( Q.shape ).toEqual( Int32Array.of(M,N) );

      if( ! Q.dtype.startsWith('float') )
            Q   = Q.mapElems('float64');
      const Q_T = Q.T,                                     I = eye( Math.min(M,N) );
      if( M >= N ) expect( matmul2(Q_T,Q) ).toBeAllCloseTo(I);
      if( M <= N ) expect( matmul2(Q,Q_T) ).toBeAllCloseTo(I);
    });


  forEachItemIn(
    function*(){
      const steps_per_binade = 3;

      const sizes = function*(){
        for(;;)
        for( let run=0*steps_per_binade; run <= 6*steps_per_binade; run++ )
          yield Math.round( 2**(run/steps_per_binade) );
      }();

      yield* zip(
        SEEDS,
        repeat(range(1,24)),
        sizes
      );
    }()
  ).it(`ortho(${dtype.map(s => `'${s}',`)}M,N,N) generates batches of orthogonal square matrices`, ([seed,M,N]) => {
    const   rng = new AleaRNG(seed);
    let Q = rng.ortho(...dtype,M,N,N);

    expect( Q.dtype ).toBe( [...dtype,'float64'][0] );
    expect( Q.shape ).toEqual( Int32Array.of(M,N,N) );

    if( ! Q.dtype.startsWith('float') )
          Q   = Q.mapElems('float64');
    const Q_T = Q.T,                        I = eye(N);
    expect( matmul2(Q,Q_T) ).toBeAllCloseTo(I);
    expect( matmul2(Q_T,Q) ).toBeAllCloseTo(I);
  });


  forEachItemIn(
    function*(){
      const steps_per_binade = 3;

      const sizes = function*(){
        for(;;)
        for( let M=0*steps_per_binade; M <= 7*steps_per_binade; M++ )
        for( let N=0*steps_per_binade; N <= 7*steps_per_binade; N++ )
          yield [
            Math.round( 2**(M/steps_per_binade) ),
            Math.round( 2**(N/steps_per_binade) )
          ];
      }();

      yield* zip(
        SEEDS,
        repeat(range(1,32)),
        sizes
      );
    }()
  ).it(`ortho(${dtype.map(s => `'${s}',`)}L,M,N) generates batches of orthogonal square matrices`, ([seed,L,[M,N]]) => {
    const   rng = new AleaRNG(seed);
    let Q = rng.ortho(...dtype,L,M,N);

    expect( Q.dtype ).toBe( [...dtype,'float64'][0] );
    expect( Q.shape ).toEqual( Int32Array.of(L,M,N) );

    if( ! Q.dtype.startsWith('float') )
          Q   = Q.mapElems('float64');
    const Q_T = Q.T,                                     I = eye( Math.min(M,N) );
    if( M >= N ) expect( matmul2(Q_T,Q) ).toBeAllCloseTo(I);
    if( M <= N ) expect( matmul2(Q,Q_T) ).toBeAllCloseTo(I);
  });


  for( const [rndo_name, rndo_method] of [
    [ 'N,N', (rng,N) => rng.ortho(...dtype,N,N) ],
    [ 'N  ', (rng,N) => rng.ortho(...dtype,N  ) ]
  ])
    forEachItemIn(
      function*(){
        const nRuns = function*(){
          for(;;)
          for( const nRuns of linspace(24*1024,48*1024, 7) )
            yield Math.round(nRuns);
        }();

        const sizes = repeat( range(1,5) );

        yield* zip(
          SEEDS.slice(0,2),
          sizes,
          nRuns
        );
      }()
    ).it(`ortho(${dtype.map(s => `'${s}',`)}  ${rndo_name}) generates evenly distributed matrices`, ([seed,N,N_RUNS]) => {
      const rng = new AleaRNG(seed),
              E = new Float64Array(N*N),
              s = new Float64Array(N*N);

      for( let run=0; run < N_RUNS; run++ )
      {
        let    Q = rndo_method(rng,N);
        expect(Q.dtype).toBe( [...dtype,'float64'][0] );
        expect(Q.shape).toEqual( Int32Array.of(N,N) );
           Q = Q.data;

        for( let i=N*N; i-- > 0; )
        {
          const    q = Q[i];
          E[i] +=  q;
          s[i] += q*q;
        }
      }

      for( let i=N*N; i-- > 0; ) {
        E[i] = E[i] / N_RUNS;
        s[i] = s[i] / N_RUNS - E[i]*E[i];
      }

      expect(E).toBeAllCloseTo(0,   {atol:1e-2, rtol:0});
      expect(s).toBeAllCloseTo(1/N, {atol:1e-2, rtol:0});
    });


  forEachItemIn(
    function*(){
      const nRuns = function*(){
        for(;;)
        for( const nRuns of linspace(24*1024,48*1024, 7) )
          yield Math.round(nRuns);
      }();

      const sizes = repeat(cartesian_prod(
        range(1,5),
        range(1,5)
      ));

      yield* zip(
        SEEDS.slice(0,2),
        sizes,
        nRuns
      );
    }()
  ).it(`ortho(${dtype.map(s => `'${s}',`)}  M,N) generates evenly distributed matrices`, ([seed, [M,N], N_RUNS]) => {
    const rng = new AleaRNG(seed),
            E = new Float64Array(M*N),
            s = new Float64Array(M*N);

    for( let run=0; run < N_RUNS; run++ )
    {
      let    Q = rng.ortho(...dtype,M,N);
      expect(Q.dtype).toBe( [...dtype,'float64'][0] );
      expect(Q.shape).toEqual( Int32Array.of(M,N) );
         Q = Q.data;

      for( let i=M*N; i-- > 0; )
      {
        const    q = Q[i];
        E[i] +=  q;
        s[i] += q*q;
      }
    }

    for( let i=M*N; i-- > 0; ) {
      E[i] = E[i] / N_RUNS;
      s[i] = s[i] / N_RUNS - E[i]*E[i];
    }

    expect(E).toBeAllCloseTo(0,               {atol:1e-2, rtol:0});
    expect(s).toBeAllCloseTo(1/Math.max(M,N), {atol:1e-2, rtol:0});
  });


  forEachItemIn(
    function*(){
      const nRuns = function*(){
        for(;;)
        for( const nRuns of linspace(64*1024,96*1024, 4) )
          yield Math.round(nRuns);
      }();

      const sizes = repeat( range(1,5) );

      yield* zip(
        SEEDS.slice(0,37),
        nRuns,
        sizes
      );
    }()
  ).it(`ortho(${dtype.map(s => `'${s}',`)}M,N,N) generates evenly distributed batches of matrices`, ([seed,M,N]) => {
    const rng = new AleaRNG(seed),
            E = new Float64Array(N*N),
            s = new Float64Array(N*N);

    let    Q = rng.ortho(...dtype,M,N,N);
    expect(Q.dtype).toBe( [...dtype,'float64'][0] );
    expect(Q.shape).toEqual( Int32Array.of(M,N,N) );
       Q = Q.data;

    for( let Q_off=Q.length; (Q_off-=N*N) >= 0; )
    for( let i=N*N; i-- > 0; )
    {
      const    q = Q[Q_off + i];
      E[i] +=  q;
      s[i] += q*q;
    }

    for( let i=N*N; i-- > 0; ) {
      E[i] = E[i] / M;
      s[i] = s[i] / M - E[i]*E[i];
    }

    expect(E).toBeAllCloseTo(0,   {atol:1e-2, rtol:0});
    expect(s).toBeAllCloseTo(1/N, {atol:1e-2, rtol:0});
  });


  forEachItemIn(
    function*(){
      const nRuns = function*(){
        for(;;)
        for( const nRuns of linspace(64*1024,96*1024, 4) )
          yield Math.round(nRuns);
      }();

      yield* zip(
        SEEDS.slice(0,37),
        nRuns,
        repeat(cartesian_prod(
          range(1,5),
          range(1,5)
        ))
      );
    }()
  ).it(`ortho(${dtype.map(s => `'${s}',`)}L,M,N) generates evenly distributed batches of matrices`, ([seed,L,[M,N]]) => {
    const rng = new AleaRNG(seed),
            E = new Float64Array(M*N),
            s = new Float64Array(M*N);

    let    Q = rng.ortho(...dtype,L,M,N);
    expect(Q.dtype).toBe( [...dtype,'float64'][0] );
    expect(Q.shape).toEqual( Int32Array.of(L,M,N) );
       Q = Q.data;

    for( let Q_off=Q.length; (Q_off-=M*N) >= 0; )
    for( let i=M*N; i-- > 0; )
    {
      const    q = Q[Q_off + i];
      E[i] +=  q;
      s[i] += q*q;
    }

    for( let i=M*N; i-- > 0; ) {
      E[i] = E[i] / L;
      s[i] = s[i] / L - E[i]*E[i];
    }

    expect(E).toBeAllCloseTo(0,               {atol:1e-2, rtol:0});
    expect(s).toBeAllCloseTo(1/Math.max(M,N), {atol:1e-2, rtol:0});
  });
});
