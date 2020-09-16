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


const rand_code_point = function(){
  const N_CODE_POINTS = ~binary_rangesearch(0,Number.MAX_SAFE_INTEGER, c => {
    try {
      String.fromCodePoint(c);
      return -1;
    }
    catch(err) {
      if( ! (err instanceof RangeError) )
        throw new Error('Assertion failed.');
      return +1;
    }
  });

  return () => String.fromCodePoint( Math.random()*N_CODE_POINTS | 0 );
}();


const SEEDS = `
  Lorem ipsum dolor sit amet, consectetur adipiscing elit. Duis bibendum massa ut ante condimentum finibus.
  Etiam ut diam arcu. Sed mattis vitae sem nec eleifend. Phasellus sed nunc felis. Nunc eleifend arcu orci,
  non aliquam ligula lacinia vitae. Etiam ultrices lacus a turpis bibendum, gravida venenatis libero egestas.
  Pellentesque euismod eu eros in elementum. Integer at auctor ligula. Nullam eu semper lorem. Etiam finibus
  nulla eu vestibulum varius. Quisque dui orci, convallis eu faucibus ut, luctus quis tellus. Aliquam in urna
  est. Sed a nunc risus. Aliquam et sodales augue.
  
  Nulla lobortis convallis quam, sed efficitur lorem. Donec est augue, laoreet sit amet elementum sit amet,
  mollis et lectus. Fusce eget elit mauris. Aliquam rutrum nibh enim, eu condimentum nisl elementum eu. In
  ac sodales lorem. Praesent sit amet risus ante. Proin bibendum sem id iaculis suscipit.
  
  Sed mauris libero, eleifend id volutpat ac, posuere at odio. Proin orci tellus, blandit et pharetra id,
  aliquet eu nunc. Mauris nibh odio, sollicitudin a varius ut, aliquet sit amet quam. Nullam at dolor ex.
  Vivamus imperdiet placerat odio. Cras quis sem dui. Sed fermentum ligula tellus, non sodales diam vulputate
  dapibus. Proin sapien urna, dictum vel lorem a, facilisis gravida justo. Ut vehicula ipsum sit amet ante
  condimentum, pulvinar pellentesque ligula finibus. Mauris iaculis, sem in commodo convallis, quam eros
  viverra arcu, id condimentum tortor arcu ut est. Morbi porta pretium mi eu malesuada. Suspendisse potenti.
  
  Nullam tincidunt augue et ante scelerisque, pretium efficitur nibh scelerisque. Nullam eget aliquet felis,
  in tristique turpis. Sed cursus rutrum nunc, nec accumsan risus. Vivamus ut nibh a ipsum sollicitudin
  suscipit sit amet et urna. Nam ut tortor a tellus lobortis euismod ut eget nisl. Class aptent taciti
  sociosqu ad litora torquent per conubia nostra, per inceptos himenaeos. In dolor nulla, laoreet nec risus
  a, tristique tincidunt ligula. Vivamus imperdiet sagittis augue eget fermentum. Nulla cursus rhoncus nibh
  quis cursus. Nam cursus ipsum id tellus vulputate molestie. Sed aliquet metus vel lacus fringilla, eu
  accumsan mi porta. Morbi tristique interdum orci, id luctus elit aliquet vitae.
  
  Integer fringilla dui sit amet dictum posuere. Aliquam vitae felis aliquet, tincidunt quam at, imperdiet sem.
  Aliquam convallis tempor lectus. Nunc molestie sit amet nunc at vestibulum. Aliquam sed dignissim turpis.
  Suspendisse sit amet tristique quam. Pellentesque non metus dapibus, pharetra libero a, aliquam tortor.
  Fusce a ornare turpis. Morbi vestibulum elit orci, vitae pharetra velit aliquam et. Nam pretium eu dui
  ac lobortis.

  Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla viverra lobortis nulla in finibus.
  Ut sit amet urna tortor. Morbi venenatis augue nisi, vitae convallis augue mattis ac. Integer eget
  sapien consequat, blandit purus ac, pellentesque orci. Donec vel velit ac est rutrum pretium in sed
  nulla. Nam nec neque quis lorem auctor elementum vitae eget libero. Suspendisse ultricies dui
  sollicitudin pretium eleifend. Duis feugiat dapibus sodales. Aliquam volutpat porttitor ante, id
  ornare ligula laoreet nec. Donec vel massa at augue euismod viverra et dignissim erat. Maecenas
  rhoncus at mi a dapibus. Vestibulum eget venenatis elit. Nunc vel facilisis mauris, ac consectetur
  orci. Orci varius natoque penatibus et magnis dis parturient montes, nascetur ridiculus mus. Curabitur
  mattis ante ut ante auctor, non mattis mi elementum.
  
  Donec blandit, risus ut viverra ornare, est felis ullamcorper libero, sed ullamcorper nisi eros non risus.
  Nullam pharetra lobortis ipsum. Cras tempus magna justo, id gravida nisi sollicitudin in. Phasellus a nunc
  nunc. Nullam eu dictum sem. Donec non magna commodo tellus laoreet hendrerit non non justo. Sed tincidunt
  rhoncus convallis. Nulla vitae arcu interdum, sollicitudin lorem sit amet, porta felis. Praesent ac nibh
  eu felis lobortis molestie ut eget massa. Aenean semper quam eget sodales efficitur. Nullam eu convallis
  eros. Suspendisse laoreet cursus pretium. Aliquam eget diam porttitor, dapibus metus quis, congue metus.
  
  Pellentesque non convallis purus. Vivamus interdum magna eu massa vehicula, a porttitor tellus volutpat.
  Aenean eget dapibus dolor. Fusce sit amet ligula nunc. Nam condimentum sit amet mauris ut pulvinar. Nullam
  tempus est urna, et dapibus dui congue nec. Quisque vitae pharetra ligula. Quisque nunc velit, condimentum
  eget nisl quis, molestie iaculis tellus. Nunc sed dignissim orci, a tempor nulla. Phasellus aliquet commodo
  sapien, fringilla tincidunt lectus eleifend ac. Fusce vitae ullamcorper mauris, non commodo quam. Curabitur
  tincidunt in enim vel faucibus. Aliquam erat volutpat. Cras luctus mi nec commodo sollicitudin. Pellentesque
  in nulla ut neque bibendum pulvinar. Morbi lacinia varius enim quis sollicitudin.
  
  Nulla facilisi. Nam vel nulla dolor. Sed placerat auctor auctor. Quisque nibh tortor, molestie vulputate
  venenatis nec, laoreet at enim. Pellentesque vitae quam tempus, posuere urna ac, tempor elit. Aenean turpis
  nisi, tempus quis eros in, eleifend volutpat libero. Sed venenatis dignissim sem, nec auctor sem convallis
  lobortis. Etiam vitae dolor laoreet, consequat felis sit amet, molestie nisl. Suspendisse lorem est, consequat id
  condimentum in, tristique in ligula. Quisque tincidunt arcu eu massa consequat, a consectetur sem eleifend. Maecenas
  turpis massa, mollis in enim nec, accumsan venenatis sapien. Donec sodales vestibulum nibh, vel vulputate dolor eleifend ac.
  
  Aliquam auctor neque ullamcorper felis rutrum, sit amet bibendum metus vehicula. Phasellus sagittis convallis ex.
  Quisque imperdiet lobortis est, sit amet tempor diam feugiat non. Nullam non tristique odio. Fusce bibendum dolor
  tellus, accumsan ultrices ipsum molestie eu. Maecenas fringilla justo a ante accumsan, vel consectetur ligula rutrum.
  Morbi consectetur libero eget feugiat viverra. Donec quis porttitor quam. Vestibulum auctor eros vel mauris bibendum,
  at consectetur ex congue. Integer sagittis tincidunt velit, faucibus imperdiet elit pellentesque in. Curabitur egestas
  facilisis imperdiet. In hac habitasse platea dictumst. Etiam in maximus mauris. Sed vel diam tellus. Nunc gravida, dolor
  at pretium scelerisque, diam tortor finibus mi, quis hendrerit nunc mauris nec arcu. Suspendisse dapibus convallis lectus
  vitae mollis.
  
  In sollicitudin vulputate nibh, ac placerat turpis tempus sed. Nullam faucibus, sapien a venenatis sollicitudin, justo
  orci ultricies dolor, ac dapibus nunc tellus eu libero. Donec interdum mi quis purus ullamcorper, nec tempor orci tempor.
  Vivamus ex lacus, finibus ac ante quis, placerat elementum quam. Duis non tristique nulla. Suspendisse faucibus, lectus sed
  fringilla malesuada, augue dolor rhoncus risus, sit amet faucibus velit mauris in enim. Duis nec convallis elit.
  
  In consequat turpis a aliquet euismod. Mauris hendrerit augue metus, non suscipit nulla egestas at. Phasellus quis venenatis
  lorem. Vestibulum facilisis a sapien in fermentum. Cras imperdiet, turpis non molestie sollicitudin, turpis nulla aliquam ex,
  quis dapibus ante felis non ipsum. Nam tempor ipsum vel quam imperdiet dictum. Pellentesque augue purus, ultricies nec velit at,
  hendrerit varius augue. Etiam elit eros, viverra ut interdum quis, posuere non dolor. Cras at mattis nisl. Pellentesque augue
  lacus, scelerisque pharetra justo non, rutrum venenatis arcu. Suspendisse at ante sollicitudin, mattis tortor eu, tempus lorem.
  Nam auctor nunc vel congue posuere. Integer finibus maximus dolor in semper.
  
  Nam facilisis ante quis massa gravida, ut semper eros molestie. Morbi id erat vel lorem interdum interdum a ut ante. Etiam
  venenatis, tortor quis tristique bibendum, ex dui porta lectus, convallis rutrum velit quam nec enim. Sed pulvinar arcu urna,
  at mollis libero tincidunt ut. Nullam quis ipsum et velit gravida cursus eu at nulla. Nunc eu libero neque. Aliquam varius,
  sapien vitae porta placerat, enim nisi fringilla diam, in porta magna justo auctor turpis. Nunc vel ipsum ut tortor luctus luctus
  in a mauris. Curabitur elementum facilisis magna et vestibulum. Pellentesque auctor ac risus sed mollis. Duis dapibus egestas gravida.
  
  Pellentesque varius pretium libero, eu blandit leo ultrices eget. Vestibulum ante ipsum primis in faucibus orci luctus et ultrices
  posuere cubilia curae; Quisque dui leo, condimentum ut lacus in, scelerisque pretium erat. Duis eu commodo nisi. Vestibulum pretium,
  eros vitae facilisis ultrices, felis ante imperdiet est, at varius enim felis a massa. Proin tincidunt feugiat nisi, quis rhoncus
  arcu ornare sit amet. Morbi varius elit et ligula viverra venenatis. Nulla tincidunt in nunc ac semper. Phasellus rutrum venenatis
  lacus, vel pulvinar sapien aliquet scelerisque. Etiam tempus neque nulla, eu dictum magna luctus id.
  
  Nulla imperdiet, arcu sit amet rutrum fringilla, nibh lorem scelerisque ligula, dapibus tristique nibh ipsum vitae enim. Curabitur
  pellentesque rhoncus ligula, non cursus urna porta vitae. Vivamus pellentesque sagittis tincidunt. Duis vel aliquet ex, in cursus
  justo. Nullam vel elit ligula. Nulla nec nulla finibus, ultrices sem ac, sagittis augue. Pellentesque porttitor tempor neque eu convallis.
  
  Nulla a accumsan urna. Nullam mattis gravida leo, a viverra odio venenatis ut. Sed egestas laoreet enim nec porttitor. Fusce sem
  justo, finibus tincidunt magna et, feugiat semper nisi. Sed nec felis at elit rhoncus luctus ut sed neque. Suspendisse non risus
  sodales, porttitor dui a, ultrices odio. Nunc at sem orci. Nunc ut lorem justo. Morbi est nibh, porta vitae blandit sit amet,
  tincidunt id augue. Curabitur porttitor nunc magna, quis consectetur mi volutpat a.
  
  In hac habitasse platea dictumst. Mauris semper accumsan mi eget aliquet. Mauris tempus justo sit amet felis placerat vulputate.
  Suspendisse potenti. Aenean dignissim quam massa, a consectetur magna rhoncus ut. Duis egestas leo at tortor efficitur pretium. In
  aliquet odio eu libero vehicula finibus. Quisque rutrum dolor nec risus dapibus, et faucibus felis malesuada. Nunc mollis purus sit
  amet pulvinar dignissim. Ut pharetra et risus eu facilisis. Phasellus ligula enim, facilisis malesuada lorem vel, facilisis consequat
  nisl. Morbi lobortis at neque ut euismod.
  
  In lectus tellus, rutrum in mi sit amet, vestibulum feugiat sem. Curabitur et ligula et nunc tristique tempus. Maecenas sapien est,
  euismod sit amet velit nec, placerat maximus erat. Pellentesque rhoncus ultricies tellus, at rutrum justo volutpat nec. Etiam non leo
  vel libero luctus condimentum sed quis nulla. Duis aliquet lacus in erat lacinia ornare. Maecenas id purus a metus gravida dapibus.
  Pellentesque in lacus in libero maximus lacinia vel at leo. Pellentesque ut elit sed lacus commodo tempor.
  
  Nullam ullamcorper nunc turpis, cursus congue orci facilisis vitae. Maecenas facilisis odio sed eros tempor pellentesque. In pretium
  turpis a diam ultrices, et dictum felis varius. Donec ultrices auctor dolor eu viverra. Interdum et malesuada fames ac ante ipsum
  primis in faucibus. Phasellus interdum, elit at luctus pellentesque, purus neque sagittis turpis, eget convallis orci sem finibus erat.
  Ut metus ipsum, dapibus sit amet sem sed, congue pharetra elit. Donec nec dapibus orci, a luctus urna. Vivamus sed pulvinar ligula.
  Donec blandit nunc felis, vel fringilla eros malesuada et. Donec non lorem facilisis, congue ante porttitor, aliquet risus. Phasellus
  maximus lacus sed mauris mollis euismod. Morbi pulvinar diam nisl, at lacinia purus cursus cursus. Donec et ante quam. Etiam tristique
  sed nisi vel vehicula.
  
  Nulla non erat vel erat commodo auctor aliquet non tortor. Quisque rhoncus tempus metus et commodo. Aenean sollicitudin placerat bibendum.
  Nulla at elit iaculis, gravida urna egestas, malesuada diam. Donec gravida erat at ipsum blandit, in porttitor enim sagittis. Suspendisse
  quis mauris odio. Quisque et facilisis arcu. Mauris eget lobortis turpis. Maecenas sed pharetra tellus. In imperdiet nulla velit, sit amet
  sollicitudin arcu euismod a. Sed sit amet dui non sapien vehicula sollicitudin. Morbi malesuada libero arcu, sit amet porttitor purus
  lobortis vitae. Vestibulum in erat id augue ullamcorper auctor.
  
  Sed posuere eros id lacus tempor tincidunt. Vestibulum ultrices laoreet massa, molestie congue neque lacinia quis. Nullam consectetur lacus
  at mauris convallis, efficitur condimentum massa porttitor. Vestibulum sed nibh id felis facilisis feugiat sit amet vitae libero. Proin varius
  leo velit, vitae suscipit ante facilisis id. Integer sodales libero accumsan erat tempus, at placerat dui tristique. Aliquam eu congue sem.
  Maecenas erat mi, ornare eu justo lacinia, dapibus luctus orci.
  
  Praesent vehicula ex tortor, nec condimentum augue euismod sit amet. Duis sagittis metus est, auctor feugiat nulla pretium in. Nunc sit amet
  viverra lacus. In hac habitasse platea dictumst. Fusce lacus nunc, cursus accumsan enim ut, ultrices semper augue. Nunc id libero efficitur,
  feugiat urna ut, vestibulum nulla. Phasellus erat purus, pharetra eget laoreet quis, laoreet eget turpis. Nunc sit amet nunc in purus ultricies
  mollis a id enim.
  
  Curabitur id nunc ac elit gravida accumsan vel ut sem. Nam non justo nisi. Cras ac placerat magna. Cras ac purus orci. Duis tempus, sem id
  volutpat pulvinar, est tortor pretium ipsum, vitae hendrerit metus lectus at mauris. Etiam ut ligula risus. Aliquam a ante sit amet erat
  viverra commodo.
  
  Vivamus iaculis dolor dui, et luctus est efficitur ut. Mauris tristique massa id magna dapibus sagittis. Curabitur ante erat, ornare eget
  efficitur nec, feugiat in ligula. Fusce egestas orci sed molestie pharetra. Donec auctor imperdiet odio, vitae efficitur mauris rhoncus blandit.
  Aenean efficitur ligula orci, quis consequat ante gravida molestie. Donec in tempor nulla.
  
  Phasellus gravida faucibus turpis, in rutrum mauris consequat sed. Aenean ex tortor, consectetur ac arcu vel, molestie feugiat mauris. Quisque
  vel mi in risus tempor tempor. Mauris quis velit et enim ultricies semper vel quis quam. Phasellus varius hendrerit nibh. Morbi scelerisque velit
  a magna tempus viverra. Donec et enim sollicitudin, auctor purus vel, elementum urna. Donec sit amet orci quis lacus aliquet feugiat in et odio.
  Donec tristique mi sit amet ex cursus, quis porttitor odio pharetra. Donec vitae justo neque. Cras risus augue, dapibus a imperdiet vitae, fermentum
  id urna. Curabitur pharetra pharetra lorem vitae placerat. Donec scelerisque, sem vitae mattis pulvinar, odio justo tempus nisi, ac lobortis nisl
  turpis nec dui. Ut et metus venenatis, cursus magna eu, sagittis ligula. Duis et sapien vitae turpis auctor aliquam vel a erat. Quisque a lectus leo.
`.split('.')
 .map( s => s.trim().replace(/\s+/g,' ') )
 .filter( s => s.length !== 0 );


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
