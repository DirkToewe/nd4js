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

import {pyon_parse} from './pyon'
import {forEachItemIn} from '../jasmine_utils'


describe('pyon', () => {

  beforeEach( () => {
    jasmine.addCustomEqualityTester(
      (x,y) => {
        if( x === 0 )
          return y === 0;
      }
    );
  })

  forEachItemIn([
    [ '-0',  -0 ],
    [ '+0',  +0 ],
    [  '1',   1 ],
    [ '-2',  -2 ],
    [ '+3',   3 ],
    [ '.4',  .4 ],
    ['-.5', -.5 ],
    ['+.6',  .6 ],
    [  '7.',  7.],
    [ '-8.', -8.],
    [ '+9.',  9.],
    ['False', false],
    ['True',  true ],
    ['{}', {}],
    [`[1,"2",3.4]`, [1,"2",3.4]],
    ['{ "x": 1, "y": [1,2] }', { "x": 1, "y": [1,2] }]
  ]).it(`pyon_parse_gen works on hand-crafted examples`, ([str,obj]) => {
    expect( pyon_parse(str) ).toEqual(obj);
  })

  forEachItemIn(
    function*(){
      function rs()
      {
        const N = Math.trunc(Math.random()*3);
        let result = '';
        for( let i=0; i < N; i++ )
          switch( Math.trunc(Math.random()*3) )
          {
            default: throw new Error('Assertion failed.');
            case 0 : result +=  ' '; break;
            case 1 : result += '\n'; break;
            case 2 : result += '\t'; break;
          }
        return result;
      }

      function* examples( depth )
      {
        switch( Math.trunc(Math.random()*5) )
        {
          default: throw new Error('Assertion failed.');
                 case 0: { const randInt  = Math.trunc(Math.random()*2e3 - 1e3);                           yield [`${rs() }${randInt  }${ rs()}`, randInt  ] }
          break; case 1: { const randFloat=            Math.random()*2e3 - 1e3 ;                           yield [`${rs() }${randFloat}${ rs()}`, randFloat] }
          break; case 2: { const randStr  = "hhellow\uAAAA,x11!?"              ; if( Math.random() < 0.5 ) yield [`${rs()}"${randStr  }"${rs()}`, randStr  ]
                                                                                 else                      yield [`${rs()}'${randStr  }'${rs()}`, randStr  ] } // <- totally random string...
          break; case 3: yield ["False", false]
          break; case 4: yield ["True",  true ]
        }

        if( depth > 0 )
        {
          for( const [str1,obj1] of examples(depth-1) )
          {
            switch( Math.trunc(Math.random()*4) )
            {
              default: throw new Error('Assertion failed.');
                     case 0: yield [`${rs()}[${str1}]${rs()}`,         [obj1]]
              break; case 1: yield [`${rs()}[${str1},${rs()}]${rs()}`, [obj1]]
              break; case 2: yield [`${rs()}(${str1})${rs()}`,         [obj1]]
              break; case 3: yield [`${rs()}(${str1},${rs()})${rs()}`, [obj1]]
            }

            for( const [str2,obj2] of examples(depth-1) )
            {
              switch( Math.trunc(Math.random()*4) )
              {
                default: throw new Error('Assertion failed.');
                       case 0: yield [`${rs()}[${str1},${str2}]${rs()}`,         [obj1, obj2]]
                break; case 1: yield [`${rs()}[${str1},${str2},${rs()}]${rs()}`, [obj1, obj2]]
                break; case 2: yield [`${rs()}(${str1},${str2})${rs()}`,         [obj1, obj2]]
                break; case 3: yield [`${rs()}(${str1},${str2},${rs()})${rs()}`, [obj1, obj2]]
              }
              if( Math.random() < 0.5 ) yield [`${rs()}{${str1}:${str2}}${rs()}`,        {[obj1]:obj2}]
              else                      yield [`${rs()}{${str1}:${str2},${rs()}}${rs()}`,{[obj1]:obj2}]

              for( const [str3,obj3] of examples(depth-1) )
              {
                switch( Math.trunc(Math.random()*4) )
                {
                  default: throw new Error('Assertion failed.');
                         case 0: yield [`${rs()}[${str1},${str2},${str3}]${rs()}`,         [obj1,obj2,obj3]]
                  break; case 1: yield [`${rs()}[${str1},${str2},${str3},${rs()}]${rs()}`, [obj1,obj2,obj3]]
                  break; case 2: yield [`${rs()}(${str1},${str2},${str3})${rs()}`,         [obj1,obj2,obj3]]
                  break; case 3: yield [`${rs()}(${str1},${str2},${str3},${rs()})${rs()}`, [obj1,obj2,obj3]]
                }

                for( const [str4,obj4] of examples(depth-1) )
                {
                  switch( Math.trunc(Math.random()*4) )
                  {
                    default: throw new Error('Assertion failed.');
                           case 0: yield [`${rs()}[${str1},${str2},${str3},${str4}]${rs()}`,         [ obj1, obj2, obj3, obj4]]
                    break; case 1: yield [`${rs()}[${str1},${str2},${str3},${str4},${rs()}]${rs()}`, [ obj1, obj2, obj3, obj4]]
                    break; case 2: yield [`${rs()}(${str1},${str2},${str3},${str4})${rs()}`,         [ obj1, obj2, obj3, obj4]]
                    break; case 3: yield [`${rs()}(${str1},${str2},${str3},${str4},${rs()})${rs()}`, [ obj1, obj2, obj3, obj4]]
                  }
                  if( Math.random() < 0.5 ) yield [`${rs()}{${str1}:${str2},${str3}:${str4}}${rs()}`,         {[obj1]:obj2,[obj3]:obj4}]
                  else                      yield [`${rs()}{${str1}:${str2},${str3}:${str4},${rs()}}${rs()}`, {[obj1]:obj2,[obj3]:obj4}]
                }
              }
            }
          }
        }
      }

      for( let i=0; i < 8; i++ )
        yield* examples(2);
    }()
  ).it(`pyon_parse works on generated examples`, ([str,obj]) => {
//    console.log(str, obj)
    expect( pyon_parse(str) ).toEqual(obj);
  })

})
