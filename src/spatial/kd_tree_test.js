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
import {_rand_int, _shuffle} from '../_test_data_generators'

import {compare as compare_arrays} from "../arrays/comparator";

import {FrobeniusNorm} from "../la/norm";

import {KDTree} from "./kd_tree";


describe('KDTree', () => {

  forEachItemIn(
    function*(){
      for( let run=0; run++ < 173; )
      {
        const ndim = _rand_int(1,4),
              size = _rand_int(1,337);

        const points =  Array.from({length: size}, () =>
                        Array.from({length: ndim}, () => _rand_int(-1337, +1337) ));

        for( let query=0; query++ < 73; )
        {
          _shuffle(points);
          const queryPoint = Array.from( {length: ndim}, () => _rand_int(-1337, +1337) );
          yield Object.freeze([ queryPoint, points.slice() ]);
        }
      }
    }()
  ).it('works given random examples', ([queryPoint,points]) => {
    const distance = pt => {
      const len = pt.length;
      if( queryPoint.length !== len ) throw new Error('Assertion failed.');

      const norm = new FrobeniusNorm();
      for( let i=len; i--> 0; )
        norm.include( queryPoint[i] - pt[i] );

      return norm.result;
    };

    const tree = new KDTree(points);

    for( let repeat=0; repeat++ < 3; )
    {
      const sorted = Object.freeze([...tree.nearest_gen(queryPoint)]);

      for( let i=sorted.length; --i > 0; )
        expect                ( distance(sorted[i-1]) )
          .toBeLessThanOrEqual( distance(sorted[i  ]) );

      expect    ( sorted.slice().sort(compare_arrays) )
        .toEqual( points.slice().sort(compare_arrays) );
    }
  })

})
