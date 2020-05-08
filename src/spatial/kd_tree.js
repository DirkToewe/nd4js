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

import {_rand_int} from '../_test_data_generators'

import {FrobeniusNorm} from "../la/norm";

import {NAryHeap} from "./_nary_heap";


// TODO:
//   - Add NDArray support


const EUCLIDEAN_DISTANCE = (x,y) =>
{
  const len = x.length;
  if(   len!==y.length )
    throw new Error('Assertion failed.');

  const norm = new FrobeniusNorm();
  for( let i=len; i-- > 0; )
    norm.include(x[i] - y[i]);

  return norm.result;
}


class Branch
{
  constructor( axis, threshold, c0, c1 )
  {
    Object.assign(this,{axis,threshold,c0,c1});
    Object.seal(this);
  }
}


class Leaf
{
  constructor( point )
  {
    this.point = point;
    Object.seal(this);
  }
}


class HeapedBranch
{
  constructor( distance, branch, nearest )
  {
    this.key = distance;
    this.nearest= nearest;
    Object.assign(this,branch);
  }
}


class HeapedLeaf
{
  constructor( distance, leaf )
  {
    this.key    = distance;
    this.point  = leaf.point;
  }
}


export class KDTree
{
  constructor( points )
  {
    points = [...points];
    const ndim = points[0].length;

    if( points.some(pt => pt.length !== ndim) )
      throw new Error('Assertion failed.');

    this.root = function build_tree( axis, start, stop )
    {
      if( !(start >= 0            ) ) throw new Error('Assertion failed.');
      if( !(start <  points.length) ) throw new Error('Assertion failed.');
      if( !(stop  >  start        ) ) throw new Error('Assertion failed.');
      if( !(stop  <= points.length) ) throw new Error('Assertion failed.');
      if( !(axis >= 0  ) ) throw new Error('Assertion failed.');
      if( !(axis < ndim) ) throw new Error('Assertion failed.');

      const swap = (i,j) => {
        if( !(i >= start) ) throw new Error('Assertion failed.');
        if( !(j >= start) ) throw new Error('Assertion failed.');
        if( !(i <  stop ) ) throw new Error('Assertion failed.');
        if( !(j <  stop ) ) throw new Error('Assertion failed.');
        const magnum_pi = points[i];
                          points[i] = points[j];
                                      points[j] = magnum_pi;
      }

      if( 1 === stop-start )
        return new Leaf(points[start]);

      for(;;)
      { // random split very similar to the quick select algorithm
        const threshold = points[ _rand_int(start,stop) ][ axis ];

        // partition using the threshold
        let l=start,
            r=start;
        for( let i=start; i < stop; i++ )
        { const pi = points[i][axis];
          if(   pi <= threshold ) { swap(i,   r);
          if(   pi <  threshold )   swap(l++, r);  r++; }
        }

        if( !(l >= start) ) throw new Error('Assertion failed.');
        if( !(l <= r    ) ) throw new Error('Assertion failed.');
        if( !(r <= stop ) ) throw new Error('Assertion failed.');

//*DEBUG*/        for( let i=start; i <  l   ; i++ ) if( ! (points[i][axis] < threshold) ) throw new Error('Assertion failed.');
//*DEBUG*/        for( let i=r    ; i <  stop; i++ ) if( ! (points[i][axis] > threshold) ) throw new Error('Assertion failed.');
//*DEBUG*/        for( let i=l    ; i <  r   ; i++ ) if( ! (points[i][axis]===threshold) ) throw new Error('Assertion failed.');

        let mid = start+stop >>> 1;
            mid = Math.max(mid, l);
            mid = Math.min(mid, r);

        const ax = (axis+1) % ndim;

        return new Branch(
          axis,
          threshold,
          build_tree(ax, start,mid),
          build_tree(ax,       mid,stop)
        );
      }
    }(0, 0,points.length);
  }

  *nearest_gen( queryPoint )
  {
    const heap = new NAryHeap();

    const enqueue = (nearest, node) => heap.add(
      node instanceof Branch
        ? new HeapedBranch(EUCLIDEAN_DISTANCE(queryPoint,nearest   ), node, nearest)
        : new HeapedLeaf  (EUCLIDEAN_DISTANCE(queryPoint,node.point), node)
    );

    enqueue(queryPoint, this.root);

    while( heap.size > 0 )
    {
      const item = heap.popMin();
      if(   item instanceof HeapedBranch )
      {
        const {nearest, axis, threshold, c0, c1} = item;
        let n0=nearest,
            n1=nearest;
        if(    nearest[axis] > threshold ) { n0 = n0.slice(); n0[axis] = threshold; }
        if(    nearest[axis] < threshold ) { n1 = n1.slice(); n1[axis] = threshold; }
        enqueue(n0,c0);
        enqueue(n1,c1);
      }
      else
        yield item.point;
    }
  }
}
