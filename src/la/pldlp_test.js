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
import {NDArray} from '../nd_array'
import {tabulate} from '../tabulate'
import {_shuffle} from '../_test_data_generators'
import {zip_elems} from '../zip_elems'

import {diag,
        diag_mat} from './diag'
import {permute_rows,
        permute_cols,
      unpermute_rows,
      unpermute_cols} from './permute'
import {pldlp_decomp,
        pldlp_solve,
        pldlp_l,
        pldlp_d,
        pldlp_p} from './pldlp'
import {matmul,
        matmul2} from './matmul'

import {pldlp_decomp_test_data} from './pldlp_test_data'


describe('PLDLᵀPᵀ Decomposition (Bunch-Kaufman LDLᵀ)', () => {
  beforeEach( () => {
    jasmine.addMatchers(CUSTOM_MATCHERS)
  })


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from

      for( let run=1024; run-- > 0; )
      {
        let ndim = randInt(0,5),
          shapes = Array.from({length: ndim}, () => randInt(1,8));
        shapes = [
          shapes,
          shapes.slice( randInt(0,ndim) )
        ];
        _shuffle(shapes);

        for( let d=ndim; d > 0; d-- )
        for( let i=randInt(0,2); i-- > 0; ) {
          const     shape = shapes[randInt(0,2)],
                j = shape.length - d
          if(0<=j)  shape[j] = 1
        }

        const N = randInt(1,24);
        shapes[0].push(N,N);
        shapes[1].push(N  );

        const P = tabulate(shapes[1],  'int32', (...idx) => idx.pop() ),
             LD = tabulate(shapes[0],'float64', (...indices) => {
               const [i,j] = indices.slice(-2);
               return i===j ?(Math.random()*1.6 + 0.4) * (Math.random() < 0.5 ? +1 : -1)
                    : i < j ? 0
                    :         Math.random()*2.4 - 1.2;
             });

        for( let P_off=0; P_off < P.data.length; P_off+=N )
        {
          const p = P.data.subarray(P_off,P_off+N);
          _shuffle(p);
          for( let i=1; i < N; i++ )
            p[i++] ^= -(Math.random() < 0.25);
        }

        yield [LD,P]
      }
    }()
  ).it('pldlp_l works on random examples', ([LD,P]) => {
    const Pi = P.reshape(...P.shape, 1);

    const L = zip_elems([LD,Pi], 'float64', (L,Pi,...idx) => {
      const j = idx.pop(),
            i = idx.pop();
      if(           j===i  ) return 1;
      if( Pi < 0 && j===i-1) return 0;
                    return L;
    });

    const l = pldlp_l(LD,P);
    
    expect(l).toBeAllCloseTo(L);
  });


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from

      for( let run=1024; run-- > 0; )
      {
        let ndim = randInt(0,5),
          shapes = Array.from({length: ndim}, () => randInt(1,8));
        shapes = [
          shapes,
          shapes.slice( randInt(0,ndim) )
        ];
        _shuffle(shapes);

        for( let d=ndim; d > 0; d-- )
        for( let i=randInt(0,2); i-- > 0; ) {
          const     shape = shapes[randInt(0,2)],
                j = shape.length - d
          if(0<=j)  shape[j] = 1
        }

        const N = randInt(1,24);
        shapes[0].push(N,N);
        shapes[1].push(N  );

        const P = tabulate(shapes[1],  'int32', (...idx) => idx.pop() ),
             LD = tabulate(shapes[0],'float64', (...indices) => {
               const [i,j] = indices.slice(-2);
               return i===j ?(Math.random()*1.6 + 0.4) * (Math.random() < 0.5 ? +1 : -1)
                    : i < j ? 0
                    :         Math.random()*2.4 - 1.2;
             });

        for( let P_off=0; P_off < P.data.length; P_off+=N )
        {
          const p = P.data.subarray(P_off,P_off+N);
          _shuffle(p);
          for( let i=1; i < N; i++ )
            p[i++] ^= -(Math.random() < 0.25);
        }

        yield [LD,P]
      }
    }()
  ).it('pldlp_d works on random examples', ([LD,P]) => {
    const Pi = P.reshape(...P.shape, 1);

    const D = zip_elems([LD,LD.T, Pi,Pi.T], 'float64', (L_ij,L_ji, P_i,P_j, ...idx) => {
      const j = idx.pop(),
            i = idx.pop();
      if(            j===i  ) return L_ij;
      if( P_i < 0 && j===i-1) return L_ij;
      if( P_j < 0 && i===j-1) return L_ji;
                              return 0;
    });

    const d = pldlp_d(LD,P);
    
    expect(d).toBeAllCloseTo(D);
  });


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from

      for( let run=4096; run-- > 0; )
      {
        let ndim = randInt(0,5),
          shapes = Array.from({length: ndim}, () => randInt(1,8));
        shapes = [
          shapes,
          shapes.slice( randInt(0,ndim) )
        ];
        _shuffle(shapes);

        for( let d=ndim; d > 0; d-- )
        for( let i=randInt(0,2); i-- > 0; ) {
          const     shape = shapes[randInt(0,2)],
                j = shape.length - d
          if(0<=j)  shape[j] = 1
        }

        const N = randInt(1,24);
        shapes[0].push(N,N);
        shapes[1].push(N  );

        const P = tabulate(shapes[1],  'int32', (...idx) => idx.pop() ),
             LD = tabulate(shapes[0],'float64', (...indices) => {
               const [i,j] = indices.slice(-2);
               return i===j ?(Math.random()*1.6 + 0.4) * (Math.random() < 0.5 ? +1 : -1)
                    : i < j ? 0
                    :         Math.random()*2.4 - 1.2;
             });

        for( let P_off=0; P_off < P.data.length; P_off+=N )
        {
          const p = P.data.subarray(P_off,P_off+N);
          _shuffle(p);
          for( let i=1; i < N; i++ )
            p[i++] ^= -(Math.random() < 0.25);
        }

        yield [LD,P]
      }
    }()
  ).it('pldlp_p works on random examples', ([LD,P]) => {
    const Q = P.mapElems( p => p ^ -(p<0) );

    const  p = pldlp_p(LD,P);
    expect(p).toBeAllCloseTo(Q);
  });


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from

      for( let run=4096; run-- > 0; )
      {
        let ndim = randInt(0,4),
          shapes = Array.from({length: ndim}, () => randInt(1,8));
        shapes = [
          shapes,
          shapes.slice( randInt(0,ndim) ),
          shapes.slice( randInt(0,ndim) )
        ];
        _shuffle(shapes);

        for( let d=ndim; d > 0; d-- )
        for( let i=randInt(0,2); i-- > 0; ) {
          const    shape = shapes[randInt(0,2)],
               j = shape.length - d
          if(0<=j) shape[j] = 1
        }

        const M = randInt(1,24),
              N = randInt(1,24);
        shapes[0].push(M,M);
        shapes[1].push(M  );
        shapes[2].push(M,N);

        const y = tabulate(shapes[2],'float64', () => Math.random()*2-1),
              P = tabulate(shapes[1],  'int32', (...idx) => idx.pop() ),
              LD= tabulate(shapes[0],'float64', (...indices) => {
                const [i,j] = indices.slice(-2);
                return i===j ?(Math.random()*1.6 + 0.4) * (Math.random() < 0.5 ? +1 : -1)
                     : i < j ? 0
                     :         Math.random()*2.4 - 1.2;
              });

        for( let P_off=0; P_off < P.data.length; P_off += M )
        {
          const p = P.data.subarray(P_off,P_off+M);
          _shuffle(p);
          for( let i=1; i < M; i++ )
            p[i++] ^= -(Math.random() < 0.25);
        }

        yield [LD,P,y]
      }
    }()
  ).it('pldlp_solve works on random examples', ([LD,P,Y]) => {
    const Pi = P.reshape(...P.shape, 1);

    const L = pldlp_l(LD,P),
          D = pldlp_d(LD,P),
          Q = pldlp_p(LD,P);
    
    const X = pldlp_solve(LD,P,Y),
          y = unpermute_rows( matmul(L,D,L.T, permute_rows(X,Q)), Q);
    
    expect(y).toBeAllCloseTo(Y, {atol: 1e-6});
  });


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from

      for( let run=512; run-- > 0; )
      {
        let ndim = randInt(2,6),
          shapes = [ Array.from({length: ndim}, () => randInt(1,8)) ]
        shapes.splice( randInt(0,2), 0, shapes[0].slice( randInt(0,ndim) ) )

        for( let d=ndim; d > 0; d-- )
        for( let i=randInt(0,2); i-- > 0; ) {
          const    shape = shapes[randInt(0,2)],
               j = shape.length - d
          if(0<=j) shape[j] = 1
        }

        const M = randInt(1,24); shapes[0].push(M,M)
        const N = randInt(1,24); shapes[1].push(M,N)

        const y = tabulate(shapes[1],'float64', () => Math.random()*8-4),
              A = tabulate(shapes[0],'float64', () => Math.random()*8-4),
              S = zip_elems([A,A.T], 'float64', (a,aT) => (a+aT) / 2);
        Object.freeze(S.data.buffer)
        Object.freeze(S)
        Object.freeze(y.data.buffer)
        Object.freeze(y)

        yield [S,y]
      }
    }()
  ).it('pldlp_decomp + pldlp_solve works on random examples', ([S,Y]) => {
    const [LD,P] = pldlp_decomp(S),
              X  = pldlp_solve(LD,P,Y),
              y  = matmul2(S,X);

    expect(y).toBeAllCloseTo(Y)
  });


  forEachItemIn(
    pldlp_decomp_test_data()
  ).it('pldlp_decomp works on generated examples', ALDP => {
    const [A,L,D,P] = ALDP.map( x => {
      Object.freeze(x);
      Object.freeze(x.data.buffer);
      return x;
    });
    if(  ! A.dtype.startsWith('float') ) throw new Error('Assertion failed.');
    if(  ! L.dtype.startsWith('float') ) throw new Error('Assertion failed.');
    if(  ! D.dtype.startsWith('float') ) throw new Error('Assertion failed.');
    expect(P.dtype).toBe('int32');

    expect(A.shape[A.ndim-2])
     .toBe(A.shape[A.ndim-1]);

    expect(L).toBeLowerTriangular();
    expect(D).toBeTridiagonal();

    const [ld,q]= pldlp_decomp(A),
              l = pldlp_l(ld,q),
              d = pldlp_d(ld,q),
              p = pldlp_p(ld,q);

    expect(q.dtype).toEqual('int32');
    expect(p).toBeAllCloseTo(P, {rtol:0, atol:0});

    expect(ld.dtype).toBe(A.dtype);
    expect(ld.shape).toEqual(A.shape);
    expect(ld).toBeLowerTriangular();

    expect(d.dtype).toBe(A.dtype);
    expect(d).toBeTridiagonal();
    expect(d).toBeAllCloseTo(D);

    expect(l.dtype).toBe(A.dtype);
    expect(l).toBeLowerTriangular();
    expect(l).toBeAllCloseTo(L);
  });


  forEachItemIn(
    function*(){
      const randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from

      for( let run=1024; run-- > 0; )
      {
        const shape = Int32Array.from({ length: randInt(2,6) }, () => randInt(1,8) );
              shape[shape.length-2] =
              shape[shape.length-1] = randInt(1,32);

        const A = tabulate(shape,'float64',(...indices) => Math.random()*8 - 4),
              S = zip_elems([A,A.T], 'float64', (a,aT) => (a+aT) / 2);
        Object.freeze(S);
        Object.freeze(S.data.buffer);

        yield S;
      }
    }()
  ).it('pldlp_decomp works on random examples', S => {
    expect(S.shape[S.ndim-2])
     .toBe(S.shape[S.ndim-1]);

    const [LD,P]= pldlp_decomp(S),
              L = pldlp_l(LD,P),
              D = pldlp_d(LD,P),
              p = pldlp_p(LD,P);

    expect(P.shape).toEqual( S.shape.slice(0,-1) );

    expect(LD.shape).toEqual(S.shape);
    expect(LD).toBeLowerTriangular();

    expect(D.shape).toEqual(S.shape);
    expect(D).toBeTridiagonal();

    expect(L.shape).toEqual(S.shape);
    expect(L).toBeLowerTriangular();

    let s = matmul(L, D, L.T);
        s = unpermute_rows(s,p);
        s = unpermute_cols(s,p);

    expect(s).toBeAllCloseTo(S);
  });


  forEachItemIn(
    function*(){
      for( let N=3; N < 1024; N = N*1.4 | 0 )
      {
        const A = new NDArray(
            Int32Array.of(N,N),
          Float64Array.from({length: N*N}, () => Math.random()*8 - 4)
        );
        const         S = zip_elems([A,A.T], 'float64', (a,aT) => (a+aT) / 2);
        Object.freeze(S.data.buffer)
        Object.freeze(S)
        yield         S
      }
    }()
  ).it('pldlp_decomp works on random large examples', S => {
    expect(S.shape[S.ndim-2])
     .toBe(S.shape[S.ndim-1]);

    const [LD,P]= pldlp_decomp(S),
              L = pldlp_l(LD,P),
              D = pldlp_d(LD,P),
              p = pldlp_p(LD,P);

    expect(P.shape).toEqual( S.shape.slice(0,-1) );

    expect(LD.shape).toEqual(S.shape);
    expect(LD).toBeLowerTriangular();

    expect(D.shape).toEqual(S.shape);
    expect(D).toBeTridiagonal();

    expect(L.shape).toEqual(S.shape);
    expect(L).toBeLowerTriangular();

    let s = matmul(L, D, L.T);
        s = unpermute_rows(s,p);
        s = unpermute_cols(s,p);

    expect(s).toBeAllCloseTo(S);
  });
})
