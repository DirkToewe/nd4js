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
{
  'use strict'

  if( undefined == nd )
    nd = require('./nd.js');

   //
  // SETUP HOMEBREW TEST FRAMEWORK
 //
  function test( name, procedure )
  {
    console.log('Testing '+name+'...')

      const assert = (assertion,err_msg) => {
        if( ! assertion )
          throw new Error( null  ==  err_msg
            ? 'Assertion failed.'
            : 'Assertion failed: ' + err_msg + '.'
          )
      }
      assert.err = (procedure,err_msg) => {
        try {
          procedure()
        }
        catch(err) { return }
        throw new Error( null == err_msg
          ? 'Error expected but none thrown.'
          : 'Error expected but none thrown: ' + err_msg + '.'
        )
      }
       //
      // EQUALITY ASSERTIONS
     //
      assert.eq = (expect, actual, err_msg) => {
        if( expect != actual )
          throw new Error( null  ==  err_msg
            ?             '"'+expect+'" expected but "'+actual+'" encountered.'
            : err_msg + '. "'+expect+'" expected but "'+actual+'" encountered.'
          )
      }
       //
      // STRICT EQUALITY ASSERTIONS
     //
      assert.strict_eq = (expect, actual, err_msg) => {
        if( expect !== actual )
          throw new Error( null  ==  err_msg
            ?             '"'+expect+'" expected but "'+actual+'" encountered.'
            : err_msg + '. "'+expect+'" expected but "'+actual+'" encountered.'
          )
      }
       //
      // ARRAY EQUALITY ASSERTIONS
     //
      assert.array_eq = (expect, actual, err_msg) => {
        if( expect.length !== actual.length || expect.some( (x,i) => x != actual[i] ) )
          throw new Error( null  ==  err_msg
            ?             '['+expect+'] expected but ['+actual+'] encountered.'
            : err_msg + '. ['+expect+'] expected but ['+actual+'] encountered.'
          )
      }
       //
      // ARRAY STRICT EQUALITY ASSERTIONS
     //
      assert.array_strict_eq = (expect, actual, err_msg) => {
        if( expect.length !== actual.length || expect.some( (x,i) => x !== actual[i] ) )
          throw new Error( null  ==  err_msg
            ?             '{'+expect+'} expected but {'+actual+'} encountered.'
            : err_msg + '. {'+expect+'} expected but {'+actual+'} encountered.'
          )
      }
       //
      // ARRAY NUMERICALLY TOLERANT EQUALITY ASSERTIONS
     //
      assert.ndarray_all_close = (expect, actual, err_msg) => {
        nd.Array.from([expect,actual], (x,y,...indices) => {
          if( ! nd.math.is_close(x,y) ) {
            let msg = '{\n'+expect+'\n} expected but {\n'+actual+'\n} encountered.\n'+x+' != '+y+' at index ['+indices+']';
            throw new Error( null  ==  err_msg ? msg : err_msg + msg );
          }
        });
      }
       //
      // OBJECTS IDENTITY ASSERTIONS
     //
      assert.is_not = (a, b, err_msg) => {
        if( Object.is(a,b) )
          throw new Error( null  ==  err_msg
            ?             'The two references must point to different objects.'
            : err_msg + '. The two references must point to different objects.'
          )
      }

    console.time('Time')
    procedure(assert)
    console.timeEnd('Time')

    console.log('Finished testing "'+name+'" successfully.\n')
  }

   //
  // LINEAR ALGEBRA TESTS
 //
//  test('nd.la.schur_decomp#BENCHMARK', assert => {
//    for( let run=8; run-- > 0; )
//    {
//      console.log(`Run ${run}`);
//      let randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from,
//          shape = [1024,1024];
//
//      const
//         A  = nd.tabulate(shape, 'float64', () => Math.random()*2 - 1 ),
//        [N] = shape.slice(-1);
//      Object.freeze(A.data.buffer);
//
//      console.time('nd.la.schur_decomp');
//      const [Q,T] = nd.la.schur_decomp(A);
//      console.timeEnd('nd.la.schur_decomp');
//
//      assert.eq('float64', Q.dtype) 
//      assert.eq('float64', T.dtype)
//
//      const a = nd.la.matmul(Q, T, Q.T);
//      assert.eq('float64', a.dtype)
//
//      assert.ndarray_all_close( 0, nd.la.tril(T,-2) )
//
//      for( const t of nd.la.diag(T,-1).reshape(-1,N-1) )
//        for( let i=1; i < N-1; i++ )
//          assert( Math.abs( t(i)*t(i-1) ) == 0 );
//
//      assert.array_eq(shape, Q.shape)
//      assert.array_eq(shape, T.shape)
//
//      assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul(Q,Q.T) );
//      assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul(Q.T,Q) );
//
//      assert.ndarray_all_close(A,a);
//    }
//  });
//  process.exit();

//  test('tf.linalg.qr#BENCHMARK', assert => {
//    const tf = require('@tensorflow/tfjs');
//
//    for( let run=1024; run-- > 0; )
//    {
//      const
//        N = 1024,
//        M = N,
//        L = Math.min(N,M),
//        shape = [N,M];
//      let A = Array.from({ length: N }, () =>
//              Array.from({ length: M }, () => Math.random()*2 - 1 ))
//      A = tf.tensor(A)
//
//      console.time('QR');
//      let [Q,R] = tf.linalg.qr(A);
//      console.timeEnd('QR');
//
//      A = nd.tabulate(A.shape, (...idx) => A.get(...idx) );
//      Q = nd.tabulate(Q.shape, (...idx) => Q.get(...idx) );
//      R = nd.tabulate(R.shape, (...idx) => R.get(...idx) );
//      const a = nd.la.matmul(Q,R);
//
//      assert.array_eq([...shape.slice(0,-2), N, L], Q.shape)
//      assert.array_eq([...shape.slice(0,-2), L, M], R.shape)
//
//      const is_close = (x,y) => {
//        const atol = 1e-5,
//              rtol = 1e-5,
//               tol = atol + rtol * nd.math.max(
//                nd.math.abs(x),
//                nd.math.abs(y)
//              );
//        return nd.math.abs( nd.math.sub(x,y) ) <= tol;
//      };
//
//      nd.Array.from([A,a], (A,a) => assert( is_close(A,a) )  )
//      if( N < M ) {
//           nd.Array.from([nd.la.eye(N), nd.la.matmul(Q,Q.T)], (x,y) => assert( is_close(x,y) ) )
//           nd.Array.from([nd.la.eye(N), nd.la.matmul(Q.T,Q)], (x,y) => assert( is_close(x,y) ) )
//      }
//      else nd.Array.from([nd.la.eye(M), nd.la.matmul(Q.T,Q)], (x,y) => assert( is_close(x,y) ) )
//
//      R.forEntries( (x,...idx) => {
//        const [i,j] = idx.slice(-2);
//        if( i > j )
//          assert( is_close(0,x), `R[${idx}] = ${x} != 0.0` )
//      });
//    }
//  })

/*
  test('nd.la.qr_decomp#BENCHMARK', assert => {
    for( let run=32; run-- > 0; )
    {
      const
        N = 1024,
        M = N,
        A = nd.tabulate( [N,M], 'float64', () => Math.random()*2-1 );
      console.time('QR');
      const [Q,R] = nd.la.qr_decomp_full(A);
      console.timeEnd('QR');
      R.forEntries( (R_ij,...indices) => {
        const [i,j] = indices.slice(-2);
        if( i > j && ! (Math.abs(R_ij) < 1e-8) )
          throw new Error(`R is not triangular:\n${R}`);
      });
    }
  })
*/
/*
  test('nd.la.bidiag_decomp#BENCHMARK', assert => {
    for( let run=32; run-- > 0; )
    {
      const
        N = 1024,
        M = 1024,
        A = nd.tabulate( [N,M], 'float64', () => Math.random()*2-1 );
      console.time('bidiag');
      const [U,B,V] = nd.la.bidiag_decomp(A);
      console.timeEnd('bidiag');
      B.forEntries( (R_ij,...indices) => {
        const [i,j] = indices.slice(-2);
        if( i > j && ! (Math.abs(R_ij) < 1e-8) )
          throw new Error(`R is not triangular:\n${R}`);
      });
    }
  })
*/
/*
  test('nd.la.lu_decomp#BENCHMARK', assert => {
    for( let run=32; run-- > 0; )
    {
      const N = 1024;
      let A = nd.tabulate( [N,N], 'float64', () => Math.random()*2-1 );
      console.time('LU');
      const [LU,P] = nd.la.lu_decomp(A);
      console.timeEnd('LU');
      A = nd.tabulate( [N,N], (i,j) => A(P(i),j) );

      const
        L = nd.la.tril(LU).mapEntries( (x,i,j) => i==j ? 1 : (i<j ? 0 : x) ),
        U = nd.la.triu(LU),
        a = nd.la.matmul(L,U);

      assert.array_eq(LU.shape,A.shape);
      assert.array_eq(P .shape,A.shape.slice(0,-1));
      assert.ndarray_all_close(A,a);
    }
  })
*/
/*
  test('nd.la.cholesky_decomp#BENCHMARK', assert => {
    for( let run=32; run-- > 0; )
    {
      const N = 1024;
      let A = nd.tabulate( [N,N], 'float64', (i,j) => i < j ? 0 : i == j ? Math.random()*0.9+0.1 : Math.random()*0.2-0.1 );
      A = nd.la.matmul2(A,A.T);
      console.time('Cholesky');
      const L = nd.la.cholesky_decomp(A);
      console.timeEnd('Cholesky');

      const a = nd.la.matmul2(L,L.T);

      assert.array_eq(L.shape,A.shape);
      assert.ndarray_all_close(A,a);
    }
  })
*/


//  test('numeric.svd#BENCHMARK', assert => {
//    nj = require('/home/dtub/Dropbox/Workspaces/JavaScript/numeric.js')
//    for( let run=1024; run-- > 0; )
//    {
//      const N = 128,
//            M = N,
//            L = Math.min(N,M);
//      let A = Array.from(
//        { length: N },
//        () =>  Float64Array.from({ length: N }, () => Math.random()*2-1 )
//      );
//  //    A = nd.la.tril(A,+4);
//  //    A = nd.la.triu(A,-4);
//  //    A = nd.la.triu(A);
//      console.time('SVD');
//      let {U,S:SV,V} = nj.svd(A);
//      console.timeEnd('SVD');
//  
//      U = nd.array(U)
//      SV= nd.array(SV)
//      V = nd.array(V).T
//  
//      A = nd.array(A);
//  
//      const
//        D = nd.la.diag_mat(SV),
//        a = nd.la.matmul(U,D,V);
//
//      for( sv of SV.reshape(-1,L) )
//      {
//        sv.forEach( x => assert(x >= 0) )
//        for( let i=1; i < sv.shape[0]; i++ )
//          assert( sv(i-1) >= sv(i) )
//      }
//
//      if( N >= M ) {
//        assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul2(U.T,U) );
//        assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul2(V.T,V) );
//        assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul2(V,V.T) );
//      }
//      else {
//        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul2(U.T,U) );
//        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul2(U,U.T) );
//        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul2(V,V.T) );
//      }
//      assert.ndarray_all_close( 0, nd.la.tril(D,-1) );
//      assert.ndarray_all_close( 0, nd.la.triu(D,+1) );
//      assert.ndarray_all_close( A, a );
//    }
//  })

//  test('complex128#EXPERIMENTS', assert => {
//    const arr = nd.tabulate([10,10], 'complex128', (i,j) => new Complex(i+1,j+1))
//    console.log( arr.toString() );
//  });


  test('nd.la.schur#EXPERIMENTS1', assert => {
    const A = nd.array([
      [  1,  -3,  3],
      [  3,  -5,  3],
      [  6,  -6,  4]
    ]);

    const [Q,T] = nd.la.schur_decomp(A);
    const EV = nd.la.schur_eigenvals(T);

    console.log(
      EV.dtype,
      EV.mapEntries( x => nd.math.mul( x.toFixed(6), 1 ) )
    );

    console.log('T:')
    console.log( T.mapEntries( x => nd.math.mul( x.toFixed(6), 1 ) ) )

    const [VALS,VECS] = nd.la.schur_eigen(Q,T);

    console.log('VALS:\n',VALS.mapEntries( x => nd.math.mul( x.toFixed(6), 1 ) ) );
    console.log('VECS:')
    console.log( VECS.mapEntries( x => nd.math.mul( x.toFixed(6), 1 ) ) );
  });


  test('nd.la.schur#EXPERIMENTS2', assert => {
    const A = nd.array([
      [  6, -13 ],
      [  1,   0 ]
    ]);

    const [Q,T] = nd.la.schur_decomp(A);
    const EV = nd.la.schur_eigenvals(T);

    console.log(
      EV.dtype,
      EV.mapEntries( x => nd.math.mul( x.toFixed(6), 1 ) )
    );

    console.log('T:')
    console.log( T.mapEntries( x => nd.math.mul( x.toFixed(6), 1 ) ) )

    const [VALS,VECS] = nd.la.schur_eigen(Q,T);

    console.log('VALS:\n',VALS.mapEntries( x => nd.math.mul( x.toFixed(6), 1 ) ) );
    console.log('VECS:')
    console.log( VECS.mapEntries( x => nd.math.mul( x.toFixed(6), 1 ) ) );
  });


  test('nd.la.schur#EXPERIMENTS3', assert => {
    const A = nd.array([
      [  5,  8,  16 ],
      [  3,  1,   8 ],
      [ -4, -4, -11 ]
    ]);

    const [Q,T] = nd.la.schur_decomp(A);
    const EV = nd.la.schur_eigenvals(T);

    console.log(
      EV.dtype,
      EV.mapEntries( x => nd.math.mul( x.toFixed(6), 1 ) )
    );

    console.log('T:')
    console.log( T.mapEntries( x => nd.math.mul( x.toFixed(6), 1 ) ) )

    const [VALS,VECS] = nd.la.schur_eigen(nd.la.eye(3),T);

    console.log('VALS:\n',VALS.mapEntries( x => nd.math.mul( x.toFixed(6), 1 ) ) );
    for( let i=0; i < 3; i++ ) {
      VECS.modify([i,1], x => x.div( VECS(2,1) ) );
      VECS.modify([i,2], x => x.div( VECS(2,2) ) );
    }
    console.log('VECS:')
    console.log( VECS.mapEntries( x => nd.math.mul( x.toFixed(6), 1 ) ) );
  });


  test('nd.la.schur#EXPERIMENTS4', assert => {
    const T = nd.array([
      [ -1.716439, -19.897788, - 0.085292 ],
      [  0.226823, - 0.283561, -12.805964 ],
      [         0,          0, - 3        ]
    ]);

    console.log('T:')
    console.log( T.mapEntries( x => nd.math.mul( x.toFixed(6), 1 ) ) )

    const [Λ,V] = nd.la.schur_eigen(nd.la.eye(3),T);

    console.log('Λ:')
    console.log( Λ.mapEntries( x => nd.math.mul( x.toFixed(6), 1 ) ) );
    for( let i=0; i < 3; i++ ) {
//      VECS.modify([i,0], x => x.div( VECS(1,0) ) );
//      VECS.modify([i,1], x => x.div( VECS(1,1) ) );
    }
    console.log('V:')
    console.log( V.mapEntries( x => nd.math.mul( x.toFixed(6), 1 ) ) );
    console.log('colNorms(V):')
    console.log(
      V.mapEntries('float64', nd.math.abs)
       .reduce(-2,nd.math.hypot)
       .mapEntries('float64', x => nd.math.mul( x.toFixed(6), 1 ) )
    );
  });


  test('nd.la.schur#EXPERIMENTS5', assert => {
    const A = nd.array([
      [  2,  1,   0 ],
      [  0,  2,   1 ],
      [  0,  0,   2 ]
    ]);

    const [Q,T] = nd.la.schur_decomp(A);
    const EV = nd.la.schur_eigenvals(T);

    console.log(
      EV.dtype,
      EV.mapEntries( x => nd.math.mul( x.toFixed(6), 1 ) )
    );

    console.log('T:')
    console.log( T.mapEntries( x => nd.math.mul( x.toFixed(6), 1 ) ) )

    const [VALS,VECS] = nd.la.schur_eigen(Q,T);

    console.log('VALS:\n',VALS.mapEntries( x => nd.math.mul( x.toFixed(6), 1 ) ) );
//    for( let i=0; i < 3; i++ ) {
//      VECS.modify([i,1], x => x.div( VECS(2,1) ) );
//      VECS.modify([i,2], x => x.div( VECS(2,2) ) );
//    }
    console.log('VECS:')
    console.log( VECS.mapEntries( x => nd.math.mul( x.toFixed(6), 1 ) ) );
  });


  test('nd.la._schur_preprocess_balance#EXPERIMENTS', assert => {
    const nRuns = 1024;
    let average = 0;
    for( let run=nRuns; run-- > 0; )
    {
      const
        randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from,
        ndim    = randInt(2,3),
        shape   = Int32Array.from({ length: ndim }, () => randInt(1,128) );
  
      shape[ndim-2] = shape[ndim-1];
  
      const
        [N] = shape.slice(-1),
         A  = nd.tabulate(shape, 'float64', () => Math.random()*2 - 1 );

      for( let i=0; i < N; i++ )
        A.modify([randInt(0,N),randInt(0,N)], x => 16*x)

      const before = Math.hypot(...A.data);
      let D; [D,s] = nd.la._schur_preprocess_balance(A,2.0)
//      let D; [D,s] = nd.la._schur_preprocess_balance(A,Infinity)
      const after = Math.hypot(...s.data);
      average += after / before;
    }
    console.log(`Average reduction: ${ average / nRuns }`)
  });


  test('indexing from (outside -> inside)', assert => {
    for( let run=1024; run-- > 0; )
    {
      const
        randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from,
        start = randInt(0,128),
        end   = randInt(1,128) + start,
        visited = Array.from({ length: end }, () => false );
  
      for( let j=1; j <= (end-start); j++ ) {
        const i = (start+end >>> 1) + ((j%2)*2-1)*(j >>> 1)
        assert( i >= start )
        assert( i < end )
        assert( ! visited[i] )
        visited[i] = true;
      }
  
      for( let i=start; i < end; i++ )
        assert( visited[i] )
    }
  });


  test('schur 2x2', assert => {
    for( let run=1024*32; run-- > 0; )
    {
      let a = Math.random()*1024 - 512, b = Math.random()*1024 - 512,
          c = Math.random()*1024 - 512, d = Math.random()*1024 - 512;
      const tr = a+d,
            det= a*d - b*c;

      // The goal is to find a givens rotation that Schur-decomposes a real-eigenvalue 2x2 matrix.
      // ┌                ┐ ┌      ┐ ┌                ┐   ┌        ┐
      // │ cos(α) -sin(α) │ │ a  b │ │ cos(α) -sin(α) │ ! │ l1   p │
      // │                │ │      │ │                │ = │        │
      // │ sin(α)  cos(α) │ │ c  d │ │ sin(α)  cos(α) │   │  0  l2 │
      // └                ┘ └      ┘ └                ┘   └        ┘
      //
      // => 0 == (a⋅sin(α) + c⋅cos(α))⋅cos(α) - (b⋅sin(α) + d⋅cos(α))⋅sin(α)
      // => 0 == (a-d)⋅sin(2⋅α) + (b+c)⋅cos(2⋅α) + (c-b) =: f(α)
      //          !
      // => |c-b| ≤ max{|a-d|,|b+c|}
      //
      // A simple and very numerically accurate solution would be Binary Search.
      // In order to do that, we have to bracket a solution. So let's determine
      // the extrema of f(α).
      //
      // f'(α_max) =!= 0 = 2*(a-d)⋅cos(2⋅α_max) - 2*(b+c)⋅sin(2⋅α_max) => α_max = atan2( a-d, b+c ) / 2 + n*π/2

      const F = α => (
          ( a*Math.sin(α) + c*Math.cos(α) ) * Math.cos(α)
        - ( b*Math.sin(α) + d*Math.cos(α) ) * Math.sin(α)
      );

      if( tr*tr >= 4*det ) {
        let α_min =         Math.atan2(a-d,b+c) / 2,
            α_max = α_min + Math.PI/2 * (α_min <= 0 ? +1 : -1);
        assert( F(α_min)*F(α_max) <= 0 )
        const α = nd.opt.root1d(F,α_min,α_max);
        assert( Math.abs( F(α) ) <= 1e-12 );
      }
    }
  })


  test('pq-formula', assert => {
    for( let run=1024*1024; run-- > 0; )
    {
      let b = Math.random()*1024 - 512,
          c = Math.random()*1024 - 512;
      [b,c] = [ -(b+c), b*c ];
      const sign = b >= 1 ? 1.0 : -1.0,
            x0 = 0.5 *     (b + sign*Math.sqrt(b*b - 4*c)), 
            x1 = 2.0 * c / (b + sign*Math.sqrt(b*b - 4*c)); 

      assert( Math.abs(x0*x0 - b*x0 + c) <= 1e-8 );
      assert( Math.abs(x1*x1 - b*x1 + c) <= 1e-8 );
    }
  });


  test('Givens Rotation Backpropagation', assert => {
    const is_close = (x,y) => {
      const atol = 1e-4,
            rtol = 1e-4,
             tol = atol + rtol * Math.max(
              Math.abs(x),
              Math.abs(y)
            );
      return Math.abs( x - y ) <= tol;
    };

    /** Numeric differentiation.
     */
    const numdiff = f => (x,y) => {
      const d = Math.sqrt(Number.EPSILON);
      x = Float64Array.from(x);
      y = Float64Array.from(y);
      const df_dx = new Float64Array(x.length),
            df_dy = new Float64Array(y.length);
      const f0 = f(x,y);
      for( let i=x.length; i-- > 0; )
      {
        const x_i = x[i];
                    x[i] += d; df_dx[i] = ( f(x,y) - f0 ) / (x[i] - x_i);
                    x[i] = x_i;
        const y_i = y[i];
                    y[i] += d; df_dy[i] = ( f(x,y) - f0 ) / (y[i] - y_i);
                    y[i] = y_i;
      }
      return [df_dx, df_dy];
    }

    /** Test functions that is some weighted sum of a Givens Rotation.
     */
    const f = (x,y) => {
      x = Float64Array.from(x);
      y = Float64Array.from(y);
      if( x.length != y.length )
        throw new Error('Assertion failed!');

      const norm = Math.hypot(x[0],y[0]),
        c = x[0] / norm,
        s = y[0] / norm;

      for( let i=x.length; i-- > 0; ) {
        const x_i = x[i],
              y_i = y[i];
        x[i] = s*y_i + c*x_i;
        y[i] = c*y_i - s*x_i;
      }
      let sum=0.0;
      for( let i=x.length; i-- > 0; ) {
        sum += (i+1)*x[i];
        sum -= (i+1)*y[i];
      }
      return sum;
    };

    /** The analytic gradients (via backpropagation)
     */
    const f_grad = (x,y) => {
      x = Float64Array.from(x);
      y = Float64Array.from(y);
      if( x.length != y.length )
        throw new Error('Assertion failed!');

      const norm = Math.hypot(x[0],y[0]),
        c = x[0] / norm,
        s = y[0] / norm;

      // there's almost certainly some underflow issues here that need to be fixed
      const dc_dx0 =   y[0]*y[0] / Math.pow(norm,3),
            dc_dy0 = - y[0]*x[0] / Math.pow(norm,3),
            ds_dx0 = - x[0]*y[0] / Math.pow(norm,3),
            ds_dy0 =   x[0]*x[0] / Math.pow(norm,3);

      // thats what usually fed to as "input" to the backpropagation
      const df_dout1 =  Float64Array.from( x, (_,i) => +i+1 ),
            df_dout2 =  Float64Array.from( y, (_,i) => -i-1 );
      const df_dx = new Float64Array(x.length),
            df_dy = new Float64Array(y.length);

      for( let i=x.length; i-- > 0; )
      {
        df_dx[i] += df_dout1[i]*c - df_dout2[i]*s;
        df_dy[i] += df_dout1[i]*s + df_dout2[i]*c;

        df_dx[0] += df_dout1[i]*(ds_dx0*y[i] + dc_dx0*x[i])  +  df_dout2[i]*(dc_dx0*y[i] - ds_dx0*x[i]);
        df_dy[0] += df_dout1[i]*(ds_dy0*y[i] + dc_dy0*x[i])  +  df_dout2[i]*(dc_dy0*y[i] - ds_dy0*x[i]);
      }
      return [df_dx, df_dy];
    };

    const f_num = numdiff(f);

    for( let run=1024; run-- > 0; )
    {
      const N = 4,
        x = Float64Array.from({ length: N }, () => Math.random()*2 - 1 ),
        y = Float64Array.from({ length: N }, () => Math.random()*2 - 1 );

      const [u,v] = f_grad(x,y);
      const [s,t] = f_num (x,y);

      for( let i=N; i-- > 0; )
      {
        assert( is_close(u[i],s[i]), `\n  [${u}]\n   =/=\n  [${s}]` );
        assert( is_close(v[i],t[i]), `\n  [${v}]\n   =/=\n  [${t}]` );
      }
    }
  });


//  test('numeric.eig#BENCHMARK', assert => {
//    nj = require('/home/dtitx/Dropbox/Workspaces/JavaScript/numeric.js')
//    for( let run=8; run-- > 0; )
//    {
//      const N = 256,
//            M = N,
//            L = Math.min(N,M);
//      let A = Array.from(
//        { length: N },
//        () =>  Float64Array.from({ length: N }, () => Math.random()*2-1 )
//      );
//    //    A = nd.la.tril(A,+4);
//    //    A = nd.la.triu(A,-4);
//    //    A = nd.la.triu(A);
//      console.time('nj.EIG');
//      let {U,S:SV,V} = nj.eig(A);
//      console.timeEnd('nj.EIG');
//    }
//  })


//  test('nd.la._svd_decomp_jac_1sided#BENCHMARK', assert => {
//    for( let run=4; run-- > 0; )
//    {
//      const
//        N = 128,
//        M = N,
//        L = Math.min(N,M),
//        A = nd.tabulate( [N,M], 'float64', (i,j) => Math.random()*2-1 );
//    //  A = nd.la.tril(A,+4);
//    //  A = nd.la.triu(A,-4);
//    //  A = nd.la.triu(A);
//      console.time('SVD_1sided');
//      const [U,SV,V] = nd.la._svd_decomp_jac_1sided(A);
//      console.timeEnd('SVD_1sided');
//    
//      const
//        D = nd.la.diag_mat(SV),
//        a = nd.la.matmul(U,D,V);
//
//      for( sv of SV.reshape(-1,L) )
//      {
//        sv.forEach( x => assert(x >= 0) )
//        for( let i=1; i < sv.shape[0]; i++ )
//          assert( sv(i-1) >= sv(i) )
//      }
//
//      if( N >= M ) {
//        assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul2(U.T,U) );
//        assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul2(V.T,V) );
//        assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul2(V,V.T) );
//      }
//      else {
//        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul2(U.T,U) );
//        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul2(U,U.T) );
//        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul2(V,V.T) );
//      }
//      assert.ndarray_all_close( 0, nd.la.tril(D,-1) );
//      assert.ndarray_all_close( 0, nd.la.triu(D,+1) );
//      assert.ndarray_all_close( A, a );
//    }
//  });


//  test('nd.la._svd_decomp_jac_classic#BENCHMARK', assert => {
//    for( let run=4; run-- > 0; )
//    {
//      const
//        N = 128,
//        M = N,
//        L = Math.min(N,M);
//
//      let A = nd.tabulate( [N,M], 'float64', (i,j) => Math.random()*2-1 );
////      A = nd.la.tril(A,+4);
////      A = nd.la.triu(A,-4);
////      A = nd.la.triu(A);
//      console.time('SVD_classic');
//      const [U,SV,V] = nd.la._svd_decomp_jac_classic(A);
//      console.timeEnd('SVD_classic');
//
//      const D = nd.la.diag_mat(SV);
//
//      for( sv of SV.reshape(-1,L) )
//      {
//        sv.forEach( x => assert(x >= 0) )
//        for( let i=1; i < sv.shape[0]; i++ )
//          assert( sv(i-1) >= sv(i) )
//      }
//  
//      const a = nd.la.matmul(U,D,V);
//  
//      if( N >= M ) {
//        assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul2(U.T,U) );
//        assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul2(V.T,V) );
//        assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul2(V,V.T) );
//      }
//      else {
//        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul2(U.T,U) );
//        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul2(U,U.T) );
//        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul2(V,V.T) );
//      }
//      assert.ndarray_all_close( 0, nd.la.tril(D,-1) );
//      assert.ndarray_all_close( 0, nd.la.triu(D,+1) );
//      assert.ndarray_all_close( A, a );
//    }
//  })


//  test('nd.la._svd_decomp_jac_round_robin#BENCHMARK', assert => {
//    for( let run=4; run-- > 0; )
//    {
//      const N = 256;
//      let A = nd.tabulate( [N,N], 'float64', (i,j) => Math.random()*2-1 );
////      A = nd.la.tril(A,+4);
////      A = nd.la.triu(A,-4);
////      A = nd.la.triu(A);
//      console.time('SVD');
//      const [U,sv,V] = nd.la._svd_decomp_jac_round_robin(A);
//      console.timeEnd('SVD');
//
//      const
//        D = nd.la.diag_mat(sv),
//        a = nd.la.matmul(U,D,V);
//
//      assert.ndarray_all_close(A,a);
//      assert.ndarray_all_close( 0, nd.la.tril(D,-1) );
//      assert.ndarray_all_close( 0, nd.la.triu(D,+1) );
//    }
//  })


//  test('nd.la.srrqr_decomp', assert => {
//    for( let run=1024; run-- > 0; )
//    {
//      let
//        randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from,
//        ndim = randInt(2,5),
//        shape = Int32Array.from({ length: ndim }, () => randInt(1,24) );
//
//      if( shape[ndim-2] < shape[ndim-1] ) {
//        const tmp = shape[ndim-2]; shape[ndim-2] = shape[ndim-1]; shape[ndim-1] = tmp;
//      }
//
//      let
//        [N,M] = shape.slice(-2),
//         L    = Math.min(N,M),
//         A    = nd.tabulate(shape, 'float64', () => Math.random()*2 - 1 ),
//        [Q,R,P] = nd.la.srrqr_decomp(A),
//         a = nd.la.matmul(Q,R);
//
//      // apply column permutations
//      A = nd.tabulate( shape, (...idx) => {
//        idx[ndim-1] = P(...idx.slice(0,-2), idx[ndim-1])
//        return A(...idx)
//      });
//
////      for( const diag of nd.la.diag(R).reshape(-1,L) ) {
////        for( let i=1; i < L; i++ )
////          assert( diag(i-1) >= diag(i) );
////      }
//
//      assert.array_eq([...shape.slice(0,-2), N, L], Q.shape);
//      assert.array_eq([...shape.slice(0,-2), L, M], R.shape);
//      if( N < M ) {
//        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul(Q,Q.T) );
//        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul(Q.T,Q) );
//      }
//      else assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul(Q.T,Q) );
//      assert.ndarray_all_close(A,a);
//      assert.ndarray_all_close( 0, nd.la.tril(R,-1) );
//    }
//  })


//  test('nd.la._svd_decomp_jac_round_robin', assert => {
//    for( let run=1024; run-- > 0; )
//    {
//      let
//        randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from,
//        ndim    = randInt(2,5),
//        shape   = Int32Array.from({ length: ndim }, () => randInt(1,24) ),
//        [N,M]   = shape.slice(-2),
//         L      = Math.min(N,M),
//         A      = nd.tabulate(shape, 'float64', () => Math.random()*2 - 1 );
//
//      // CREATE SOME RANK DEFICIENCIES
//      if( Math.random() < 0.5 )
//      {
//        const a = A.reshape(-1,N,M);
//        for( let k=a.shape[0]; k-- > 0; )
//        {
//          for( let i=0; i < N; i++ )
//            if( Math.random() < 0.1 )
//            {
//              const
//                l = randInt(0,N),
//                scale = Math.random()*4 - 2;
//              for( let j=0; j < M; j++ ) a.set( [k,i,j], scale*a(k,l,j) );
//            }
//        }
//      }
//
//      const
//        [U,SV,V]= nd.la._svd_decomp_jac_round_robin(A),
//         D      = nd.la.diag_mat(SV);
//
//      for( sv of SV.reshape(-1,L) )
//      {
//        sv.forEach( x => assert(x >= 0) )
//        for( let i=1; i < sv.shape[0]; i++ )
//          assert( sv(i-1) >= sv(i) )
//      }
//
//      const a = nd.la.matmul(U,D,V);
//
//      if( N >= M ) {
//        assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul2(U.T,U) );
//        assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul2(V.T,V) );
//        assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul2(V,V.T) );
//      }
//      else {
//        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul2(U.T,U) );
//        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul2(U,U.T) );
//        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul2(V,V.T) );
//      }
//      assert.ndarray_all_close( 0, nd.la.tril(D,-1) );
//      assert.ndarray_all_close( 0, nd.la.triu(D,+1) );
//      assert.ndarray_all_close( A, a );
//    }
//  });


//  test('nd.la.hessenberg_decomp#BENCHMARK', assert => {
//    for( let run=1024; run-- > 0; )
//    {   
//      const
//        N = 1024,
//        shape = [N,N],
//        A = nd.tabulate(shape, 'float64', () => Math.random()*2 - 1 );
//      Object.freeze(A.data.buffer);
//
//      console.time('HESS');
//      const [U,H] = nd.la.hessenberg_decomp(A);
//      console.timeEnd('HESS');
//
//      const a = nd.la.matmul(U,H,U.T);
//
//      assert.array_eq(shape, U.shape)
//      assert.array_eq(shape, H.shape)
//
//      assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul(U,U.T) );
//      assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul(U.T,U) );
//      assert.ndarray_all_close(A,a);
//
//      assert.ndarray_all_close( 0, nd.la.tril(H,-2) )
//    }
//  });


  test('nd.la.eigen', assert => {
    for( let run=1024; run-- > 0; )
    {
      let
        randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from,
        ndim = randInt(2,5),
        shape = Int32Array.from({ length: ndim }, () => randInt(1,24) );
      shape[ndim-2] = shape[ndim-1];

      const A  = nd.tabulate(shape, 'float64', () => Math.random()*2 - 1 ),
           [N] = shape.slice(-1);
      Object.freeze(A.data.buffer);

      const [Λ,V] = nd.la.eigen(A);

      assert.array_eq( shape            , V.shape )
      assert.array_eq( shape.slice(0,-1), Λ.shape )

      // ASSERT THAT THE EIGENVECTORS ARE NORMALIZED
      assert.ndarray_all_close( 1, V.mapEntries('float64',nd.math.abs).reduce(-2,'float64',nd.math.hypot) );

      const ΛV = nd.la.matmul2(A,V);

      const λ = nd.Array.from([ΛV,V], (x,y) => nd.math.mul( x, nd.math.conj(y) ) ).reduce( -2, nd.math.add );
      assert.ndarray_all_close(Λ,λ);


      const λv = nd.Array.from([V,Λ.sliceEntries('...','new',[])], 'complex128', nd.math.mul );
      assert.ndarray_all_close(ΛV,λv);
    }
  });


  test('nd.la.schur_decomp', assert => {
    for( let run=1024; run-- > 0; )
    {
      let
        randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from,
        ndim = randInt(2,5),
        shape = Int32Array.from({ length: ndim }, () => randInt(1,24) );
      shape[ndim-2] = shape[ndim-1];

      const
         A  = nd.tabulate(shape, 'float64', () => Math.random()*2 - 1 ),
        [N] = shape.slice(-1);
      Object.freeze(A.data.buffer);

      const        
        [Q,T] = nd.la.schur_decomp(A),
         a    = nd.la.matmul(Q, T, Q.T);

      assert.ndarray_all_close( 0, nd.la.tril(T,-2) )

      for( const t of T.reshape(-1,N,N) )
        for( let i=2; i < N; i++ )
          assert( Math.abs( t(i,i-1)*t(i-1,i-2) ) == 0 );

      assert.array_eq(shape, Q.shape)
      assert.array_eq(shape, T.shape)
      assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul(Q,Q.T) );
      assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul(Q.T,Q) );

      assert.ndarray_all_close(A,a);
    }
  });


  test('nd.la._schur_preprocess_balance#1', assert => {
    for( let run=1024; run-- > 0; )
    { const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from,
            ndim    = randInt(2,5),
            shape   = Int32Array.from({ length: ndim }, () => randInt(1,24) );
  
      shape[ndim-2] = shape[ndim-1];
  
      let [N] = shape.slice(-1),
           A  = nd.tabulate(shape, 'float64', () => {
             const result = Math.random()*2 - 1;
             return Math.random() < 0.1 ? 8*result : result;
           });

      const before = A.reduce( (x,y) => Math.hypot(x,y) );
      let D; [D,S] = nd.la._schur_preprocess_balance(A,2)
      const after = S.reduce( (x,y) => Math.hypot(x,y) );
      assert( after <= before );

      const a = nd.Array.from(
        [D.sliceEntries('...',[],'new'   ), S,
         D.sliceEntries('...',   'new',[])],
        (D_i,S_ij,D_j) => D_i * S_ij / D_j
      );

      nd.Array.from([A,a], (A_ij,a_ij,...idx) => assert(A_ij == a_ij, `${A_ij} != ${a_ij} ad index [${idx}]`) );
    }
  });


  test('nd.la._schur_preprocess_balance#2', assert => {
    const absMax = (x,y) => Math.max( Math.abs(x), Math.abs(y) );

    for( let run=1024; run-- > 0; )
    { const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from,
            ndim    = randInt(2,5),
            shape   = Int32Array.from({ length: ndim }, () => randInt(1,24) );
  
      shape[ndim-2] = shape[ndim-1];
  
      let [N] = shape.slice(-1),
           A  = nd.tabulate(shape, 'float64', () => {
             const result = Math.random()*2 - 1;
             return Math.random() < 0.1 ? 8*result : result;
           });

      const before = A.reduce(absMax);
      let D; [D,S] = nd.la._schur_preprocess_balance(A,Infinity)
      const after = S.reduce(absMax);
      assert( after <= before );

      const a = nd.Array.from(
        [D.sliceEntries('...',[],'new'   ), S,
         D.sliceEntries('...',   'new',[])],
        (D_i,S_ij,D_j) => D_i * S_ij / D_j
      );

      nd.Array.from([A,a], (A_ij,a_ij,...idx) => assert(A_ij == a_ij, `${A_ij} != ${a_ij} ad index [${idx}]`) );
    }
  });


  test('nd.la._svd_decomp_jac_1sided', assert => {
    for( let run=1024; run-- > 0; )
    {
//      if( run % 16 == 0 )
//        console.log(`Run ${run}`)
      let
        randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from,
        ndim    = randInt(2,5),
        shape   = Int32Array.from({ length: ndim }, () => randInt(1,24) ),
        [N,M]   = shape.slice(-2),
         L      = Math.min(N,M),
         A      = nd.tabulate(shape, 'float64', () => Math.random()*2 - 1 );

      // CREATE SOME RANK DEFICIENCIES
//      if( Math.random() < 0.5 )
//      {
//        const a = A.reshape(-1,N,M);
//        for( let k=a.shape[0]; k-- > 0; )
//          for( let i=0; i < N; i++ )
//            if( Math.random() < 0.1 )
//            { const l = randInt(0,N),
//                scale = Math.random()*4 - 2;
//              for( let j=0; j < M; j++ ) a.set( [k,i,j], scale*a(k,l,j) );
//            }
//      }

      const
        [U,SV,V] = nd.la._svd_decomp_jac_1sided(A),
        D        = nd.la.diag_mat(SV);

      for( sv of SV.reshape(-1,L) )
      {
        sv.forEach( x => assert(x >= 0) )
        for( let i=1; i < sv.shape[0]; i++ )
          assert( sv(i-1) >= sv(i) )
      }

      const a = nd.la.matmul(U,D,V);

      if( N >= M ) {
        assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul2(U.T,U) );
        assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul2(V.T,V) );
        assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul2(V,V.T) );
      }
      else {
        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul2(U.T,U) );
        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul2(U,U.T) );
        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul2(V,V.T) );
      }
      assert.ndarray_all_close( 0, nd.la.tril(D,-1) );
      assert.ndarray_all_close( 0, nd.la.triu(D,+1) );
      assert.ndarray_all_close( A, a );
    }
  });


  test('nd.la._svd_decomp_jac_classic', assert => {
    for( let run=1024; run-- > 0; )
    {
//      console.log(`Test ${run}`)
      let
        randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from,
        ndim    = randInt(2,5),
        shape   = Int32Array.from({ length: ndim }, () => randInt(1,24) ),
        [N,M]   = shape.slice(-2),
         L      = Math.min(N,M),
         A      = nd.tabulate(shape, 'float64', () => Math.max(0, Math.random() - 0.2)*2.5 - 1 );

      // CREATE SOME RANK DEFICIENCIES
      if( Math.random() < 0.5 )
      {
        const a = A.reshape(-1,N,M);
        for( let k=a.shape[0]; k-- > 0; )
        {
          for( let i=0; i < N; i++ )
            if( Math.random() < 0.1 )
            {
              const
                l = randInt(0,N),
                scale = Math.random()*4 - 2;
              for( let j=0; j < M; j++ ) a.set( [k,i,j], scale*a(k,l,j) );
            }
        }
      }

      const
        [U,SV,V]= nd.la._svd_decomp_jac_classic(A),
         D      = nd.la.diag_mat(SV);

      for( sv of SV.reshape(-1,L) )
      {
        sv.forEntries( x => assert(x >= -0.0) )
        for( let i=1; i < sv.shape[0]; i++ )
          assert( sv(i-1) >= sv(i) )
      }

      const a = nd.la.matmul(U,D,V);

      if( N >= M ) {
        assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul2(U.T,U) );
        assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul2(V.T,V) );
        assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul2(V,V.T) );
      }
      else {
        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul2(U.T,U) );
        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul2(U,U.T) );
        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul2(V,V.T) );
      }
      assert.ndarray_all_close( 0, nd.la.tril(D,-1) );
      assert.ndarray_all_close( 0, nd.la.triu(D,+1) );
      assert.ndarray_all_close( A, a );
    }
  });


  test('nd.la.rrqr_decomp', assert => {
    for( let run=1024; run-- > 0; )
    {
      let
        randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from,
        ndim = randInt(2,5),
        shape = Int32Array.from({ length: ndim }, () => randInt(1,24) );

      if( shape[ndim-2] < shape[ndim-1] ) {
        const tmp = shape[ndim-2]; shape[ndim-2] = shape[ndim-1]; shape[ndim-1] = tmp;
      }

      let
        [N,M] = shape.slice(-2),
         L    = Math.min(N,M),
         A    = nd.tabulate(shape, 'float64', () => Math.random()*2 - 1 ),
        [Q,R,P] = nd.la.rrqr_decomp(A),
         a = nd.la.matmul(Q,R);

      // apply column permutations
      A = nd.tabulate( shape, (...idx) => {
        idx[ndim-1] = P(...idx.slice(0,-2), idx[ndim-1])
        return A(...idx)
      });

      for( const diag of nd.la.diag(R).reshape(-1,L) ) {
        for( let i=1; i < L; i++ )
          assert( diag(i-1) >= diag(i) );
      }

      assert.array_eq([...shape.slice(0,-2), N, L], Q.shape);
      assert.array_eq([...shape.slice(0,-2), L, M], R.shape);
      if( N < M ) {
        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul(Q,Q.T) );
        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul(Q.T,Q) );
      }
      else assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul(Q.T,Q) );
      assert.ndarray_all_close(A,a);
      assert.ndarray_all_close( 0, nd.la.tril(R,-1) );

      R.forEntries( (x,...idx) => {
        const [i,j] = idx.slice(-2);
        if( i > j )
          assert( 0.0 == x, `R[${idx}] = ${x} != 0.0` )
      });
    }
  });


  test('nd.la.qr_decomp', assert => {
    for( let run=1024; run-- > 0; )
    {
      let
        randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from,
        ndim = randInt(2,5),
        shape = Int32Array.from({ length: ndim }, () => randInt(1,24) );
        
      const
        [N,M] = shape.slice(-2),
         L    = Math.min(N,M),
         A    = nd.tabulate(shape, 'float64', () => Math.random()*2 - 1 ),
        [Q,R] = nd.la.qr_decomp(A),
         a = nd.la.matmul(Q,R);

      assert.array_eq([...shape.slice(0,-2), N, L], Q.shape)
      assert.array_eq([...shape.slice(0,-2), L, M], R.shape)
      assert.ndarray_all_close( 0, nd.la.tril(R,-1) )

      assert.ndarray_all_close(A,a);
      if( N < M ) {
        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul(Q,Q.T) );
        assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul(Q.T,Q) );
      }
      else assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul(Q.T,Q) );

      R.forEntries( (x,...idx) => {
        const [i,j] = idx.slice(-2);
        if( i > j )
          assert( 0.0 == x, `R[${idx}] = ${x} != 0.0` )
      });
    }
  });


  test('nd.la.hessenberg_decomp', assert => {
    for( let run=1024; run-- > 0; )
    {
      let
        randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from,
        ndim = randInt(2,5),
        shape = Int32Array.from({ length: ndim }, () => randInt(1,24) );
      shape[ndim-2] = shape[ndim-1];
        
      const
        [N] = shape.slice(-1),
         A  = nd.tabulate(shape, 'float64', () => Math.random()*2 - 1 );
      Object.freeze(A.data.buffer);
      const
        [U,H] = nd.la.hessenberg_decomp(A),
         a    = nd.la.matmul(U,H,U.T);

      assert.array_eq(shape, U.shape)
      assert.array_eq(shape, H.shape)

      assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul(U,U.T) );
      assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul(U.T,U) );
      assert.ndarray_all_close(A,a);

      assert.ndarray_all_close( 0, nd.la.tril(H,-2) )
    }
  });


  test('Round Robin Loop', assert => {
    const randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from;
    for( let run=1024; run-- > 0; )
    {
      const
        N = randInt(1,128),
        visits = nd.tabulate([N,N], 'int32', () => 0);
      // https://en.wikipedia.org/wiki/Round-robin_tournament#Scheduling_algorithm
      const
        nRounds= N - (N+1)%2, 
        nGames = N+1 >>> 1; // <- number of games per round (including 1 Bye)
      for( let r=0; r < nRounds; r++ )
      for( let g=0; g < nGames ; g++ ) { // <- games can be executed in parallel (which is goal of the exercise)
        const
          i =  g==0 ? 0 : 1 +          (g+r-1) % (2*nGames-1),
          j =             1 + (2*nGames-g+r-2) % (2*nGames-1);
        if( i < N && j < N ) {
          visits.modify([i,j], v => v+1);
          visits.modify([j,i], v => v+1);
        }
      }
      visits.forEntries( (v,i,j) => {
        if(i==j) assert.eq(v,0);
        else     assert.eq(v,1);
      });
    }
  });


  test('nd.la.diag', assert => {
    for( let run=1024; run-- > 0; )
    {
      const
        randInt= (from,until) => Math.floor( Math.random() * (until-from) ) + from,
        ndim   = randInt(2,5),
        shape  = Int32Array.from({ length: ndim }, () => randInt(1,8) ),
        [N,M]  = shape.slice(-2),
        off    = randInt(-N+1,+M),
        D      = nd.tabulate(shape, 'float64', () => Math.random()*2 - 1 ),
        d      = nd.la.diag(D,off);

      d.forEntries( (d_k, ...idx) => {
        const [k] = idx.slice(-1);
        assert.eq( d_k, D(
          ...idx.slice(0,-1),
          Math.max(k,k-off),
          Math.max(k,k+off)
        ) );
      });
    }
  });


  test('nd.la.diag_mat', assert => {
    for( let run=1024; run-- > 0; )
    {
      const
        randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from,
        ndim    = randInt(2,5),
        shape   = Int32Array.from({ length: ndim }, () => randInt(1,4) ),
        [N]     = shape.slice(-1),
         d      = nd.tabulate(shape, 'float64', () => Math.random()*2 - 1 ),
         D      = nd.la.diag_mat(d);
      assert.array_eq([...shape,N], D.shape);
      D.forEntries( (D_ij,...idx) => {
        const [i,j] = idx.slice(-2);
        if( i != j ) assert( D_ij == 0 );
        else         assert( D_ij == d(...idx.slice(0,-1)) )
      });
    }
  });


  test('Jacobi SVD Rotation', assert => {
    for( let i=0; i < 1024; i++ )
    {
      const [a,b,c,d] = Array.from({ length: 4 }, () => Math.max( 0, Math.random()-0.01 )/0.99*2 - 1 );
      const
        m = Math.atan2(c-b, a+d),
        p = Math.atan2(b+c, d-a);
      const
        α = (p+m)/2,
        β = (p-m)/2;
      const
        U = [
          [Math.cos(β), -Math.sin(β)],
          [Math.sin(β),  Math.cos(β)]
        ],
        A = [
          [a,b],
          [c,d]
        ],
        V = [
          [ Math.cos(α), Math.sin(α)],
          [-Math.sin(α), Math.cos(α)]
        ],
        D = nd.la.matmul(U,A,V);
      assert.ndarray_all_close( 0, D(0,1) )
      assert.ndarray_all_close( 0, D(1,0) )
    }
  })


  test('nd.la.bidiag_decomp', assert => {
    for( let run=1024; run-- > 0; )
    {
      let
        randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from,
        ndim = randInt(2,5),
        shape = Int32Array.from({ length: ndim }, () => randInt(1,24) );

      const
        [N,M] = shape.slice(-2),
         A    = nd.tabulate(shape, 'float64', () => Math.random()*2 - 1 );

      const [U,B,V]= nd.la.bidiag_decomp(A);

      const a = nd.la.matmul(U,B,V);

      const absMax = (x,y) => nd.math.max(
        nd.math.abs(x),
        nd.math.abs(y)
      );
      assert.eq( 0, nd.la.tril(B,-1).reduce(absMax) );
      assert.eq( 0, nd.la.triu(B,+2).reduce(absMax) );
      if( N >= M ) {
        assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul2(U.T,U) );
        assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul2(V.T,V) );
        assert.ndarray_all_close( nd.la.eye(M), nd.la.matmul2(V,V.T) );
      }
      else {
        assert.ndarray_all_close( nd.la.eye(N)  , nd.la.matmul2(U.T,U) );
        assert.ndarray_all_close( nd.la.eye(N)  , nd.la.matmul2(U,U.T) );
        assert.ndarray_all_close( nd.la.eye(N+1), nd.la.matmul2(V,V.T) );
      }
      assert.ndarray_all_close( a, A );
      B.forEntries( (B,...idx) => {
        const [i,j] = idx.slice(-2);
        if( i > j   ) assert( B == 0.0 );
        if( i < j-1 ) assert( B == 0.0 );
      });
    }
  })


  test('nd.la.qr_decomp_full', assert => {
    for( let run=1024; run-- > 0; )
    {
      const
        randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from,
        ndim = randInt(2,5),
        shape = Int32Array.from({ length: ndim }, () => randInt(1,24) );

      const
        [N,M] = shape.slice(-2),
         A = nd.tabulate(shape, 'float64', () => Math.random()*2 - 1 ),
        [Q,R] = nd.la.qr_decomp_full(A),
         a = nd.la.matmul(Q,R);

      assert.array_eq(shape,R.shape)
      assert.array_eq([...shape.slice(0,-1), shape[ndim-2]], Q.shape)
      assert.ndarray_all_close( 0, nd.la.tril(R,-1) )
      assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul(Q,Q.T) );
      assert.ndarray_all_close( nd.la.eye(N), nd.la.matmul(Q.T,Q) );
      assert.ndarray_all_close(a,A);
    }
  })


  test('nd.la.qr_solve', assert => {
    for( let run=1024; run-- > 0; )
    {
      const
        randInt = (from,until) => Math.floor( Math.random() * (until-from) ) + from,
        A_shape = Int32Array.from({ length: randInt(2,5) }, () => randInt(1,24) ),
        y_shape = A_shape.slice( randInt(0,A_shape.length-1) );

      A_shape[A_shape.length-2] = A_shape[A_shape.length-1];
      y_shape[y_shape.length-2] = A_shape[A_shape.length-1];
      y_shape[y_shape.length-1] = randInt(1,16);

      const
        A = nd.tabulate(A_shape, 'float64', () => Math.random()*2-1),
        y = nd.tabulate(y_shape, 'float64', () => Math.random()*2-1),
        [Q,R] = nd.la.qr_decomp(A),
        x = nd.la.qr_solve( Q, R, y );

      assert.ndarray_all_close( y, nd.la.matmul(A,x) );
    }
  })


  test('nd.la.lu_solve', assert => {
    for( let run=1024; run-- > 0; )
    {
      let
        randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from,
        LU_shape = Int32Array.from({ length: randInt(2,8) }, () => randInt(1,8) ),
         P_shape = LU_shape.slice( randInt(0,LU_shape.length-1), -1 );
         y_shape = LU_shape.slice( randInt(0,LU_shape.length-2)     );
      switch( randInt(0,3) )
      {
        case 0: break;
        case 1:
           y_shape = Array.from({ length: randInt(2,8) }, () => randInt(1,8) );
           P_shape = y_shape.slice( randInt(0,y_shape.length-1), -1 );
          LU_shape = y_shape.slice( randInt(0,y_shape.length-2)     );
          break;
        case 2:
          P_shape  = Array.from({ length: randInt(1,7) }, () => randInt(1,8) );
          y_shape  = P_shape.slice( randInt(0,P_shape.length-1) );
          LU_shape = P_shape.slice( randInt(0,P_shape.length-1) );
          LU_shape.push( randInt(1,8) );
           y_shape.push( randInt(1,8) );
          break;
      }
      LU_shape[LU_shape.length-2] = LU_shape[LU_shape.length-1];
       y_shape[ y_shape.length-2] = LU_shape[LU_shape.length-1];
       P_shape[ P_shape.length-1] = LU_shape[LU_shape.length-1];

      for(
        let LU=LU_shape.length-2,
             P= P_shape.length-1,
             y= y_shape.length-2;
        LU-- > 0 && P-- > 0 && y-- > 0;
      )
        switch( randInt(0,4) ) {
          default: throw new Error();
          case 0: break;
          case 1: LU_shape[LU] = 1; break;
          case 2:  P_shape[P]  = 1; break;
          case 3:  y_shape[y]  = 1; break;
        }

      let
        y = nd.tabulate(y_shape,  'float64', () => Math.random()*2-1),
        P = nd.tabulate(P_shape.slice(0,-1), () => {
          const idx = Int32Array.from({ length: P_shape[P_shape.length-1] }, (_,i) => i );
          // SHUFFLE
          for( let i=idx.length-1; --i > 0; ) {
            const
              j = randInt(0,i+1),
              idx_i = idx[i]; idx[i] = idx[j]; idx[j] = idx_i;
          }
          return idx;
        });

      P = nd.tabulate(P_shape, 'int32', (...idx) => P(...idx.slice(0,-1))[idx[idx.length-1]])

      const
        LU = nd.tabulate(LU_shape,'float64', (...indices) => {
          const [i,j] = indices.slice(-2);
          if( i==j ) return Math.random()*1 + 0.5;
          if( i< j ) return 0;
          return Math.random()*2e-1 - 1e-1;
        }),
        L = LU.mapEntries( (LU_ij,...indices) => {
          [i,j] = indices.slice(-2);
          if( i == j ) return 1;
          if( i >  j ) return LU_ij;
          return 0;
        }),
        U = LU.mapEntries( (LU_ij,...indices) => {
          [i,j] = indices.slice(-2);
          if( i <= j ) return LU_ij;
          return 0;
        }),
        x = nd.la.lu_solve(LU,P,y),
        Y = nd.la.matmul(L,U,x);

      y = nd.Array.from([y,P.sliceEntries('...','new')], (y) => y)
      P = nd.Array.from([P,y.sliceEntries('...', 0   )], (P) => P)
      y = nd.tabulate( y.shape, 'float64', (...idx) => {
        idx[idx.length-2] = P( ...idx.slice(0,-1) );
        return y(...idx);
      })

      assert.ndarray_all_close(Y,y);
    }
  })


  test('nd.la.lu_decomp#1', assert => {
    let
      A = nd.array([
        [2, -1,  6],
        [1, -2,  3],
        [0,  2,-10]
      ]),
      [LU,P] = nd.la.lu_decomp(A),
       L     = nd.la.tril(LU).mapEntries( (x, ...idx) => idx[0] == idx[1] ? 1 : x ),
        U    = nd.la.triu(LU);
    A = nd.tabulate( A.shape, (i,j) => A(P(i),j) );

    assert.array_eq( A.shape, LU.shape );
    assert.ndarray_all_close( nd.la.matmul(L,U), A );
  })


  test('nd.la.lu_decomp#2', assert => {
    for( let run=1024; run-- > 0; )
    {
      let
        randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from,
        ndim = randInt(2,5),
        A_shape = Int32Array.from({ length: ndim }, () => randInt(1,24) );
      A_shape[A_shape.length-2] = A_shape[A_shape.length-1];

      let A = nd.tabulate( A_shape, () => Math.random()*2 - 1 );

      const [LU,P] = nd.la.lu_decomp(A);
      A = nd.tabulate( A_shape, (...idx) => {
        idx[ndim-2] = P( ...idx.slice(0,-1) );
        return A(...idx);
      });

      const
        L = LU.mapEntries( (x, ...idx) => { const [i,j] = idx.slice(-2); return i<j ? 0 : (i==j ? 1 : x) } ),
        U = nd.la.triu(LU),
        a = nd.la.matmul2(L,U);

      assert.array_eq(LU.shape,A.shape);
      assert.array_eq(P .shape,A.shape.slice(0,-1));
      assert.ndarray_all_close(A,a);
    }
  })


  test('nd.la.triu_solve', assert => {
    for( let run=1024; run-- > 0; )
    {
      let
        randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from,
        U_shape = Int32Array.from({ length: randInt(2,8) }, () => randInt(1,8) ),
        y_shape = U_shape.slice( randInt(0,U_shape.length-2) );
      if( Math.random() < 0.5 )
      {
        y_shape = Array.from({ length: randInt(2,8) }, () => randInt(1,8) );
        U_shape = y_shape.slice( randInt(0,y_shape.length-2) );
      }
      U_shape[U_shape.length-2] = U_shape[U_shape.length-1];
      y_shape[y_shape.length-2] = U_shape[U_shape.length-1];

      for( let U=U_shape.length-2, y=y_shape.length-2; U-- > 0 && y-- > 0; )
        switch( randInt(0,3) )
        {
          case 0: break;
          case 1: U_shape[U] = 1; break;
          case 2: y_shape[y] = 1; break;
        }

      const
        y = nd.tabulate(y_shape,'float64', () => Math.random()*2-1),
        U = nd.tabulate(U_shape,'float64', (...indices) => {
          const [i,j] = indices.slice(-2);
          if( i==j ) return Math.random()*1 + 0.5;
          if( i> j ) return 0;
          return Math.random()*2e-1 - 1e-1;
        }),
        x = nd.la.triu_solve(U,y),
        Y = nd.la.matmul2(U,x);

      assert.ndarray_all_close(y,Y);
    }
  })


  test('nd.la.tril_solve', assert => {
    for( let run=1024; run-- > 0; )
    {
      let
        randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from,
        L_shape = Int32Array.from({ length: randInt(2,8) }, () => randInt(1,8) ),
        y_shape = L_shape.slice( randInt(0,L_shape.length-2) );
      if( Math.random() < 0.5 )
      {
        y_shape = Array.from({ length: randInt(2,8) }, () => randInt(1,8) );
        L_shape = y_shape.slice( randInt(0,y_shape.length-2) );
      }
      L_shape[L_shape.length-2] = L_shape[L_shape.length-1];
      y_shape[y_shape.length-2] = L_shape[L_shape.length-1];

      for( let L=L_shape.length-2, y=y_shape.length-2; L-- > 0 && y-- > 0; )
        switch( randInt(0,3) )
        {
          case 0: break;
          case 1: L_shape[L] = 1; break;
          case 2: y_shape[y] = 1; break;
        }

      const
        y = nd.tabulate(y_shape,'float64', () => Math.random()*2-1),
        L = nd.tabulate(L_shape,'float64', (...indices) => {
          const [i,j] = indices.slice(-2);
          if( i==j ) return Math.random()*1 + 0.5;
          if( i< j ) return 0;
          return Math.random()*2e-1 - 1e-1;
        }),
        x = nd.la.tril_solve(L,y),
        Y = nd.la.matmul2(L,x);

      assert.ndarray_all_close( 0, nd.la.triu(L,1) );
      assert.ndarray_all_close(y,Y);
    }
  })


  test('nd.la.cholesky_solve', assert => {
    for( let run=1024; run-- > 0; )
    {
      let
        randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from,
        L_shape = Int32Array.from({ length: randInt(2,8) }, () => randInt(1,8) ),
        y_shape = L_shape.slice( randInt(0,L_shape.length-2) );
      if( Math.random() < 0.5 )
      {
        y_shape = Array.from({ length: randInt(2,8) }, () => randInt(1,8) );
        L_shape = y_shape.slice( randInt(0,y_shape.length-2) );
      }
      L_shape[L_shape.length-2] = L_shape[L_shape.length-1];
      y_shape[y_shape.length-2] = L_shape[L_shape.length-1];

      for( let L=L_shape.length-2, y=y_shape.length-2; L-- > 0 && y-- > 0; )
        switch( randInt(0,3) )
        {
          case 0: break;
          case 1: L_shape[L] = 1; break;
          case 2: y_shape[y] = 1; break;
        }

      const
        y = nd.tabulate(y_shape,'float64', () => Math.random()*2-1),
        L = nd.tabulate(L_shape,'float64', (...indices) => {
          const [i,j] = indices.slice(-2);
          if( i==j ) return Math.random()*1 + 0.5;
          if( i< j ) return 0;
          return Math.random()*2e-1 - 1e-1;
        }),
        x = nd.la.cholesky_solve(L,y),
        Y = nd.la.matmul(L,L.T,x);

      assert.ndarray_all_close(y,Y);
    }
  })


  test('nd.la.cholesky_decomp', assert => {
    for( let run=1024; run-- > 0; )
    {
      let
        randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from,
        shape = Int32Array.from({ length: randInt(2,5) }, () => randInt(1,24) );
      shape[shape.length-2] = shape[shape.length-1];
      const L = nd.tabulate(shape,'float64',(...indices) => {
        const [i,j] = indices.slice(-2);
        if( i==j ) return Math.random()*1 + 0.5;
        if( i< j ) return 0;
        return Math.random()*2e-1 - 1e-1;
      });
      const
        LLT = nd.la.matmul(L,L.T),
        l = nd.la.cholesky_decomp(LLT);

      assert.array_eq(l.shape,L.shape);
      assert.ndarray_all_close(l,L);
    }
  })


  test('nd.la.matmul#1', assert => {
    const
      a = nd.array([[1],[2]]),
      b = nd.array([[30,40,50]]),
      c = nd.la.matmul(a,b)
    assert.array_eq(c.shape,[2,3])
    assert.array_eq(c.data,[1*30,1*40,1*50,2*30,2*40,2*50])
  })


  test('nd.la.matmul#2', assert => {
    const
      a = nd.array([
        [1,2,3],
        [4,5,6]
      ]),
      b = nd.array([
        [ 70, 80],
        [ 90,100],
        [110,120]
      ]),
      c = nd.la.matmul(a,b)
    assert.array_eq(c.shape,[2,2])
    assert.array_eq(c.data,[
      (1*70+2*90+3*110), (1*80+2*100+3*120),
      (4*70+5*90+6*110), (4*80+5*100+6*120)
    ])
  })


  test('nd.la.matmul#3', assert => {
    for( let run=1024; run-- > 0; )
    {
      let
        randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from,
        cShape = Int32Array.from({ length: randInt(2,8) }, () => randInt(1,8) ),
        aShape = cShape.slice(),
        bShape = cShape.slice( randInt(0,cShape.length-2) );

      for( let a=aShape.length-2, b=bShape.length-2; a-- > 0 && b-- > 0; )
        switch( randInt(0,3) )
        {
          case 0: break;
          case 1: aShape[a] = 1; break;
          case 2: bShape[b] = 1; break;
        }

      if( Math.random() < 0.5 ) {
        const tmp = aShape; aShape = bShape; bShape = tmp;
      }
      aShape[aShape.length-2] = randInt(1,5)
      aShape[aShape.length-1] = bShape[bShape.length-2]
      bShape[bShape.length-1] = randInt(1,5)

      const
        a = nd.tabulate(aShape, 'float64', () => Math.random()*2 - 1 ),
        b = nd.tabulate(bShape, 'float64', () => Math.random()*2 - 1 ),
        c = nd.la.matmul(a,b),
        C = nd.Array.from([
          a.sliceEntries('...','new'),
          b.sliceEntries('...','new',[],[])
        ], 'float64', (x,y) => x*y ).reduce([-2], (x,y) => x+y );

      assert.ndarray_all_close(C,c)
    }
  })

   //
  // TEST TRANSFORMATIONS
 //
  test('nd.Array.transpose#1', assert => { 
    const arr = nd.array([1,2,3])
    for( const axes of [ [], [0] ] )
    {
      const arr_T = arr.transpose(...axes)
      assert.array_eq(arr_T.shape, arr.shape)
      assert.array_eq(arr_T.data,  arr.data)
      assert.is_not  (arr_T.data,  arr.data)
    }
  })

  test('nd.Array.transpose#2', assert => { 
    const arr = nd.array([
      [1,2,3],
      [4,5,6]
    ])
    for( const axes of [ [], [1], [1,0] ] )
    {
      const arr_T = arr.transpose(...axes)
      assert.array_eq(arr_T.shape, [3,2])
      assert.array_eq(arr_T.data, [1,4, 2,5, 3,6])
    }
    for( const axes of [ [0], [0,1] ] )
    {
      const arr_T = arr.transpose(...axes)
      assert.array_eq(arr_T.shape, arr.shape)
      assert.array_eq(arr_T.data,  arr.data )
    }
  })

  test('nd.Array.transpose#3', assert => { 
    const arr = nd.array([
      [[ 1, 2, 3, 4],
       [ 5, 6, 7, 8],
       [ 9,10,11,12]],
      [[13,14,15,16],
       [17,18,19,20],
       [21,22,23,24]],
    ])
    for( const axes of [ [], [0,2], [0,2,1] ] )
    {
      const arr_T = arr.transpose(...axes)
      assert.array_eq(arr_T.shape, [2,4,3])
      assert.array_eq(arr_T.data, [
         1, 5, 9,  2, 6,10,  3, 7,11,  4, 8,12,
        13,17,21, 14,18,22, 15,19,23, 16,20,24
      ])
    }
    for( const axes of [ [1], [1,0], [1,0,2] ] )
    {
      const arr_T = arr.transpose(...axes)
      assert.array_eq(arr_T.shape, [3,2,4])
      assert.array_eq(arr_T.data, [
         1, 2, 3, 4, 13,14,15,16,  5, 6, 7, 8,
        17,18,19,20,  9,10,11,12, 21,22,23,24
      ])
    }
    for( const axes of [ [2], [2,0], [2,0,1] ] )
    {
      const arr_T = arr.transpose(...axes)
      assert.array_eq(arr_T.shape, [4,2,3])
      assert.array_eq(arr_T.data, [
         1, 5, 9, 13,17,21,  2, 6,10, 14,18,22,
         3, 7,11, 15,19,23,  4, 8,12, 16,20,24
      ])
    }
  })

   //
  // FACTORY METHODS
 //
  test('nd.array', assert => {
    function test(...shape)
    {
      function array(d, ...indices)
      {
        if( shape.length === d )
          return indices.reduce( (a,b) => 10*a+b )
        else
          return Array.from({ length: shape[d] }, (_,i) => array(d+1, ...indices, i) )
      }
      array = nd.array( array(0) )
      assert.array_eq(shape, array.shape)
      function test(d, ...indices)
      {
        if( shape.length === d )
          assert.strict_eq(
            indices.reduce( (a,b) => 10*a+b ),
            array(...indices)
          )
        else for( let i=0; i < shape[d]; i++ )
          test(d+1, ...indices, i)
      }
      test(0)
    }
    for( let i=1; i < 8; i++ ) {test(i)
    for( let j=1; j < 8; j++ ) {test(i,j)
    for( let k=1; k < 8; k++ ) {test(i,j,k)
    for( let l=1; l < 8; l++ ) {test(i,j,k,l)}}}}
  })

  test('nd.Array.from#1', assert => {
    const
      a = new nd.Array(Int32Array.of(2,2), [1,2,3,4]),
      b = new nd.Array(Int32Array.of(2,2), [5,6,7,8])
      c = nd.Array.from([a,b], (a_ij,b_ij,i,j) => {
        assert.strict_eq(a(i,j), a_ij)
        assert.strict_eq(b(i,j), b_ij)
        return 10*a_ij + b_ij
      })
    assert.array_strict_eq([2,2], c.shape)
    assert.array_strict_eq([
      15, 26,
      37, 48
    ], c.data)
  })

  test('nd.Array.from#2', assert => {
    const
      a = new nd.Array(Int32Array.of(3,1), [1,2,3]),
      b = new nd.Array(Int32Array.of(1,4), [4,5,6,7])
      c = nd.Array.from([a,b], (a_i0,b_j,i,j) => {
        assert.strict_eq(a(i,0), a_i0)
        assert.strict_eq(b(0,j),   b_j)
        return 10*a_i0 + b_j
      })
    assert.array_strict_eq([3,4], c.shape)
    assert.array_strict_eq([
      14, 15, 16, 17,
      24, 25, 26, 27,
      34, 35, 36, 37,
    ], c.data)
  })

  test('nd.Array.from#3', assert => {
    function test(...shape)
    {
      function test(d, shapeA, shapeB)
      {
        if( shape.length === d )
        {
          shapeA = Int32Array.from(shapeA);
          shapeB = Int32Array.from(shapeB);
          function toA(indices){ return shapeA.map( (shp,i) => indices[i+shape.length-shapeA.length] % shp ) }
          function toB(indices){ return shapeB.map( (shp,i) => indices[i+shape.length-shapeB.length] % shp ) }
          const
            a = new nd.Array(shapeA, Int32Array.from({ length: shapeA.reduce((a,b) => a*b, 1) }, (_,i) => (i+1)*1000 ) ),
            b = new nd.Array(shapeB, Int32Array.from({ length: shapeB.reduce((a,b) => a*b, 1) }, (_,i) => (i+1)      ) ),
            c =     nd.Array.from([a,b], (aVal,bVal,...indices) => {
              assert.strict_eq(a(...toA(indices)), aVal)
              assert.strict_eq(b(...toB(indices)), bVal)
              return a+b
            })
          assert.err( () => nd.Array.from([a,b]) )
          assert.err( () => nd.Array.from([a,b], 'int32') )
          assert.err( () => nd.Array.from([a,b], 'not_a_type', (aVal,bVal,...indices) => {}) )
        }
        else {
          test(d+1, [...shapeA,shape[d]], [...shapeB,shape[d]])
          test(d+1, [...shapeA,      1 ], [...shapeB,shape[d]])
          test(d+1, [...shapeA,shape[d]], [...shapeB,      1 ])
          test(d+1, [...shapeA,      1 ], [...shapeB,      1 ])
          if( shapeA.length === 0 ) test(d+1, [], [...shapeB,      1 ])
          if( shapeB.length === 0 ) test(d+1, [...shapeA,      1 ], [])
        }
      }
      test(0,[],[])
    }
    for( let i=2; i <= 4; i++ ) {test(i)
    for( let j=2; j <= 4; j++ ) {test(i,j)
    for( let k=2; k <= 4; k++ ) {test(i,j,k)
    for( let l=2; l <= 4; l++ ) {test(i,j,k,l)}}}}
  })

  test('nd.Array.from#4', assert => {
    function test(shapeA, shapeB)
    {
      let broadcastCompatible = true
      for( let i=shapeA.length, j=shapeB.length; i-- > 0 && j-- > 0; )
        broadcastCompatible &= shapeA[i] === shapeB[j] || shapeA[i] === 1 || shapeB[j] === 1
      if( ! broadcastCompatible ) {
        shapeA = Int32Array.from(shapeA);
        shapeB = Int32Array.from(shapeB);
        const
          a = new nd.Array(shapeA, Int32Array.from({ length: shapeA.reduce((a,b) => a*b, 1) }, (_,i) => (i+1)*1000 ) ),
          b = new nd.Array(shapeB, Int32Array.from({ length: shapeB.reduce((a,b) => a*b, 1) }, (_,i) => (i+1)      ) )
        assert.err( () => nd.Array.from([a,b]                                           ) )
        assert.err( () => nd.Array.from([a,b],            (aVal,bVal,...indices) => null) )
        assert.err( () => nd.Array.from([a,b], 'float32', (aVal,bVal,...indices) => null) )
        assert.err( () => nd.Array.from([a,b], 'float32'                                ) )
      }
    }
    for( let i=1; i <= 6; i++ ) {
    for( let j=1; j <= 6; j++ ) {test([i],[j])
    for( let k=1; k <= 6; k++ ) {test([i,j],[k]);     test([i],[j,k]);
    for( let l=1; l <= 6; l++ ) {test([i],[j,k,l]);   test([i,j],[k,l]);   test([i,j,k],[l])
    for( let m=1; m <= 6; m++ ) {test([i],[j,k,l,m]); test([i,j],[k,l,m]); test([i,j,k],[l,m]); test([i,j,k,l],[m])}}}}}
  })

  test('nd.tabulate', assert => {
    function test(...shape)
    {
      shape = Int32Array.from(shape);
      arr = nd.tabulate(shape, (...indices) => indices.reduce((a,b) => 10*a+b) )
      assert.array_eq(shape, arr.shape)
      function test(d,...indices)
      {
        if( d == shape.length )
          assert.strict_eq( indices.reduce((a,b) => 10*a+b), arr(...indices) )
        else for( let i=0; i < shape[d]; i++ )
          test(d+1, ...indices, i)
      }
      test(0)
    }
    for( let i=1; i < 10; i++ ) {test(i)
    for( let j=1; j < 10; j++ ) {test(i,j)
    for( let k=1; k < 10; k++ ) {test(i,j,k)
    for( let l=1; l < 10; l++ ) {test(i,j,k,l)}}}}
  })

   //
  // GENERAL
 //
  test('nd.Array.toString#1', assert => {
    const a = nd.array([
      [1,2,3],
      [4,5,6]
    ])
    assert.strict_eq('[[1,2,3],[4,5,6]]', a.toString().replace(/\s/g,'') )
  })

  test('nd.Array.toString#2', assert => {
    const a = nd.array([
      [[1,2],[3,4]],
      [[5,6],[7,8]]
    ])
    assert.strict_eq('[[[1,2],[3,4]],[[5,6],[7,8]]]', a.toString().replace(/\s/g,'') )
  })

  test('nd.Array(...indices)#1', assert => {
    const arr = new nd.Array(Int32Array.of(2,3), [
      1,2,3,
      4,5,6
    ])
    function test()
    {
      for( let [i,j,a_ij] of [
        [0,0, 1], [0,1, 2], [0,2, 3],
        [1,0, 4], [1,1, 5], [1,2, 6]
      ])
      {
        assert.strict_eq( a_ij, arr(i,  j  ) )
        assert.strict_eq( a_ij, arr(i-2,j  ) )
        assert.strict_eq( a_ij, arr(i,  j-3) )
        assert.strict_eq( a_ij, arr(i-2,j-3) )
      }
      for( let j=-100; j <= 100; j++ ) {
      for( let i=   2; i <= 100; i++ ) assert.err( () => arr(i,j) )
      for( let i=-100; i <   -2; i++ ) assert.err( () => arr(i,j) )
      }
      for( let i=-100; i <= 100; i++ ) {
      for( let j=   3; j <= 100; j++ ) assert.err( () => arr(i,j) )
      for( let j=-100; j <   -3; j++ ) assert.err( () => arr(i,j) )
      }
    }
    test(); arr.data =   Int32Array.from(arr.data)
    test(); arr.data = Float64Array.from(arr.data)
    test(); arr.data = Float32Array.from(arr.data)
    test()
  })

  test('nd.Array(...indices)#2', assert => {
    const arr = new nd.Array(Int32Array.of(4), [1,2,3,4])
    function test()
    {
      for( let [i,a_i] of [ [0,1], [1,2], [2,3], [3,4] ])
      {
        assert.strict_eq( a_i, arr(i  ) )
        assert.strict_eq( a_i, arr(i-4) )
      }
      for( let i=   4; i <= 100; i++ ) assert.err( () => arr(i) )
      for( let i=-100; i <   -4; i++ ) assert.err( () => arr(i) )
    }
    test(); arr.data =   Int32Array.from(arr.data)
    test(); arr.data = Float64Array.from(arr.data)
    test(); arr.data = Float32Array.from(arr.data)
    test()
  })

  test('nd.Array(...indices)#3', assert => {
    function test(...shape) {
      assert( shape.every( d => 0 < d && d < 10 ) )
      const
        len = shape.length,
        arr = new nd.Array( Int32Array.from(shape), Array.from(
          { length: shape.reduce( (len,d) => len*d, 1 ) },
          (_,i) => i
        )),
        visited = Uint8Array.from(arr.data, () => 0 )
      function test(d, ...multi_idx)
      {
        if( shape.length == d )
        {
          let
            idx = 0,
            stride = 1,
            out_of_bounds = false
          for( let i=len; i-- > 0; stride *= shape[i] )
          {
            idx += multi_idx[i] * stride
            if(0 > multi_idx[i])
              idx  +=  shape[i] * stride
            out_of_bounds |= (
                 multi_idx[i] < -shape[i] 
              || multi_idx[i] >= shape[i]
            )
          }
          if( out_of_bounds )
            assert.err( () => arr(...multi_idx) )
          else {
            assert( 4*(1<<len) >= ++visited[idx] )
            assert.strict_eq( idx, arr(...multi_idx) )
          }
        }
        else for( let i=-8; i < +8; i++ )
          test(d+1, ...multi_idx, i)
      }
      test(0); arr.data =   Int32Array.from(arr.data)
      test(0); arr.data = Float64Array.from(arr.data)
      test(0); arr.data = Float32Array.from(arr.data)
      test(0)
      assert( visited.every( x => x === 4*(1<<len) ) )
    }
    test()
    for( let i=1; i <= 4; i++ ) {test(i)
    for( let j=1; j <= 4; j++ ) {test(i,j)
    for( let k=1; k <= 4; k++ ) {test(i,j,k)}}}
  })

   //
  // TRANSFORMATIONS
 //

  test('nd.Array.sliceEntries#1', assert => {
    const a = new nd.Array(Int32Array.of(3,4), [
      11, 12, 13, 14,
      21, 22, 23, 24,
      31, 32, 33, 34
    ])
    function test()
    {
      let b
      
      for( row of [0,1,2] )
      {
        b = a.sliceEntries(row)
        assert.strict_eq(a.dtype, b.dtype)
        assert.array_eq([4], b.shape)
        assert.array_eq([a(row,0),a(row,1),a(row,2),a(row,3)], b.data)
      }
  
      for( rest of ['...', [], [,,], [,,1], [,4,], [,4,1], [0,,], [0,,1], [0,4,], [0,4,1]])
        for( row of [0,1,2] )
        {
          b = a.sliceEntries(row,rest)
          assert.strict_eq(a.dtype, b.dtype)
          assert.array_eq([4], b.shape)
          assert.array_eq([a(row,0),a(row,1),a(row,2),a(row,3)], b.data)
        }
  
      for( rest of ['...',[], [,,], [,,1], [,3,], [,3,1], [0,,], [0,,1], [0,3,], [0,3,1]])
        for( col of [0,1,2,3] )
        {
          b = a.sliceEntries(rest,col)
          assert.strict_eq(a.dtype, b.dtype)
          assert.array_eq([3], b.shape)
          assert.array_eq([a(0,col),a(1,col),a(2,col)], b.data)
        }
  
      b = a.sliceEntries([1,,],[1,3])
      assert.strict_eq(a.dtype, b.dtype)
      assert.array_eq([2,2], b.shape)
      assert.array_eq([22,23,32,33], b.data)

      b = a.sliceEntries([,,2],[1,3])
      assert.strict_eq(a.dtype, b.dtype)
      assert.array_eq([2,2], b.shape)
      assert.array_eq([12,13,32,33], b.data)

      b = a.sliceEntries([,,2],[1,,2])
      assert.strict_eq(a.dtype, b.dtype)
      assert.array_eq([2,2], b.shape)
      assert.array_eq([12,14,32,34], b.data)
    }
    test(); a.data =  Int32Array.from(a.data)
    test(); a.data =Float64Array.from(a.data)
    test(); a.data =Float32Array.from(a.data)
    test()
  })

  test('nd.Array.reduce#1', assert => {
    const a = nd.array([
      [[1,2], [3,4]],
      [[5,6], [7,8]]
    ])
    assert.strict_eq(36, a.reduce( (a,b) => a+b ) )
    let
    b = a.reduce([0,1,2], (a,b) => a+b)
    assert.array_strict_eq([], b.shape)
    assert.array_strict_eq([36], b.data)

    b = a.reduce([], (a,b) => a+b )
    assert.array_strict_eq(a.shape, b.shape)
    assert.array_strict_eq(a.data,  b.data )

    b = a.reduce(0, (a,b) => a+b)
    assert.array_strict_eq([2,2], b.shape)
    assert.array_strict_eq([6,8,10,12], b.data)
    b = a.reduce(1, (a,b) => a+b)
    assert.array_strict_eq([2,2], b.shape)
    assert.array_strict_eq([4,6,12,14], b.data)
    b = a.reduce(2, (a,b) => a+b)
    assert.array_strict_eq([2,2], b.shape)
    assert.array_strict_eq([3,7,11,15], b.data)

    b = a.reduce([0,1], (a,b) => a+b)
    assert.array_strict_eq([2], b.shape)
    assert.array_strict_eq([16,20], b.data)
    b = a.reduce([0,2], (a,b) => a+b)
    assert.array_strict_eq([2], b.shape)
    assert.array_strict_eq([14,22], b.data)
    b = a.reduce([1,2], (a,b) => a+b)
    assert.array_strict_eq([2], b.shape)
    assert.array_strict_eq([10,26], b.data)
  })

  test('nd.Array.reduce#2', assert => {
    const a = nd.array([
      [1,2,3],
      [4,5,6]
    ])
    assert.strict_eq(21, a.reduce( (a,b) => a+b ) )
    let
    b = a.reduce([0,1], (a,b) => a+b)
    assert.array_strict_eq([], b.shape)
    assert.array_strict_eq([21], b.data)

    b = a.reduce([], (a,b) => a+b )
    assert.array_strict_eq(a.shape, b.shape)
    assert.array_strict_eq(a.data,  b.data )

    b = a.reduce(0, (a,b) => a+b)
    assert.array_strict_eq([3], b.shape)
    assert.array_strict_eq([5,7,9], b.data)
    b = a.reduce(1, (a,b) => a+b)
    assert.array_strict_eq([2], b.shape)
    assert.array_strict_eq([6,15], b.data)
  })

  test('nd.concat#1', assert => {
    const
      a = nd.array([
        [11,12,13],
        [21,22,23]
      ]),
      b = nd.array([
        [31,32,33]
      ]),
      c = nd.array([
        [41,42,43],
        [51,52,53],
        [61,62,63],
      ])
    let
    d = nd.concat(0,[a,b,c])
    assert.array_strict_eq([6,3], d.shape)
    assert.array_strict_eq([11,12,13, 21,22,23, 31,32,33, 41,42,43, 51,52,53, 61,62,63], d.data)
  })

  test('nd.stack#1', assert => {
    const
      a = nd.array([
        [11,12,13],
        [21,22,23]
      ]),
      b = nd.array([
        [31,32,33],
        [41,42,43]
      ]),
      c = nd.array([
        [51,52,53],
        [61,62,63]
      ])
    let
    d = nd.stack(0,[a,b,c])
    assert.array_strict_eq([3,2,3], d.shape)
    assert.array_strict_eq([11,12,13, 21,22,23, 31,32,33, 41,42,43, 51,52,53, 61,62,63], d.data)

    d = nd.stack(1,[a,b,c])
    assert.array_strict_eq([2,3,3], d.shape)
    assert.array_strict_eq([11,12,13, 31,32,33, 51,52,53, 21,22,23, 41,42,43, 61,62,63], d.data)

    d = nd.stack(2,[a,b,c])
    assert.array_strict_eq([2,3,3], d.shape)
    assert.array_strict_eq([11,31,51, 12,32,52, 13,33,53, 21,41,61, 22,42,62, 23,43,63], d.data)
  })

   //
  // DATA TYPE OPERATIONS
 //
  test('nd.Array.dtype', assert => {
    assert.strict_eq('float64', new nd.Array(Int32Array.of(3  ),Float64Array.of(1,2,3)       ).dtype )
    assert.strict_eq('float32', new nd.Array(Int32Array.of(2,3),Float32Array.of(1,2,3,4,5,6) ).dtype )
    assert.strict_eq(  'int32', new nd.Array(Int32Array.of(1,3),  Int32Array.of(1,2,3)       ).dtype )
    assert.strict_eq( 'object', new nd.Array(Int32Array.of(3,3),         [1,2,3,4,5,6,7,8,9] ).dtype )
  })

  test('nd._check_dtype', assert => {
    nd._check_dtype(  'int32')
    nd._check_dtype('float32')
    nd._check_dtype('float64')
    nd._check_dtype( 'object')
    assert.err( () => nd._check_dtype(12) )
    assert.err( () => nd._check_dtype('not_a_type') )
  })

  test('nd.dtypes', assert => {
    assert(    'object' in nd.dtypes )
    assert('complex128' in nd.dtypes )
    assert(   'float64' in nd.dtypes )
    assert(   'float32' in nd.dtypes )
    assert(     'int32' in nd.dtypes )
    assert( Object.keys(nd.dtypes).length === 5 )
    assert.strict_eq(       Array, nd.dtypes.object )
    assert.strict_eq(Float64Array, nd.dtypes.float64)
    assert.strict_eq(Float32Array, nd.dtypes.float32)
    assert.strict_eq(  Int32Array, nd.dtypes.  int32)
  })

  test('nd.dtypeof', assert => {
    assert.strict_eq(     'int32', nd.dtypeof( 1337) )
    assert.strict_eq(     'int32', nd.dtypeof(new Complex(1337)) )
    assert.strict_eq(   'float64', nd.dtypeof(1.337) )
    assert.strict_eq(   'float64', nd.dtypeof(new Complex(1.337)) )
    assert.strict_eq('complex128', nd.dtypeof(new Complex(1,2)) )
    assert.strict_eq(   'object' , nd.dtypeof('str') )
    assert.strict_eq(   'object' , nd.dtypeof({x:2}) )
    assert.strict_eq(   'object' , nd.dtypeof([1,2]) )
  })

  test('nd.super_dtype', assert => {
    assert.strict_eq(  'int32', nd.super_dtype('int32', 'int32') )

    assert.strict_eq('float32', nd.super_dtype('int32',  'float32') )
    assert.strict_eq('float32', nd.super_dtype('float32','float32') )
    assert.strict_eq('float32', nd.super_dtype('float32',  'int32') )

    assert.strict_eq('float64', nd.super_dtype('int32',  'float64') )
    assert.strict_eq('float64', nd.super_dtype('float32','float64') )
    assert.strict_eq('float64', nd.super_dtype('float64','float64') )
    assert.strict_eq('float64', nd.super_dtype('float64',  'int32') )
    assert.strict_eq('float64', nd.super_dtype('float64','float32') )

    assert.strict_eq('complex128', nd.super_dtype(     'int32','complex128') )
    assert.strict_eq('complex128', nd.super_dtype(   'float32','complex128') )
    assert.strict_eq('complex128', nd.super_dtype(   'float64','complex128') )
    assert.strict_eq('complex128', nd.super_dtype('complex128','complex128') )
    assert.strict_eq('complex128', nd.super_dtype('complex128',     'int32') )
    assert.strict_eq('complex128', nd.super_dtype('complex128',   'float32') )
    assert.strict_eq('complex128', nd.super_dtype('complex128',   'float64') )

    assert.strict_eq('object', nd.super_dtype(     'int32','object') )
    assert.strict_eq('object', nd.super_dtype(   'float32','object') )
    assert.strict_eq('object', nd.super_dtype(   'float64','object') )
    assert.strict_eq('object', nd.super_dtype('complex128','object') )
    assert.strict_eq('object', nd.super_dtype('object',    'object') )
    assert.strict_eq('object', nd.super_dtype('object',     'int32') )
    assert.strict_eq('object', nd.super_dtype('object',   'float32') )
    assert.strict_eq('object', nd.super_dtype('object',   'float64') )
    assert.strict_eq('object', nd.super_dtype('object','complex128') )
  })

  test('nd.is_subdtype', assert => {
    assert( nd.is_subdtype(     'int32','object') )
    assert( nd.is_subdtype(   'float32','object') )
    assert( nd.is_subdtype(   'float64','object') )
    assert( nd.is_subdtype('complex128','object') )
    assert( nd.is_subdtype(    'object','object') )

    assert( nd.is_subdtype(     'int32','complex128') )
    assert( nd.is_subdtype(   'float32','complex128') )
    assert( nd.is_subdtype(   'float64','complex128') )
    assert( nd.is_subdtype('complex128','complex128') )
    assert(!nd.is_subdtype(    'object','complex128') )

    assert( nd.is_subdtype(     'int32','float64') )
    assert( nd.is_subdtype(   'float32','float64') )
    assert( nd.is_subdtype(   'float64','float64') )
    assert(!nd.is_subdtype('complex128','float64') )
    assert(!nd.is_subdtype(    'object','float64') )

    assert( nd.is_subdtype(     'int32','float32') )
    assert( nd.is_subdtype(   'float32','float32') )
    assert(!nd.is_subdtype(   'float64','float32') )
    assert(!nd.is_subdtype('complex128','float32') )
    assert(!nd.is_subdtype(    'object','float32') )

    assert( nd.is_subdtype(     'int32','int32') )
    assert(!nd.is_subdtype(   'float32','int32') )
    assert(!nd.is_subdtype(   'float64','int32') )
    assert(!nd.is_subdtype('complex128','int32') )
    assert(!nd.is_subdtype(    'object','int32') )
  });

  console.log('All tests finished!');
}
