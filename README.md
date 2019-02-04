# Introduction
`nd.js` is a lightweight JavaScript library for ND-Arrays including some optimization functionality and one of the most complete linear algebra modules for the Web. It is strongly inspired by [NumPy](http://www.numpy.org/). There are, however, some key differences. Broadcasting, slicing and reshape work in a similar way as in NumPy. Instead of the predefined operations (+, -, *, /, sin, cos, ...), `nd.js` relies on functional-style map- and zip-like methods.

A function reference can be found [here](https://dirktoewe.github.io/ndjs/doc.html).

  * [Building and Testing](#building-and-testing)
  * [Array Instantiation](#array-instantiation)
  * [Random Access](#random-access)
  * [Unary Operations (sin, cos, exp, ...)](#unary-operations-sin-cos-exp-)
  * [Binary Operations (+, -, *, /, ...)](#binary-operations------)
  * [Ternary Operations (?:, ...)](#ternary-operations--)
  * [Linear Algebra](#linear-algebra)

# Building and Testing
`nd.js` is built and tested using [NPM](https://www.npmjs.com/). To initialize the project open the command line, navigate to the project directory and call:

```
npm i
```

To cover as many test cases as possible with little effort, `nd.js` mostly uses randomized testing instead of hand crafted test cases. As a result, testing the entire project takes rather long (~50 minutes). If You want to run only a subset of the tests during development, change the [glob](https://en.wikipedia.org/wiki/Glob_(programming)) pattern of the `file` setting inside of `karma.conf.js`, e.g. use `'src/la/**/*_test.js'` instead of `'src/**/*_test.js'` to test linear algebra methods only.

In order to run the tests call:
```
npm run test
```

To build/bundle the library, call:
```
npm run build
```

`nd.js` has some development dependencies, most notably [Babel](https://babeljs.io/) and [Webpack](https://webpack.js.org/) for bundling and [Jasmine](https://jasmine.github.io/) for testing. There are however no deployment dependencies as of yet.

# Array Instantiation
[nd.array](https://dirktoewe.github.io/ndjs/doc.html#nd.array) allows to create nd.Arrays in a well-readable and intuitive way, using nested JavaScript Arrays.

```js
const a = nd.array([
  [1,0,0],
  [0,2,0],
  [0,0,3]
])
```

An nd.Array can also be created from its entries' indices using [nd.tabulate](https://dirktoewe.github.io/ndjs/doc.html#nd.tabulate).

*Input:*
```js
const a = nd.tabulate([1000,1000], (i,j) => i==j ? 1 : 0 )
console.log( a.toString() )
```
*Output:*
```
[[ 1,  0,  0,  0,  0,  ...990 more...,  0,  0,  0,  0,  0],
 [ 0,  1,  0,  0,  0,  ...990 more...,  0,  0,  0,  0,  0],
 [ 0,  0,  1,  0,  0,  ...990 more...,  0,  0,  0,  0,  0],
 [ 0,  0,  0,  1,  0,  ...990 more...,  0,  0,  0,  0,  0],
 [ 0,  0,  0,  0,  1,  ...990 more...,  0,  0,  0,  0,  0],
  ...990 more...,
 [ 0,  0,  0,  0,  0,  ...990 more...,  1,  0,  0,  0,  0],
 [ 0,  0,  0,  0,  0,  ...990 more...,  0,  1,  0,  0,  0],
 [ 0,  0,  0,  0,  0,  ...990 more...,  0,  0,  1,  0,  0],
 [ 0,  0,  0,  0,  0,  ...990 more...,  0,  0,  0,  1,  0],
 [ 0,  0,  0,  0,  0,  ...990 more...,  0,  0,  0,  0,  1]]
```

# Random Access
Since nd.Array extends Function, elements can be read by calling the array object.

*Input:*
```js
const a = nd.array([
  [1,2],
  [3,4]
])

console.log( a(0,0) )
console.log( a(0,1) )
console.log( a(1,0) )
console.log( a(1,1) )
```
*Output:*
```js
1
2
3
4
```

[nd.NDArray.set](https://dirktoewe.github.io/ndjs/doc.html#nd.NDArray.set) allows writing array entries.

*Input:*
```js
const a = nd.tabulate([3,3], () => 0 )

for( let i=3; i-- > 0; )
for( let j=3; j-- > 0; )
  a.set( [i,j], 10*(i+1)+(j+1) );

console.log( a.toString() );
```

*Output:*
```js
[[ 11, 12, 13 ],
 [ 21, 22, 23 ],
 [ 31, 32, 33 ]]
```

If array elements are to me modified, [nd.NDArray.modify](https://dirktoewe.github.io/ndjs/doc.html#nd.NDArray.modify) is a concise
alternative to using [nd.NDArray.get](https://dirktoewe.github.io/ndjs/doc.html#nd.NDArray.get) and
[nd.NDArray.set](https://dirktoewe.github.io/ndjs/doc.html#nd.NDArray.set).

*Input:*
```js
const a = nd.array([[1, 2, 3],
                    [4, 5, 6]])

a.modify([0,1], x => 7*x)

console.log( a.toString() );
```

*Output:*
```js
[[  1, 14,  3 ],
 [  4,  5,  6 ]]
```


# Unary Operations (sin, cos, exp, ...)
[nd.NDArray.mapElems](https://dirktoewe.github.io/ndjs/doc.html#nd.NDArray.mapElems) is used to apply unary operations on an nd.Array.

*Input:*
```js
const a = nd.array([
  [1,2],
  [3,4]
])

const b = a.mapElems( (a_ij, i,j) => i==j ? a_ij : a_ij*a_ij )
console.log( b.toString() )
```
*Output:*
```js
[[         1,          4],
 [         9,          4]]
```

# Binary Operations (+, -, *, /, ...)
[nd.zip_elems](https://dirktoewe.github.io/ndjs/doc.html#nd.zip_elems) can be used to apply binary operations on two nd.Arrays.

*Input:*
```js
const
  a = nd.array([
    [11,12,13],
    [21,22,23],
    [31,32,33]
  ]),
  b = nd.array([
    [1,2,3],
    [4,5,6],
    [7,8,9]
  ])

const c = nd.zip_elems([a,b], (a_ij,b_ij, i,j) => i==j ? a_ij : b_ij )
console.log( c.toString() )
```
*Output:*
```js
[[        11,          2,          3],
 [         4,         22,          6],
 [         7,          8,         33]]
```

[nd.zip_elems](https://dirktoewe.github.io/ndjs/doc.html#nd.zip_elems) supports [NumPy-style broadcasting](https://docs.scipy.org/doc/numpy/user/basics.broadcasting.html).

*Input:*
```js
const
  a = nd.array([
    [1],
    [2],
    [3]
  ]),
  b = nd.array([1,2,3])

const c = nd.zip_elems([a,b], (x,y) => 10*x + y )
console.log( c.toString() )
```
*Output:*
```js
[[        11,         12,         13],
 [        21,         22,         23],
 [        31,         32,         33]]
```

# Ternary Operations (?:, ...)
[nd.zip_elems](https://dirktoewe.github.io/ndjs/doc.html#nd.zip_elems) can also be used for any n-ary operation, such as ternary conditional operator in JavaScript.

*Input:*
```js
const
  flags = nd.array([true, false, false, true]),
  a = nd.array([
    [10],
    [20],
    [30],
    [40]
  ]),
  b = nd.array([1,2,3,4])

const c = nd.zip_elems([flags,a,b], (f,x,y) => f ? x : y )
console.log( c.toString() )
```
*Output:*
```js
[[        10,          2,          3,         10],
 [        20,          2,          3,         20],
 [        30,          2,          3,         30],
 [        40,          2,          3,         40]]
```

# Linear Algebra
`nd.js` now offers a fairly wide variety of Linear Algebra operations in the [nd.la](https://dirktoewe.github.io/ndjs/doc.html#nd.la) subpackage.
[la.matmul](https://dirktoewe.github.io/ndjs/doc.html#nd.la.matmul) computes the matrix product of two or more matrices. The order of multiplication
is automatically optimized to minimize the number of floating point operations.

*Input:*
```js
const v = nd.array([[1,2,-3]]).T,
      A = nd.array([
        [1,2,3],
        [4,5,6],
        [7,8,9]
      ]);

console.log( nd.la.matmul(v.T, A, v) );
```
*Output:*
```js
[[ 144 ]]
```

Available operations and decompositions:
  * [Linear Equation Solution](https://dirktoewe.github.io/ndjs/doc.html#nd.la.solve)
  * [Lower Triangular Solution](https://dirktoewe.github.io/ndjs/doc.html#nd.la.tril_solve)
  * [Upper Triangular Solution](https://dirktoewe.github.io/ndjs/doc.html#nd.la.triu_solve)
  * [Linear Least Squares](https://dirktoewe.github.io/ndjs/doc.html#nd.la.lstsq)
  * [Eigen Solution](https://dirktoewe.github.io/ndjs/doc.html#nd.la.eigen)
  * [Schur Decomposition](https://dirktoewe.github.io/ndjs/doc.html#nd.la.schur_decomp)
  * [Singular Value Decomposition (SVD)](https://dirktoewe.github.io/ndjs/doc.html#nd.la.svd_decomp)
  * [QR Decomposition](https://dirktoewe.github.io/ndjs/doc.html#nd.la.qr_decomp)
  * [Rank-Revaling QR Decomposition](https://dirktoewe.github.io/ndjs/doc.html#nd.la.rrqr_decomp)
  * [LU Decomposition](https://dirktoewe.github.io/ndjs/doc.html#nd.la.lu_decomp)
  * [Cholesky Decomposition](https://dirktoewe.github.io/ndjs/doc.html#nd.la.cholesky_decomp)
  * [Hessenberg Decomposition](https://dirktoewe.github.io/ndjs/doc.html#nd.la.hessenberg_decomp)
  * [Bidiagonal Decomposition](https://dirktoewe.github.io/ndjs/doc.html#nd.la.bidiag_decomp)

