# Introduction
`nd4js` is a lightweight JavaScript library for ND-Arrays including some optimization functionality and one of the most complete linear algebra modules for the web. It is strongly inspired by [NumPy](http://www.numpy.org/). There are, however, some key differences. Broadcasting, slicing and reshape work in a similar way as in NumPy. Instead of the predefined operations (+, -, *, /, sin, cos, ...), `nd4js` relies on functional-style map- and zip-like methods.

The goal of this library is to explore and improve the author's understanding of linear algebra, optimization and computational numerics in general. The author is therefore constantly working on the edge of his knowledge. Despite his best efforts, no guarantees can be made about performance, numeric accuracy, reproducability or anything at all. There will be bugs, underflows, overflows, NaNs, memory issues and the like. If You find a bug - and don't want to keep it - feel free file an issue or a PR.

A function reference can be found [here](https://dirktoewe.github.io/nd4js/doc.html).

  * [Installation](#installation)
  * [Building and Testing](#building-and-testing)
  * [Array Instantiation](#array-instantiation)
  * [Random Access](#random-access)
  * [Unary Operations (sin, cos, exp, ...)](#unary-operations-sin-cos-exp-)
  * [Binary Operations (+, -, *, /, ...)](#binary-operations------)
  * [Ternary Operations (?:, ...)](#ternary-operations--)
  * [Linear Algebra](#linear-algebra)
  * [Optimization](#optimization)

# Installation

```
npm i nd4js
```

# Building and Testing
`nd4js` is built and tested using [NPM](https://www.npmjs.com/). To initialize the project open the command line, navigate to the project directory and call:

```
npm i
```

To cover as many test cases as possible with little effort, `nd4js` mostly uses randomized testing instead of hand crafted test cases. As a result, testing the entire project takes rather long (~50 minutes). If You want to run only a subset of the tests during development, change the [glob](https://en.wikipedia.org/wiki/Glob_(programming)) pattern of the `file` setting inside of `karma.conf.js`, e.g. use `'src/la/**/*_test.js'` instead of `'src/**/*_test.js'` to test linear algebra methods only.

In order to run the tests call:
```
npm run test
```

To build/bundle the library, call:
```
npm run build
```

`nd4js` has some development dependencies, most notably [Babel](https://babeljs.io/) and [Webpack](https://webpack.js.org/) for bundling and [Jasmine](https://jasmine.github.io/) for testing. There are however no deployment dependencies as of yet.

# Array Instantiation
[nd.array](https://dirktoewe.github.io/nd4js/doc.html#nd.array) allows to create NDArray instances in a well-readable and intuitive way, using nested JavaScript Arrays.

```js
const nd = require('nd4js')

const a = nd.array([
  [1,0,0],
  [0,2,0],
  [0,0,3]
])
```

An NDArray can also be created from its entries' indices using [nd.tabulate](https://dirktoewe.github.io/nd4js/doc.html#nd.tabulate).

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

[nd.NDArray.set](https://dirktoewe.github.io/nd4js/doc.html#nd.NDArray.set) allows writing array entries.

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

If array elements are to be modified, [nd.NDArray.modify](https://dirktoewe.github.io/nd4js/doc.html#nd.NDArray.modify) is a concise
alternative to using [nd.NDArray.get](https://dirktoewe.github.io/nd4js/doc.html#nd.NDArray.get) and
[nd.NDArray.set](https://dirktoewe.github.io/nd4js/doc.html#nd.NDArray.set).

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
[nd.NDArray.mapElems](https://dirktoewe.github.io/nd4js/doc.html#nd.NDArray.mapElems) is used to apply unary operations on an NDArray.

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
[nd.zip_elems](https://dirktoewe.github.io/nd4js/doc.html#nd.zip_elems) can be used to apply binary operations on two NDArrays.

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

[nd.zip_elems](https://dirktoewe.github.io/nd4js/doc.html#nd.zip_elems) supports [NumPy-style broadcasting](https://docs.scipy.org/doc/numpy/user/basics.broadcasting.html).

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
[nd.zip_elems](https://dirktoewe.github.io/nd4js/doc.html#nd.zip_elems) can also be used for any n-ary operation, such as ternary conditional operator in JavaScript.

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
`nd4js` now offers a fairly wide variety of Linear Algebra operations in the [nd.la](https://dirktoewe.github.io/nd4js/doc.html#nd.la) subpackage.
[la.matmul](https://dirktoewe.github.io/nd4js/doc.html#nd.la.matmul) computes the matrix product of two or more matrices. The order of multiplication
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
  * [Linear Equation Systems](https://dirktoewe.github.io/nd4js/doc.html#nd.la.solve) (incl. [Lower](https://dirktoewe.github.io/nd4js/doc.html#nd.la.tril_solve) and [Upper](https://dirktoewe.github.io/nd4js/doc.html#nd.la.triu_solve) Triangular Systems).
  * [Linear Least Squares](https://dirktoewe.github.io/nd4js/doc.html#nd.la.lstsq)
  * [Eigen Solution](https://dirktoewe.github.io/nd4js/doc.html#nd.la.eigen)
  * [Schur Decomposition](https://dirktoewe.github.io/nd4js/doc.html#nd.la.schur_decomp)
  * [Singular Value Decomposition (SVD)](https://dirktoewe.github.io/nd4js/doc.html#nd.la.svd_decomp)
  * [QR Decomposition](https://dirktoewe.github.io/nd4js/doc.html#nd.la.qr_decomp)
  * [Rank-Revaling QR Decomposition](https://dirktoewe.github.io/nd4js/doc.html#nd.la.rrqr_decomp)
  * [Strong Rank-Revaling QR Decomposition](https://dirktoewe.github.io/nd4js/doc.html#nd.la.drrqr_decomp)
  * [URV Decomposition](https://dirktoewe.github.io/nd4js/doc.html#nd.la.urv_decomp_full)
  * [Cholesky Decomposition](https://dirktoewe.github.io/nd4js/doc.html#nd.la.cholesky_decomp)
  * [LDLᵀ Decomposition](https://dirktoewe.github.io/nd4js/doc.html#nd.la.ldl_decomp)
  * [PLDLᵀPᵀ Decomposition (Bunch-Kaufmann)](https://dirktoewe.github.io/nd4js/doc.html#nd.la.pldlp_decomp)
  * [Hessenberg Decomposition](https://dirktoewe.github.io/nd4js/doc.html#nd.la.hessenberg_decomp)
  * [Bidiagonal Decomposition](https://dirktoewe.github.io/nd4js/doc.html#nd.la.bidiag_decomp)
  * [LU Decomposition](https://dirktoewe.github.io/nd4js/doc.html#nd.la.lu_decomp)

# Optimization

`nd4js` offers the following linear and nonlinear optimization methods and solvers:
  * [Parameter-Linear Curve Fitting (Linear Regression)](https://dirktoewe.github.io/nd4js/doc.html#nd.opt.fit_lin)
  * Nonlinear Curve Fitting and Least Square Optimization:
    * [Levenberg-Marquardt](https://dirktoewe.github.io/nd4js/doc.html#nd.opt.lsq_lm_gen)
    * [Dogleg Method](https://dirktoewe.github.io/nd4js/doc.html#nd.opt.lsq_dogleg_gen)
    * [L-BFGS](https://dirktoewe.github.io/nd4js/doc.html#nd.opt.lsq_lbfgs_gen)
  * Orthogonal Least-Squares Regression:
    * [Levenberg-Marquardt](https://dirktoewe.github.io/nd4js/doc.html#nd.opt.odr_lm_gen)
    * [Dogleg Method](https://dirktoewe.github.io/nd4js/doc.html#nd.opt.odr_dogleg_gen)
  * General Nonlinear Optimization
    * [L-BFGS Optimizer](https://dirktoewe.github.io/nd4js/doc.html#nd.opt.min_lbfgs_gen)
    * [L-BFGS-B Optimizer (Bound Constrained)](https://dirktoewe.github.io/nd4js/doc.html#nd.opt.min_lbfgsb_gen)
  * 1D Derivative-Free Root Finding
    * [Bisection](https://dirktoewe.github.io/nd4js/doc.html#nd.opt.root1d_bisect)
    * [Brent's method](https://dirktoewe.github.io/nd4js/doc.html#nd.opt.root1d_brent)
    * [Ford's (Illinois-type) method](https://dirktoewe.github.io/nd4js/doc.html#nd.opt.root1d_illinois)
  * 1D Derivative-Free Optimization
    * [Golden Section Search 1D Minimizer](https://dirktoewe.github.io/nd4js/doc.html#nd.opt.min1d_gss)

