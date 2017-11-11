# nd.js
nd.js is a small library for ND-Arrays in JavaScript. It is strongly inspired by [NumPy](http://www.numpy.org/). There are, however, some key differences. Broadcasting, slicing and reshape work in a similar way as in NumPy. Instead of the predefined operations (+, -, * /, sin, cos, ...), nd.js relies on a functional-style map- and zip-like methods.

A function reference can be found [here](https://dirktoewe.github.io/ndjs/doc.html).

# Array Instantiation
[nd.array](https://dirktoewe.github.io/ndjs/doc.html#array) allows to create nd.Arrays in a well-readable and intuitive way, using nested JavaScript Arrays.

```js
const a = nd.array([
  [1,0,0],
  [0,2,0],
  [0,0,3]
])
```

An nd.Array can also be created from its entries' indices using [nd.tabulate](https://dirktoewe.github.io/ndjs/doc.html#tabulate).

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

# Unary Operations (sin, cos, exp, ...)
[nd.Array.map](https://dirktoewe.github.io/ndjs/doc.html#Array.map) is used to apply unary operations on an nd.Array.

*Input:*
```js
const a = nd.array([
  [1,2],
  [3,4]
])

const b = a.map( (a_ij, i,j) => i==j ? a_ij : a_ij*a_ij )
console.log( b.toString() )
```
*Output:*
```js
[[         1,          4],
 [         9,          4]]
```

# Binary Operations (+, -, *, /, ...)
[nd.Array.from](https://dirktoewe.github.io/ndjs/doc.html#Array.from_STATIC) can be used to apply binary operations on two nd.Arrays.

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

const c = nd.Array.from([a,b], (a_ij,b_ij, i,j) => i==j ? a_ij : b_ij )
console.log( c.toString() )
```
*Output:*
```js
[[        11,          2,          3],
 [         4,         22,          6],
 [         7,          8,         33]]
```

[nd.Array.from](https://dirktoewe.github.io/ndjs/doc.html#Array.from_STATIC) supports [NumPy-style broadcasting](https://docs.scipy.org/doc/numpy/user/basics.broadcasting.html).

*Input:*
```js
const
  a = nd.array([
    [1],
    [2],
    [3]
  ]),
  b = nd.array([1,2,3])

const c = nd.Array.from([a,b], (x,y) => 10*x + y )
console.log( c.toString() )
```
*Output:*
```js
[[        11,         12,         13],
 [        21,         22,         23],
 [        31,         32,         33]]
```

# Ternary Operations (?:, ...)
[nd.Array.from](https://dirktoewe.github.io/ndjs/doc.html#Array.from_STATIC) can also be used for any n-ary operation, such as ternary conditional operator in JavaScript.

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

const c = nd.Array.from([flags,a,b], (f,x,y) => f ? x : y )
console.log( c.toString() )
```
*Output:*
```js
[[        10,          2,          3,         10],
 [        20,          2,          3,         20],
 [        30,          2,          3,         30],
 [        40,          2,          3,         40]]
```
