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

import * as nd from '.'


nd.help_str = func_or_obj => func_or_obj.__doc__ || 'No documentation available.'
nd.help     = func_or_obj => console.log( nd.help_str(func_or_obj) )


export * from '.'
export const help     = nd.help
export const help_str = nd.help_str



nd.dt.Complex.__doc__ = `\
Rudimentary implementation of the Complex number type, used
mainly for the calculation of eigenvalues. Other than that,
complex dtype is not yet supported in NDJS.
`

  //
 // FACTORY METHODS
//
nd.NDArray.constructor.__doc__ = `\
Creates a new NDArray with given shape and data array. This
constructor may not create any protection copies of the given
data but use it directly instead. This method is intended
for internal use mostly. More user-friendly and convenient
factory methods are available like nd.array, nd.tabulate or
NDArray.from.

WARNING: This method freezes the buffer that's underlying the shape it's given.

Parameters
----------
shape: Int32Array
  The shape of new NDArray. May be used directly (no protection
  copy may be perfomed). The buffer behind this Int32Array will be frozen.
data: dtype[]
  The flat data array to be used in this NDArray. May be used
  directly (no protection copy may be performed).

Returns
-------
ndarray: NDArray
  A new NDArray.

Examples
--------
>>> a = new NDArray([2,3,2], [111,112,121,122,131,132,211,212,221,222,231,232])
>>> console.log( a.toString() )
  [[[111,112],
    [121,122],
    [131,132]],
   [[211,212],
    [221,222],
    [231,232]]]
`



nd.zip_elems.__doc__ = `\
Creates a new NDArray from one or more NDArrays, using a function to
merge the values of the NDArray(s) to the values of the new NDArray.

This can be used to perform unary (sin, cos, exp, ...) and binary
(+, -, *, ...) operations.

This method supports NumPy-style broadcasting of the NDArray arguments.

Parameters
----------
ndarray: NDArray or NDArray[]
  The NDArray or list of NDArrays from which the new NDArray is built.
dtype: String
  [OPTIONAL] The data type of the new NDArray.
mapper: (...values, ...indices) => dtype
  A function used to determine the entry values of the new NDArray from
  the entry values of the old NDArray(s). May be omitted if only a single
  NDArray is given by ndarray, in which case the entry values are copied.

Returns
-------
ndarray: NDArray
  The newly created NDArray, where:
  ndarray(i0,i1,...,i[n]) = mapper(ndarray[0](i0,i1,...,i[n]), ndarray[1](i0,i1,...,i[n]), ..., ndarray[m](i0,i1,...,i[n]), i0,i1,...,i[n]).

Examples
--------
>>> let a = nd.array([1,2,3,4])
>>> let b = nd.zip_elems(a, x => x*x)
>>> console.log( b.toString() )
  [1, 4, 9, 16]


>>> let c = nd.array([[1],[2]])
>>> let d = nd.array([1,2,3,4])
>>> let e = zip_elems([c,d], (x,y) => 10*x+y )
>>> console.log( e.toString() )
  [[11,12,13,14],
   [21,22,23,24]]
`



nd.array.__doc__ = `\
Tries to heuristically create an NDArray from the input. This
method is usually used to create NDArrays nested arrays, which
allows for a creation of NDArray which is well readable by
humans.

  * If the input is an object containing a 'shape' and 'data'
    property, a new NDArray is created using said shape and data
    directly (creating a protection copy of data).
  * If the input is a (possibly nested) JavaScript array-like, it
    is converted to an NDArray, copying the array data in the
    process. The nesting defines the shape of the resulting
    NDArray. 
  * Otherwise the data is interpreted as a scalar value and
    is put into an NDArray of shape [].

Parameters
----------
dtype: String
  [Optional] The dtype of the returned NDArray.
content: { length: int[], data: dtype[] } or dtype or dtype[] or dtype[][] or ...
  The content of the returned NDArray.

Returns
-------
ndarray: NDArray
  An NDArray representation of the content.

Examples
--------
>>> let content = { shape: [2,3], data: [1,2,3,4,5,6] }
>>> let a = nd.array(content)
>>> console.log( a.toString() )
  [[1,2,3],
   [4,5,6]]
>>> content.data[3] = 7
>>> console.log( a.toString() )
  [[1,2,3],
   [7,5,6]]


>>> let content = { shape: [2,3], data: [1,2,3,4,5,6] }
>>> let a = nd.array('int32',content)
>>> console.log( a.toString() )
  [[1,2,3],
   [4,5,6]]
>>> content.data[3] = 7
>>> console.log( a.toString() )
  [[1,2,3],
   [4,5,6]]


>>> let a = nd.array([
...   [1,2,3,4],
...   [5,6,7,8]
... ])
>>> console.log( a.toString() )
  [[1,2,3,4],
   [5,6,7,8]]


>>> console.log( nd.array(12) )
  { [NDArray: self] shape: [], data: Int32Array[12] }
`


nd.asarray.__doc__ = `\
Similar to \`nd.array(content)\` except that input is not copied if
it already is an \`NDArray\`.

Parameters
----------
content: { length: int[], data: dtype[] } or dtype or dtype[] or dtype[][] or ...
  The content of the returned NDArray.

Returns
-------
ndarray: NDArray
  \`ndarray = content instanceof NDArray ? content : nd.array(content);\`
`



nd.tabulate.__doc__ = `\
Creates a new NDArray of the given shape, calling a function for
each entry index and using the returned value as entry value of
the new NDArray.

Parameters
----------
shape: int[]
  The shape of the new NDArray.
dtype: String
  [OPTIONAL] The data type of the newly created NDArray.
idx2val: (...int) => dtype
  The function returning the entry value for a given entry index as input.

Returns
-------
ndarray: NDArray
  The newly tabulated NDArray, where:
  ndarray(i0,i1,...,i[n]) = idx2val(i0,i1,...i[n])

Examples
--------
>>> let a = nd.tabulate([3,2], (i,j) => 10*(i+1) + (j+1) )
>>> console.log( a.toString() )
  [[11,12],
   [21,22],
   [31,32]]
`



  //
 // GENERAL
//
nd.NDArray.__doc__ = `\
An n-dimensional, axis-aligned grid of entries, very similar
to NumPy ndarrays in Python.

Attributes
----------
shape: int[]
  The shape of the NDArray.
data : dtype[]
  The flat data array backing this NDArray. May be either
  a standard JavaScript Array, or one of JavaScripts primitive
  array types. See nd.dtypes for a list of supported data types.
ndim: int
  The number of dimensions in this NDArray, i.e. the number of
  indices used to address the entries. Equivalent to shape.length.
  
dtype: String
  The data type of this NDArray. Specifies which type of values
  may be stored in the NDArray. See nd.dtypes for the available
  data types and their respective (primitive) array types.
`


nd.NDArray.prototype.call.__doc__ = `\
Returns the value of the entry specified by the given index.

Parameters
----------
this: NDArray
  The array whose entry is being read.
indices: ...int
  The multi-index of the entry whose value is to be returned.

Returns
-------
value: dtype

Examples
--------
>>> let a = nd.array([[1,2,3],[4,5,6]])
>>> console.log( a(0,0), a(0,1), a(0,2), a(1,0), a(1,1), a(1,2) )
  1 2 3 4 5 6
`
nd.NDArray.prototype.apply.__doc__ = nd.NDArray.prototype.call.__doc__



nd.NDArray.prototype.set.__doc__ = `\
Sets the value of the specified entry to the given value.

Parameters
----------
indices: int[]
  The multi index of the entry that is to be changed.
value: dtype
  The new value for the entry at indices.

Examples
--------
>>> let a = nd.array([ [0,0,0], [0,0,0], [0,0,0] )
>>> a.set([0,0], 1)
>>> a.set([0,1], 2)
>>> a.set([1,1], 3)
>>> a.set([2,2], 4)
>>> console.log( a.toString() )
  [[1,2,0],
   [0,3,0],
   [0,0,4]]
`



nd.NDArray.prototype.modify.__doc__ = `\
Modifies an array entriy by applying the given functions to it.
This method call is equivalent to:
\`ndarray.set( indices, modifier(ndarray(...indices), ...indices) )\`

Parameters
----------
indices: int[]
  The multi index of the entry that is to modified.
modifier: (dtype, ...int) => dtype
  The function that is applied to the NDArray entry.
`



nd.NDArray.prototype.toString.__doc__ = `\
Returns a readable string representation of this NDArray.

Parameters
----------
max_len: int
  [OPTIONAL] The maximum number of elements along each to be represented
  in the String. If there is more elements along an axis than that, an
  ellipsis (...) is inserted instead of the central elements.

Returns
-------
repr: String
  A readable string representation of this nd.Arrray.

Examples
--------
>>> let a = nd.tabulate([6,6], (i,j) => 10*(i+1) + (j+1) )
>>> console.log( a.toString() )
  [[11,12,13,14,15,16],
   [21,22,23,24,25,26],
   [31,32,33,34,35,36],
   [41,42,43,44,45,46],
   [51,52,53,54,55,56],
   [61,62,63,64,65,66]]

>>> console.log( a.toString(4) )
  [[11, 12, ...2 more..., 15, 16],
   [21, 22, ...2 more..., 25, 26],
    ...2 more...,
   [51, 52, ...2 more..., 55, 56],
   [61, 62, ...2 more..., 65, 66]]
`



  //
 // ITERATIONS
//
nd.NDArray.prototype[Symbol.iterator].__doc__ =`\
Returns an iterator that iterates through the slices of this NDArray along
the first axis. Equivalent to:

function*() {
  for( let i=0; i < this.shape[0]; i++ )
    yield this.slice(i)
}

Returns
-------
iter: Iterator
  An iterator that yields the slices of this NDArray along its first axis.

Examples
--------
>>> a = nd.array([[1,2],[3,4],[5,6]])
>>> for( let row of a )
...   console.log( row.toString() )
  [1,2]
  [3,4]
  [5,6]
`



nd.NDArray.prototype.forEach.__doc__ = `\
Calls the given callback for each slice of this NDArray along the
first axis. Equivalent to:

for( let i=0; i < this.shape[0]; i++ )
  consumer( this.slice(i), i )

Parameters
----------
consumer: (slice: NDArray, index: int) => ()
  The callback to be called for each slice of this NDArray along the first axis.

Examples
--------
>>> let a = nd.array([[1,2,3],[4,5,6]])
>>> a.forEach( (row,i) => console.log(\`\${i} -> \${row}\`) )
  0 -> [1,2,3]
  1 -> [4,5,6]
`



nd.NDArray.prototype.elems.__doc__ = `\
Returns an iterator of all multiindex-value pairs of entries in
this NDArray.

Returns
-------
iter: *[[...int], dtype]
  An iterator of multiindex-value pairs, one for each entry in this
  NDArray.

Examples
--------
>>> a = nd.array([[1,2,3],[4,5,6]])
>>> for( let [[i,j],a_ij] of a.entries() )
...   console.log(\`a[\${i},\${j}] = \${val}\`)
  a[0,0] = 1
  a[0,1] = 2
  a[0,2] = 3
  a[1,0] = 4
  a[1,1] = 5
  a[1,2] = 6
`



nd.NDArray.prototype.forElems.__doc__ = `\
Calls the given callback for each entry in this NDArray.

Parameters
----------
consumer: (val: dtype, indices: ...int) => ()

Examples
--------
>>> let a = nd.array([[1,2,3],[4,5,6]])
>>> a.forEntries( (a_ij,i,j) => console.log(\`a[\${i},\${j}] = \${a_ij}\`) )
  a[0,0] = 1
  a[0,1] = 2
  a[0,2] = 3
  a[1,0] = 4
  a[1,1] = 5
  a[1,2] = 6
`



  //
 // TRANSFORMATION
//
nd.NDArray.prototype.valueOf.__doc__ = `
If this NDArray is scalar (shape == []), the only entry's value is returned,
otherwise this NDArray itself is returned.

This method rarely has to be called explicitly but is used by JavaScript
whenever it deems it appropriate.

Returns
-------
value: dtype or NDArray
  Returns the only entry if this NDArray is scalar. Returns this otherwise.

Example
-------
>>> a = nd.array(12)
>>> console.log( a+1 )
  13
`



nd.NDArray.prototype.mapElems.__doc__ = `\
Creates a new NDArray by applying the given mapping function to each entry
of this array and writing the results back into a new array.

Equivalent to NDArray.from(this, dtype, mapper).

Parameters
----------
dtype: String
  [OPTIONAL] The data type of the newly created NDArray.
mapper: (value, ...indices) => dtype
  The function used to map the old NDArray's entries to the values of the
  new NDArray.

Returns
-------
mapped: NDArray
  The newly created NDArray, where:
  mapped(i0,i1,...,i[n]) = mapper(this(i0,i1,...,i[n]), i0,i1,...,i[n]).

Examples
--------
>>> let a = nd.array([[1,2],[3,4]])
>>> let b = a.mapElems( x => x*x )
>>> console.log( b.toString() )
  [[1, 4],
   [9,16]]
`



nd.NDArray.prototype.transpose.__doc__ = `\
Reorders the axes of this NDArray. By default the last two axes are
swapped.

Parameters
----------
axes: ...int
  The order in which the axes of the this NDArray appear in the
  transposed array. If indices are missing from \`axes\`, they are
  appended to \`axes\` in order.

Returns
-------
transposed: NDArray
  A transposed copy A of this NDArray, where:
  \`A[i[0], i[1], ...]  =  this[i[axes[0]], i[axes[1], ...]\`

Examples
--------
>>> let a = nd.array([[1,2,3]])
>>> console.log( a.transpose() )
  [[1],
   [2],
   [3]]

>>> let a = nd.array([ [[1,2,3], [4,5,6]], [[7,8,9], [10,11,12]] ])
>>> console.log( a.transpose(2).toString() )
  [[[ 1,  4],
    [ 7, 10]],
   [[ 2,  5],
    [ 8, 11]],
   [[ 3,  6],
    [ 9, 12]]]
`



nd.NDArray.prototype.reshape.__doc__ = `\
Returns view of this NDArray with a different shape. This is similar
to NumPy's reshape with a 'C'-order.

Parameters
----------
shape: ...int
  The shape of the view. May contain a single -1 entry, in which case
  the respective axis' size is inferred.

Returns
-------
reshaped: NDArray
  A reshaped view of this NDArray.

Examples
--------
>>> let a = nd.array([1,2,3,4,5,6,7,8,9])
>>> let b = a.reshape(3,-1)
>>> console.log( b.toString() )
  [[1,2,3],
   [4,5,6],
   [7,8,9]]
`



nd.NDArray.prototype.reduceElems.__doc__ = `\
Uses the given binary operator to reduce the entries of of this NDArray
along the specified axes. If no axes are specified all entries are reduce
to a single value and said value is returned instead of an NDArray.

Parameters
----------
axes: int or int[]
  [OPTIONAL] The axes along which the array is to be reduced. If not defined
  the NDArray is reduced along all axes and a scalar value is returned instead
  of an NDArray.
dtype: String
  The data type of the reduced NDArray. Has to be a super-dtype of the original
  NDArray.
reducer: (dtype,dtype) => dtype
  The function used to reduce the entries along the axes.

Returns
-------
reduced: NDArray or dtype
  The reduced NDArray if axes were specified or the reduced value id no
  axes were specified.

Examples
--------
>>> let a = nd.array([
...   [1,2,3],
...   [4,5,6]
... ])
>>> console.log( a.reduce( (x,y) => x+y ) )
  21


>>> console.log( a.reduce( [0,1], (x,y) => x+y ) )
  { [NDArray: self] shape: [], data: [21] }

>>> console.log( a.reduce( 0, (x,y) => x+y ).toString() )
  [5,7,9]

>>> console.log( a.reduce( [1], (x,y) => x+y ).toString() )
  [6,15]
`



nd.NDArray.prototype.sliceElems.__doc__ = `\
Extracts a sub-region specified by combination of indices, ranges, newaxis symbols.

Parameters
----------
slices: int or nd.newaxis or nd.ellipsis or [start,stop,step]
  The slices to be taken along each axis, e.g 3 would only take the fourth entries
  along an axis, [,,] would take all elements along and axis, [1,,] would take all
  but the first entries along an axis, [,,3] would take ever third element along
  an axis, [,-1,] would take all but the last entries along an axis, nd.newaxis
  would insert a new axis of size 1.

  nd.ellipsis can be used to fill up with [,,] for the remaining axes.

Returns
-------
sliced: NDArray
  A sliced subregion of this NDArray.

Examples
--------
>>> let a = nd.array([
...   [11,12,13,14],
...   [21,22,23,24],
...   [31,32,33,34]
... ])
>>> console.log( a.sliceElems('...', -2).toString() )
  [13, 23, 33]

>>> console.log( a.sliceElems([,,2]).toString() )
  [[11,12,13,14],
   [31,32,33,34]]

>>> console.log( a.sliceElems(2,[1,3,]).toString() )
  [32,33]
`



nd.stack.__doc__ = `\
Arranges a list of NDArrays into a new NDArray. The NDArrays are stacked
along a newly axis, inserted at the specified index. All stacked NDArrays
must have the same shape.

Parameters
----------
axis: int
  [OPTIONAL] The index at which the new index is to be inserted. Default value: 0.
dtype: String
  [OPTIONAL] The type of the stacked NDArray. Should be a super-dtype of all NDArrays'
  dtypes.
ndarrays: NDArray[]
  The NDArrays that are to be stacked. All NDArrays must have the same shape.

Returns
-------
stacked: NDArray

Examples
--------
>>> let a = nd.array([1,2,3])
>>> let b = nd.array([4,5,6])
>>> console.log( nd.stack([a,b]).toString )
  [[1,2,3],
   [4,5,6]]

>>> console.log( nd.stack(1,'float64',[a,b]).toString )
  [[1.0, 4.0],
   [2.0, 5.0],
   [3.0, 6.0]]
`



nd.concat.__doc__ = `\
Arranges a list of NDArrays into a new NDArray. The NDArrays are concatenated
along an existing specified axis. Aside from said axis, the shape of all NDArrays
must be the same.

Parameters
----------
axis: int
  [OPTIONAL] The axis along which the NDArrays are concatenated. Default value: 0.
dtype: String
  [OPTIONAL] The type of the concatenated NDArray. Should be a super-dtype of all NDArrays'
  dtypes.
ndarrays: NDArray[]
  The NDArrays that are to be concatenated. All NDArrays must have the same shape.

Returns
-------
concatenated: NDArray


Examples
--------
>>> let a = nd.array([[1,2],[3,4]])
>>> let b = nd.array([[5,6],[7,8]])
>>> console.log( nd.concat([a,b]).toString() )
  [[1,2],
   [3,4],
   [5,6],
   [7,8]]
   
>>> console.log( nd.concat(1,[a,b]).toString() )
  [[1,2,5,6],
   [3,4,7,8]]
`



nd.help.__doc__ = `\
Outputs a documentation string for the given method ND.JS method to the console. This method
is intended to be used in interactive mode to explore and learn the ND.JS API.

This method is only available in the non-minfied version of ND.JS.

Parameters
----------
fun: Function
  The function or package ND.JS API method for which the documentation is to be retrieved.

Examples
--------
>>> nd.help(nd.la)
  The Linear Algebra subpackage of ND.JS. Contains methods for the handling and decomposition
  of matrices as well as solvers for linear equations and linear least square systems.

  As convention the last two axes of an \`NDArray\` are considered as the matrix dimensions.
  Higher dimensional arrays are considered "arrays of matrices" and the broadcasting rules
  apply for the leading dimenions, i.e. for all dimensions of the shape except the last two.
`



nd.help_str.__doc__ = `\
Returns a documentation string for the given method ND.JS method to the console. This method
is used by \`create_doc_jsdom()\` to create the documentation and by \`nd.help()\` to
print help to the console.

This method is only available in the non-minfied version of ND.JS.

Parameters
----------
fun: Function
  The function or package ND.JS API method for which the documentation is to be retrieved.

Returns
-------
help: String
  A documentation String for the given method or package.

Examples
--------
>>> console.log( nd.help(nd.la) )
  "The Linear Algebra subpackage of ND.JS. Contains methods for the handling and decomposition
   of matrices as well as solvers for linear equations and linear least square systems.

   As convention the last two axes of an \`NDArray\` are considered as the matrix dimensions.
   Higher dimensional arrays are considered "arrays of matrices" and the broadcasting rules
   apply for the leading dimenions, i.e. for all dimensions of the shape except the last two."
`



  //
 // OPTIMIZATION
//
nd.opt.min_lbfgs_gen.__doc__ = `\
Iteratively minimizes a function using the L-BFGS method. An indefinite number of solutions is
returned unless the line search does not make any further progress, in which case a
\`LineSearchNoProgressError\` is thrown. The user must check for a proper stopping condition
her-/himself.

Parameters
----------
fg: (x: NDArray[N]) => [f: NDArray[], g: NDArray[N]]
  A method return both the function value and gradients of the optimized function for
  the given input \`x\`. 
x0: NDArray[N]
  The starting point for the minimization.
options: {
  historySize=8: int
    The number of past value-gradient pairs that are memoized in order to approximate the Hessian.
  lineSearch=nd.opt.line_search.strong_wolfe()
    The line search method used to advance the solution along the current search direction.
    Must at least satisfy the Wolfe Condition.
  negDir0 = g=>g: (g: NDArray[N]) => (negDir: NDArray[N])
    Returns the gradient-descent-like initial search direction. Can be used to control the initial
    search step.
}

Returns
-------
approximations: Iterator<[x: NDArray[N], f: NDArray[] g: NDArray[N]]>
  An iterator over the approximations/iterations made by the L-BFGS method. Will return an indefinite
  amount of approximations unless a \`LineSearchNoProgressError\` is thrown.

Throws
------
noProgress: nd.opt.LineSearchNoProgressError
  If the optimization is not making any more progress. This is usally happens when the optimizer
  is already very close to the minimum.

References
----------
.. [1] https://en.wikipedia.org/wiki/Limited-memory_BFGS
.. [2] https://en.wikipedia.org/wiki/Wolfe_conditions

Example
-------
>>> const fg = ([x,y]) => [
...   nd.array(  (x-1)**2 + 100*(y-x*x)**2 ),
...   nd.array([ (x-1) *2 - 400*(y-x*x) *x,
...                         200*(y-x*x) ])
... ];
... let x,f,g, nIter = -1;
... try {
...   for( [x,f,g] of nd.opt.min_lbfgs_gen(fg, /*x0=*/[0,0], {negDir: g => g.mapElems(x => x*0.1)}) )
...   {
...     const gNorm = g.reduceElems(nd.math.hypot);
...     if( ++nIter > 1e3 )
...       throw new Error('Too many iterations.');
...     if( gNorm <= 1e-8 )
...       break;
...   }
... }
... catch(err) {
...   if( ! (err instanceof nd.opt.LineSearchNoProgressError) )
...     throw err;
...   console.log('No progress.');
... }
... console.log('Solution:', x);
  Solution: [ 0.9999999998879853, 0.9999999997700713 ]
`


nd.opt.root1d_bisect.__doc__ = `\
Finds a single root of a continuous univariate function using the bisection method.

Parameters
----------
F: float => float
  Continuous function for which a root is to be found.
x1: float
  One side of the range containing a root. F(x1) and F(x2) must have opposite signs (or be zero).
x2: float
  Other side of the range containing a root. F(x1) and F(x2) must have opposite signs (or be zero).

Returns
-------
x0: float
  A root of F, i.e. F(x0) ≈ 0.

References
----------
.. [1] https://en.wikipedia.org/wiki/Bisection_method

Example
-------
>>> const sqrt9 = nd.opt.root1d_bisect(x => x*x-9, 0, 9);
... console.log('sqrt(9) =', sqrt9);
  sqrt(9) = 3
`


nd.opt.min1d_gss.__doc__ = `\
Finds a single minimum of a continuous univariate function using Golden Section Search.

Parameters
----------
F: float => float
  Continuous function for which a minimum is to be found.
x1: float
  One side of the range containing a root. Must be left of the minimum.
x2: float
  Other side of the range containing a root. Must be right of the minimum.

Returns
-------
xMin: float
  A local minimum of F.

References
----------
.. [1] https://en.wikipedia.org/wiki/Golden-section_search

Example
-------
>>> const sqrt9 = nd.opt.min1d_gss(x => Math.abs(x*x-9), -9, +9);
... console.log('sqrt(9) =', sqrt9);
  sqrt(9) = 3
`


nd.opt.fit_lin.__doc__ = `\
Fits a parameter-linear linear function using Linear Regression or Ridge Regression.

Parameters
----------
x: float[n_samples(,n_inputs)]
  Function inputs of the sample points that the function is fit through.
y: float[n_samples]
  Function outputs of the sample points that the function is fit through.
regularization: float
  [optional] A regularization factor that penalizes the norm of
  the function coefficients. The larger the regularization, the
  smaller the coefficients tend to be, i.e. the function becomes
  smoother and more approximative than interpolative.
funcs: (float[n_inputs] => float)[n_coeffs]
  Components of the function that is to be fit through (x,y).
  Each component has it's own linear parameter/coefficient, i.e.
  \`fit(x) = c[0]*funcs[0](x) + c[1]*funcs[1](x) + ...\`.

Returns
-------
fit {coeffs: float[n]}: float[n_coeffs]
  A function that fits the points given by x and y. The coefficients
  and the function components can be retrieved as \`fit.coeffs\` and
  \`fit.funcs\` properties of the function.

References
----------
.. [1] https://en.wikipedia.org/wiki/Linear_regression

Example
-------
>>> const [a,b,c,d] = [3,1,2,4];
... 
... // Underlying function
... const fun = ([x,y]) => a + b*x + c*x*y + d * Math.sin(y);
... 
... // Generate sample points to fit to
... let X = [],
...     Y = [];
... for( let xi = -2; xi <= +2; xi += 0.5 ) {
... for( let yi = -2; yi <= +2; yi += 0.5 ) {
...   X.push(     [xi,yi]  );
...   Y.push( fun([xi,yi]) );
... }}
...
... // Fitted function
... const fit = nd.opt.fit_lin(
...   X,Y,
...   [([x,y]) => Math.sin(y),
...    ([x,y]) => x*y,
...    ([x,y]) => x,
...    ([x,y]) => 1]
... );
... 
... console.log(
...   'Fitted coefficients: [d,c,b,a] =',
...   [...fit.coeffs.map(x => x.toFixed(12))]
... );
... 
... console.log( 'fun([-1,-1]) =',  fun([-1,-1]) );
... console.log( 'fit([-1,-1]) =',  fit([-1,-1]) );
  Fitted coefficients: [d,c,b,a] = [ 4, 2, 1, 3 ]
  fun([-1,-1]) = 0.634116060768414
  fit([-1,-1]) = 0.6341160607684131
`


nd.opt.fit_lm_gen.__doc__ = `\
Nonlinear least squares fit of a function to the given sample points using
the Levenberg-Marquardt trust region method.

Parameters
----------
x: float[n_samples,n_inputs]
  Function inputs of the sample points that the function is fit through.
y: float[n_samples]
  Function outputs of the sample points that the function is fit through.
fg: (params: float[nParam]) => (x: float[nDim]) => [f: float, g: float[nParam]]
  Functions whose parameters are fitted. Returns the function value and its
  gradients with respect to its parameters.
p0: float[nParam]
  Starting values for the parameters for fitting.
opt : {
  r0: float
    [optional] Starting value for trust region radius.
  rMin: float
    [optional] Lower trust region radius bounds.
  rMax: float
    [optional] Upper trust region radius bounds.
  rTol: float
    [optional] Relative tolerance by which the iteration adheres to the trust region radius.
  lmLower: float ∈ (0,1)
    [optional] A small number that avoids that the Levenberg-Marquardt parameter becomes
    too small too quickly during iteration.
  shrinkLower: float ∈ (0,shrinkUpper]
    [optional] During a single iteration, the trust region radius may at most decrease by this factor.
  shrinkUpper: float ∈ [shrinkLower,1)
    [optional] During a single iteration, the trust region radius may at least decrease by this factor.
  grow: float
    [optional] During a single iteration, the trust region radius may increase by this factor.
  expectGainMin: float
    [optional] The iteration is expected to provide at least \`expectGainMin\` of the predicted improvement.
    Otherwise the trust region is shrunk.
  expectGainMax: float
    [optional] The iteration is expected to provide at most \`expectGainMax\` of the predicted improvement.
    Otherwise the trust region is increased.
}

Returns
-------
Iterator<[
  p: float[nParam]
    Current parameter values of the iteration.
  mse: float
    Mean squared error.
  mse_grad: float[nParam]
    Gradient of the mean squared error with respect to the parameters.
  res: float[]
    Residuals of the fit.
  res_jac: float[]
    Jaciobian of the residuals with respect to the parameters.
]>

References
----------
.. [1] https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm

Example
-------
>>> // Determine half-life of beer using dataset from:
... // "Multivariate Analyses of Beer Foam Stand" by James J. Hackbart
... const time = nd.array([[ 0,    15,    30,    45,    60,   90,   120,   150,   180,   210,   240,   270,   300   ]]).T,
...     height = nd.array( [17.40, 15.10, 13.10, 11.60, 10.60, 8.70,  7.40,  6.35,  5.40,  4.50,  3.80,  3.30,  2.90] );
... 
... // Exponential decay function.
... const fG = ([H0, c]) => ([t]) => [
...    H0 * 2**(-t*c),
...   [     2**(-t*c),                    // <- d(fG) / d(H0)
...    H0 * 2**(-t*c) * Math.log(2) * -t] // <- d(fG) / d(c)
... ];
... 
... for( const [[H0, c], mse, mse_grad] of nd.opt.fit_lm_gen(time, height, fG, /*p0=*/[1,1]) )
...   if( nd.la.norm(mse_grad) <= 1e-6 ) {
...     console.log({H0, half_life: 1/c});
...     break;
...   }
  { H0: 16.386111806695105, half_life: 108.02447348677644 }
`


  //
 // I/O
//
nd.io.istr_parse.__doc__ = `\
Parses an inlineable string representation (ISTR) of an \`NDArray\`. ISTR is intended
for inlining large \`NDArray\`s into JavaScript code or JSON files. The string representation
consists of a data type string, followed by a JSON array for the shape and finally a Base64
encoded, little-endian representation of the \`NDArray\`'s data/content. Except for within
the type string or in between the digits of a shape integer, ISTR is whitespace-agnostic.

Parameters
----------
chars: Iterable<char>
  The character sequence that is to be parsed.

Returns
-------
ndarray: NDArray
  The \`NDArray\` parsed from \`chars\`.

Examples
--------
>>> const a = nd.io.istr_parse(\`int32[7,6]
...   CwAAAAwAAAANAAAADgAAAA8AAAAQAAAAFQAAABYAAAAXAAAAGAAAABkAAAAaAAAAHwAAACAAAAAhAAAAIgAAACMAAAAkAAAAKQAAACoAAAArAAAALAAAAC0AAAAuAAAA
...   MwAAADQAAAA1AAAANgAAADcAAAA4AAAAPQAAAD4AAAA/AAAAQAAAAEEAAABCAAAARwAAAEgAAABJAAAASgAAAEsAAABMAAAA\`
... );
>>> console.log( a.toString() )
  [[ 11, 12, 13, 14, 15, 16 ],
   [ 21, 22, 23, 24, 25, 26 ],
   [ 31, 32, 33, 34, 35, 36 ],
   [ 41, 42, 43, 44, 45, 46 ],
   [ 51, 52, 53, 54, 55, 56 ],
   [ 61, 62, 63, 64, 65, 66 ],
   [ 71, 72, 73, 74, 75, 76 ]]
`



nd.io.istr_stringify.__doc__ = `\
Converts an \`NDArray\` into an inlineable string representation (ISTR). ISTR is intended
for inlining large \`NDArray\`s into JavaScript code or JSON files. The string representation
consists of a data type string, followed by a JSON array for the shape and finally a Base64
encoded, little-endian representation of the \`NDArray\`'s data/content. Except for within
the type string or in between the digits of a shape integer, ISTR is whitespace-agnostic.

Parameters
----------
ndarray: NDArray
  The array that is to be converted to an inlineable string representation.
options: {linewidth=128}
  [optional] The maximum number of Base64 characters per line before a newline
  character is inserted.

Returns
-------
istr: string
  An inlineable string representation of \`ndarray\`.

Examples
--------
>>> const a = nd.tabulate([7,6], 'int32', (i,j) => 10*i+j+11)
>>> console.log('\`' + nd.io.istr_stringify(a) + '\`')
  \`int32[7,6]
  CwAAAAwAAAANAAAADgAAAA8AAAAQAAAAFQAAABYAAAAXAAAAGAAAABkAAAAaAAAAHwAAACAAAAAhAAAAIgAAACMAAAAkAAAAKQAAACoAAAArAAAALAAAAC0AAAAuAAAA
  MwAAADQAAAA1AAAANgAAADcAAAA4AAAAPQAAAD4AAAA/AAAAQAAAAEEAAABCAAAARwAAAEgAAABJAAAASgAAAEsAAABMAAAA\`
`



nd.io.npy_serialize.__doc__ = `\
Serializes an \`NDArray\` into an \`Uint8Array\` using version 1.0 of the NumPy NPY-format.
NPY is popular matrix exchange format that is implemented in several programming languages.
It is well-suited for network streaming or file system storage.

In order to reduce memory overhead, consider using \`nd.io.npy_serialize_gen()\` which returns
an \`Iterable<uint8>\` instead of an \`Uint8Array\`.

Parameters
----------
ndarray: NDArray
  The array that is to be serialized.

Returns
-------
bytes: Uint8Array
  The bytes of the NPY representation of \`ndarray\`.
`



  //
 // LINEAR ALGEBRA
//
nd.la.__doc__ = `\
The Linear Algebra subpackage of ND.JS. Contains methods for the handling and decomposition
of matrices as well as solvers for linear equations and linear least square systems.

As convention the last two axes of an \`NDArray\` are considered as the matrix dimensions.
Higher dimensional arrays are considered "arrays of matrices" and the broadcasting rules
apply for the leading dimenions, i.e. for all dimensions of the shape except the last two.
`



nd.la.eye.__doc__ = `\
Returns an NDArray of identity matrices, where the last two axes are the matrix dimensions.

Parameters
----------
shape: ...int
  The shape of the resulting array. If only one value N is given, an N*N identity
  matrix is returned.

Returns
-------
I: NDArray
  The array of identity matrices, where the last two axes are the matrix dimensions.

Examples
--------
>>> const I = nd.la.eye(2);
>>> console.log( I.toString() );
  [[1,0],
   [0,1]]

>>> const J = nd.la.eye(2,3,4);
>>> console.log( J.toString() );
  [[[ 1, 0, 0, 0],
    [ 0, 1, 0, 0],
    [ 0, 0, 1, 0]],
   [[ 1, 0, 0, 0],
    [ 0, 1, 0, 0],
    [ 0, 0, 1, 0]]]
`



nd.la.tril.__doc__ = `\
Returns a copy of an array with all entries outside of the lower triangular part set to 0.
The last two axes are considered the matrix dimensions.

Parameters
----------
A: NDArray[...,N,M]
  The array of matrices whose lower triangular parts are to be returned.
offset: int
  [Optional] The offset from the main diagonal of entries included in the  result.
  If \`offset = +1\` the entries above/right of the main diagonal are not set to zero.
  If \`offset = -2\` the entries on the main diagonal and one below/left of if are set to zero.
  The default value is 0.

Returns
-------
L: NDArray[...,N,M]
  The lower triangular part of A, where \`L(i,j) == i >= j-offset ? A(i,j) : 0\`.

Examples
--------
>>> const A = nd.tabulate([2,3,3], (k,i,j) => 100*(k+1) + 10*(i+1) + (j+i) );
>>> console.log( A.toString() );
  [[[ 110, 111, 112 ],
    [ 121, 122, 123 ],
    [ 132, 133, 134 ]],
   [[ 210, 211, 212 ],
    [ 221, 222, 223 ],
    [ 232, 233, 234 ]]]

>>> const L1 = nd.la.tril(A);
>>> console.log( L1.toString() );
  [[[ 110,   0,   0 ],
    [ 121, 122,   0 ],
    [ 132, 133, 134 ]],
   [[ 210,   0,   0 ],
    [ 221, 222,   0 ],
    [ 232, 233, 234 ]]]
>>> const L2 = nd.la.tril(A,-1);
>>> console.log( L2.toString() );
  [[[   0,   0,  0 ],
    [ 121,   0,  0 ],
    [ 132, 133,  0 ]],
   [[   0,   0,  0 ],
    [ 221,   0,  0 ],
    [ 232, 233,  0 ]]]
`



nd.la.triu.__doc__ = `\
Returns a copy of an array with all entries outside of the upper triangular part set to 0.
The last two axes are considered the matrix dimensions.

Parameters
----------
A: NDArray[...,N,M]
  The array of matrices whose upper triangular parts are to be returned.
offset: int
  [Optional] The offset from the main diagonal of entries included in the result.
  If \`offset = +1\` the main diagonal is set to zero.
  If \`offset = -2\` the entries on the main diagonal and one below/left of if are not set to zero.
  The default value is 0.

Returns
-------
U: NDArray[...,N,M]
  The upper triangular part of A, where \`L(i,j) == i <= j-offset ? A(i,j) : 0\`.

Examples
--------
>>> const A = nd.tabulate([2,3,3], (k,i,j) => 100*(k+1) + 10*(i+1) + (j+i) );
>>> console.log( A.toString() );
  [[[ 110, 111, 112 ],
    [ 121, 122, 123 ],
    [ 132, 133, 134 ]],
   [[ 210, 211, 212 ],
    [ 221, 222, 223 ],
    [ 232, 233, 234 ]]]

>>> const L1 = nd.la.triu(A);
>>> console.log( L1.toString() );
[[[ 110, 111, 112 ],
  [   0, 122, 123 ],
  [   0,   0, 134 ]],
 [[ 210, 211, 212 ],
  [   0, 222, 223 ],
  [   0,   0, 234 ]]]
>>> const L2 = nd.la.triu(A,+1);
>>> console.log( L2.toString() );
[[[ 0, 111, 112 ],
  [ 0,   0, 123 ],
  [ 0,   0,   0 ]],
 [[ 0, 211, 212 ],
  [ 0,   0, 223 ],
  [ 0,   0,   0 ]]]
`



nd.la.diag_mat.__doc__ = `\
Creates a (main) diagonal matrix with the given diagonal values.

Parameters
----------
diag: NDArray[...,N]
  The diagonal values.

Returns
-------
D: NDArray[...,N,N]
  Such that \`D(i,j) == i==j ? diag(i) : 0\`.

Examples
--------
>>> const D = nd.la.diag_mat([[1,2,3],[4,5,6]]);
>>> console.log( D.toString() );
  [[[ 1, 0, 0 ],
    [ 0, 2, 0 ],
    [ 0, 0, 3 ]],
   [[ 4, 0, 0 ],
    [ 0, 5, 0 ],
    [ 0, 0, 6 ]]]
`



nd.la.matmul.__doc__ = `\
Computes the matrix product of a series of matrices. The order
of multiplications is optimized to minimize the number of floating
point operations required. The last two axes are considered to be
the matrix dimensions. For the leading axes, broadcasting rules apply.
Along with \`nd.la.matmul2\`, this is one of the only two methods
that currently allows complex inputs.

Parameters
----------
matrices: ...NDArray
  The matrices that are to be multiplied.

Returns
-------
product: NDArray
  Where \`product = matrices[0] @ matrices[1] @ ... @ matrices[n-1]\`.

Examples
--------
>>> const v = nd.tabulate([1e6,1], () => 1);
>>> console.log( v.toString() );
  [[ 1 ],
   [ 1 ],
   [ 1 ],
   [ 1 ],
   [ 1 ],
   ...999990 more...,
   [ 1 ],
   [ 1 ],
   [ 1 ],
   [ 1 ],
   [ 1 ]]
>>> const u = nd.la.matmul(v,v.T,v,v.T,v);
>>> console.log( u.toString() );
  [[1000000000000],
   [1000000000000],
   [1000000000000],
   [1000000000000],
   [1000000000000],
   ...999990 more...,
   [1000000000000],
   [1000000000000],
   [1000000000000],
   [1000000000000],
   [1000000000000]]
`


nd.la.matmul2.__doc__ = `\
Computes the matrix product of exactly two matrices. The last two axes
are considered to be the matrix dimensions. For the leading axes,
broadcasting rules apply. Along with \`nd.la.matmul\`, this is one
of the only two methods that currently allows complex inputs.

Parameters
----------
a: NDArray[...,N,K]
  The left factor matrices.
b: NDArray[...,K,M]
  The right factor matrices.

Returns
-------
product: NDArray[...,N,M]
  Where \`product = a @ b\`.

Examples
--------
>>> const A = [
...   [2, 0],
...   [0,-1]
... ];
>>> const x = nd.array([
...   [[1,2]],
...   [[3,4]],
...   [[5,6]]
... ]).T;
>>> console.log( x.toString() );
  [ [[1],
     [2]],

    [[3],
     [4]],

    [[5],
     [6]] ]
>>> const y = nd.la.matmul2(A,x);
>>> console.log( y.toString() );
  [ [[ 2 ],
     [-2 ]],
  
    [[ 6 ],
     [-4 ]],
  
    [[10 ],
     [-6 ]] ]
`



nd.la.cholesky_decomp.__doc__ = `\
Computes the Cholesky Decomposition of an \`NDArray\` of symmetric, positive
definite (square) matrices. The implementation assume the matrices to be symmetric
and only looks at the lower triangular values. The last two axes of the \`NDArray\`
are considered to be the matrix dimensions.

Parameters
----------
S: NDArray[...,N,N]
  An \`NDArray\` of symmetric, positive definite (square) matrices. In other words
  S is symmetric and has only positive eigenvalues. 

Returns
-------
L: NDArray[...,N,N]
  A lower triangular matrix, such that \`S = L @ L.T\`.

Examples
--------
>>> const S = [
...   [ 25, -50],
...   [-50, 101]
... ];
>>> const L = nd.la.cholesky_decomp(S);
>>> console.log( nd.la.matmul2(L,L.T).toString() );
  [[ 25, -50],
   [-50, 101]]
`



nd.la.cholesky_solve.__doc__ = `\
Given a cholesky decomposition and the right hand side of a linear equation system,
this method computes the result of said system.

Parameters
----------
L: NDArray[...,N,N]
  An array of lower triangular (square) matrices. L is assumed to be lower triangular
  without looking at the upper triangular values.
y: NDArray[...,N,M]
  The right hand side of the linear equations system.

Returns
-------
x: NDArray[...,N,M]
  Such that \`L @ L.T @ x == y\`.

Examples
--------
>>> const S = [
...   [ 25, -50],
...   [-50, 101]
... ];
>>> const y = [
...   [[1],
...    [2]],
...   [[3],
...    [4]]
... ];
>>> const L = nd.la.cholesky_decomp(S);
>>> const x = nd.la.cholesky_solve(L,y);
>>> console.log( nd.la.matmul2(S,x) );
 [[[1.0],
   [2.0]],
  [[3.0],
   [4.0]]]
`



nd.la.tril_solve.__doc__ = `\
Given the lower triangular (square) matrix and the right hand side of a linear
equation system, this method computes the result.

Parameters
----------
L: NDArray[...,N,N]
  The lower triangular matrix of the linear equations system. L is assumed to
  be lower triangular without ever looking at the values in the upper triangular
  region.
y: NDArray[...,N,M]
  The right hand side matrix of the linear equation system.

Returns
-------
x: NDArray[...,N,M]
  Such that \`L @ x = y\`.
`



nd.la.triu_solve.__doc__ = `\
Given the upper triangular (square) matrix and the right hand side of a linear
equation system, this method computes the result.

Parameters
----------
U: NDArray[...,N,N]
  The upper triangular matrix of the linear equations system. U is assumed to
  be upper triangular without ever looking at the values in the lower triangular
  region.
y: NDArray[...,N,M]
  The right hand side matrix of the linear equation system.

Returns
-------
x: NDArray[...,N,M]
  Such that \`U @ x = y\`.
`



nd.la.lu_decomp.__doc__ = `\
Given an \`NDArray\` of square matrices, this method Computes the
LU(P) decomposition with Column Pivotization.

Parameters
----------
A: NDArray[...,N,N]

Returns
-------
LU: NDArray[...,N,N]
  A matrix containing both the lower and upper triangular part of the
  LU decomposition. The diagonal of ther lower triangular matrix L
  is not contained in LU. It only containes ones and is implied.
P : NDArray[...,N]
  The order in which the row indices of A appear in the LU decomposition, i.e.:
  \`A(P[i],j) == (L @ U)[i,j]\`.
`



nd.la.lu_solve.__doc__ = `\
Given the LU(P) decomposition and the right hand side of a
Linear Equations System, this method computes the result of
said system.

Parameters
----------
LU: NDArray[...,N,N]
P : NDArray[...,N]
y : NDArray[...,N,M]

Returns
-------
x : NDArray[...,N,M]
  The solution of the Linear Equation System, such that:
  \`y[P[i],:] == (L @ U @ x)[i,:] \`
`



nd.la.qr_decomp_full.__doc__ = `\
Computes the (full) QR Decomposition of a matrix. The QR
Decomposition can be used to solve both Linear Equations
and Linear Least Square problems with high numeric accuracy.
Under normal circumstances, the incomplete QR Decomposition
(\`nd.la.qr_decomp\`) is to be preferred over this method
as it may be significantly more memory efficient.

Parameters
----------
A: NDArray[...,N,M]
  The matrix for which the QR Decomposition is computed.

Returns
-------
Q: NDArray[...,N,N]
  An orthogonal square matrix, i.e. \`Q @ Q.T == Q.T @ Q == nd.la.eye(N)\`.
R: NDArray[...,N,M]
  An upper triangular matrix, such that: \`A = Q @ R\`
`



nd.la.hessenberg_decomp.__doc__ = `\
Computes the (Upper) Hessenberg Decomposition of a matrix. It
is worth mentioning that the Hessenberg Decomposition of a
symmetric matrix is tridiagonal.

Parameters
----------
A: NDArray[...,N,N]
  The matrix for which the Hessemberg Decomposition is computed.

Returns
-------
U: NDArray[...,N,N]
  An orthogonal matrix, i.e. \`U @ U.T == U.T @ U == nd.la.eye(N)\`.
H: NDArray[...,N,N]
  An upper Hessemberg matrix, i.e. \`nd.la.tril(H,-2) == 0\`, such that \`U @ H @ U.T == A\`.
`


nd.la.eigen.__doc__ = `\
Returns the eigenpairs of a real square matrix.

Parameters
----------
A: NDArray[...,N,N]
  The non-symmetric real squre matrix for which the eigenvalues are computed.

Returns
-------
Λ: NDArray[...,N]
  A matrix containing all eigenvalues.
V: NDArray[...,N,N]
  A matrix containing the eigenvectors corresponding to Λ as columns, i.e.:
  \`Λ[i]*V[j,i] == (A @ V)[j,i]\`.
  The columns are normalized using the 2-norm. 
`



nd.la.schur_decomp.__doc__ = `\
Computes the (real) Schur Decomposition of a matrix. The
Schur Decomposition is used to compute the Eigenvalues and
Eigenvectors.

Parameters
----------
A: NDArray[...,N,N]

Returns
-------
Q: NDArray[...,N,N]
  An orthogonal square matrix.
T: NDArray[...,N,N]
  An upper quasi-triangular matrix, such that \`Q @ T @ Q.T == A\`.
  The diagonal of T consists of 2x2 blocks that represent complex
  eigenpairs and 1x1 matrix that represent a real eigenvalue.
`



nd.la.qr_decomp.__doc__ = `\
Computes the (economic) QR Decomposition of a matrix. The
QR Decomposition can be used to solve both Linear Equations
and Linear Least Square problems with high numeric accuracy.
Under normal circumstances, this method is to
be preferred over the full QR Decomposition (\`nd.la.qr_decomp_full\`)
as it may be significantly more memory efficient.

Parameters
----------
A: NDArray[...,N,M]
  The matrix for which the QR Decomposition is computed.

Returns
-------
Q: NDArray[...,N,min(N,M)]
  An orthogonal rectangular matrix, i.e. \`Q.T @ Q == nd.la.eye(min(N,M))\`.
R: NDArray[...,min(N,M),M]
  An upper triangular matrix, such that: \`A = Q @ R\`
`



nd.la.rrqr_decomp.__doc__ = `\
Computes the economic Rank-Revealing QR Decomposition of a matrix.

Parameters
----------
A: NDArray[...,N,M]
  The matrix for which the QR Decomposition is computed.

Returns
-------
Q: NDArray[...,N,min(N,M)]
  An orthogonal rectangular matrix, i.e. \`Q.T @ Q == nd.la.eye(min(N,M))\`.
R: NDArray[...,min(N,M),M]
  An upper triangular matrix, where \`R[i,i] >= R[j,j]\` if and only if \`i <= j\`.
P: NDArray[...,M]
  The permuted column indices, such that:
  \`(Q @ R)[:,j] == A[:,P[j]]\`
`



nd.la.qr_lstsq.__doc__ = `\
Given the QR Decompostion and the right hand side of a full-rank
linear equation system, this method solves it. If the Linear
Equation System is under-determined, the Linear Least Square
solution is computed. This method is no suited for rank-deficient
or under-determined systems, in which case \`rrqr_lstsq\` or
\`svd_lstsq\` can be used instead.

Parameters
----------
Q: NDArray[...,N,K]
  The orthognal rectangular matrix of the QR Decomposition.
R: NDArray[...,K,M]
  The upper triangular matrix of the QR Decomposition.
y: NDArray[...,N,L]
  The right hand side of the Linear Equation System.

Returns
-------
x: NDArray[...,N,L]
  Such that \`‖(Q @ R @ x) - y‖₂\` is minimal.
`



nd.la.bidiag_decomp.__doc__ = `\
Computes the Upper Bidiagonal Decomposition of a matrix. Bidiagonal
Decomposition is a common preprocessing step for the Singular Value
Decomposition.

Parameters
----------
A: NDArray[...,N,M]
  The matrix for which the Upper Bidiagonal Decomposition is
  computed.

Returns
-------
U: NDArray[..., N, min(N,M)]
  An orthognal rectangular matrix.
B: NDArray[..., min(N,M), N >= M ? M : N+1 ]
  An upper bidiagonal rectangular matrix.
V: NDArray[..., N >= M ? M : N+1, M]
  An orthogonal rectangular matrix. Such that \`A = U @ B @ V @\`
`



nd.la.solve.__doc__ = `\
Solves a full-rank Linear Equations System (LES). The method tries to
detect rank-deficient systems and throw a \`SingularMatrixSolveError\`,
which contains some least squares solution to the rank-deficient
system as property \`x\`.

Parameters
----------
A: NDArray[...,N,K]
  The matrix of the LES.
y: NDArray[...,K,M]
  The right hand side of the LES.

Returns
-------
x: NDArray[...,N,M]
  The solution of \`(A @ x) = y\`.
`



nd.la.svd_decomp.__doc__ = `\
Computes the Singular Value Decomposition (SVD) of a matrix.

Parameters
----------
A: NDArray[...,N,M]
  The matrix for which the SVD is computed.

Returns
-------
U : NDArray[...,N,min(N,M)]
  An orthogonal, rectangular matrix.
sv: NDArray[..., min(N,M) ]
  The singular values of A sorted in descending order.
V : NDArray[...,min(N,M),M]
  An orthogonal, rectangular matrix, such that \`A == U @ diag(sv) @ V\`.
`
