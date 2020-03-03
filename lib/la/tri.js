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

var _interopRequireDefault = require("@babel/runtime/helpers/interopRequireDefault");

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.tril = tril;
exports.triu = triu;
exports._tril_solve = _tril_solve;
exports.tril_solve = tril_solve;
exports._triu_solve = _triu_solve;
exports.triu_solve = triu_solve;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _dt = require("../dt");

var _nd_array = require("../nd_array");

function tril(m) {
  var k = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0;
  m = (0, _nd_array.asarray)(m);
  if (m.ndim < 2) throw new Error('Input must be at least 2D.');
  return m.mapElems(m.dtype, function (x) {
    for (var _len = arguments.length, indices = new Array(_len > 1 ? _len - 1 : 0), _key = 1; _key < _len; _key++) {
      indices[_key - 1] = arguments[_key];
    }

    var _indices$slice = indices.slice(-2),
        _indices$slice2 = (0, _slicedToArray2["default"])(_indices$slice, 2),
        i = _indices$slice2[0],
        j = _indices$slice2[1];

    return i < j - k ? 0 : x;
  });
}

function triu(m) {
  var k = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0;
  m = (0, _nd_array.asarray)(m);
  if (m.ndim < 2) throw new Error('Input must be at least 2D.');
  return m.mapElems(m.dtype, function (x) {
    for (var _len2 = arguments.length, indices = new Array(_len2 > 1 ? _len2 - 1 : 0), _key2 = 1; _key2 < _len2; _key2++) {
      indices[_key2 - 1] = arguments[_key2];
    }

    var _indices$slice3 = indices.slice(-2),
        _indices$slice4 = (0, _slicedToArray2["default"])(_indices$slice3, 2),
        i = _indices$slice4[0],
        j = _indices$slice4[1];

    return i > j - k ? 0 : x;
  });
}

function _tril_solve(M, N, O, L, L_off, X, X_off) {
  if (!(M <= N)) throw new Error('Assertion failed.');
  M |= 0;
  N |= 0;
  O |= 0;
  L_off |= 0;
  X_off |= 0; // FORWARD SUBSTITUTION

  for (var i = 0; i < M; i++) {
    for (var k = 0; k < i; k++) {
      for (var j = 0; j < O; j++) {
        X[X_off + O * i + j] -= L[L_off + N * i + k] * X[X_off + O * k + j];
      }
    }

    for (var _j = 0; _j < O; _j++) {
      X[X_off + O * i + _j] /= L[L_off + N * i + i];
    }
  }
}

function tril_solve(L, Y) {
  L = (0, _nd_array.asarray)(L);
  if (L.ndim < 2) throw new Error('tril_solve(L,Y): L.ndim must be at least 2.');
  Y = (0, _nd_array.asarray)(Y);
  if (Y.ndim < 2) throw new Error('tril_solve(L,Y): Y.ndim must be at least 2.');

  var _Y$shape$slice = Y.shape.slice(-2),
      _Y$shape$slice2 = (0, _slicedToArray2["default"])(_Y$shape$slice, 2),
      M = _Y$shape$slice2[0],
      N = _Y$shape$slice2[1];

  if (L.shape[L.ndim - 2] !== M) throw new Error("tril_solve(L,Y): L and Y don't match.");
  if (L.shape[L.ndim - 1] !== M) throw new Error('tril_solve(L,Y): Last two dimensions of L must be quadratic.');
  var ndim = Math.max(L.ndim, Y.ndim),
      L_shape = L.shape,
      Y_shape = Y.shape,
      X_shape = new Int32Array(ndim);
  X_shape.fill(1, 0, -2);
  X_shape[ndim - 2] = M;
  X_shape[ndim - 1] = N; // FIND COMMON (BROADCASTED) SHAPE

  for (var _i = 0, _arr = [L_shape, Y_shape]; _i < _arr.length; _i++) {
    var shape = _arr[_i];

    for (var i = ndim - 2, j = shape.length - 2; i-- > 0 && j-- > 0;) {
      if (1 === X_shape[i]) X_shape[i] = shape[j];else if (X_shape[i] != shape[j] && shape[j] != 1) throw new Error('tril_solve(L,Y): L and Y not broadcast-compatible.');
    }
  } // GENERATE RESULT DATA


  var dtype = [L, Y].every(function (A) {
    return A.dtype === 'float32';
  }) ? 'float32' : 'float64';
  L = L.data;
  Y = Y.data;
  var X = new _dt.ARRAY_TYPES[dtype](X_shape.reduce(function (a, b) {
    return a * b;
  }));
  var L_off = 0,
      L_stride = 1,
      Y_off = 0,
      Y_stride = 1,
      X_off = 0;

  function solv(d) {
    if (d === ndim - 2) {
      L_stride = M * M;
      Y_stride = M * N; // COPYING y

      for (var i = Y_stride; i-- > 0;) {
        X[X_off + i] = Y[Y_off + i];
      }

      _tril_solve(M, M, N, L, L_off, X, X_off);

      L_off += L_stride;
      Y_off += Y_stride;
      X_off += Y_stride;
    } else {
      for (var l = X_shape[d];; l--) {
        solv(d + 1);
        if (l == 1) break;
        if (!(L_shape[d - ndim + L_shape.length] > 1)) L_off -= L_stride;
        if (!(Y_shape[d - ndim + Y_shape.length] > 1)) Y_off -= Y_stride;
      }

      L_stride *= L_shape[d - ndim + L_shape.length] || 1;
      Y_stride *= Y_shape[d - ndim + Y_shape.length] || 1;
    }
  }

  solv(0);
  return new _nd_array.NDArray(X_shape, X);
}

function _triu_solve(M, N, O, U, U_off, X, X_off) {
  if (!(M <= N)) throw new Error('Assertion failed.');
  M |= 0;
  N |= 0;
  O |= 0;
  U_off |= 0;
  X_off |= 0;
  if (!(0 <= U_off)) throw new Error('Assertion failed.');
  if (!(0 <= X_off)) throw new Error('Assertion failed.');
  if (!(M * N <= U.length - U_off)) throw new Error('Assertion failed.');
  if (!(M * O <= X.length - X_off)) throw new Error('Assertion failed.'); // BACKWARD SUBSTITUTION

  for (var i = M; i-- > 0;) {
    for (var j = O; j-- > 0;) {
      for (var k = M; --k > i;) {
        X[X_off + O * i + j] -= U[U_off + N * i + k] * X[X_off + O * k + j];
      }

      X[X_off + O * i + j] /= U[U_off + N * i + i];
    }
  }
}

function triu_solve(U, Y) {
  U = (0, _nd_array.asarray)(U);
  if (U.ndim < 2) throw new Error('triu_solve(U,Y): U.ndim must be at least 2.');
  Y = (0, _nd_array.asarray)(Y);
  if (Y.ndim < 2) throw new Error('triu_solve(U,Y): Y.ndim must be at least 2.');

  var _Y$shape$slice3 = Y.shape.slice(-2),
      _Y$shape$slice4 = (0, _slicedToArray2["default"])(_Y$shape$slice3, 2),
      M = _Y$shape$slice4[0],
      N = _Y$shape$slice4[1];

  if (U.shape[U.ndim - 2] !== M) throw new Error("triu_solve(U,Y): U and Y don't match.");
  if (U.shape[U.ndim - 1] !== M) throw new Error('triu_solve(U,Y): Last two dimensions of U must be quadratic.');
  var ndim = Math.max(U.ndim, Y.ndim),
      U_shape = U.shape,
      Y_shape = Y.shape,
      X_shape = new Int32Array(ndim);
  X_shape.fill(1, 0, -2);
  X_shape[ndim - 2] = M;
  X_shape[ndim - 1] = N; // FIND COMMON (BROADCASTED) SHAPE

  for (var _i2 = 0, _arr2 = [U_shape, Y_shape]; _i2 < _arr2.length; _i2++) {
    var shape = _arr2[_i2];

    for (var i = ndim - 2, j = shape.length - 2; i-- > 0 && j-- > 0;) {
      if (1 === X_shape[i]) X_shape[i] = shape[j];else if (X_shape[i] != shape[j] && shape[j] != 1) throw new Error('triu_solve(U,Y): U and Y not broadcast-compatible.');
    }
  } // GENERATE RESULT DATA


  var dtype = [U, Y].every(function (A) {
    return A.dtype === 'float32';
  }) ? 'float32' : 'float64';
  U = U.data;
  Y = Y.data;
  var X = new _dt.ARRAY_TYPES[dtype](X_shape.reduce(function (a, b) {
    return a * b;
  }));
  var U_off = U.length,
      U_stride = 1,
      Y_off = Y.length,
      Y_stride = 1,
      X_off = X.length;

  function solv(d) {
    if (d === ndim - 2) {
      U_stride = M * M;
      Y_stride = M * N;
      U_off -= U_stride;
      Y_off -= Y_stride;
      X_off -= Y_stride; // COPYING y

      for (var i = Y_stride; i-- > 0;) {
        X[X_off + i] = Y[Y_off + i];
      }

      _triu_solve(M, M, N, U, U_off, X, X_off);
    } else {
      for (var l = X_shape[d];; l--) {
        solv(d + 1);
        if (l == 1) break;
        if (!(U_shape[d - ndim + U_shape.length] > 1)) U_off += U_stride;
        if (!(Y_shape[d - ndim + Y_shape.length] > 1)) Y_off += Y_stride;
      }

      U_stride *= U_shape[d - ndim + U_shape.length] || 1;
      Y_stride *= Y_shape[d - ndim + Y_shape.length] || 1;
    }
  }

  solv(0);
  return new _nd_array.NDArray(X_shape, X);
}