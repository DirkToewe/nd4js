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
exports._ldl_decomp = _ldl_decomp;
exports.ldl_decomp = ldl_decomp;
exports._ldl_solve = _ldl_solve;
exports.ldl_solve = ldl_solve;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _dt = require("../dt");

var _kahan_sum = require("../kahan_sum");

var _nd_array = require("../nd_array");

// THE FOLLOWING IMPLEMENTATION HAS BETTER CACHE ALIGNMENT BUT REQUIRES AN EXTRA ARRAY OF LENGTH M-1
//export function _ldl_decomp(M,N, LD,LD_off, V)
//{
//  if( ! (M   <= N       ) ) throw new Error('Assertion failed.');
//  if( ! (M-1 <= V.length) ) throw new Error('Assertion failed.');
//
//  // https://arxiv.org/abs/1111.4144
//  for( let j=0; j < M; j++ )
//  {
//    for( let k=0; k < j; k++ )
//      V[k] = LD[LD_off + N*j+k] * LD[LD_off + N*k+k];
//
//    for( let i=j; i < M; i++ )
//    {
//      for( let k=0; k < j; k++ )
//        LD[LD_off + N*i+j] -= LD[LD_off + N*i+k] * V[k];
//
//      if( j < i )
//        LD[LD_off + N*i+j] /= LD[LD_off + N*j+j];
//    }
//  }
//}
function _ldl_decomp(M, N, LD, LD_off) {
  if (!(M <= N)) throw new Error('Assertion failed.'); // https://arxiv.org/abs/1111.4144

  for (var j = 0; j < M; j++) {
    for (var k = 0; k < j; k++) {
      var V_k = LD[LD_off + N * j + k] * LD[LD_off + N * k + k];

      for (var i = j; i < M; i++) {
        LD[LD_off + N * i + j] -= LD[LD_off + N * i + k] * V_k;
      }
    }

    for (var _i = j; ++_i < M;) {
      LD[LD_off + N * _i + j] /= LD[LD_off + N * j + j];
    }
  }
}

function ldl_decomp(S) {
  S = (0, _nd_array.asarray)(S);

  var DTypeArray = _dt.ARRAY_TYPES[S.dtype === 'float32' ? 'float32' : 'float64'],
      shape = S.shape,
      _shape$slice = shape.slice(-2),
      _shape$slice2 = (0, _slicedToArray2["default"])(_shape$slice, 2),
      N = _shape$slice2[0],
      M = _shape$slice2[1];

  S = S.data;
  if (N !== M) throw new Error('Last two dimensions must be quadratic.');
  var LD = new DTypeArray(S.length);

  for (var LD_off = 0; LD_off < LD.length; LD_off += N * N) {
    for (var i = 0; i < N; i++) {
      for (var j = 0; j <= i; j++) {
        LD[LD_off + N * i + j] = S[LD_off + N * i + j];
      }
    }

    _ldl_decomp(N, N, LD, LD_off);
  }

  return new _nd_array.NDArray(shape, LD);
}

function _ldl_solve(M, N, O, LD, LD_off, X, X_off) {
  if (0 !== LD.length % 1) throw new Error('Assertion failed.');
  if (0 !== X.length % 1) throw new Error('Assertion failed.');
  if (0 !== LD_off % 1) throw new Error('Assertion failed.');
  if (0 !== X_off % 1) throw new Error('Assertion failed.');
  if (0 !== M % 1) throw new Error('Assertion failed.');
  if (0 !== N % 1) throw new Error('Assertion failed.');
  if (0 !== O % 1) throw new Error('Assertion failed.');
  if (LD.length - LD_off < M * N) throw new Error('Assertion failed.');
  if (X.length - X_off < M * O) throw new Error('Assertion failed.');
  if (!(0 < M)) throw new Error('Assertion failed.');
  if (!(0 < N)) throw new Error('Assertion failed.');
  if (!(0 < O)) throw new Error('Assertion failed.');
  if (!(M <= N)) throw new Error('Assertion failed.'); // FORWARD SUBSTITUTION

  for (var i = 1; i < M; i++) {
    for (var k = 0; k < i; k++) {
      for (var j = 0; j < O; j++) {
        X[X_off + O * i + j] -= LD[LD_off + N * i + k] * X[X_off + O * k + j];
      }
    }
  } // SCALING


  for (var _i2 = 0; _i2 < M; _i2++) {
    for (var _j = 0; _j < O; _j++) {
      X[X_off + O * _i2 + _j] /= LD[LD_off + N * _i2 + _i2];
    }
  } // BACKWARD SUBSTITUTION


  for (var _k = M; _k-- > 1;) {
    for (var _i3 = _k; _i3-- > 0;) {
      for (var _j2 = O; _j2-- > 0;) {
        X[X_off + O * _i3 + _j2] -= LD[LD_off + N * _k + _i3] * X[X_off + O * _k + _j2];
      }
    }
  }
}

function ldl_solve(LD, y) {
  LD = (0, _nd_array.asarray)(LD);
  y = (0, _nd_array.asarray)(y);
  if (LD.ndim < 2) throw new Error('ldl_solve(LD,y): LD must be at least 2D.');
  if (y.ndim < 2) throw new Error('ldl_solve(LD,y): y must be at least 2D.');

  var _LD$shape$slice = LD.shape.slice(-2),
      _LD$shape$slice2 = (0, _slicedToArray2["default"])(_LD$shape$slice, 2),
      N = _LD$shape$slice2[0],
      M = _LD$shape$slice2[1],
      _y$shape$slice = y.shape.slice(-2),
      _y$shape$slice2 = (0, _slicedToArray2["default"])(_y$shape$slice, 2),
      I = _y$shape$slice2[0],
      J = _y$shape$slice2[1];

  if (N != M) throw new Error('ldl_solve(LD,y): Last two dimensions of LD must be quadratic.');
  if (I != M) throw new Error("ldl_solve(LD,y): LD and y don't match.");
  var ndim = Math.max(LD.ndim, y.ndim),
      shape = Int32Array.from({
    length: ndim
  }, function () {
    return 1;
  });
  shape[ndim - 2] = I;
  shape[ndim - 1] = J; // FIND COMMON (BROADCASTED) SHAPE

  for (var _i4 = 0, _arr = [LD, y]; _i4 < _arr.length; _i4++) {
    var arr = _arr[_i4];

    for (var i = ndim - 2, j = arr.ndim - 2; i-- > 0 && j-- > 0;) {
      if (1 === shape[i]) shape[i] = arr.shape[j];else if (shape[i] != arr.shape[j] && arr.shape[j] != 1) throw new Error('Shapes are not broadcast-compatible.');
    }
  } // GENERATE RESULT DATA


  var dtype = LD.dtype === 'float32' && y.dtype === 'float32' ? 'float32' : 'float64',
      x_dat = new _dt.ARRAY_TYPES[dtype](shape.reduce(function (a, b) {
    return a * b;
  }, 1)),
      LD_dat = LD.data,
      y_dat = y.data;
  var LD_off = 0,
      LD_stride = 1,
      y_off = 0,
      y_stride = 1,
      x_off = 0;

  function solv(d) {
    if (d === ndim - 2) {
      LD_stride = N * N;
      y_stride = N * J; // COPY y

      for (var _i5 = 0; _i5 < y_stride; _i5++) {
        x_dat[x_off + _i5] = y_dat[y_off + _i5];
      }

      _ldl_solve(N, N, J, LD_dat, LD_off, x_dat, x_off);

      LD_off += LD_stride;
      y_off += y_stride;
      x_off += y_stride;
      return;
    }

    for (var l = shape[d];; l--) {
      solv(d + 1);
      if (l == 1) break;
      if (!(LD.shape[d - ndim + LD.ndim] > 1)) LD_off -= LD_stride;
      if (!(y.shape[d - ndim + y.ndim] > 1)) y_off -= y_stride;
    }

    LD_stride *= LD.shape[d - ndim + LD.ndim] || 1;
    y_stride *= y.shape[d - ndim + y.ndim] || 1;
  }

  solv(0);
  return new _nd_array.NDArray(shape, x_dat);
}