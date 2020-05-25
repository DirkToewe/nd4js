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
exports._pldlp_decomp = _pldlp_decomp;
exports.pldlp_decomp = pldlp_decomp;
exports.pldlp_l = pldlp_l;
exports.pldlp_d = pldlp_d;
exports.pldlp_p = pldlp_p;
exports._pldlp_solve = _pldlp_solve;
exports.pldlp_solve = pldlp_solve;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _dt = require("../dt");

var _kahan_sum = require("../kahan_sum");

var _nd_array = require("../nd_array");

// REFERENCES
// ----------
// .. [1] "DSYTF2.f"
//         Reference-LAPACK v3.9.0
//         https://github.com/Reference-LAPACK/lapack/blob/v3.9.0/SRC/dsytf2.f
// .. [2] "Matrix Computations" 4th Edition
//         Chapter 4   "Special Linear Systems"
//         Section 4.4 "Symmetric Indefinite Systems"
//         pp. 186ff
//         Hindustan Book Agency, 2015
var _pldlp_decomp_1x1 = function _pldlp_decomp_1x1(M, N, LD, LD_off, k) {
  if (0 !== M % 1) throw new Error('Assertion failed.');
  if (0 !== N % 1) throw new Error('Assertion failed.');
  if (0 !== k % 1) throw new Error('Assertion failed.');
  if (0 !== LD_off % 1) throw new Error('Assertion failed.');
  if (0 !== LD.length % 1) throw new Error('Assertion failed.');
  if (!(0 <= M)) throw new Error('Assertion failed.');
  if (!(M <= N)) throw new Error('Assertion failed.');
  if (!(k >= 0)) throw new Error('Assertion failed.');
  if (!(k < M)) throw new Error('Assertion failed.');
  var D_kk = LD[LD_off + N * k + k];

  for (var i = k; ++i < M;) {
    var LD_ik = LD[LD_off + N * i + k] / D_kk;

    for (var j = k; j++ < i;) {
      LD[LD_off + N * i + j] -= LD_ik * LD[LD_off + N * j + k];
    }
  }

  for (var _i = k; ++_i < M;) {
    LD[LD_off + N * _i + k] /= D_kk;
  }
};

var _pldlp_decomp_2x2 = function _pldlp_decomp_2x2(M, N, LD, LD_off, k) {
  if (0 !== M % 1) throw new Error('Assertion failed.');
  if (0 !== N % 1) throw new Error('Assertion failed.');
  if (0 !== k % 1) throw new Error('Assertion failed.');
  if (0 !== LD_off % 1) throw new Error('Assertion failed.');
  if (0 !== LD.length % 1) throw new Error('Assertion failed.');
  if (!(0 <= M)) throw new Error('Assertion failed.');
  if (!(M <= N)) throw new Error('Assertion failed.');
  if (!(k >= 0)) throw new Error('Assertion failed.');
  if (!(k < M - 1)) throw new Error('Assertion failed.'); // INVERSE OF DIAGONAL 2X2 BLOCK

  var tmp = LD[LD_off + N * (k + 1) + k],
      D00 = LD[LD_off + N * (k + 1) + (k + 1)] / tmp,
      D11 = LD[LD_off + N * k + k] / tmp,
      D10 = 1 / (D00 * D11 - 1) / tmp; // <- TODO: understand this underflow-safe Vodoo magic

  for (var j = k + 2; j < M; j++) {
    var s0 = D10 * (D00 * LD[LD_off + N * j + k] - LD[LD_off + N * j + k + 1]),
        s1 = D10 * (D11 * LD[LD_off + N * j + k + 1] - LD[LD_off + N * j + k]);

    for (var i = j; i < M; i++) {
      LD[LD_off + N * i + j] -= LD[LD_off + N * i + k] * s0 + LD[LD_off + N * i + k + 1] * s1;
    }

    LD[LD_off + N * j + k] = s0;
    LD[LD_off + N * j + k + 1] = s1;
  }
};

var _pldlp_swap = function _pldlp_swap(M, N, LD, LD_off, P, P_off, k, l) {
  if (!(k <= l)) throw new Error('Assertion failed: ' + JSON.stringify({
    k: k,
    l: l
  }));
  if (!(M <= N)) throw new Error('Assertion failed.');
  var P_k = P[P_off + k];
  P[P_off + k] = P[P_off + l];
  P[P_off + l] = P_k;

  for (var i = 0; i < k; i++) {
    var LD_ki = LD[LD_off + N * k + i];
    LD[LD_off + N * k + i] = LD[LD_off + N * l + i];
    LD[LD_off + N * l + i] = LD_ki;
  }

  for (var _i2 = k; ++_i2 < l;) {
    var LD_ik = LD[LD_off + N * _i2 + k];
    LD[LD_off + N * _i2 + k] = LD[LD_off + N * l + _i2];
    LD[LD_off + N * l + _i2] = LD_ik;
  }

  var LD_kk = LD[LD_off + N * k + k];
  LD[LD_off + N * k + k] = LD[LD_off + N * l + l];
  LD[LD_off + N * l + l] = LD_kk;

  for (var _i3 = l; ++_i3 < M;) {
    var _LD_ik = LD[LD_off + N * _i3 + k];
    LD[LD_off + N * _i3 + k] = LD[LD_off + N * _i3 + l];
    LD[LD_off + N * _i3 + l] = _LD_ik;
  }
};

function _pldlp_decomp(M, N, LD, LD_off, P, P_off) {
  if (!(M <= N)) throw new Error('Assertion failed.');

  for (var i = 0; i < M; i++) {
    P[P_off + i] = i;
  }

  var α = (Math.sqrt(17) + 1) / 8;

  for (var k = 0; k < M; k++) {
    var is1x1 = true;

    if (k < M - 1) {
      var r = void 0,
          A_rk = -Infinity,
          A_kk = Math.abs(LD[LD_off + N * k + k]);

      for (var _i4 = k; ++_i4 < M;) {
        var A_ik = Math.abs(LD[LD_off + N * _i4 + k]);

        if (!(A_rk >= A_ik)) {
          // <- handles NaN
          A_rk = A_ik;
          r = _i4;
        }
      }

      if (!(0 < Math.max(A_rk, A_kk))) throw new Error('_pldlp_decomp(M,N, LD,LD_off, P,P_off): Zero column or NaN encountered.');

      if (A_kk < α * A_rk) {
        var s = void 0,
            A_sr = -Infinity;

        for (var _i5 = k; _i5 < r; _i5++) {
          var A_ri = Math.abs(LD[LD_off + N * r + _i5]);

          if (!(A_sr >= A_ri)) {
            // <- handles NaN
            A_sr = A_ri;
            s = _i5;
          }
        }

        for (var _i6 = r; ++_i6 < M;) {
          var A_ir = Math.abs(LD[LD_off + N * _i6 + r]);

          if (!(A_sr >= A_ir)) {
            // <- handles NaN
            A_sr = A_ir;
            s = _i6;
          }
        }

        if (r === s) throw new Error('Assertion failed.');

        if (A_kk < α * A_rk * (A_rk / A_sr)) {
          var A_rr = Math.abs(LD[LD_off + N * r + r]);

          if (A_rr < α * A_sr) {
            is1x1 = false;
            P[P_off + r] ^= -1;
            ++k;
          }

          if (k !== r) _pldlp_swap(M, N, LD, LD_off, P, P_off, k, r);
        }
      }
    }

    if (is1x1) _pldlp_decomp_1x1(M, N, LD, LD_off, k);else _pldlp_decomp_2x2(M, N, LD, LD_off, k - 1);
  }
}

function pldlp_decomp(S) {
  S = (0, _nd_array.asarray)(S);

  var DTypeArray = _dt.ARRAY_TYPES[S.dtype === 'float32' ? 'float32' : 'float64'],
      LD_shape = S.shape,
      P_shape = LD_shape.slice(0, -1),
      _LD_shape$slice = LD_shape.slice(-2),
      _LD_shape$slice2 = (0, _slicedToArray2["default"])(_LD_shape$slice, 2),
      N = _LD_shape$slice2[0],
      M = _LD_shape$slice2[1];

  S = S.data;
  if (N !== M) throw new Error('Last two dimensions must be quadratic.');
  var LD = new DTypeArray(S.length),
      P = new Int32Array(S.length / N);

  for (var LD_off = 0; LD_off < LD.length; LD_off += N * N) {
    for (var i = 0; i < N; i++) {
      for (var j = 0; j <= i; j++) {
        LD[LD_off + N * i + j] = S[LD_off + N * i + j];
      }
    }

    var P_off = LD_off / N;

    _pldlp_decomp(N, N, LD, LD_off, P, P_off);
  }

  return [new _nd_array.NDArray(LD_shape, LD), new _nd_array.NDArray(P_shape, P)];
}

function pldlp_l(LD, P) {
  if (null == P) {
    var _LD = LD;

    var _LD2 = (0, _slicedToArray2["default"])(_LD, 2);

    LD = _LD2[0];
    P = _LD2[1];
  }

  LD = (0, _nd_array.asarray)(LD);
  P = (0, _nd_array.asarray)(P);
  if (LD.ndim < 2) throw new Error('pldlp_solve(LD,P,y): LD must be at least 2D.');
  if (P.ndim < 1) throw new Error('pldlp_solve(LD,P,y): P must be at least 1D.');
  var N = LD.shape[LD.ndim - 2];
  if (N !== LD.shape[LD.ndim - 1]) throw new Error('pldlp_solve(LD,P,y): LD must be square.');
  if (N !== P.shape[P.ndim - 1]) throw new Error('pldlp_solve(LD,P,y): LD and P shape mismatch.');
  var ndim = Math.max(LD.ndim, P.ndim + 1),
      shape = Int32Array.from({
    length: ndim
  }, function () {
    return 1;
  });
  shape[ndim - 2] = N;
  shape[ndim - 1] = N; // FIND COMMON (BROADCASTED) SHAPE

  for (var _i7 = 0, _arr = [LD.shape.slice(0, -2), P.shape.slice(0, -1)]; _i7 < _arr.length; _i7++) {
    var shp = _arr[_i7];

    for (var i = ndim - 2, j = shp.length; i-- > 0 && j-- > 0;) {
      if (1 === shape[i]) shape[i] = shp[j];else if (shape[i] !== shp[j] && shp[j] !== 1) throw new Error('Shapes are not broadcast-compatible.');
    }
  } // GENERATE RESULT DATA


  var DType = [LD, P].every(function (x) {
    return x.dtype === 'float32';
  }) ? 'float32' : 'float64',
      DTypeArray = _dt.ARRAY_TYPES[DType],
      L_dat = new DTypeArray(shape.reduce(function (a, b) {
    return a * b;
  }, 1)),
      LD_dat = LD.data,
      P_dat = P.data;
  var LD_off = 0,
      LD_stride = 1,
      P_off = 0,
      P_stride = 1,
      L_off = 0;

  function solv(d) {
    if (d === ndim - 2) {
      LD_stride = N * N;
      P_stride = N;

      for (var _i8 = 0; _i8 < N; _i8++) {
        var I = _i8 - (P_dat[P_off + _i8] < 0);

        for (var _j = 0; _j < I; _j++) {
          L_dat[L_off + N * _i8 + _j] = LD_dat[LD_off + N * _i8 + _j];
        }

        L_dat[L_off + N * _i8 + _i8] = 1; //        L_dat[L_off + N*i+i] = LD_dat[LD_off + N*i+i];
        //        if( P[P_off + i] < 0 )
        //          L_dat[ L_off + N* i   + i-1] =
        //          L_dat[ L_off + N*(i-1)+ i  ] =
        //         LD_dat[LD_off + N* i   + i-1];
      }

      LD_off += LD_stride;
      L_off += LD_stride;
      P_off += P_stride;
      return;
    }

    for (var l = shape[d];; l--) {
      solv(d + 1);
      if (l === 1) break;
      if (!(LD.shape[d - ndim + LD.ndim] > 1)) LD_off -= LD_stride;
      if (!(P.shape[d - ndim + P.ndim + 1] > 1)) P_off -= P_stride;
    }

    LD_stride *= LD.shape[d - ndim + LD.ndim] || 1;
    P_stride *= P.shape[d - ndim + P.ndim + 1] || 1;
  }

  solv(0);
  return new _nd_array.NDArray(shape, L_dat);
}

function pldlp_d(LD, P) {
  if (null == P) {
    var _LD3 = LD;

    var _LD4 = (0, _slicedToArray2["default"])(_LD3, 2);

    LD = _LD4[0];
    P = _LD4[1];
  }

  LD = (0, _nd_array.asarray)(LD);
  P = (0, _nd_array.asarray)(P);
  if (LD.ndim < 2) throw new Error('pldlp_solve(LD,P,y): LD must be at least 2D.');
  if (P.ndim < 1) throw new Error('pldlp_solve(LD,P,y): P must be at least 1D.');
  var N = LD.shape[LD.ndim - 2];
  if (N !== LD.shape[LD.ndim - 1]) throw new Error('pldlp_solve(LD,P,y): LD must be square.');
  if (N !== P.shape[P.ndim - 1]) throw new Error('pldlp_solve(LD,P,y): LD and P shape mismatch.');
  var ndim = Math.max(LD.ndim, P.ndim + 1),
      shape = Int32Array.from({
    length: ndim
  }, function () {
    return 1;
  });
  shape[ndim - 2] = N;
  shape[ndim - 1] = N; // FIND COMMON (BROADCASTED) SHAPE

  for (var _i9 = 0, _arr2 = [LD.shape.slice(0, -2), P.shape.slice(0, -1)]; _i9 < _arr2.length; _i9++) {
    var shp = _arr2[_i9];

    for (var i = ndim - 2, j = shp.length; i-- > 0 && j-- > 0;) {
      if (1 === shape[i]) shape[i] = shp[j];else if (shape[i] !== shp[j] && shp[j] !== 1) throw new Error('Shapes are not broadcast-compatible.');
    }
  } // GENERATE RESULT DATA


  var DType = [LD, P].every(function (x) {
    return x.dtype === 'float32';
  }) ? 'float32' : 'float64',
      DTypeArray = _dt.ARRAY_TYPES[DType],
      D_dat = new DTypeArray(shape.reduce(function (a, b) {
    return a * b;
  }, 1)),
      LD_dat = LD.data,
      P_dat = P.data;
  var LD_off = 0,
      LD_stride = 1,
      P_off = 0,
      P_stride = 1,
      D_off = 0;

  function solv(d) {
    if (d === ndim - 2) {
      LD_stride = N * N;
      P_stride = N;

      for (var _i10 = 0; _i10 < N; _i10++) {
        D_dat[D_off + N * _i10 + _i10] = LD_dat[LD_off + N * _i10 + _i10];
        if (P_dat[P_off + _i10] < 0) D_dat[D_off + N * _i10 + _i10 - 1] = D_dat[D_off + N * (_i10 - 1) + _i10] = LD_dat[LD_off + N * _i10 + _i10 - 1];
      }

      LD_off += LD_stride;
      D_off += LD_stride;
      P_off += P_stride;
      return;
    }

    for (var l = shape[d];; l--) {
      solv(d + 1);
      if (l === 1) break;
      if (!(LD.shape[d - ndim + LD.ndim] > 1)) LD_off -= LD_stride;
      if (!(P.shape[d - ndim + P.ndim + 1] > 1)) P_off -= P_stride;
    }

    LD_stride *= LD.shape[d - ndim + LD.ndim] || 1;
    P_stride *= P.shape[d - ndim + P.ndim + 1] || 1;
  }

  solv(0);
  return new _nd_array.NDArray(shape, D_dat);
}

function pldlp_p(LD, P) {
  if (null == P) {
    var _LD5 = LD;

    var _LD6 = (0, _slicedToArray2["default"])(_LD5, 2);

    LD = _LD6[0];
    P = _LD6[1];
  }

  LD = (0, _nd_array.asarray)(LD);
  P = (0, _nd_array.asarray)(P);
  if (LD.ndim < 2) throw new Error('pldlp_solve(LD,P,y): LD must be at least 2D.');
  if (P.ndim < 1) throw new Error('pldlp_solve(LD,P,y): P must be at least 1D.');
  var N = LD.shape[LD.ndim - 2];
  if (N !== LD.shape[LD.ndim - 1]) throw new Error('pldlp_solve(LD,P,y): LD must be square.');
  if (N !== P.shape[P.ndim - 1]) throw new Error('pldlp_solve(LD,P,y): LD and P shape mismatch.');

  for (var i = LD.ndim - 2, j = P.ndim - 1; i-- > 0 && j-- > 0;) {
    if (P.shape[j] !== 1 && LD.shape[i] !== 1 && LD.shape[i] !== P.shape[j]) throw new Error('Shapes are not broadcast-compatible.');
  }

  var P_shape = P.shape;
  P = Int32Array.from(P.data);
  LD = undefined;
  var test = new Uint8Array(N);

  for (var P_off = P.length; (P_off -= N) >= 0;) {
    test.fill(false);

    for (var _i11 = N; _i11-- > 0;) {
      var p = P[P_off + _i11];

      if (p < 0) {
        if (!(0 < _i11) || P[P_off + _i11 - 1] < 0) throw new Error('pldlp_solve(LD,P,y): P is invalid.');
        p ^= -1;
      }

      if (!(p >= 0)) throw new Error('pldlp_solve(LD,P,y): P contains out-of-bounds indices.');
      if (!(p < N)) throw new Error('pldlp_solve(LD,P,y): P contains out-of-bounds indices.');
      if (test[p]) throw new Error('pldlp_solve(LD,P,y): P contains duplicate indices.');
      test[p] = true;
      P[P_off + _i11] = p;
    }

    if (test.some(function (x) {
      return !x;
    })) throw new Error('pldlp_solve(LD,P,y): P is invalid.');
  }

  return new _nd_array.NDArray(P_shape, P);
}

function _pldlp_solve(M, N, O, LD, LD_off, P, P_off, X, X_off, tmp) {
  if (0 !== LD.length % 1) throw new Error('Assertion failed.');
  if (0 !== P.length % 1) throw new Error('Assertion failed.');
  if (0 !== X.length % 1) throw new Error('Assertion failed.');
  if (0 !== tmp.length % 1) throw new Error('Assertion failed.');
  if (0 !== LD_off % 1) throw new Error('Assertion failed.');
  if (0 !== P_off % 1) throw new Error('Assertion failed.');
  if (0 !== X_off % 1) throw new Error('Assertion failed.');
  if (0 !== M % 1) throw new Error('Assertion failed.');
  if (0 !== N % 1) throw new Error('Assertion failed.');
  if (0 !== O % 1) throw new Error('Assertion failed.');
  if (LD.length - LD_off < M * N) throw new Error('Assertion failed.');
  if (P.length - P_off < M) throw new Error('Assertion failed.');
  if (X.length - X_off < M * O) throw new Error('Assertion failed.');
  if (tmp.length < M * O) throw new Error('Assertion failed.');
  if (!(0 < M)) throw new Error('Assertion failed.');
  if (!(0 < N)) throw new Error('Assertion failed.');
  if (!(0 < O)) throw new Error('Assertion failed.');
  if (!(M <= N)) throw new Error('Assertion failed.'); // PERMUTE INPUT

  for (var i = 0; i < M; i++) {
    var I = P[P_off + i];
    I ^= -(I < 0);

    for (var j = 0; j < O; j++) {
      tmp[O * i + j] = X[X_off + O * I + j];
    }
  } // FORWARD SUBSTITUTION


  for (var _i12 = 1; _i12 < M; _i12++) {
    var K = _i12 - (P[P_off + _i12] < 0);

    for (var k = 0; k < K; k++) {
      for (var _j2 = 0; _j2 < O; _j2++) {
        tmp[O * _i12 + _j2] -= LD[LD_off + N * _i12 + k] * tmp[O * k + _j2];
      }
    }
  } // SCALING


  for (var _i13 = 0; _i13 < M; _i13++) {
    if (_i13 < M - 1 && P[P_off + _i13 + 1] < 0) {
      // 2x2 BLOCK
      var s = LD[LD_off + N * (_i13 + 1) + _i13],
          D00 = LD[LD_off + N * (_i13 + 1) + _i13 + 1] / s,
          D11 = LD[LD_off + N * _i13 + _i13] / s,
          D10 = 1 / (D00 * D11 - 1) / s; // <- TODO: understand this underflow-safe Vodoo magic (see [1])

      for (var _j3 = 0; _j3 < O; _j3++) {
        var y0 = tmp[O * _i13 + _j3],
            y1 = tmp[O * (_i13 + 1) + _j3];
        tmp[O * _i13 + _j3] = D10 * (D00 * y0 - y1);
        tmp[O * (_i13 + 1) + _j3] = D10 * (D11 * y1 - y0);
      }

      ++_i13;
    } else {
      // 1x1 BLOCK
      for (var _j4 = 0; _j4 < O; _j4++) {
        tmp[O * _i13 + _j4] /= LD[LD_off + N * _i13 + _i13];
      }
    }
  } // BACKWARD SUBSTITUTION


  for (var _k = M; _k-- > 1;) {
    var _I = _k - (P[P_off + _k] < 0);

    for (var _i14 = _I; _i14-- > 0;) {
      for (var _j5 = O; _j5-- > 0;) {
        tmp[O * _i14 + _j5] -= LD[LD_off + N * _k + _i14] * tmp[O * _k + _j5];
      }
    }
  } // (UN)PERMUTE OUTPUT


  for (var _i15 = 0; _i15 < M; _i15++) {
    var _I2 = P[P_off + _i15];
    _I2 ^= -(_I2 < 0);

    for (var _j6 = 0; _j6 < O; _j6++) {
      X[X_off + O * _I2 + _j6] = tmp[O * _i15 + _j6];
    }
  }
}

function pldlp_solve(LD, P, y) {
  if (null == y) {
    if (null == P) throw new Error('Assertion failed.');
    y = P;
    var _LD7 = LD;

    var _LD8 = (0, _slicedToArray2["default"])(_LD7, 2);

    LD = _LD8[0];
    P = _LD8[1];
  }

  LD = (0, _nd_array.asarray)(LD);
  P = (0, _nd_array.asarray)(P);
  y = (0, _nd_array.asarray)(y);
  if (LD.ndim < 2) throw new Error('pldlp_solve(LD,P,y): LD must be at least 2D.');
  if (P.ndim < 1) throw new Error('pldlp_solve(LD,P,y): P must be at least 1D.');
  if (y.ndim < 2) throw new Error('pldlp_solve(LD,P,y): y must be at least 2D.');
  var N = LD.shape[LD.ndim - 2],
      O = y.shape[y.ndim - 1];
  if (N !== LD.shape[LD.ndim - 1]) throw new Error('pldlp_solve(LD,P,y): LD must be square.');
  if (N !== P.shape[P.ndim - 1]) throw new Error('pldlp_solve(LD,P,y): LD and P shape mismatch.');
  if (N !== y.shape[y.ndim - 2]) throw new Error('pldlp_solve(LD,P,y): LD and y shape mismatch.');
  var ndim = Math.max(LD.ndim, P.ndim + 1, y.ndim),
      shape = Int32Array.from({
    length: ndim
  }, function () {
    return 1;
  });
  shape[ndim - 2] = N;
  shape[ndim - 1] = O; // FIND COMMON (BROADCASTED) SHAPE

  for (var _i16 = 0, _arr3 = [LD.shape.slice(0, -2), P.shape.slice(0, -1), y.shape.slice(0, -2)]; _i16 < _arr3.length; _i16++) {
    var shp = _arr3[_i16];

    for (var i = ndim - 2, j = shp.length; i-- > 0 && j-- > 0;) {
      if (1 === shape[i]) shape[i] = shp[j];else if (shape[i] !== shp[j] && shp[j] !== 1) throw new Error('Shapes are not broadcast-compatible.');
    }
  } // GENERATE RESULT DATA


  var DType = [LD, P, y].every(function (x) {
    return x.dtype === 'float32';
  }) ? 'float32' : 'float64',
      DTypeArray = _dt.ARRAY_TYPES[DType],
      tmp = new DTypeArray(N * O),
      x_dat = new DTypeArray(shape.reduce(function (a, b) {
    return a * b;
  }, 1)),
      LD_dat = LD.data,
      P_dat = P.data,
      y_dat = y.data;
  var LD_off = 0,
      LD_stride = 1,
      P_off = 0,
      P_stride = 1,
      y_off = 0,
      y_stride = 1,
      x_off = 0;

  function solv(d) {
    if (d === ndim - 2) {
      LD_stride = N * N;
      P_stride = N;
      y_stride = N * O; // COPY y

      for (var _i17 = 0; _i17 < y_stride; _i17++) {
        x_dat[x_off + _i17] = y_dat[y_off + _i17];
      }

      _pldlp_solve(N, N, O, LD_dat, LD_off, P_dat, P_off, x_dat, x_off, tmp);

      LD_off += LD_stride;
      P_off += P_stride;
      y_off += y_stride;
      x_off += y_stride;
      return;
    }

    for (var l = shape[d];; l--) {
      solv(d + 1);
      if (l === 1) break;
      if (!(LD.shape[d - ndim + LD.ndim] > 1)) LD_off -= LD_stride;
      if (!(P.shape[d - ndim + P.ndim + 1] > 1)) P_off -= P_stride;
      if (!(y.shape[d - ndim + y.ndim] > 1)) y_off -= y_stride;
    }

    LD_stride *= LD.shape[d - ndim + LD.ndim] || 1;
    P_stride *= P.shape[d - ndim + P.ndim + 1] || 1;
    y_stride *= y.shape[d - ndim + y.ndim] || 1;
  }

  solv(0);
  return new _nd_array.NDArray(shape, x_dat);
}