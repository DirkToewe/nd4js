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
exports._bidiag_decomp_horiz = _bidiag_decomp_horiz;
exports.bidiag_decomp = bidiag_decomp;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _nd_array = require("../nd_array");

var _dt = require("../dt");

var _giv_rot = require("./_giv_rot");

var _norm4 = require("./norm");

var _transpose_inplace2 = require("./transpose_inplace");

function _bidiag_decomp_vert(M, N, U, U_off, B, V, BV_off) {
  M |= 0;
  N |= 0;
  U_off |= 0;
  BV_off |= 0;
  var NORM = new _norm4.FrobeniusNorm();
  if (M < N) throw new Error('Assertion failed.'); // init V to identity

  for (var i = N; i-- > 0;) {
    V[BV_off + N * i + i] = 1;
  }

  for (var _i = 0; _i < N; _i++) {
    // ELIMINATE ELEMENTS BELOW (i,i) USING GIVENS
    var ii = U_off + N * _i + _i;

    for (var j = _i; ++j < M;) {
      var ji = U_off + N * j + _i;
      var B_ji = U[ji];
      if (0 === B_ji) continue;
      var B_ii = U[ii];

      var _giv_rot_qr2 = (0, _giv_rot._giv_rot_qr)(B_ii, B_ji),
          _giv_rot_qr3 = (0, _slicedToArray2["default"])(_giv_rot_qr2, 3),
          c = _giv_rot_qr3[0],
          s = _giv_rot_qr3[1],
          norm = _giv_rot_qr3[2];

      if (s !== 0) {
        // <- might happen in case of underflow
        if (c < 0) {
          c *= -1;
          s *= -1;
          norm *= -1;
        }

        (0, _giv_rot._giv_rot_rows)(U, N - 1 - _i, ii + 1, ji + 1, c, s);
        U[ii] = norm;
      }

      U[ji] = s;
    }

    if (_i < N - 2) {
      // ELIMINATE ELEMENTS RIGHT OF (i,i+1) USING HOUSEHOLDER
      NORM.reset();

      for (var _j = N; --_j > _i + 1;) {
        NORM.include(U[U_off + N * _i + _j]);
      }

      if (NORM.max === 0) continue;

      var _norm = NORM.resultIncl(U[ii + 1]) * (U[ii + 1] > 0 ? -1 : +1);

      NORM.include(U[ii + 1] -= _norm);
      var max = NORM.max,
          div = Math.sqrt(NORM.sum);

      for (var _j2 = _i; ++_j2 < N;) {
        U[U_off + N * _i + _j2] = U[U_off + N * _i + _j2] / max / div;
      } // apply householder to right of V


      for (var _j3 = 0; _j3 < N; _j3++) {
        var sum = 0;

        for (var k = _i; ++k < N;) {
          sum += V[BV_off + N * _j3 + k] * U[U_off + N * _i + k];
        }

        sum *= 2;

        for (var _k = _i; ++_k < N;) {
          V[BV_off + N * _j3 + _k] -= U[U_off + N * _i + _k] * sum;
        }
      } // apply householder to right of B


      for (var _j4 = _i; ++_j4 < M;) {
        var _sum = 0;

        for (var _k2 = _i; ++_k2 < N;) {
          _sum += U[U_off + N * _j4 + _k2] * U[U_off + N * _i + _k2];
        }

        _sum *= 2;

        for (var _k3 = _i; ++_k3 < N;) {
          U[U_off + N * _j4 + _k3] -= U[U_off + N * _i + _k3] * _sum;
        }
      }

      U[ii + 1] = _norm;
      U.fill(0.0, ii + 2, U_off + N * (_i + 1));
    }
  }

  (0, _transpose_inplace2._transpose_inplace)(N, V, BV_off); // COPY U -> B

  for (var _i2 = N; --_i2 > 0;) {
    B[BV_off + N * _i2 + _i2] = U[U_off + N * _i2 + _i2];
    B[BV_off + N * (_i2 - 1) + _i2] = U[U_off + N * (_i2 - 1) + _i2];
  }

  B[BV_off] = U[U_off]; // COMPUTE U
  // init upper right triangle of U

  for (var _i3 = 0; _i3 < N; _i3++) {
    for (var _j5 = _i3; _j5 < N; _j5++) {
      U[U_off + N * _i3 + _j5] = +(_i3 === _j5);
    }
  }

  for (var _i4 = N; _i4-- > 0;) {
    for (var _j6 = M; --_j6 > _i4;) {
      var _s = U[U_off + N * _j6 + _i4];
      if (0 === _s) continue;
      U[U_off + N * _j6 + _i4] = 0;

      var _c = Math.sqrt((1 + _s) * (1 - _s));

      (0, _giv_rot._giv_rot_rows)(U, N - _i4, U_off + N * _j6 + _i4, U_off + N * _i4 + _i4, _c, _s);
    }
  }
}

function _bidiag_decomp_square(N, U, B, V, off) {
  N |= 0;
  off |= 0;
  var NORM = new _norm4.FrobeniusNorm(); // INIT U & V TO IDENTITY

  for (var i = N; i-- > 0;) {
    U[off + N * i + i] = 1;
  }

  for (var _i5 = N; _i5-- > 0;) {
    V[off + N * _i5 + _i5] = 1;
  }

  for (var _i6 = 0; _i6 < N - 1; _i6++) {
    // ELIMINATE ELEMENTS BELOW (i,i) USING GIVENS
    var ii = off + N * _i6 + _i6;

    for (var j = _i6; ++j < N;) {
      var ji = off + N * j + _i6;
      var B_ji = B[ji];
      if (0 === B_ji) continue;

      var B_ii = B[ii],
          _giv_rot_qr4 = (0, _giv_rot._giv_rot_qr)(B_ii, B_ji),
          _giv_rot_qr5 = (0, _slicedToArray2["default"])(_giv_rot_qr4, 3),
          c = _giv_rot_qr5[0],
          s = _giv_rot_qr5[1],
          _norm2 = _giv_rot_qr5[2];

      B[ji] = 0;
      if (0 === s) continue; // <- can happen due to underflow

      B[ii] = _norm2;
      (0, _giv_rot._giv_rot_rows)(B, N - 1 - _i6, ii + 1, ji + 1, c, s);
      (0, _giv_rot._giv_rot_rows)(U, 1 + j, off + N * _i6, off + N * j, c, s);
    } // ELIMINATE ELEMENTS RIGHT OF (i,i+1) USING HOUSEHOLDER


    NORM.reset();

    for (var _j7 = N; --_j7 > _i6 + 1;) {
      NORM.include(B[off + N * _i6 + _j7]);
    }

    if (NORM.max === 0) continue;
    var norm = NORM.resultIncl(B[ii + 1]) * (B[ii + 1] > 0 ? -1 : +1);
    NORM.include(B[ii + 1] -= norm);
    var max = NORM.max,
        div = Math.sqrt(NORM.sum);

    for (var _j8 = _i6; ++_j8 < N;) {
      B[off + N * _i6 + _j8] = B[off + N * _i6 + _j8] / max / div;
    } // apply householder to right of V


    for (var _j9 = 0; _j9 < N; _j9++) {
      var sum = 0;

      for (var k = _i6; ++k < N;) {
        sum += V[off + N * _j9 + k] * B[off + N * _i6 + k];
      }

      sum *= 2;

      for (var _k4 = _i6; ++_k4 < N;) {
        V[off + N * _j9 + _k4] -= B[off + N * _i6 + _k4] * sum;
      }
    } // apply householder to right of B


    for (var _j10 = _i6; ++_j10 < N;) {
      var _sum2 = 0;

      for (var _k5 = _i6; ++_k5 < N;) {
        _sum2 += B[off + N * _j10 + _k5] * B[off + N * _i6 + _k5];
      }

      _sum2 *= 2;

      for (var _k6 = _i6; ++_k6 < N;) {
        B[off + N * _j10 + _k6] -= B[off + N * _i6 + _k6] * _sum2;
      }
    }

    B[ii + 1] = norm;
    B.fill(0.0, ii + 2, off + N * (_i6 + 1));
  }

  (0, _transpose_inplace2._transpose_inplace)(N, V, off);
  (0, _transpose_inplace2._transpose_inplace)(N, U, off);
}

function _bidiag_decomp_horiz(M, N, U, U_off, B, B_off, V, V_off) {
  M |= 0;
  N |= 0;
  U_off |= 0;
  B_off |= 0;
  V_off = V_off + N | 0;
  if (!(M <= N)) throw Error('Assertion failed.');
  var NORM = new _norm4.FrobeniusNorm(); // INIT U TO IDENTITY

  for (var i = M; i-- > 0;) {
    U[U_off + M * i + i] = 1;
  }

  outer_loop: for (var _i7 = 0; _i7 < M && _i7 < N - 1; _i7++) {
    // ELIMINATE ELEMENTS BELOW (i,i) USING GIVENS
    var ii = V_off + N * _i7 + _i7;

    inner_loop: for (var j = _i7; ++j < M;) {
      var ji = V_off + N * j + _i7;
      var B_ji = V[ji];
      if (0 === B_ji) continue inner_loop;

      var B_ii = V[ii],
          _giv_rot_qr6 = (0, _giv_rot._giv_rot_qr)(B_ii, B_ji),
          _giv_rot_qr7 = (0, _slicedToArray2["default"])(_giv_rot_qr6, 3),
          c = _giv_rot_qr7[0],
          s = _giv_rot_qr7[1],
          _norm3 = _giv_rot_qr7[2];

      V[ji] = 0;
      if (0 === s) continue inner_loop;
      V[ii] = _norm3;
      (0, _giv_rot._giv_rot_rows)(V, N - 1 - _i7, ii + 1, ji + 1, c, s);
      (0, _giv_rot._giv_rot_rows)(U, 1 + j, U_off + M * _i7, U_off + M * j, c, s);
    } // ELIMINATE ELEMENTS RIGHT OF (i,i+1) USING HOUSEHOLDER
    // compute householder vector


    var above = V_off + N * (_i7 - 1); // <- write householder to row (i-1)

    NORM.reset();

    for (var _j11 = N - 1;;) {
      var V_j = V[above + _j11] = V[V_off + N * _i7 + _j11];
      if (--_j11 <= _i7) break;
      NORM.include(V_j);
    }

    if (0 === NORM.max) {
      V[above + (_i7 + 1)] = 0;
      continue outer_loop;
    }

    var norm = NORM.resultIncl(V[above + (_i7 + 1)]) * (V[above + (_i7 + 1)] > 0 ? -1 : +1);
    NORM.include(V[above + (_i7 + 1)] -= norm);
    var max = NORM.max,
        div = Math.sqrt(NORM.sum);

    for (var _j12 = _i7; ++_j12 < N;) {
      V[above + _j12] = V[above + _j12] / max / div;
    } // apply householder to right of B (stored in V)


    for (var _j13 = _i7; ++_j13 < M;) {
      var sum = 0;

      for (var k = _i7; ++k < N;) {
        sum += V[V_off + N * _j13 + k] * V[above + k];
      }

      sum *= 2;

      for (var _k7 = _i7; ++_k7 < N;) {
        V[V_off + N * _j13 + _k7] -= V[above + _k7] * sum;
      }
    }

    V[ii + 1] = norm;
  } // TRANSPOSE U


  (0, _transpose_inplace2._transpose_inplace)(M, U, U_off); // COPY V -> B

  for (var _i8 = M; _i8-- > 0;) {
    B[B_off + (M + 1) * _i8 + (_i8 + 1)] = _i8 < N - 1 ? V[V_off + N * _i8 + (_i8 + 1)] : 0;
    B[B_off + (M + 1) * _i8 + _i8] = V[V_off + N * _i8 + _i8];
  } // COMPUTE V


  for (var _i9 = Math.min(M, N - 1);;) {
    --_i9;
    V.fill(0.0, V_off + N * _i9, V_off + N * (_i9 + 1));
    V[V_off + N * _i9 + (_i9 + 1)] = 1;
    if (_i9 < 0) break;

    var _above = V_off + N * (_i9 - 1); // <- read householder from row (i-1)
    // apply householder to right of V


    for (var _j14 = _i9; _j14 < M; _j14++) {
      var _sum3 = 0;

      for (var _k8 = _i9; ++_k8 < N;) {
        _sum3 += V[V_off + N * _j14 + _k8] * V[_above + _k8];
      }

      _sum3 *= 2;

      for (var _k9 = _i9; ++_k9 < N;) {
        V[V_off + N * _j14 + _k9] -= V[_above + _k9] * _sum3;
      }
    }
  }
}

function bidiag_decomp(A) {
  A = (0, _nd_array.asarray)(A);
  if (A.ndim < 2) throw new Error('bidiag_decomp(A): A must be at least 2D.');
  if (A.dtype.startsWith('complex')) throw new Error('bidiag_decomp(A): complex A not yet supported.');

  var DType = A.dtype === 'float32' ? 'float32' : 'float64',
      DTypeArray = _dt.ARRAY_TYPES[DType],
      U_shape = Int32Array.from(A.shape),
      B_shape = Int32Array.from(U_shape),
      V_shape = Int32Array.from(U_shape),
      _A$shape$slice = A.shape.slice(-2),
      _A$shape$slice2 = (0, _slicedToArray2["default"])(_A$shape$slice, 2),
      M = _A$shape$slice2[0],
      N = _A$shape$slice2[1],
      len = A.data.length / (M * N) | 0,
      I = Math.min(M, N) | 0,
      J = (M >= N ? I : I + 1) | 0; // #columns of bidiag. matrix B


  U_shape[U_shape.length - 1] = I;
  B_shape[B_shape.length - 2] = I;
  B_shape[B_shape.length - 1] = J;
  V_shape[V_shape.length - 2] = J;

  var _ref = function () {
    if (M === N) {
      var _B = DTypeArray.from(A.data),
          _U = new DTypeArray(_B.length),
          _V = new DTypeArray(_B.length);

      for (var off = 0; off < _B.length; off += N * N) {
        _bidiag_decomp_square(N, _U, _B, _V, off);
      }

      return [_U, _B, _V];
    } else if (M > N) {
      var _U2 = DTypeArray.from(A.data);

      A = undefined;

      var _B2 = new DTypeArray(len * N * N),
          _V2 = new DTypeArray(len * N * N);

      for (var U_off = 0, BV_off = 0; BV_off < _B2.length; U_off += M * N, BV_off += N * N) {
        _bidiag_decomp_vert(M, N, _U2, U_off, _B2, _V2, BV_off);
      }

      return [_U2, _B2, _V2];
    } else {
      /*M < N*/
      var _V3 = new DTypeArray(len * (M + 1) * N);
      /*DEBUG*/


      if (!(M < N)) throw new Error('Assertion failed.');
      A = A.data;

      for (var j = 0, i = N; i < _V3.length; i += N) {
        for (var _I = i + M * N; i < _I; j++, i++) {
          _V3[i] = A[j];
        }
      }

      A = undefined;

      var _U3 = new DTypeArray(len * M * M),
          _B3 = new DTypeArray(len * M * (M + 1));

      for (var _U_off = 0, B_off = 0, V_off = 0; B_off < _B3.length; _U_off += M * M, B_off += M * (M + 1), V_off += (M + 1) * N) {
        _bidiag_decomp_horiz(M, N, _U3, _U_off, _B3, B_off, _V3, V_off);
      }

      return [_U3, _B3, _V3];
    }
  }(),
      _ref2 = (0, _slicedToArray2["default"])(_ref, 3),
      U = _ref2[0],
      B = _ref2[1],
      V = _ref2[2];

  return [new _nd_array.NDArray(U_shape, U), new _nd_array.NDArray(B_shape, B), new _nd_array.NDArray(V_shape, V)];
} // TODO: bidiag_decomp_full
// TODO: bidiag_solve
// TODO: bidiag_lstsq?