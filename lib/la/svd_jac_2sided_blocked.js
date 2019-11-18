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
exports.svd_jac_2sided_blocked = svd_jac_2sided_blocked;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _nd_array = require("../nd_array");

var _dt = require("../dt");

var _matmul = require("./matmul");

var _qr = require("./qr");

var _transpose_inplace = require("./transpose_inplace");

var _svd_jac_utils = require("./_svd_jac_utils");

var _giv_rot = require("./_giv_rot");

var B = 4,
    // <- block size should match cache line size
BB = B * B; // The memory of `S` is organized in a tiled order with a tile size
// of [B,B]. If for example BR=2 and BC=3, a matrix of shape [4,6]
// would have the following memory order.
// ┌                       ┐
// │  0  1    4  5    8  9 │
// │  2  3    6  7   10 11 │
// │                       │
// │ 12 13   16 17   20 21 │
// │ 14 15   18 19   22 23 │
// └                       ┘
// The purpose of this memory order is to make column (Givens) rotations
// faster.
//
// Let's say the cache line size is L. If we had a fortran-style,
// column major memory order that would correspond to a tile size of
// [L,1]. Say we have an [L,L] block and we want to access one row and
// one column of it. Accessing the column requires L cache lines to be
// read/written. Accessing the row only takes one cache line for a
// total of (L+1) cache lines.
//
// If we, on the other hand, used a tile size of [√(L),√(L)], both
// row and column access would require √(L) cache lines to be
// read/written each for a total of only 2*√(L) cache lines.
//
// ('B' as in Block)
//
// A power of two square block size is chosen such that it barely fits
// at least on cache line (64 bytes).
//
// Some quick and dirty tests have shown that is this method is faster
// for matrix sizes 800 and above.

function svd_jac_2sided_blocked(A) {
  A = (0, _nd_array.asarray)(A);
  if (A.dtype.startsWith('complex')) throw new Error('svd_jac_1sided(A): A.dtype must be float.');
  var shape = A.shape,
      N = shape[shape.length - 2] | 0; // QR DECOMPOSE RECTANGULAR MATRICES AND WORK ON (QUADRATIC) R

  {
    var M = shape[shape.length - 1]; // if A is not square use QR Decomposition

    if (N > M) {
      var _qr_decomp = (0, _qr.qr_decomp)(A),
          _qr_decomp2 = (0, _slicedToArray2["default"])(_qr_decomp, 2),
          Q = _qr_decomp2[0],
          R = _qr_decomp2[1],
          _svd_jac_2sided_block = svd_jac_2sided_blocked(R),
          _svd_jac_2sided_block2 = (0, _slicedToArray2["default"])(_svd_jac_2sided_block, 3),
          _U = _svd_jac_2sided_block2[0],
          _sv = _svd_jac_2sided_block2[1],
          _V = _svd_jac_2sided_block2[2];

      return [(0, _matmul.matmul2)(Q, _U), _sv, _V];
    }

    if (N < M) {
      var _qr_decomp3 = (0, _qr.qr_decomp)(A.T),
          _qr_decomp4 = (0, _slicedToArray2["default"])(_qr_decomp3, 2),
          _Q = _qr_decomp4[0],
          _R = _qr_decomp4[1],
          _svd_jac_2sided_block3 = svd_jac_2sided_blocked(_R),
          _svd_jac_2sided_block4 = (0, _slicedToArray2["default"])(_svd_jac_2sided_block3, 3),
          _U2 = _svd_jac_2sided_block4[0],
          _sv2 = _svd_jac_2sided_block4[1],
          _V2 = _svd_jac_2sided_block4[2];

      (0, _transpose_inplace.transpose_inplace)(_V2);
      return [_V2, _sv2, (0, _matmul.matmul2)(_Q, _U2).T];
    }
  } // ALLOCATE RESULT DATA

  var DType = A.dtype === 'float32' ? 'float32' : 'float64',
      DTypeArray = _dt.ARRAY_TYPES[DType],
      TOL = Math.pow(N * (0, _dt.eps)(DType), 2),
      n = Math.ceil(N / B) | 0,
      U = DTypeArray.from(A.data);
  A = undefined; // <- potentially allow GC

  var S = new DTypeArray(n * n * B * B),
      // <- tempory storage for decomposition
  V = new DTypeArray(U.length),
      sv = new DTypeArray(U.length / N),
      ord = Int32Array.from({
    length: N
  }, function (_, i) {
    return i;
  });
  if (1 > N) throw new Error('Assertion failed.');

  if (1 == N) {
    for (var i = U.length; i-- > 0;) {
      if (U[i] < +0.0) {
        U[i] *= -1.0;
        sv[i] = -1.0;
      } else sv[i] = +1.0;
    }

    return [new _nd_array.NDArray(shape, sv), new _nd_array.NDArray(shape.slice(0, -1), U), new _nd_array.NDArray(shape, V.fill(1))];
  }

  for (var UV_off = 0, sv_off = 0; sv_off < sv.length; UV_off += N * N, sv_off += N) {
    // MOVE FROM U TO S
    for (var _i = 0; _i < N; _i++) {
      for (var j = 0; j < N; j++) {
        var k = BB * (n * (_i / B | 0) + (j / B | 0)) + B * (_i % B) + j % B | 0;
        S[k] = U[UV_off + N * _i + j];
        U[UV_off + N * _i + j] = +(_i === j);
      }
    }

    ; // INIT V TO IDENTITY

    for (var _i2 = 0; _i2 < N; _i2++) {
      for (var _j = 0; _j < N; _j++) {
        V[UV_off + N * _i2 + _j] = +(_i2 === _j);
      }
    } //
    // JACOBI SVD ITERATIONS
    //


    for (var finished = false; !finished;) {
      finished = true;

      for (var _Q2 = n; _Q2-- > 0;) {
        for (var P = _Q2; P >= 0; P--) {
          for (var q = B; q-- > 0;) {
            for (var p = P === _Q2 ? q : B; p-- > 0;) {
              var S_pp = S[BB * (n * P + P) + B * p + p],
                  S_pq = S[BB * (n * P + _Q2) + B * p + q],
                  S_qp = S[BB * (n * _Q2 + P) + B * q + p],
                  S_qq = S[BB * (n * _Q2 + _Q2) + B * q + q]; // stopping criterion inspiredy by:
              //  "Jacobi's Method is More Accurate than QR"
              //   by James Demmel
              //   SIAM J. Matrix Anal. Appl, vol. 13, pp. 1204-1245, 1992

              if (!(S_pq * S_pq + S_qp * S_qp > Math.abs(S_pp * S_qq) * TOL)) continue;
              finished = false;

              var _svd_jac_angles2 = (0, _svd_jac_utils._svd_jac_angles)(S_pp, S_pq, S_qp, S_qq),
                  _svd_jac_angles3 = (0, _slicedToArray2["default"])(_svd_jac_angles2, 4),
                  cα = _svd_jac_angles3[0],
                  sα = _svd_jac_angles3[1],
                  cβ = _svd_jac_angles3[2],
                  sβ = _svd_jac_angles3[3]; // ROTATE ROWS IN S


              for (var _i3 = P * n * BB + p * B, _j2 = _Q2 * n * BB + q * B, I = _i3 + n * BB; _i3 < I; _i3 += BB, _j2 += BB) {
                for (var _k = 0; _k < B; _k++) {
                  var S_ki = cα * S[_k + _i3] + sα * S[_k + _j2];
                  S[_k + _j2] = -sα * S[_k + _i3] + cα * S[_k + _j2];
                  S[_k + _i3] = S_ki;
                }
              } // ROTATE COLUMNS IN S


              for (var _i4 = P * BB + p + n * n * B * B, _j3 = _Q2 * BB + q + n * n * B * B; (_i4 -= n * BB) >= 0;) {
                _j3 -= n * BB;

                for (var _k2 = BB; (_k2 -= B) >= 0;) {
                  var S_ik = S[_i4 + _k2] * cβ + S[_j3 + _k2] * -sβ;
                  S[_j3 + _k2] = S[_i4 + _k2] * sβ + S[_j3 + _k2] * cβ;
                  S[_i4 + _k2] = S_ik;
                }
              } // ENTRIES (k,l) AND (l,k) ARE REMAINDERS (CANCELLATION ERROR) FROM ELIMINATION => SHOULD BE SAFELY ZEROABLE


              S[BB * (n * P + _Q2) + B * p + q] = S[BB * (n * _Q2 + P) + B * q + p] = 0.0;
              var rowP = UV_off + N * (B * P + p),
                  rowQ = UV_off + N * (B * _Q2 + q); // ROTATE U & V

              (0, _giv_rot._giv_rot_rows)(U, N, rowP, rowQ, cα, sα);
              (0, _giv_rot._giv_rot_rows)(V, N, rowP, rowQ, cβ, -sβ);
            }
          }
        }
      }
    }

    for (var _i5 = 0; _i5 < N; _i5++) {
      var _k3 = BB * (n + 1) * (_i5 / B | 0) + (B + 1) * (_i5 % B) | 0;

      sv[sv_off + _i5] = S[_k3];
    }

    (0, _svd_jac_utils._svd_jac_post_skip1)(N, U, V, UV_off, sv, sv_off, ord);
  }

  return [new _nd_array.NDArray(shape, U), new _nd_array.NDArray(shape.slice(0, -1), sv), new _nd_array.NDArray(shape, V)];
}