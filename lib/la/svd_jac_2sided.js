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
exports.svd_jac_2sided = svd_jac_2sided;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _dt = require("../dt");

var _nd_array = require("../nd_array");

var _matmul = require("./matmul");

var _qr = require("./qr");

var _transpose_inplace = require("./transpose_inplace");

var _svd_jac_utils = require("./_svd_jac_utils");

var _giv_rot = require("./_giv_rot");

function svd_jac_2sided(A) {
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
          _svd_jac_2sided = svd_jac_2sided(R),
          _svd_jac_2sided2 = (0, _slicedToArray2["default"])(_svd_jac_2sided, 3),
          _U = _svd_jac_2sided2[0],
          _sv = _svd_jac_2sided2[1],
          _V = _svd_jac_2sided2[2];

      return [(0, _matmul.matmul2)(Q, _U), _sv, _V];
    }

    if (N < M) {
      var _qr_decomp3 = (0, _qr.qr_decomp)(A.T),
          _qr_decomp4 = (0, _slicedToArray2["default"])(_qr_decomp3, 2),
          _Q = _qr_decomp4[0],
          _R = _qr_decomp4[1],
          _svd_jac_2sided3 = svd_jac_2sided(_R),
          _svd_jac_2sided4 = (0, _slicedToArray2["default"])(_svd_jac_2sided3, 3),
          _U2 = _svd_jac_2sided4[0],
          _sv2 = _svd_jac_2sided4[1],
          _V2 = _svd_jac_2sided4[2];

      (0, _transpose_inplace.transpose_inplace)(_V2);
      return [_V2, _sv2, (0, _matmul.matmul2)(_Q, _U2).T];
    }
  } // ALLOCATE RESULT DATA

  var DType = A.dtype === 'float32' ? 'float32' : 'float64',
      DTypeArray = _dt.ARRAY_TYPES[DType],
      TOL = Math.pow(N * (0, _dt.eps)(DType), 2),
      B = DType === 'float32' ? 64 / 4 : 64 / 8,
      // <- block size should match cache line size
  U = DTypeArray.from(A.data);
  A = undefined; // <- potentially allow GC

  var S = new DTypeArray(N * N),
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
        S[N * _i + j] = U[UV_off + N * _i + j];
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
      finished = true; // SWEEP

      for (var _Q2 = 0; _Q2 < N; _Q2 += B) {
        for (var P = 0; P <= _Q2; P += B) {
          for (var q = _Q2; q < _Q2 + B && q < N; q++) {
            for (var p = P; p < P + B && p < q; p++) {
              var S_pp = S[N * p + p],
                  S_pq = S[N * p + q],
                  S_qp = S[N * q + p],
                  S_qq = S[N * q + q]; // stopping criterion inspiredy by:
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
                  sβ = _svd_jac_angles3[3]; // ROTATE S


              (0, _giv_rot._giv_rot_rows)(S, N, N * p, N * q, cα, sα);
              (0, _giv_rot._giv_rot_cols)(S, N, p, q, cβ, sβ); // ENTRIES (k,l) AND (l,k) ARE REMAINDERS (CANCELLATION ERROR) FROM ELIMINATION => SHOULD BE SAFELY ZEROABLE

              S[N * p + q] = S[N * q + p] = 0.0; // ROTATE U & V

              (0, _giv_rot._giv_rot_rows)(U, N, UV_off + N * p, UV_off + N * q, cα, sα);
              (0, _giv_rot._giv_rot_rows)(V, N, UV_off + N * p, UV_off + N * q, cβ, -sβ);
            }
          }
        }
      }
    }

    (0, _svd_jac_utils._svd_jac_post)(N, U, S, V, UV_off, sv, sv_off, ord);
  }

  return [new _nd_array.NDArray(shape, U), new _nd_array.NDArray(shape.slice(0, -1), sv), new _nd_array.NDArray(shape, V)];
}