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
exports._svd_dc_bidiag = _svd_dc_bidiag;
exports._svd_dc = _svd_dc;
exports.svd_dc = svd_dc;
exports._svd_dc_neves = exports._svd_dc_2x3 = exports._svd_dc_1x2 = void 0;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _dt = require("../dt");

var _nd_array = require("../nd_array");

var _transpose_inplace = require("./transpose_inplace");

var _giv_rot = require("./_giv_rot");

var _svd_jac_utils = require("./_svd_jac_utils");

var _norm3 = require("./norm");

var _bidiag = require("./bidiag");

var _svd_dc_1x2 = function _svd_dc_1x2(N, U, U_off, F, B_off, V_off) {
  N |= 0;
  U_off |= 0;
  B_off |= 0;
  V_off |= 0;
  var c = F[B_off],
      s = F[B_off + 1];
  var norm = Math.hypot(c, s);

  if (0 !== norm) {
    c /= norm;
    s /= norm;
    F[B_off] = norm;
    F[B_off + 1] = NaN;
    F[V_off] = c;
    F[V_off + 1] = -s;
    F[V_off + N * 1] = s;
    F[V_off + N * 1 + 1] = c;
  } else {
    F[V_off] = 1;
    F[V_off + N * 1 + 1] = 1;
  }

  U[U_off] = 1;
};
/* Computes the SVD of the 2x3 matrix B
 *     ┌            ┐
 *     │ b1  b2   0 │
 * B = │            │
 *     │ 0   b3  b4 │
 *     └            ┘
 */


exports._svd_dc_1x2 = _svd_dc_1x2;

var _svd_dc_2x3 = function _svd_dc_2x3(N, U, U_off, F, B_off, V_off) {
  U_off |= 0;
  B_off |= 0;
  V_off |= 0;
  N |= 0;
  var M = N - 1 | 0;
  var b1 = F[B_off],
      b2 = F[B_off + 1],
      b3 = F[B_off + 2],
      b4 = F[B_off + 3]; // STEP 1:
  //   RQ-decompose down to 2x2 problem

  if (0 !== b4) {
    // eliminate B[1,2] (stored in b4)
    var norm = Math.hypot(b3, b4),
        _ca = b3 / norm,
        _sa = b4 / norm;

    F[V_off + N * 1 + 1] = _ca;
    F[V_off + N * 2 + 1] = _sa;
    b3 = norm;
    b4 = -_sa * b2;
    b2 = _ca * b2;

    if (0 !== b4) {
      // eliminate B[0,2] (stored in b4)
      var _norm = Math.hypot(b1, b4),
          _cb = b1 / _norm,
          _sb = b4 / _norm;

      b1 = _norm;
      F[V_off] = _cb;
      F[V_off + N * 1] = _sb * -_sa;
      F[V_off + N * 2] = _sb * _ca;
      F[V_off + 2] = -_sb;
      F[V_off + N * 1 + 2] = _cb * -_sa;
      F[V_off + N * 2 + 2] = _cb * _ca;
    } else {
      F[V_off] = 1;
      F[V_off + N * 1 + 2] = -_sa;
      F[V_off + N * 2 + 2] = _ca;
    }
  } else {
    F[V_off] = 1;
    F[V_off + N * 1 + 1] = 1;
    F[V_off + N * 2 + 2] = 1;
  } // STEP 2;
  //   USE SINGLE JACOBI ROTATION TO COMPUTE SVD


  var _svd_jac_angles2 = (0, _svd_jac_utils._svd_jac_angles)(b1, b2, 0, b3),
      _svd_jac_angles3 = (0, _slicedToArray2["default"])(_svd_jac_angles2, 4),
      ca = _svd_jac_angles3[0],
      sa = _svd_jac_angles3[1],
      cb = _svd_jac_angles3[2],
      sb = _svd_jac_angles3[3]; // WRITE SINGULAR VALUES TO B  


  F[B_off + 1] = NaN; // <- TODO remove after testing

  F[B_off + 3] = NaN; // <- TODO remove after testing

  F[B_off] = ca * b1 * cb - (ca * b2 + sa * b3) * sb;
  b3 = (-sa * b2 + ca * b3) * cb - sa * b1 * sb;
  var s = b3 < 0 ? -1 : +1;
  F[B_off + 2] = s * b3; // WRITE U

  U[U_off] = ca;
  U[U_off + 1] = -sa * s;
  U[U_off + M * 1] = sa;
  U[U_off + M * 1 + 1] = ca * s; // ROTATE V

  F[V_off + 1] = sb * F[V_off];
  F[V_off] *= cb;
  var V1 = F[V_off + N * 1] * cb - F[V_off + N * 1 + 1] * sb;
  F[V_off + N * 1 + 1] = F[V_off + N * 1] * sb + F[V_off + N * 1 + 1] * cb;
  F[V_off + N * 1] = V1;
  var V2 = F[V_off + N * 2] * cb - F[V_off + N * 2 + 1] * sb;
  F[V_off + N * 2 + 1] = F[V_off + N * 2] * sb + F[V_off + N * 2 + 1] * cb;
  F[V_off + N * 2] = V2;
};
/* Computes the singular values of a matrix shaped like a backwards seven.
 *     ┌                        ┐
 *     │ d1                     │
 *     │                        │
 *     │     d2                 │
 *     │       .                │
 * B = │         .              │
 *     │           .            │
 *     │            d[n-1]      │
 *     │                        │
 *     │ z1  . . .  z[n-1] z[n] │
 *     └                        ┘
 * See:
 *   - "INTRODUCTION OF DOUBLE DIVIDE AND CONQUER AND THE RECENT PROGRESS", TARO KONDA & YOSHIMASA NAKAMURA
 *
 * Keep in mind that B is stored in a sparse memory layout.
 * 
 * B = Array(d1, z1, d2, z2, ..., d[n]=0, z[n])
 * 
 * where d[i] >= d[i+1]
 * 
 * Keep in mind that V is stored in a column-major order in this method.
 */


exports._svd_dc_2x3 = _svd_dc_2x3;

var _svd_dc_neves = function _svd_dc_neves(N, n, U, U_off, F, B_off, V_off, I) {
  U_off |= 0;
  B_off |= 0;
  V_off |= 0;
  N |= 0;
  var M = N - 1 | 0;
  n |= 0;
  var m = n - 1 | 0;
  if (I.length < M * 3) throw new Error('Assertion failed: Integer work matrix I too small.');
  if (F.length < M * (M + 2) + M * 2
  /*B*/
  + (M + 1) * (M + 1)
  /*V*/
  ) throw new Error('Assertion failed: Scalar work matrix F too small.'); // Amount of temp. float memory:
  //   - m*m entries to store the matrix to update U and V
  //   -   m entries to store the (shifted) singular values
  //   -   m entries to compute the matrix multiplication (row by row)

  var σ_off = M * (M + 2) - m,
      mm_off = M * (M + 2) - m * 2,
      W_off = M * (M + 2) - m * (m + 2); // <- mm as in "matrix multiplication"
  // Amount of temp. int memory:
  //   - m entries for the outer (row   ) order of U & V update matrices
  //   - m entries for the inner (column) order of U & V update matrices
  //   - m entries for the rotation pairings from step 2

  var rot_off = m * 2,
      inn_off = m,
      // <- inner order
  out_off = 0; // <- outer order

  if (n < 2) throw new Error('Assertion failed.'); //     DIAGONAL ELEMENTS: d[i] = F[B_off + 2*i  ]
  // OFF-DIAGONAL ELEMENTS: z[i] = F[B_off + 2*i+1]
  // TODO remove checks after testing
  // ensure d[-1] === 1

  if (F[B_off + 2 * m - 2] !== 0) throw new Error('Assertion failed.'); // ensure d[i] >= d[i+1]

  for (var i = 1; i < m; i++) {
    if (F[B_off + 2 * (i - 1)] < F[B_off + 2 * i]) throw new Error('Assertion failed.');
  }

  var NORM = new _norm3.FrobeniusNorm();

  var _ref = function () {
    for (var _i = 0; _i < m; _i++) {
      NORM.include(F[B_off + 2 * _i + 1]);
    }

    var zNorm = NORM.result;

    for (var _i2 = 0; _i2 < m; _i2++) {
      NORM.include(F[B_off + 2 * _i2]);
    }

    var scale = NORM.result;
    if (scale === 0) scale = 1;
    return [zNorm / scale, scale];
  }(),
      _ref2 = (0, _slicedToArray2["default"])(_ref, 2),
      zNorm = _ref2[0],
      scale = _ref2[1]; // normalize


  for (var _i3 = 0; _i3 < 2 * m; _i3++) {
    F[B_off + _i3] /= scale;
  }

  var TOL = (0, _dt.eps)(U instanceof Float32Array ? 'float32' : 'float64'); // STEP 1: DEFLATION
  //   TEMPORARILY MOVE VALUES WITH z[j] ≈ 0 TO THE LEFT SUCH THAT
  //
  //       ┌                                   ┐
  //       │ d1                                │
  //       │    .                              │
  //       │      .                            │
  //       │        .                          │
  //       │        d[i]                       │
  //       │                                   │
  //   B = │            d[i+1]                 │
  //       │                 .                 │
  //       │                   .               │
  //       │                     .             │
  //       │                       d[n-1]      │
  //       │                                   │
  //       │ 0 . . . 0  z[i+1] ... z[n-1] z[n] │
  //       └                                   ┘
  //
  //   This way only the lower right quadrant needs to be solved
  //   iteratively (using the secular equations). The actual
  //   implementation moves the deflated values to σ and not to
  //   the left.

  var n0 = function () {
    var n0 = 0;

    for (var j = m - 1, _i4 = m - 1; _i4-- > 0;) {
      var di = F[B_off + 2 * _i4],
          zi = F[B_off + 2 * _i4 + 1],
          oi = I[out_off + _i4];

      if (Math.abs(zi) / TOL <= di) {
        // <- FIXME this could be estimated more accurately
        F[σ_off + n0] = di;
        I[inn_off + n0] = oi; // <- used as temp. for outerOrder

        ++n0;
      } else {
        --j;
        F[B_off + 2 * j] = di;
        F[B_off + 2 * j + 1] = F[B_off + 2 * _i4 + 1];
        I[inn_off + j] = oi; // <- used as temp. for outerOrder
      }
    }

    return n0;
  }(); // innerOrder just used as temp. memory, move to outerOrder


  for (var _i5 = 0; _i5 < n0; _i5++) {
    I[out_off + _i5] = I[inn_off + _i5];
  } // STEP 2:
  //   HANDLE DUPLICATE VALUES ON THE DIAGONAL
  //
  // Lets say d[i] == d[i+1] and z[i] != 0 with 0 < i < n-1
  //     ┌                                    ┐
  //     │ d1                                 │
  //     │   .                                │
  //     │     .                              │
  //     │       .                            │
  //     │        d[i]                        │
  //     │                                    │
  // B = │             d[i+1]                 │
  //     │                  .                 │
  //     │                    .               │
  //     │                      .             │
  //     │                        d[n-1]      │
  //     │                                    │
  //     │ z1 ... z[i] z[i+1] ... z[n-1] z[n] │
  //     └                                    ┘
  //
  // We can now use a left and right Givens rotation
  // to eliminate z[i+1].
  // n := √(z²ᵢ + z²ᵢ₊₁)
  // c := zᵢ   / n
  // s := zᵢ₊₁ / n
  //
  //     ┌                            ┐
  //     │ 1                          │
  //     │   .                        │
  //     │     .                      │
  //     │       .                    │
  //     │         1                  │
  //     │                            │
  //     │            c -s            │
  // G = │                            │
  //     │            s  c            │
  //     │                            │
  //     │                  1         │
  //     │                    .       │
  //     │                      .     │
  //     │                        .   │
  //     │                          1 │
  //     └                            ┘
  //
  //                           ┌                                   ┐
  //                           │ d1                                │
  //                           │   .                               │
  //                           │     .                             │
  //                           │       .                           │
  //                           │        d[i]                       │
  //                           │                                   │
  // B = GᵀG⋅B⋅GᵀG = Gᵀ⋅B'⋅G = │            d[i+1]                 │
  //                           │                 .                 │
  //                           │                   .               │
  //                           │                     .             │
  //                           │                       d[n-1]      │
  //                           │                                   │
  //                           │ z1 ...  0    n  . . . z[n-1] z[n] │
  //                           └                                   ┘
  //
  // Thus z'[i] is now 0 and can be moved to the left side (similar to Step 1)


  var n1 = function () {
    var n1 = n0;

    for (var j = m - 1, _i6 = m - 1; _i6-- > n0;) {
      var di = F[B_off + 2 * _i6],
          dj = F[B_off + 2 * j],
          oi = I[inn_off + _i6],
          zi = F[B_off + 2 * _i6 + 1],
          zj = F[B_off + 2 * j + 1],
          z = Math.hypot(zi, zj);

      if ((di - dj) / TOL <= di || !isFinite(Math.sqrt(m) / (di - dj))) // <- TODO find better criteria
        {
          var c = zj / z,
              s = zi / z;
          F[W_off + 2 * n1] = c;
          F[W_off + 2 * n1 + 1] = s;
          F[B_off + 2 * j + 1] = z;
          F[B_off + 2 * j] = F[σ_off + n1] = di;
          I[out_off + n1] = oi; // <- used as temp. for outerOrder

          I[rot_off + n1] = j;
          ++n1;
        } else {
        --j;
        F[B_off + 2 * j] = di;
        F[B_off + 2 * j + 1] = F[B_off + 2 * _i6 + 1];
        I[out_off + j] = oi;
      }
    }

    return n1;
  }();

  for (var _i7 = 2 * n0; _i7 < 2 * n1; _i7++) {
    F[B_off + _i7] = F[W_off + _i7];
    F[W_off + _i7] = 0;
  } // STEP 3:
  //   SOLVE THE SECULAR EQUATIONS
  //
  //      !            zₗ²
  //   0  =  1 + Σₗ ────────      k = 0, ..., m-1
  //                dₗ² - σₖ²
  //


  var _loop = function _loop(_i8) {
    var sLo = F[B_off + 2 * _i8],
        sHi = n1 < _i8 ? F[B_off + 2 * _i8 - 2] : sLo + zNorm;
    if (sLo > sHi) throw new Error('Assertion failed.');

    var shift = function () {
      if (_i8 > n1) {
        var mid = (sLo + sHi) / 2;
        var sum = 1;

        for (var _i32 = n1; _i32 < m; _i32++) {
          var _di = F[B_off + 2 * _i32],
              _zi = F[B_off + 2 * _i32 + 1];
          sum += _zi / (_di - mid) * (_zi / (_di + mid));
        }

        if (!isFinite(sum)) throw new Error('Assertion failed.');

        if (sum < 0) {
          var _s2 = sHi;
          sLo = sLo - sHi;
          sHi = -Number.MIN_VALUE;
          return _s2;
        }
      }

      var s = sLo;
      sHi = sHi - sLo;
      sLo = +Number.MIN_VALUE;
      return s;
    }(); // TODO: Make bisection more accurate (see nd.opt.root1d_bisect)


    for (;;) {
      // bisect
      var _s3 = (sLo + sHi) / 2;

      if (_s3 === sLo || _s3 === sHi) {
        F[σ_off + _i8] = _s3; // FIXME at this point sLo and sHi still have to be compared

        break;
      } // evalue the secular equation


      var sum = 1;

      for (var _i33 = n1; _i33 < m; _i33++) {
        var _di2 = F[B_off + 2 * _i33],
            _zi2 = F[B_off + 2 * _i33 + 1];
        sum += _zi2 / (_di2 - shift - _s3) * (_zi2 / (_di2 + shift + _s3));
      }

      if (!isFinite(sum)) throw new Error('Assertion failed.'); // adjust bisection bounds

      if (sum <= 0) sLo = _s3;
      if (sum >= 0) sHi = _s3;
    }
  };

  for (var _i8 = n1; _i8 < m; _i8++) {
    _loop(_i8);
  }

  if (Math.abs(F[B_off + 2 * m - 1]) === 0) {
    F[B_off + 2 * m - 2] = 0;
    F[B_off + 2 * m - 1] = 0;
    F[σ_off + m - 1] = 0;
  } // STEP 4:
  //   RECOMPUTE z TO ORTHOGONALIZE U & V
  //   (as originally suggested by Gu and Eisenstat)


  {
    var σn = F[σ_off + m - 1],
        sn = F[B_off + 2 * (m - 1 - (σn < 0))]; // <- shift

    for (var _i9 = n1; _i9 < m; _i9++) {
      var di = F[B_off + 2 * _i9];
      var zi = (sn - di + σn) * (sn + di + σn);

      for (var j = n1; j < _i9; j++) {
        var σj = F[σ_off + j],
            sj = F[B_off + 2 * (j - (σj < 0))],
            // <- shift
        dj = F[B_off + 2 * j];
        zi *= (sj - di + σj) / (dj - di) * ((sj + di + σj) / (dj + di));
      }

      for (var _j = _i9; _j < m - 1; _j++) {
        var _j2 = F[σ_off + _j],
            _sj = F[B_off + 2 * (_j - (_j2 < 0))],
            // <- shift
        _dj = F[B_off + 2 * _j + 2];
        zi *= (_sj - di + _j2) / (_dj - di) * ((_sj + di + _j2) / (_dj + di));
      }

      F[B_off + 2 * _i9 + 1] = Math.sign(F[B_off + 2 * _i9 + 1]) * Math.sqrt(zi);
    }
  } // triple-merge-sort singular values from deflation
  // and the ones from the secular equation solutions

  for (var h = n0 - 1, _i10 = n1 - 1, _j3 = n1, k = 0; k < m; k++) {
    var val = -Infinity,
        best = 3;

    if (_j3 < m) {
      var _j4 = F[σ_off + _j3];
      best = 2;
      val = _j4 + F[B_off + 2 * (_j3 - (_j4 < 0))];
    }

    if (_i10 >= n0) {
      var σi = F[σ_off + _i10];

      if (!(σi < val)) {
        best = 1;
        val = σi;
      }
    }

    if (h >= 0) {
      var σh = F[σ_off + h];

      if (!(σh < val)) {
        best = 0;
        val = σh;
      }
    }

    switch (best) {
      case 0:
        I[inn_off + h--] = k;
        continue;

      case 1:
        I[inn_off + _i10--] = k;
        continue;

      case 2:
        I[inn_off + _j3++] = k;
        continue;

      default:
        throw new Error('Assertion failed.');
    }
  } // STEP 5:
  //   UPDATE U


  for (var _i11 = n1; _i11 < m; _i11++) {
    var _i12 = F[σ_off + _i11],
        si = F[B_off + 2 * (_i11 - (_i12 < 0))]; // <- shift

    NORM.reset();

    for (var _j5 = n1; _j5 < m - 1; _j5++) {
      var _dj2 = F[B_off + 2 * _j5],
          zj = F[B_off + 2 * _j5 + 1],
          W_ij = zj / (_dj2 - si - _i12) * (_dj2 / (_dj2 + si + _i12));
      NORM.include(F[W_off + m * _i11 + _j5] = W_ij);
    }

    NORM.include(F[W_off + m * _i11 + m - 1] = -1);
    var norm = NORM.result;
    if (!(0 < norm)) throw new Error('Assertion failed.');

    for (var _j6 = n1; _j6 < m; _j6++) {
      F[W_off + m * _i11 + _j6] /= norm;
    }
  } // transpose dense part of W


  for (var _i13 = n1; _i13 < m; _i13++) {
    for (var _j7 = _i13; ++_j7 < m;) {
      var _W_ij = F[W_off + m * _i13 + _j7];
      F[W_off + m * _i13 + _j7] = F[W_off + m * _j7 + _i13];
      F[W_off + m * _j7 + _i13] = _W_ij;
    }
  }

  if (n0 < n1) {
    // init deflated region in W
    for (var _i14 = n0; _i14 < n1; _i14++) {
      F.fill(0.0, W_off + m * _i14 + n0, W_off + m * _i14 + m);
      F[W_off + m * _i14 + _i14] = 1;
    }

    for (var _i15 = n1; _i15 < m; _i15++) {
      F.fill(0.0, W_off + m * _i15 + n0, W_off + m * _i15 + n1);
    } // rotate W


    for (var _i16 = n1; _i16-- > n0;) {
      var _j8 = I[rot_off + _i16];

      if (_j8 < m - 1) {
        var c = F[B_off + 2 * _i16],
            s = F[B_off + 2 * _i16 + 1];
        (0, _giv_rot._giv_rot_rows)(F, m - _i16, W_off + m * _i16 + _i16, W_off + m * _j8 + _i16, c, s);
      }
    }
  } // UPDATE U: U = U⋅Wᵀ


  for (var r = 0; r < m; r++) {
    // compute dense part of row (matrix multiplication)
    // U is fairly sparse so the matrix multiplication
    // is designed to benefit from that fact
    F.fill(0.0, mm_off + n0, mm_off + m);

    for (var _i17 = n0; _i17 < m; _i17++) {
      var U_ri = U[U_off + M * r + I[out_off + _i17]];
      if (0 !== U_ri) for (var _j9 = n0; _j9 < m; _j9++) {
        F[mm_off + _j9] += U_ri * F[W_off + m * _i17 + _j9];
      }
    } // compute deflated part of row


    for (var _i18 = 0; _i18 < n0; _i18++) {
      var _c = I[out_off + _i18];
      F[mm_off + _i18] = U[U_off + M * r + _c];
    } // write back row


    for (var _i19 = 0; _i19 < m; _i19++) {
      var _c2 = I[inn_off + _i19];
      U[U_off + M * r + _c2] = F[mm_off + _i19];
    }
  } // STEP 6:
  //   UPDATE V


  for (var _i20 = n1; _i20 < m; _i20++) {
    var _i21 = F[σ_off + _i20],
        _si = F[B_off + 2 * (_i20 - (_i21 < 0))];
    NORM.reset();

    for (var _j10 = n1; _j10 < m; _j10++) {
      var _dj3 = F[B_off + 2 * _j10],
          _zj = F[B_off + 2 * _j10 + 1],
          _W_ij2 = _zj / (_dj3 - _si - _i21) / (_dj3 + _si + _i21);

      NORM.include(F[W_off + m * _i20 + _j10] = _W_ij2);
    }

    var _norm2 = NORM.result;
    if (!(0 < _norm2 || _i20 === m - 1)) throw new Error('Assertion failed.');

    for (var _j11 = n1; _j11 < m; _j11++) {
      F[W_off + m * _i20 + _j11] /= _norm2;
    }
  }

  if (0 === F[B_off + 2 * m - 1]) {
    for (var _i22 = n1; _i22 < m - 1; _i22++) {
      F[W_off + m * (m - 1) + _i22] = 0;
    }

    F[W_off + m * (m - 1) + (m - 1)] = 1;
  } // transpose dense part of W


  for (var _i23 = n1; _i23 < m; _i23++) {
    for (var _j12 = _i23; ++_j12 < m;) {
      var _W_ij3 = F[W_off + m * _i23 + _j12];
      F[W_off + m * _i23 + _j12] = F[W_off + m * _j12 + _i23];
      F[W_off + m * _j12 + _i23] = _W_ij3;
    }
  }

  if (n0 < n1) {
    // init deflated region in W
    for (var _i24 = n0; _i24 < n1; _i24++) {
      F.fill(0.0, W_off + m * _i24 + n0, W_off + m * _i24 + m);
      F[W_off + m * _i24 + _i24] = 1;
    }

    for (var _i25 = n1; _i25 < m; _i25++) {
      F.fill(0.0, W_off + m * _i25 + n0, W_off + m * _i25 + n1);
    } // rotate W


    for (var _i26 = n1; _i26-- > n0;) {
      var _j13 = I[rot_off + _i26],
          _c3 = F[B_off + 2 * _i26],
          _s = F[B_off + 2 * _i26 + 1];
      (0, _giv_rot._giv_rot_rows)(F, m - _i26, W_off + m * _i26 + _i26, W_off + m * _j13 + _i26, _c3, _s);
    }
  } // UPDATE V: V = V⋅Wᵀ


  for (var _r = 0; _r < n; _r++) {
    // compute dense part of row (matrix multiplication)
    // U is fairly sparse so the matrix multiplication
    // is designed to benefit from that fact
    F.fill(0.0, mm_off + n0, mm_off + m);

    for (var _i27 = n0; _i27 < m; _i27++) {
      var V_ri = F[V_off + N * _r + I[out_off + _i27]];
      if (0 !== V_ri) for (var _j14 = n0; _j14 < m; _j14++) {
        F[mm_off + _j14] += V_ri * F[W_off + m * _i27 + _j14];
      }
    } // compute deflated part of row


    for (var _i28 = 0; _i28 < n0; _i28++) {
      var _c4 = I[out_off + _i28];
      F[mm_off + _i28] = F[V_off + N * _r + _c4];
    } // write back row


    for (var _i29 = 0; _i29 < m; _i29++) {
      var _c5 = I[inn_off + _i29];
      F[V_off + N * _r + _c5] = F[mm_off + _i29];
    }
  } // STEP 7:
  //   GENERAL POSTPROCESSING


  for (var _i30 = n1; _i30 < m; _i30++) {
    F[σ_off + _i30] += F[B_off + 2 * (_i30 - (F[σ_off + _i30] < 0))];
  }

  for (var _i31 = 0; _i31 < m; _i31++) {
    var _j15 = I[inn_off + _i31];
    F[B_off + 2 * _j15] = F[σ_off + _i31] * scale;
    F[B_off + 2 * _j15 + 1] = NaN; // <- TEST ONLY
  }
};
/* Computes the svd of an upper bidiagonal matrix inplace.
 * The bidiagonal matrix is stored in a memory efficient
 * banded memory format where B[B_off+2*i] is the diagonal
 * entry B(i,i) and B[B_off+2*i+1] is the off-diagonal
 * element B(i,i+1).
 */


exports._svd_dc_neves = _svd_dc_neves;

function _svd_dc_bidiag(N, n, U, U_off, F, B_off, V_off, I) {
  // TODO:
  //   Replacing this recursion with a loop would likely
  //   increase performance by a significant amount
  // TODO:
  //   For small enough bidiagonal matrices (around m <= 64),
  //   the QR method could be used to compute the SVD not just
  //   fast but also more accurately, see:
  //     * "Accurate Singular Values of Bidiagonal Matrices"
  //       by James Demmel & W. Kahan
  //       Lapack Working Note (LAWN) No. 3
  //       http://www.netlib.org/lapack/lawnspdf/lawn03.pdf
  //     * "Computing Small Singular Values of Bidiagonal Matrices with Guaranteed High Relative Accuracy"
  //       by James Demmel & W. Kahan
  if (n > N) throw new Error('Assertion failed.');
  if (1 >= n) throw new Error('Assertion failed.');
  if (2 === n) return _svd_dc_1x2(N, U, U_off, F, B_off, V_off);
  if (3 === n) return _svd_dc_2x3(N, U, U_off, F, B_off, V_off);
  var M = N - 1,
      m = n - 1,
      n0 = n >>> 1,
      m0 = n0 - 1;
  if (I.length < M) throw new Error('Assertion failed: Integer work matrix I too small.');
  if (F.length < M * 2) throw new Error('Assertion failed: Scalar work matrix F too small.'); // STEP 1: DIVIDE
  // --------------
  // The idea is to divide the bidiagonal matrix into an upper left block B1 and
  // a lower right block B2 with one row of the bidiagonal with the non-zero entries
  // b1 and b2 separating the two blocks.
  //       ┏                     ┓     ┏                     ┓
  //       ┃          ┆          ┃     ┃          ┆          ┃
  //       ┃    B1    ┆          ┃     ┃ U1⋅S1⋅V1 ┆          ┃
  //       ┃          ┆          ┃     ┃          ┆          ┃
  //       ┃┄┄┄┄┄┄┄┬┄┄┼┄┄┐       ┃(SVD)┃┄┄┄┄┄┄┄┬┄┄┼┄┄┐       ┃
  //   B = ┃       ┆b1┆b2┆       ┃  =  ┃       ┆b1┆b2┆       ┃
  //       ┃       └┄┄┼┄┄┴┄┄┄┄┄┄┄┃     ┃       └┄┄┼┄┄┴┄┄┄┄┄┄┄┃
  //       ┃          ┆          ┃     ┃          ┆          ┃
  //       ┃          ┆    B2    ┃     ┃          ┆ U2⋅S2⋅V2 ┃
  //       ┃          ┆          ┃     ┃          ┆          ┃
  //       ┗                     ┛     ┗                     ┛

  _svd_dc_bidiag(N, n0, U, U_off, F, B_off, V_off, I);

  _svd_dc_bidiag(N, n - n0, U, U_off + M * n0 + n0, F, B_off + 2 * n0, V_off + N * n0 + n0, I);

  U[U_off + M * m0 + m0] = 1;
  var b1 = F[B_off + 2 * m0],
      b2 = F[B_off + 2 * m0 + 1]; //  if( ! isFinite(b1) ) throw new Error('Assertion failed.');
  //  if( ! isFinite(b2) ) throw new Error('Assertion failed.');
  // The orthogonal matrix U1,U2,V1,V2 can moved to separate left and right orthogonal
  // matrices. V1ᵀ[-1] is the last row of the transpose of V1. V2ᵀ[0] is the first row of
  // the transpose of V2.                                                   ┏                     ┓
  //     ┏                     ┓   ┏               ┓ ┏                    ┓ ┃          ┆          ┃
  //     ┃          ┆          ┃   ┃      ┆        ┃ ┃        ┆ ┆         ┃ ┃          ┆          ┃
  //     ┃ U1⋅S1⋅V1 ┆          ┃   ┃  U1  ┆        ┃ ┃   S1'  ┆0┆         ┃ ┃    V1    ┆          ┃
  //     ┃          ┆          ┃   ┃      ┆        ┃ ┃        ┆ ┆         ┃ ┃          ┆          ┃
  //     ┃┄┄┄┄┄┄┄┬┄┄┼┄┄┐       ┃   ┃┄┄┄┄┄┄┼┄┐      ┃ ┃┄┄┄┄┄┄┄┄┴┄┼┄┄┄┄┄┄┄┄┄┃ ┃          ┆          ┃
  // B = ┃       ┆b1┆b2┆       ┃ = ┃      ┆1┆      ┃ ┃b1⋅V1ᵀ[-1]┆b2⋅V2ᵀ[0]┃ ┃┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┄┄┄┄┄┃
  //     ┃       └┄┄┼┄┄┴┄┄┄┄┄┄┄┃   ┃      └┄┼┄┄┄┄┄┄┃ ┃┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┄┄┬┄┃ ┃          ┆          ┃
  //     ┃          ┆          ┃   ┃        ┆      ┃ ┃          ┆       ┆ ┃ ┃          ┆          ┃
  //     ┃          ┆ U2⋅S2⋅V2 ┃   ┃        ┆  U2  ┃ ┃          ┆   S2' ┆0┃ ┃          ┆    V2    ┃
  //     ┃          ┆          ┃   ┃        ┆      ┃ ┃          ┆       ┆ ┃ ┃          ┆          ┃
  //     ┗                     ┛   ┗               ┛ ┗                    ┛ ┃          ┆          ┃
  //                                                                        ┗                     ┛

  for (var i = 0; i < m0; i++) {
    F[B_off + 2 * i + 1] = b1 * F[V_off + N * m0 + i];
  }

  for (var _i34 = n0; _i34 < m; _i34++) {
    F[B_off + 2 * _i34 + 1] = b2 * F[V_off + N * n0 + _i34];
  } // With the following definitions:
  //                 ┏          ╷    ┓
  //   b1⋅V1ᵀ[-1] =: ┃    R1    ┆ r1 ┃
  //                 ┗          ╵    ┛
  //                 ┏          ╷    ┓
  //   b2⋅V2ᵀ[ 0] =: ┃    R2    ┆ r2 ┃
  //                 ┗          ╵    ┛
  //   h := √(r1² + r2²)
  //   c = r1 / h
  //   s = r2 / h
  //         ┏             ┓
  //         ┃             ┃
  //         ┃      W1     ┃
  //   V1 =: ┃             ┃
  //         ┃┄┄┄┄┄┄┄┄┄┄┄┄┄┃
  //         ┃      w1     ┃
  //         ┗             ┛
  //         ┏             ┓
  //         ┃             ┃
  //         ┃      W2     ┃
  //   V2 =: ┃             ┃
  //         ┃┄┄┄┄┄┄┄┄┄┄┄┄┄┃
  //         ┃      w2     ┃
  //         ┗             ┛
  // And a single Givens rotation, we can make the right-most column of the middle matrix zero.
  //                                                  ┏                     ┓                                            ┏                     ┓
  //       ┏               ┓ ┏                      ┓ ┃          ┆          ┃   ┏               ┓ ┏                    ┓ ┃          ┆          ┃
  //       ┃      ┆        ┃ ┃        ┆  ┆          ┃ ┃          ┆          ┃   ┃      ┆        ┃ ┃        ┆ ┆         ┃ ┃     W1   ┆          ┃
  //       ┃  U1  ┆        ┃ ┃   S1'  ┆ 0┆          ┃ ┃    V1    ┆          ┃   ┃  U1  ┆        ┃ ┃   S1'  ┆0┆         ┃ ┃          ┆          ┃
  //       ┃      ┆        ┃ ┃        ┆  ┆          ┃ ┃          ┆          ┃   ┃      ┆        ┃ ┃        ┆ ┆         ┃ ┃┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┄┄┄┄┄┃
  //       ┃┄┄┄┄┄┄┼┄┐      ┃ ┃┄┄┄┄┄┄┄┄┼┄┄┼┄┄┄┄┄┄┄┬┄┄┃ ┃          ┆          ┃   ┃┄┄┄┄┄┄┼┄┐      ┃ ┃┄┄┄┄┄┄┄┄┼┄┼┄┄┄┄┄┄┄┬┄┃ ┃   c*w1   ┆  s*w2    ┃
  //   B = ┃      ┆1┆      ┃ ┃   R1   ┆r1┆  R2   ┆r2┃ ┃┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┄┄┄┄┄┃ = ┃      ┆1┆      ┃ ┃   R1   ┆h┆  R2   ┆ ┃ ┃┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┄┄┄┄┄┃
  //       ┃      └┄┼┄┄┄┄┄┄┃ ┃┄┄┄┄┄┄┄┄┴┄┄┼┄┄┄┄┄┄┄┼┄┄┃ ┃          ┆          ┃   ┃      └┄┼┄┄┄┄┄┄┃ ┃┄┄┄┄┄┄┄┄┴┄┼┄┄┄┄┄┄┄┤ ┃ ┃          ┆          ┃
  //       ┃        ┆      ┃ ┃           ┆       ┆  ┃ ┃          ┆          ┃   ┃        ┆      ┃ ┃          ┆       ┆0┃ ┃          ┆    W2    ┃
  //       ┃        ┆  U2  ┃ ┃           ┆   S2' ┆ 0┃ ┃          ┆    V2    ┃   ┃        ┆  U2  ┃ ┃          ┆   S2' ┆ ┃ ┃          ┆          ┃
  //       ┃        ┆      ┃ ┃           ┆       ┆  ┃ ┃          ┆          ┃   ┃        ┆      ┃ ┃          ┆       ┆ ┃ ┃┄┄┄┄┄┄┄┄┄┄┼┄┄┄┄┄┄┄┄┄┄┃
  //       ┗               ┛ ┗                      ┛ ┃          ┆          ┃   ┗               ┛ ┗                    ┛ ┃  -s*w1   ┆  c*w2    ┃
  //                                                  ┗                     ┛                                            ┗                     ┛


  var _ref3 = function () {
    var c = b1 * F[V_off + N * m0 + m0],
        s = b2 * F[V_off + N * n0 + m];
    var h = Math.hypot(c, s);
    return [c / h, s / h, h];
  }(),
      _ref4 = (0, _slicedToArray2["default"])(_ref3, 3),
      c = _ref4[0],
      s = _ref4[1],
      h = _ref4[2]; //  if( ! isFinite(h) ) throw new Error('Assertion failed.'); 


  F[B_off + 2 * m0] = 0;
  F[B_off + 2 * m0 + 1] = h;

  if (0 !== h) {
    for (var _i35 = 0; _i35 < n0; _i35++) {
      F[V_off + N * _i35 + m] = F[V_off + N * _i35 + m0] * -s;
      F[V_off + N * _i35 + m0] *= c;
    }

    for (var _i36 = n0; _i36 < n; _i36++) {
      F[V_off + N * _i36 + m0] = F[V_off + N * _i36 + m] * s;
      F[V_off + N * _i36 + m] *= c;
    }
  } // STEP 2: CONQUER
  // ---------------
  // After step 1, row and column swaps can be used to turn the middle matrix into a Neves matrix.
  // See _svd_dc_neves() for more information about the shape of the Neves matrix. An integer array
  // called the "outer order" is used to keep track of the row and column swaps. The outer order
  // indicates where the rows and column of the middle matrix have originaally been.


  I[m - 1] = m0; // merge sort diagonal entries

  for (var _i37 = 0, j = n0, k = 0; k < m - 1; k++) {
    I[k] = j >= m || _i37 < m0 && F[B_off + 2 * _i37] >= F[B_off + 2 * j] ? _i37++ : j++;
  } //  if( I.slice(0,m).sort().some((i,j) => i !== j) ) throw new Error('Assertion failed.');


  for (var _i38 = 0; _i38 < m; _i38++) {
    var _j16 = I[_i38];
    F[2 * _i38] = F[B_off + 2 * _j16];
    F[2 * _i38 + 1] = F[B_off + 2 * _j16 + 1];
  }

  for (var _i39 = 0; _i39 < m; _i39++) {
    F[B_off + 2 * _i39] = F[2 * _i39];
    F[B_off + 2 * _i39 + 1] = F[2 * _i39 + 1];
  }

  _svd_dc_neves(N, n, U, U_off, F, B_off, V_off, I);
}

function _svd_dc(M, N, U, U_off, sv, sv_off, V, V_off, I, F) {
  if (M > N) throw new Error('Assertion failed.');
  if (I.length < M * 3) throw new Error('Assertion failed: Integer work matrix I too small.');
  if (F.length < M * (M + 2) + M * 2 + (M + 1) * (M + 1) + (M + 1) * N) throw new Error('Assertion failed: Scalar work matrix F too small.');
  var B_off = M * (M + 2),
      V1_off = B_off + M * 2,
      V2_off = V1_off + (M + 1) * (M + 1);
  F.fill(0.0, V1_off, V2_off + N); // COPY A FROM V TO V2

  for (var i = 0; i < M; i++) {
    for (var j = 0; j < N; j++) {
      F[V2_off + N + N * i + j] = V[V_off + N * i + j];
    }
  }

  V.fill(0.0, V_off, V_off + M * N);
  (0, _bidiag._bidiag_decomp_horiz)(M, N, U, U_off, F, 0, F, V2_off);

  for (var _i40 = 0; _i40 < M; _i40++) {
    F[B_off + 2 * _i40] = F[(M + 1) * _i40 + _i40];
    F[B_off + 2 * _i40 + 1] = F[(M + 1) * _i40 + _i40 + 1];
  }

  _svd_dc_bidiag(M + 1, M + 1, V, V_off, F, B_off, V1_off, I);

  for (var _i41 = 0; _i41 < M; _i41++) {
    sv[sv_off + _i41] = F[B_off + 2 * _i41];
  } // update U = U @ U2 (U2 is stored in V)


  for (var _i42 = 0; _i42 < M; _i42++) {
    F.fill(0.0, 0, M);

    for (var k = 0; k < M; k++) {
      for (var _j17 = 0; _j17 < M; _j17++) {
        F[_j17] += U[U_off + M * _i42 + k] * V[V_off + M * k + _j17];
      }
    }

    for (var _j18 = 0; _j18 < M; _j18++) {
      U[U_off + M * _i42 + _j18] = F[_j18];
    }
  } // update V = V1 @ V2 (keep in mind that VV was computed in a transposed/column-major fashion)


  V.fill(0.0, V_off, V_off + M * M);

  for (var _k = 0; _k < M + 1; _k++) {
    for (var _i43 = 0; _i43 < M; _i43++) {
      for (var _j19 = 0; _j19 < N; _j19++) {
        V[V_off + N * _i43 + _j19] += F[V1_off + (M + 1) * _k + _i43] * F[V2_off + N * _k + _j19];
      }
    }
  }
}

function svd_dc(A) {
  // SEE:
  //  - "INTRODUCTION OF DOUBLE DIVIDE AND CONQUER AND THE RECENT PROGRESS", TARO KONDA & YOSHIMASA NAKAMURA
  //  - "A Divide-and-Conquer Approach for Solving Singular Value Decomposition on a Heterogeneous System", Ding Liu & Ruixuan Li & David J. Lilja & Weijun Xiao
  A = (0, _nd_array.asarray)(A);
  if (A.dtype.startsWith('complex')) throw new Error('svd_dc(A): A.dtype must be float.');

  var _A$shape$slice = A.shape.slice(-2),
      _A$shape$slice2 = (0, _slicedToArray2["default"])(_A$shape$slice, 2),
      M = _A$shape$slice2[0],
      N = _A$shape$slice2[1];

  if (M > N) {
    var _svd_dc2 = svd_dc(A.T),
        _svd_dc3 = (0, _slicedToArray2["default"])(_svd_dc2, 3),
        _U = _svd_dc3[0],
        _sv = _svd_dc3[1],
        _V = _svd_dc3[2];

    (0, _transpose_inplace.transpose_inplace)(_U);
    return [_V.T, _sv, _U];
  }

  var V_shape = A.shape,
      U_shape = V_shape.slice(),
      sv_shape = V_shape.slice(0, -1);
  U_shape[U_shape.length - 1] = M;
  var DType = A.dtype === 'float32' ? 'float32' : 'float64',
      DTypeArray = _dt.ARRAY_TYPES[DType];
  A = A.data;
  var len = A.length / (M * N),
      V = DTypeArray.from(A);
  A = undefined;
  var U = new DTypeArray(len * M * M),
      sv = new DTypeArray(len * M),
      F = new DTypeArray(M * (M + 2)
  /*F*/
  + M * 2
  /*B*/
  + (M + 1) * (M + 1)
  /*V1*/
  + (M + 1) * Math.max(1 + M, N)
  /*V2*/
  ),
      I = new Int32Array(M * 3);

  for (var U_off = 0, sv_off = 0, V_off = 0; sv_off < sv.length; U_off += M * M, sv_off += M, V_off += M * N) {
    _svd_dc(M, N, U, U_off, sv, sv_off, V, V_off, I, F);
  }

  return [new _nd_array.NDArray(U_shape, U), new _nd_array.NDArray(sv_shape, sv), new _nd_array.NDArray(V_shape, V)];
}