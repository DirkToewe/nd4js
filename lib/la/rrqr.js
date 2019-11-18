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
exports._rrqr_rank = _rrqr_rank;
exports.rrqr_decomp_full = rrqr_decomp_full;
exports.rrqr_decomp = rrqr_decomp;
exports.rrqr_rank = rrqr_rank;
exports.rrqr_solve = rrqr_solve;
exports.rrqr_lstsq = rrqr_lstsq;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _singular_matrix_solve_error = require("./singular_matrix_solve_error");

var _nd_array = require("../nd_array");

var _dt = require("../dt");

var _permute = require("./permute");

var _norm4 = require("./norm");

var _giv_rot = require("./_giv_rot");

var _transpose_inplace2 = require("./transpose_inplace");

function _norm_update(norm, Q, i, j) {
  i |= 0;
  j <<= 1;

  for (; j < norm.length; j = j + 2 | 0) {
    var s = Math.abs(Q[i + (j >>> 1)]);

    if (s !== 0) {
      if (norm[j] < s) {
        var r = norm[j] / s;
        norm[j + 1] *= r * r;
        norm[j] = s;
      }

      s /= norm[j];
      norm[j + 1] += s * s;
    }
  }
}

function _norm(norm, i) {
  i <<= 1;
  var max = norm[i];
  return isFinite(max) ? Math.sqrt(norm[i + 1]) * max : max;
}

function _rrqr_rank(M, N, R, R_off, tmp) {
  // SEE: Gene H. Golub, Charles F. Van Golub
  //      "Matrix Computations", 4th edition
  //      Page 276f, Chap. 5.4.2 (QR with Column Pivoting) & 5.4.3 (Numerical Rank and AΠ=QR)
  var L = Math.min(M, N);
  if (!(tmp.length >= L)) throw new Error('Assertion failed.');
  var norm = new _norm4.FrobeniusNorm();

  for (var i = L; i-- > 0;) {
    for (var j = N; j-- > i;) {
      norm.include(R[R_off + N * i + j]);
    }

    tmp[i] = norm.result;
    if (!isFinite(tmp[i])) throw new Error('Infinity or NaN encountered during rank estimation.');
  }

  var dtype = tmp instanceof Float32Array ? 'float32' : 'float64',
      T = Math.sqrt((0, _dt.eps)(dtype)) * tmp[0]; // <- threshold

  var r = L;

  while (r > 0 && tmp[r - 1] <= T) {
    // <- TODO use binary search here for slightly improved performance
    --r;
  }

  return r;
}

function rrqr_decomp_full(A) {
  A = (0, _nd_array.asarray)(A);
  if (A.ndim < 2) throw new Error('A must be at least 2D.');

  var DTypeArray = _dt.ARRAY_TYPES[A.dtype === 'float32' ? 'float32' : 'float64'],
      R_shape = A.shape,
      Q_shape = Int32Array.from(R_shape),
      P_shape = Q_shape.slice(0, -1),
      _R_shape$slice = R_shape.slice(-2),
      _R_shape$slice2 = (0, _slicedToArray2["default"])(_R_shape$slice, 2),
      M = _R_shape$slice2[0],
      N = _R_shape$slice2[1];

  Q_shape[Q_shape.length - 1] = M;
  P_shape[P_shape.length - 1] = N;
  var R = DTypeArray.from(A.data); // <- we could encourage GC by setting `A = undefined` after this line

  A = undefined;
  var P = new Int32Array(R.length / M),
      // <- tracks column permutations
  Q = new DTypeArray(R.length / N * M),
      norm = new DTypeArray(N << 1); // <─ underflow-safe representation of the column norm (https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Math/hypot#Polyfill)

  for (var Q_off = 0, R_off = 0, P_off = 0; Q_off < Q.length; Q_off += M * M, R_off += M * N, P_off += N) {
    // INIT P
    for (var i = 0; i < N; i++) {
      P[P_off + i] = i;
    } // INIT Q (TO IDENTITY)


    for (var _i = 0; _i < M; _i++) {
      Q[Q_off + M * _i + _i] = 1;
    } // COMPUTE COLUMN NORM
    // (https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Math/hypot#Polyfill)


    for (var j = 0; j < M; j++) {
      _norm_update(norm, R, R_off + N * j, 0);
    } // ELIMINATE COLUMN BY COLUMN OF R (WHICH IS CURRENTLY STORED IN Q)


    for (var _i2 = 0; _i2 < N; _i2++) {
      // FIND PIVOT COLUMN THAT (HOPEFULLY) ENSURES RANK REVEAL (RRQR is inherently not guaranteed to do that)
      var p = -1,
          norm_p = -Infinity;

      for (var _j = _i2; _j < N; _j++) {
        var norm_j = _norm(norm, _j);

        if (norm_p < norm_j) {
          p = _j;
          norm_p = norm_j;
        }
      } // swap pivot to column i


      if (p !== _i2) {
        for (var _j2 = 0; _j2 < M; _j2++) {
          var ji = R_off + N * _j2 + _i2,
              jp = R_off + N * _j2 + p,
              R_ji = R[ji];
          R[ji] = R[jp];
          R[jp] = R_ji;
        }

        var P_i = P[P_off + _i2];
        P[P_off + _i2] = P[P_off + p];
        P[P_off + p] = P_i;
      } // RESET COLUMN NORM (INDEX i IS SET TO ZERO FOR THE NEXT RRQR)


      norm.fill(0.0, _i2 << 1); // ELIMINATE COLUMN BELOW DIAGONAL

      var ii = R_off + N * _i2 + _i2;

      for (var _j3 = _i2; ++_j3 < M;) {
        var _ji = R_off + N * _j3 + _i2;

        var _R_ji = R[_ji];

        if (0 !== _R_ji) {
          var R_ii = R[ii],
              _norm2 = Math.hypot(R_ii, _R_ji),
              c = R_ii / _norm2,
              s = _R_ji / _norm2;

          R[_ji] = 0;

          if (s !== 0) {
            R[ii] = _norm2;
            (0, _giv_rot._giv_rot_rows)(R, N - 1 - _i2, ii + 1, _ji + 1, c, s);
            (0, _giv_rot._giv_rot_rows)(Q, 1 + _j3, Q_off + M * _i2, Q_off + M * _j3, c, s);
          }
        }

        _norm_update(norm, R, R_off + N * _j3, _i2 + 1);
      }

      if (_i2 < M) {
        var _R_ii = R[R_off + N * _i2 + _i2];

        if (0 > _R_ii || Object.is(-0, _R_ii)) {
          for (var _j4 = _i2; _j4 < N; _j4++) {
            R[R_off + N * _i2 + _j4] *= -1;
          }

          for (var _j5 = 0; _j5 <= _i2; _j5++) {
            Q[Q_off + M * _i2 + _j5] *= -1;
          }
        }

        R[R_off + N * _i2 + _i2] = norm_p;
      }
    }

    (0, _transpose_inplace2._transpose_inplace)(M, Q, Q_off);
  }

  return [new _nd_array.NDArray(Q_shape, Q), new _nd_array.NDArray(R_shape, R), new _nd_array.NDArray(P_shape, P)];
}

function rrqr_decomp(A) {
  // TODO: implement Strong RRQR as well, e.g. as la.rrqr_decomp_strong 
  // SEE: Ming Gu, Stanley C. Eisenstat,
  //     "EFFICIENT ALGORITHMS FOR COMPUTING A STRONG RANK-REVEALING QR FACTORIZATION"
  //      https://math.berkeley.edu/~mgu/MA273/Strong_RRQR.pdf
  A = (0, _nd_array.asarray)(A);
  if (A.ndim < 2) throw new Error('A must be at least 2D.');

  var DTypeArray = _dt.ARRAY_TYPES[A.dtype === 'float32' ? 'float32' : 'float64'],
      Q_shape = A.shape,
      R_shape = Int32Array.from(Q_shape),
      P_shape = Q_shape.slice(0, -1),
      _Q_shape$slice = Q_shape.slice(-2),
      _Q_shape$slice2 = (0, _slicedToArray2["default"])(_Q_shape$slice, 2),
      N = _Q_shape$slice2[0],
      M = _Q_shape$slice2[1];

  R_shape[R_shape.length - 2] = M;
  P_shape[P_shape.length - 1] = M;
  if (N <= M) return rrqr_decomp_full(A);
  var Q = DTypeArray.from(A.data); // <- we could encourage GC by setting `A = undefined` after this line

  A = undefined;
  var P = new Int32Array(Q.length / N),
      // <- tracks column permutations
  R = new DTypeArray(Q.length / N * M),
      // <- cache cos() and sin() values to apply M column rotations to Q at once
  norm = new DTypeArray(M << 1); // <─ underflow-safe representation of the column norm (https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Math/hypot#Polyfill)

  for (var Q_off = 0, R_off = 0, P_off = 0; Q_off < Q.length; Q_off += N * M, R_off += M * M, P_off += M) {
    // INIT P
    for (var i = 0; i < M; i++) {
      P[P_off + i] = i;
    } // COMPUTE COLUMN NORM
    // (https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Math/hypot#Polyfill)


    norm.fill(0);

    for (var j = 0; j < N; j++) {
      _norm_update(norm, Q, Q_off + M * j, 0);
    } // ELIMINATE COLUMN BY COLUMN OF R (WHICH IS CURRENTLY STORED IN Q)


    for (var _i3 = 0; _i3 < M; _i3++) {
      // FIND PIVOT COLUMN THAT (HOPEFULLY) ENSURES RANK REVEAL (RRQR is inherently not guaranteed to do that)
      var p = -1,
          norm_p = -Infinity;

      for (var _j6 = _i3; _j6 < M; _j6++) {
        var norm_j = _norm(norm, _j6);

        if (norm_p < norm_j) {
          p = _j6;
          norm_p = norm_j;
        }
      } // swap pivot to column i


      if (p !== _i3) {
        for (var _j7 = 0; _j7 < N; _j7++) {
          var ji = Q_off + M * _j7 + _i3,
              jp = Q_off + M * _j7 + p,
              A_ji = Q[ji];
          Q[ji] = Q[jp];
          Q[jp] = A_ji;
        }

        var P_i = P[P_off + _i3];
        P[P_off + _i3] = P[P_off + p];
        P[P_off + p] = P_i;
      } // RESET COLUMN NORM (INDEX i IS SET TO ZERO FOR THE NEXT RRQR)


      norm.fill(0.0, _i3 + 1 << 1); // ELIMINATE COLUMN BELOW DIAGONAL

      var ii = Q_off + M * _i3 + _i3;

      for (var _j8 = _i3; ++_j8 < N;) {
        var _ji2 = Q_off + M * _j8 + _i3;

        var _A_ji = Q[_ji2];

        if (0 !== _A_ji) {
          var A_ii = Q[ii];

          var _norm3 = Math.hypot(A_ii, _A_ji),
              c = A_ii / _norm3,
              s = _A_ji / _norm3;

          if (s !== 0) {
            if (c < 0) {
              c *= -1;
              s *= -1;
              _norm3 *= -1;
            } // rotate i and j


            (0, _giv_rot._giv_rot_rows)(Q, M - 1 - _i3, ii + 1, _ji2 + 1, c, s);
            Q[ii] = _norm3;
          }

          Q[_ji2] = s;
        }

        _norm_update(norm, Q, Q_off + M * _j8, _i3 + 1);
      }

      Q[ii] = (Q[ii] < 0 || Object.is(Q[ii], -0) ? -1 : +1) * norm_p;
    } // MOVE R FROM Q -> R


    for (var _i4 = 0; _i4 < M; _i4++) {
      var R_ii = R[R_off + M * _i4 + _i4],
          _s = R_ii < 0 || Object.is(R_ii, -0) ? -1 : +1;

      for (var _j9 = _i4; _j9 < M; _j9++) {
        R[R_off + M * _i4 + _j9] = Q[Q_off + M * _i4 + _j9] * _s;
        Q[Q_off + M * _i4 + _j9] = _s * (_i4 === _j9);
      }
    } // COMPUTE Q


    for (var _i5 = M; _i5-- > 0;) {
      for (var _j10 = N; --_j10 > _i5;) {
        var _s2 = Q[Q_off + M * _j10 + _i5];
        if (0 === _s2) continue;
        Q[Q_off + M * _j10 + _i5] = 0;

        var _c = Math.sqrt((1 - _s2) * (1 + _s2));

        (0, _giv_rot._giv_rot_rows)(Q, M - _i5, Q_off + M * _j10 + _i5, Q_off + M * _i5 + _i5, _c, _s2);
      }
    }
  }

  return [new _nd_array.NDArray(Q_shape, Q), new _nd_array.NDArray(R_shape, R), new _nd_array.NDArray(P_shape, P)];
}

function rrqr_rank(R) {
  R = (0, _nd_array.asarray)(R);

  var _R$shape$slice = R.shape.slice(-2),
      _R$shape$slice2 = (0, _slicedToArray2["default"])(_R$shape$slice, 2),
      M = _R$shape$slice2[0],
      N = _R$shape$slice2[1],
      r_shape = R.shape.slice(0, -2);

  R = R.data;
  var r = new _dt.ARRAY_TYPES['int32'](R.length / M / N),
      tmp = new _dt.ARRAY_TYPES[R.dtype === 'float32' ? 'float32' : 'float64'](Math.min(M, N));

  for (var r_off = 0, R_off = 0; r_off < r.length; r_off++, R_off += M * N) {
    r[r_off] = _rrqr_rank(M, N, R, R_off, tmp);
  }

  return new _nd_array.NDArray(r_shape, r);
}

function rrqr_solve(Q, R, P, y) {
  Q = (0, _nd_array.asarray)(Q);
  R = (0, _nd_array.asarray)(R);
  var N = Q.shape[Q.ndim - 2];
  if (N !== R.shape[R.ndim - 1]) throw new Error('rrqr_solve(Q,R,P, y): Q @ R not square.');
  var x = rrqr_lstsq(Q, R, P, y),
      tmp = new _dt.ARRAY_TYPES[R.dtype === 'float32' ? 'float32' : 'float64'](N);

  for (var R_off = 0; R_off < R.data.length; R_off += N * N) {
    var rank = _rrqr_rank(N, N, R.data, R_off, tmp);

    if (rank < N) throw new _singular_matrix_solve_error.SingularMatrixSolveError(x);
  }

  return x;
}

function rrqr_lstsq(Q, R, P, y) {
  if (y == undefined) {
    var _Q, _Q2;

    if (P != undefined) throw new Error('rrqr_lstsq(Q,R,P, y): Either 2 ([Q,R,P], y) or 4 arguments (Q,R,P, y) expected.');
    y = R((_Q = Q, _Q2 = (0, _slicedToArray2["default"])(_Q, 3), Q = _Q2[0], R = _Q2[1], P = _Q2[2], _Q));
  }

  Q = (0, _nd_array.asarray)(Q);
  if (Q.ndim < 2) throw new Error('rrqr_lstsq(Q,R,P, y): Q.ndim must be at least 2.');
  R = (0, _nd_array.asarray)(R);
  if (R.ndim < 2) throw new Error('rrqr_lstsq(Q,R,P, y): R.ndim must be at least 2.');
  P = (0, _nd_array.asarray)(P);
  if (P.ndim < 1) throw new Error('rrqr_lstsq(Q,R,P, y): P.ndim must be at least 1.');
  y = (0, _nd_array.asarray)(y);
  if (y.ndim < 2) throw new Error('rrqr_lstsq(Q,R,P, y): y.ndim must be at least 2.');
  if (P.dtype !== 'int32') throw new Error('rrqr_lstsq(Q,R,P, y): P.dtype must be "int32".'); //  ________________   ______                   ___________
  // |                | |\(MxI)|                 |           |
  // |                | | \ R  |  ___________    |           |
  // |     (NxM)      | |  \   | |   (IxJ)   |   |   (NxJ)   |
  // |       Q        | |   \  | |     X     | = |     Y     |
  // |                | |    \ | |           |   |           |
  // |                | |     \| |___________|   |           |
  // |                | |  0   |                 |           |
  // |________________| |______|                 |___________|

  var _Q$shape$slice = Q.shape.slice(-2),
      _Q$shape$slice2 = (0, _slicedToArray2["default"])(_Q$shape$slice, 2),
      N = _Q$shape$slice2[0],
      M = _Q$shape$slice2[1],
      _R$shape$slice3 = R.shape.slice(-1),
      _R$shape$slice4 = (0, _slicedToArray2["default"])(_R$shape$slice3, 1),
      I = _R$shape$slice4[0],
      _y$shape$slice = y.shape.slice(-1),
      _y$shape$slice2 = (0, _slicedToArray2["default"])(_y$shape$slice, 1),
      J = _y$shape$slice2[0];

  if (N != y.shape[y.ndim - 2]) throw new Error("rrqr_lstsq(Q,R,P,y): Q and y don't match.");
  if (M != R.shape[R.ndim - 2]) throw new Error("rrqr_lstsq(Q,R,P,y): Q and R don't match.");
  if (I != P.shape[P.ndim - 1]) throw new Error("rrqr_lstsq(Q,R,P,y): R and P don't match.");
  var ndim = Math.max(Q.ndim, R.ndim, P.ndim + 1, y.ndim),
      shape = Int32Array.from({
    length: ndim
  }, function () {
    return 1;
  });
  shape[ndim - 2] = I;
  shape[ndim - 1] = J; // FIND COMMON (BROADCASTED) SHAPE

  for (var _i6 = 0, _arr = [Q, R, y]; _i6 < _arr.length; _i6++) {
    var arr = _arr[_i6];

    for (var _i11 = ndim - 2, _j14 = arr.ndim - 2; _i11-- > 0 && _j14-- > 0;) {
      if (1 === shape[_i11]) shape[_i11] = arr.shape[_j14];else if (shape[_i11] != arr.shape[_j14] && arr.shape[_j14] != 1) throw new Error('rrqr_lstsq(Q,R,P,y): Q,R,P,y not broadcast-compatible.');
    }
  }

  for (var i = ndim - 2, j = P.ndim - 1; i-- > 0 && j-- > 0;) {
    if (1 === shape[i]) shape[i] = P.shape[j];else if (shape[i] != P.shape[j] && P.shape[j] != 1) throw new Error('rrqr_lstsq(Q,R,P,y): Q,R,P,y not broadcast-compatible.');
  } // GENERATE RESULT DATA


  var DTypeArray = _dt.ARRAY_TYPES[[Q, R, y].every(function (a) {
    return a.dtype === 'float32';
  }) ? 'float32' : 'float64'],
      tmp_rank = new DTypeArray(Math.min(M, I)),
      tmp_perm = new Int32Array(I),
      x_dat = new DTypeArray(shape.reduce(function (a, b) {
    return a * b;
  }, 1)),
      Q_dat = Q.data,
      R_dat = R.data,
      P_dat = P.data,
      y_dat = y.data;
  var Q_off = 0,
      Q_stride = 1,
      R_off = 0,
      R_stride = 1,
      P_off = 0,
      P_stride = 1,
      y_off = 0,
      y_stride = 1,
      x_off = 0;

  function solv(d) {
    if (d === ndim - 2) {
      Q_stride = N * M;
      R_stride = M * I;
      P_stride = I;
      y_stride = N * J;

      var _R = _rrqr_rank(M, I, R_dat, R_off, tmp_rank); // Q.T @ y


      for (var k = 0; k < N; k++) {
        for (var _i7 = 0; _i7 < _R; _i7++) {
          for (var _j11 = 0; _j11 < J; _j11++) {
            x_dat[x_off + _i7 * J + _j11] += Q_dat[Q_off + k * M + _i7] * y_dat[y_off + k * J + _j11];
          }
        }
      } // BACKWARD SUBSTITUTION


      for (var _i8 = _R; _i8-- > 0;) {
        for (var _j12 = J; _j12-- > 0;) {
          for (var _k = _R; --_k > _i8;) {
            x_dat[x_off + _i8 * J + _j12] -= R_dat[R_off + I * _i8 + _k] * x_dat[x_off + _k * J + _j12];
          }

          x_dat[x_off + _i8 * J + _j12] /= R_dat[R_off + I * _i8 + _i8];
        }
      } // APPLY P TO X (PERMUTE ROWS)
      // https://www.geeksforgeeks.org/reorder-a-array-according-to-given-indexes/


      for (var _i9 = I; _i9-- > 0;) {
        tmp_perm[_i9] = P_dat[P_off + _i9];
      }

      for (var _i10 = I; _i10-- > 0;) {
        var _k2 = void 0;

        while ((_k2 = tmp_perm[_i10]) !== _i10) {
          if (tmp_perm[_k2] === _k2) throw new Error("rrqr_lstsq(Q,R,P,y): Invalid indices in P.");
          tmp_perm[_i10] = tmp_perm[_k2];
          tmp_perm[_k2] = _k2;

          for (var _j13 = J; _j13-- > 0;) {
            var x_ij = x_dat[x_off + J * _i10 + _j13];
            x_dat[x_off + J * _i10 + _j13] = x_dat[x_off + J * _k2 + _j13];
            x_dat[x_off + J * _k2 + _j13] = x_ij;
          }
        }
      }

      Q_off += Q_stride;
      R_off += R_stride;
      P_off += P_stride;
      y_off += y_stride;
      x_off += I * J;
      return;
    }

    for (var l = shape[d];; l--) {
      solv(d + 1);
      if (l == 1) break;
      if (!(Q.shape[d - ndim + Q.ndim] > 1)) Q_off -= Q_stride;
      if (!(R.shape[d - ndim + R.ndim] > 1)) R_off -= R_stride;
      if (!(P.shape[d - ndim + 1 + P.ndim] > 1)) P_off -= P_stride;
      if (!(y.shape[d - ndim + y.ndim] > 1)) y_off -= y_stride;
    }

    Q_stride *= Q.shape[d - ndim + Q.ndim] || 1;
    R_stride *= R.shape[d - ndim + R.ndim] || 1;
    P_stride *= P.shape[d - ndim + 1 + P.ndim] || 1;
    y_stride *= y.shape[d - ndim + y.ndim] || 1;
  }

  solv(0);
  return new _nd_array.NDArray(shape, x_dat);
}