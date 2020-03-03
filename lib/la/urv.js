'use strict';
/* This file is part of ND4JS.
 *
 * ND4JS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ND4JS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ND4JS. If not, see <http://www.gnu.org/licenses/>.
 */

var _interopRequireDefault = require("@babel/runtime/helpers/interopRequireDefault");

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports._urv_decomp_full = _urv_decomp_full;
exports.urv_decomp_full = urv_decomp_full;
exports._urv_lstsq = _urv_lstsq;
exports.urv_lstsq = urv_lstsq;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _dt = require("../dt");

var _nd_array = require("../nd_array");

var _srrqr = require("./srrqr");

var _giv_rot = require("./_giv_rot");

var _tri = require("./tri");

// TODO: add economic URV decomposition
function _urv_decomp_full(M, N, R, R_off, V, V_off, P, P_off) {
  if (0 !== M % 1) throw new Error('Assertion failed.');
  if (0 !== N % 1) throw new Error('Assertion failed.');
  if (0 !== R_off % 1) throw new Error('Assertion failed.');
  if (0 !== V_off % 1) throw new Error('Assertion failed.');
  /*DEBUG*/

  for (var i = N; i-- > 0;) {
    /*DEBUG*/
    for (var j = N; j-- > 0;) {
      /*DEBUG*/
      if (0 !== V[V_off + N * i + j]) throw new Error('Assertion failed.');
    }
  }

  M |= 0;
  N |= 0;
  R_off |= 0;
  V_off |= 0;
  if (!(M <= N)) throw new Error('Assertion failed.');

  if (M === N) {
    for (var _i = N; _i-- > 0;) {
      var _j = P[P_off + _i];
      V[V_off + N * _i + _j] = 1;
    }
  } else {
    // INIT V TO IDENTITY
    for (var _i2 = N; _i2-- > 0;) {
      V[V_off + N * _i2 + _i2] = 1;
    }

    for (var _i3 = M; _i3-- > 0;) {
      for (var _j2 = M; _j2 < N; _j2++) {
        var ij = R_off + N * _i3 + _j2,
            R_ij = R[ij];
        if (0 === R_ij) continue;
        var ii = R_off + N * _i3 + _i3,
            R_ii = R[ii];

        var _giv_rot_qr2 = (0, _giv_rot._giv_rot_qr)(R_ii, R_ij),
            _giv_rot_qr3 = (0, _slicedToArray2["default"])(_giv_rot_qr2, 3),
            c = _giv_rot_qr3[0],
            s = _giv_rot_qr3[1],
            norm = _giv_rot_qr3[2];

        if (s !== 0) {
          // ROT COLUMNS IN R
          for (var k = _i3; k-- > 0;) // <- TODO can this be made cache friendlier?
          {
            var ki = R_off + N * k + _i3,
                R_ki = R[ki],
                kj = R_off + N * k + _j2,
                R_kj = R[kj];
            R[ki] = R_kj * s + R_ki * c;
            R[kj] = R_kj * c - R_ki * s;
          } // ROT ROWS IN V


          (0, _giv_rot._giv_rot_rows)(V, 1 + _j2 - _i3, V_off + N * _i3 + _i3, V_off + N * _j2 + _i3, c, s);
          R[ii] = norm;
        }

        R[ij] = 0;
      }
    } // apply colum permutations
    // nested loop but actually only O(N) swaps


    for (var _i4 = 0; _i4 < N; ++_i4) {
      for (;;) // <- start swap cycle
      {
        var _j3 = P[P_off + _i4];
        P[P_off + _i4] = P[P_off + _j3];
        P[P_off + _j3] = _j3;
        if (_j3 <= _i4) break; // SWAP COLUMNS IN V

        for (var _k = N; _k-- > 0;) {
          var tmp = V[V_off + N * _k + _i4];
          V[V_off + N * _k + _i4] = V[V_off + N * _k + _j3];
          V[V_off + N * _k + _j3] = tmp;
        }
      }
    }
  }
}

function urv_decomp_full(A) // <- TODO add tolerance parameters to pass on to SRRQR
{
  var _srrqr_decomp_full = (0, _srrqr.srrqr_decomp_full)(A),
      _srrqr_decomp_full2 = (0, _slicedToArray2["default"])(_srrqr_decomp_full, 4),
      U = _srrqr_decomp_full2[0],
      RR = _srrqr_decomp_full2[1],
      P = _srrqr_decomp_full2[2],
      rnks = _srrqr_decomp_full2[3];

  A = undefined;

  var _RR$shape$slice = RR.shape.slice(-2),
      _RR$shape$slice2 = (0, _slicedToArray2["default"])(_RR$shape$slice, 2),
      M = _RR$shape$slice2[0],
      N = _RR$shape$slice2[1],
      V_shape = RR.shape.slice();

  V_shape[V_shape.length - 2] = N;
  P = P.data;
  var r = rnks.data,
      R = RR.data,
      V = new _dt.ARRAY_TYPES[RR.dtype](R.length / M * N);

  for (var R_off = 0, V_off = 0, P_off = 0, r_off = 0; r_off < r.length; R_off += M * N, V_off += N * N, P_off += N, r_off += 1) {
    var rnk = r[r_off];

    if (rnk < N) {
      // set close-to-zero lower right region to zero
      for (var i = rnk; i < M; i++) {
        for (var j = rnk; j < N; j++) {
          R[R_off + N * i + j] = 0;
        }
      } // <- TODO consider doing this in SRRQR already

    }

    _urv_decomp_full(rnk, N, R, R_off, V, V_off, P, P_off); // <- TODO if we init V to identity then we can do faster rotations of V. P would then be applied after computing V.

  }

  return [U, RR, new _nd_array.NDArray(V_shape, V), rnks];
}

function _urv_lstsq(rnk, I, J, K, L, M, U, U_off, R, R_off, V, V_off, X, X_off, Y, Y_off, tmp) {
  //   ┏        ┓ ┏        ┓ ┏        ┓ ┏        ┓   ┏        ┓
  //   ┃        ┃ ┃        ┃ ┃        ┃ ┃        ┃   ┃        ┃
  //   ┃ U[I,J] ┃ ┃ R[J,K] ┃ ┃ V[K,L] ┃ ┃ X[L,M] ┃ = ┃ Y[I,M] ┃
  //   ┃        ┃ ┃        ┃ ┃        ┃ ┃        ┃   ┃        ┃
  //   ┗        ┛ ┗        ┛ ┗        ┛ ┗        ┛   ┗        ┛
  if (rnk % 1 !== 0) throw new Error('Assertion failed.');
  if (I % 1 !== 0) throw new Error('Assertion failed.');
  if (J % 1 !== 0) throw new Error('Assertion failed.');
  if (K % 1 !== 0) throw new Error('Assertion failed.');
  if (L % 1 !== 0) throw new Error('Assertion failed.');
  if (M % 1 !== 0) throw new Error('Assertion failed.');
  if (U_off % 1 !== 0) throw new Error('Assertion failed.');
  if (R_off % 1 !== 0) throw new Error('Assertion failed.');
  if (V_off % 1 !== 0) throw new Error('Assertion failed.');
  if (X_off % 1 !== 0) throw new Error('Assertion failed.');
  if (Y_off % 1 !== 0) throw new Error('Assertion failed.');
  if (tmp.length < I * M) throw new Error('Assertion failed.');
  tmp.fill(0.0, 0, I * M); // TMP = U.T @ Y

  for (var k = 0; k < I; k++) {
    for (var i = 0; i < rnk; i++) {
      for (var j = 0; j < M; j++) {
        tmp[M * i + j] += U[U_off + J * k + i] * Y[Y_off + M * k + j];
      }
    }
  } // TMP = R \ TMP


  (0, _tri._triu_solve)(rnk, K, M, R, R_off, tmp, 0); // TMP = V.T @ TMP

  for (var _k2 = 0; _k2 < rnk; _k2++) {
    for (var _i5 = 0; _i5 < L; _i5++) {
      for (var _j4 = 0; _j4 < M; _j4++) {
        X[X_off + M * _i5 + _j4] += V[V_off + L * _k2 + _i5] * tmp[M * _k2 + _j4];
      }
    }
  }
}

function urv_lstsq(U, R, V, ranks, Y) {
  if (Y == null) {
    if (null != ranks || null != V) throw new Error('urv_lstsq( U,R,V,ranks, Y ): Either 2 ([U,R,V,ranks], Y) or 5 arguments (U,R,V,ranks, Y) expected.');
    Y = R;
    var _U = U;

    var _U2 = (0, _slicedToArray2["default"])(_U, 4);

    U = _U2[0];
    R = _U2[1];
    V = _U2[2];
    ranks = _U2[3];
  }

  U = (0, _nd_array.asarray)(U);
  if (!(U.ndim >= 2)) throw new Error('urv_lstsq(U,R,V, Y): U.ndim must be at least 2.');
  R = (0, _nd_array.asarray)(R);
  if (!(R.ndim >= 2)) throw new Error('urv_lstsq(U,R,V, Y): R.ndim must be at least 2.');
  V = (0, _nd_array.asarray)(V);
  if (!(V.ndim >= 2)) throw new Error('urv_lstsq(U,R,V, Y): V.ndim must be at least 2.');
  Y = (0, _nd_array.asarray)(Y);
  if (!(Y.ndim >= 2)) throw new Error('urv_lstsq(U,R,V, Y): Y.ndim must be at least 2.');
  ranks = (0, _nd_array.asarray)(ranks);
  var U_shape = U.shape,
      R_shape = R.shape,
      V_shape = V.shape,
      Y_shape = Y.shape,
      ranks_shape = ranks.shape;
  var dtype = [U, R, V, Y].every(function (x) {
    return x.dtype === 'float32';
  }) ? 'float32' : 'float64';
  var DTypeArray = _dt.ARRAY_TYPES[dtype];
  var ndim = [U, R, V, Y].reduce(function (max, _ref) {
    var ndim = _ref.ndim;
    return Math.max(max, ndim);
  }, ranks.ndim + 2);
  var X_shape = new Int32Array(ndim);
  X_shape.fill(1); // FIND COMMON (BROADCASTED) SHAPE

  for (var _i6 = 0, _arr = [U, R, V, Y]; _i6 < _arr.length; _i6++) {
    var arr = _arr[_i6];

    for (var _i7 = ndim - 2, _j5 = arr.ndim - 2; _i7-- > 0 && _j5-- > 0;) {
      if (1 === X_shape[_i7]) X_shape[_i7] = arr.shape[_j5];else if (X_shape[_i7] != arr.shape[_j5] && arr.shape[_j5] != 1) throw new Error('urv_lstsq( U,R,V,ranks, Y ): U,R,V,ranks, Y not broadcast-compatible.');
    }
  }

  for (var i = ndim - 2, j = ranks.ndim; i-- > 0 && j-- > 0;) {
    if (1 === X_shape[i]) X_shape[i] = ranks.shape[j];else if (X_shape[i] != ranks.shape[j] && ranks.shape[j] != 1) throw new Error('urv_lstsq( U,R,V,ranks, Y ): U,R,V,ranks, Y not broadcast-compatible.');
  } // CHECK MATRIX SHAPES
  //   ┏        ┓ ┏        ┓ ┏        ┓ ┏        ┓   ┏        ┓
  //   ┃        ┃ ┃        ┃ ┃        ┃ ┃        ┃   ┃        ┃
  //   ┃ U[I,J] ┃ ┃ R[J,K] ┃ ┃ V[K,L] ┃ ┃ X[L,M] ┃ = ┃ Y[I,M] ┃
  //   ┃        ┃ ┃        ┃ ┃        ┃ ┃        ┃   ┃        ┃
  //   ┗        ┛ ┗        ┛ ┗        ┛ ┗        ┛   ┗        ┛


  var _U$shape$slice = U.shape.slice(-2),
      _U$shape$slice2 = (0, _slicedToArray2["default"])(_U$shape$slice, 2),
      I = _U$shape$slice2[0],
      J = _U$shape$slice2[1],
      _V$shape$slice = V.shape.slice(-2),
      _V$shape$slice2 = (0, _slicedToArray2["default"])(_V$shape$slice, 2),
      K = _V$shape$slice2[0],
      L = _V$shape$slice2[1],
      M = Y.shape[Y.ndim - 1];

  if (R.shape[R.ndim - 2] !== J) throw new Error('urv_lstsq( U,R,V,ranks, Y ): Matrix dimensions incompatible.');
  if (R.shape[R.ndim - 1] !== K) throw new Error('urv_lstsq( U,R,V,ranks, Y ): Matrix dimensions incompatible.');

  if (J !== K) {
    if (I < L) if (I !== J) throw new Error('Assertion failed.');else if (K !== L) throw new Error('Assertion failed.');
  }

  if (!(I >= J)) throw new Error('Assertion failed.');
  if (!(K <= L)) throw new Error('Assertion failed.');
  X_shape[ndim - 2] = L;
  X_shape[ndim - 1] = M;
  ranks = ranks.data;
  U = U.data;
  R = R.data;
  V = V.data;
  Y = Y.data;
  var X = new DTypeArray(X_shape.reduce(function (m, n) {
    return m * n;
  })),
      tmp = new DTypeArray(I * M);
  var U_off = U.length,
      U_stride = 1,
      R_off = R.length,
      R_stride = 1,
      V_off = V.length,
      V_stride = 1,
      Y_off = Y.length,
      Y_stride = 1,
      ranks_off = ranks.length,
      ranks_stride = 1,
      X_off = X.length;

  function solv(d) {
    if (d === ndim - 2) {
      U_stride = I * J;
      R_stride = J * K;
      V_stride = K * L;
      Y_stride = I * M;
      ranks_stride = 1;
      U_off -= U_stride;
      R_off -= R_stride;
      V_off -= V_stride;
      Y_off -= Y_stride;
      ranks_off -= ranks_stride;
      X_off -= L * M;
      var rnk = ranks[ranks_off];

      _urv_lstsq(rnk, I, J, K, L, M, U, U_off, R, R_off, V, V_off, X, X_off, Y, Y_off, tmp);
    } else {
      for (var l = X_shape[d];; l--) {
        solv(d + 1);
        if (l == 1) break;
        if (!(U_shape[d - ndim + U_shape.length] > 1)) U_off += U_stride;
        if (!(R_shape[d - ndim + R_shape.length] > 1)) R_off += R_stride;
        if (!(V_shape[d - ndim + V_shape.length] > 1)) V_off += V_stride;
        if (!(Y_shape[d - ndim + Y_shape.length] > 1)) Y_off += Y_stride;
        if (!(ranks_shape[d - ndim + 2 + ranks_shape.length] > 1)) ranks_off += ranks_stride;
      }

      U_stride *= U_shape[d - ndim + U_shape.length] || 1;
      R_stride *= R_shape[d - ndim + R_shape.length] || 1;
      V_stride *= V_shape[d - ndim + V_shape.length] || 1;
      Y_stride *= Y_shape[d - ndim + Y_shape.length] || 1;
      ranks_stride *= ranks_shape[d - ndim + 2 + ranks_shape.length] || 1;
    }
  }

  solv(0);
  return new _nd_array.NDArray(X_shape, X);
}