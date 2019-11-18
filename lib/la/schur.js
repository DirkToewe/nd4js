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
exports.schur_eigenvals = schur_eigenvals;
exports.schur_eigen = schur_eigen;
exports.schur_decomp = schur_decomp;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _dt = require("../dt");

var _nd_array = require("../nd_array");

var _math = require("../math");

var _hessenberg = require("./hessenberg");

var _root1d_bisect = require("../opt/root1d_bisect");

var _matmul = require("./matmul");

var _mutable_complex = require("../dt/mutable_complex");

function schur_eigenvals(T) {
  T = (0, _nd_array.asarray)(T);

  var _T$shape$slice = T.shape.slice(-1),
      _T$shape$slice2 = (0, _slicedToArray2["default"])(_T$shape$slice, 1),
      N = _T$shape$slice2[0],
      Λ_shape = T.shape.slice(0, -1),
      Λ = new _dt.ARRAY_TYPES['complex128'](T.data.length / N);

  if (T.shape[T.ndim - 2] != N) throw new Error('T is not square.');
  T = T.data;
  var T_off = 0;

  var t = function t(i, j) {
    return T[T_off + N * i + j];
  };

  for (var Λ_off = 0; T_off < T.length; T_off += N * N, Λ_off += N) {
    // COMPUTE EIGENVECTORS (right -> left)
    for (var j = N - 1; j >= 0; j--) {
      var i = j - 1;

      if (0 === j || 0 == t(j, i)) {
        //
        // 1x1 BLOCK
        //
        // the eigenvalue is the diagonal value
        Λ[Λ_off + j] = t(j, j);
      } else {
        //
        // 2x2 BLOCK
        //
        // STEP1: compute eigenpairs of the 2x2 matrix
        var T_ii = t(i, i),
            T_ij = t(i, j),
            T_ji = t(j, i),
            T_jj = t(j, j),
            tr = _math.math.add(T_ii, T_jj),
            det = _math.math.sub(_math.math.mul(T_ii, T_jj), _math.math.mul(T_ij, T_ji));

        var λ1 = void 0,
            λ2 = void 0;

        var sqrt = _math.math.sqrt(_math.math.sub(_math.math.mul(tr, tr), _math.math.mul(det, 4)));

        if (tr * tr >= 4 * det) {
          throw new Error("schur_eigenvals(T): T must not contain real eigenvalued 2x2 blocks.");
          var sign = tr < 0 ? -1.0 : +1.0;
          λ1 = 0.5 * (tr + sign * sqrt);
          λ2 = det * 2.0 / (tr + sign * sqrt); // <- Citardauq Formula
        } else {
          λ1 = _math.math.mul(_math.math.add(tr, sqrt), 0.5);
          λ2 = _math.math.mul(_math.math.sub(tr, sqrt), 0.5);
        }

        Λ[Λ_off + i] = λ1;
        Λ[Λ_off + j] = λ2;
        j--;
      }
    }
  }

  return new _nd_array.NDArray(Λ_shape, Λ);
}

function schur_eigen(Q, T) {
  Q = (0, _nd_array.asarray)(Q);
  T = (0, _nd_array.asarray)(T);
  if (Q.ndim != T.ndim) throw new Error('Q.ndim != T.ndim.');

  for (var i = T.ndim; i-- > 0;) {
    if (Q.shape[i] != T.shape[i]) throw new Error('Q.shape != T.shape.');
  }

  var _T$shape$slice3 = T.shape.slice(-1),
      _T$shape$slice4 = (0, _slicedToArray2["default"])(_T$shape$slice3, 1),
      N = _T$shape$slice4[0],
      V_shape = T.shape,
      Λ_shape = V_shape.slice(0, -1); // T is quasi-triangular. For simplification let's assume it's upper triangular:
  // ┌                        ┐
  // │ λ₁ ...                 │
  // │                        │
  // │ 0   λ₂ ...             │
  // │       .                │
  // │ 0   0   .              │
  // │           .            │
  // │ .       .   λₖ ...     │
  // │ .         .   .        │
  // │ .           .   .      │
  // │                   .    │
  // │ 0    ...        0   λₙ │
  // └                        ┘
  //
  // Let's say we want to find an eigenvector for λₖ, we can just solve:
  // ┌                                   ┐ ┌     ┐
  // │ λ₁-λₖ ...                         │ │ x₁  │
  // │                                   │ │     │
  // │ 0   λ₂-λₖ ...                     │ │ x₂  │
  // │          .                        │ │ .   │
  // │ 0     0    .                      │ │ .   │
  // │              .                    │ │ .   │ ! ⇀
  // │ .         .  λₖ₋₁-λₖ ...          │ │ xₖ₋₁│ = 0
  // │                                   │ │     │
  // │ .             .    λₖ-λₖ ...      │ │ 1   │
  // │                        .          │ │     │
  // │ .                 .      .        │ │ 0   │
  // │                            .      │ │ :   │
  // │ 0     .   .   .       0    λₙ-λₖ  │ │ 0   │
  // └                                   ┘ └     ┘
  //
  // Since T is quasi-triangular this is solvable via a modified Backward Substition.
  // For multiple (equivalent) eigenvalues, there is not guaranteed to be more than
  // one linearily independent eigenvector. If that is the case, we at some point
  // arrive at an row in the backward substitution where:
  //
  // (λₛ-λₖ)xₛ = 0·xₛ = 0 - T[s,s+1]·xₛ₊₁ ... - T[s,k-1]*xₖ₋₁ - T[s,k] = -xₛ ≠ 0.
  //
  // We can resolve this by setting:
  //
  // xₜ := 0; for t > s
  // xₛ := 1
  //
  // This will of course than mean that for s and k, we have the same eigenvector.
  //


  if (T.shape[T.ndim - 2] != N) throw new Error('Q is not square.');
  var ComplexArray = _dt.ARRAY_TYPES['complex128'],
      V = ComplexArray.from(T.data);
  T = undefined;

  var Λ = new ComplexArray(V.length / N),
      V_arr = V._array,
      Λ_arr = Λ._array,
      // temporary vectors for the eigenvalue computation
  v1 = new ComplexArray(N),
      v1_arr = v1._array,
      v2 = new ComplexArray(N),
      v2_arr = v2._array,
      norm_sum = v1._array.subarray(0, N),
      norm_max = v2._array.subarray(N);

  var v_i = new _mutable_complex.MutableComplex(NaN, NaN),
      v_j = new _mutable_complex.MutableComplex(NaN, NaN),
      det = new _mutable_complex.MutableComplex(NaN, NaN);
  var TOL = NaN;
  /** Computes indices j < J of an eigenvector v using backward substition (see amazing UTF-8 art above).
   */

  function computeVec(λ, v, J) {
    var K = Math.min(N, J + 2),
        v_arr = v._array;

    for (var j = J; j-- > 0;) {
      for (var k = K; --k > j;) {
        var re0 = v_arr[2 * k + 0],
            re1 = V_arr[2 * (V_off + N * j + k) + 0],
            im0 = v_arr[2 * k + 1],
            im1 = V_arr[2 * (V_off + N * j + k) + 1];
        v_arr[2 * j + 0] -= re0 * re1 - im0 * im1;
        v_arr[2 * j + 1] -= re0 * im1 + im0 * re1;
      }

      if (0 == j || t(j, j - 1) == 0) {
        //
        // 1x1 BLOCK
        //
        var T_jj_re = V_arr[2 * (V_off + N * j + j) + 0] - λ.re,
            T_jj_im = V_arr[2 * (V_off + N * j + j) + 1] - λ.im;
        v_j.re = v_arr[2 * j + 0];
        v_j.im = v_arr[2 * j + 1];

        if (Math.hypot(T_jj_re, T_jj_im) <= TOL) {
          // <- TODO
          if (v_j.abs() <= TOL) {
            // <- TODO add test case for this zeroness test
            // v is already a valid eigenvalue, let's return it
            v_arr[2 * j + 0] = 0;
            v_arr[2 * j + 1] = 0;
          } else {
            // v is invalid, let's reset
            v_arr[2 * j] = 1.0;
            v_arr.fill(0.0, 2 * j + 1, 2 * K);
          }
        } else {
          v_j['/='](T_jj_re, T_jj_im);
          v_arr[2 * j + 0] = v_j.re;
          v_arr[2 * j + 1] = v_j.im;
        }
      } else {
        //
        // 2x2 BLOCK
        //
        var _i = j - 1;

        for (var _k = K; --_k > j;) {
          var _re = v_arr[2 * _k + 0],
              _re2 = V_arr[2 * (V_off + N * _i + _k) + 0],
              _im = v_arr[2 * _k + 1],
              _im2 = V_arr[2 * (V_off + N * _i + _k) + 1];
          v_arr[2 * _i + 0] -= _re * _re2 - _im * _im2;
          v_arr[2 * _i + 1] -= _re * _im2 + _im * _re2;
        }

        var T_ii_re = V_arr[2 * (V_off + N * _i + _i) + 0] - λ.re,
            T_ii_im = V_arr[2 * (V_off + N * _i + _i) + 1] - λ.im,
            _T_jj_re = V_arr[2 * (V_off + N * j + j) + 0] - λ.re,
            _T_jj_im = V_arr[2 * (V_off + N * j + j) + 1] - λ.im,
            T_ij_re = V_arr[2 * (V_off + N * _i + j) + 0],
            T_ij_im = V_arr[2 * (V_off + N * _i + j) + 1],
            T_ji_re = V_arr[2 * (V_off + N * j + _i) + 0],
            T_ji_im = V_arr[2 * (V_off + N * j + _i) + 1];

        det['= c0*c1'](
        /*c0=*/
        T_ii_re, T_ii_im,
        /*c1=*/
        _T_jj_re, _T_jj_im);
        det['-= c0*c1'](
        /*c0=*/
        T_ij_re, T_ij_im,
        /*c1=*/
        T_ji_re, T_ji_im);
        if (det.re === 0 && det.im === 0) throw new Error('Assertion failed.');
        v_j['= c0*c1'](
        /*c0=*/
        T_ii_re, T_ii_im,
        /*c1=*/
        v_arr[2 * j + 0], v_arr[2 * j + 1]);
        v_j['-= c0*c1'](
        /*c0=*/
        T_ji_re, T_ji_im,
        /*c1=*/
        v_arr[2 * _i + 0], v_arr[2 * _i + 1]);
        v_j['/='](det.re, det.im);
        v_i['= c0*c1'](
        /*c0=*/
        _T_jj_re, _T_jj_im,
        /*c1=*/
        v_arr[2 * _i + 0], v_arr[2 * _i + 1]);
        v_i['-= c0*c1'](
        /*c0=*/
        T_ij_re, T_ij_im,
        /*c1=*/
        v_arr[2 * j + 0], v_arr[2 * j + 1]);
        v_i['/='](det.re, det.im);
        v_arr[2 * _i + 0] = v_i.re;
        v_arr[2 * _i + 1] = v_i.im;
        v_arr[2 * j + 0] = v_j.re;
        v_arr[2 * j + 1] = v_j.im;
        j--;
      }
    }
  }

  var V_off = 0;

  var t = function t(i, j) {
    return V[V_off + N * i + j];
  };

  for (var Λ_off = 0; V_off < V.length; V_off += N * N, Λ_off += N) {
    TOL = Math.sqrt(Number.EPSILON) * function () {
      // compute Frobenius norm
      var sum = 0,
          max = 0;
      var iEnd = V_off + N * N;

      for (var _i2 = V_off; _i2 < iEnd; _i2++) {
        var elem = _math.math.abs(V[_i2]);

        if (elem != 0) {
          // <- handles NaN (by making the result NaN)
          if (elem > max) {
            sum *= Math.pow(max / elem, 2);
            max = elem;
          }

          sum += Math.pow(elem / max, 2);
        }
      }

      return Math.sqrt(sum) * max;
    }();

    if (!(TOL >= 0)) throw new Error('Assertion failed.'); // COMPUTE EIGENVECTORS (right -> left)

    for (var j = N - 1; j >= 0; j--) {
      var _i3 = j - 1;

      if (0 == j || t(j, _i3) == 0) {
        //
        // 1x1 BLOCK
        //
        // the eigenvalue is the diagonal value
        var λ = t(j, j);
        Λ[Λ_off + j] = λ; // since 0*1 is zero, the eigenequation should be solvable for vec1[j] = 1
        // (unless there is a duplicate eigenvalue with linarily non-independent eigenvectors, but that will be resolved by computeVec)

        v1.fill(0.0, 0, j + 2);
        v1[j] = 1;
        computeVec(λ, v1, j); // write the solution in the the (j+1)-th column

        for (var k = Math.min(N, j + 2); k-- > 0;) {
          V_arr[2 * (V_off + N * k + j) + 0] = v1_arr[2 * k + 0];
          V_arr[2 * (V_off + N * k + j) + 1] = v1_arr[2 * k + 1];
        }
      } else {
        //
        // 2x2 BLOCK
        //
        // STEP1: compute eigenpairs of the 2x2 matrix
        var T_ii = t(_i3, _i3),
            T_ij = t(_i3, j),
            T_ji = t(j, _i3),
            T_jj = t(j, j),
            tr = _math.math.add(T_ii, T_jj),
            _det = _math.math.sub(_math.math.mul(T_ii, T_jj), _math.math.mul(T_ij, T_ji));

        var λ1 = void 0,
            λ2 = void 0;

        var sqrt = _math.math.sqrt(_math.math.sub(_math.math.mul(tr, tr), _math.math.mul(_det, 4)));

        if (tr * tr >= 4 * _det) {
          throw new Error('The given Schur decomposition must not contain 2x2 diagonal blocks with real eigenvalues.');
          var sign = tr < 0 ? -1.0 : +1.0;
          λ1 = 0.5 * (tr + sign * sqrt);
          λ2 = _det * 2.0 / (tr + sign * sqrt); // <- Citardauq Formula
        } else {
          λ1 = _math.math.mul(_math.math.add(tr, sqrt), 0.5);
          λ2 = _math.math.mul(_math.math.sub(tr, sqrt), 0.5);
        } // TODO: the whole following section should be feasible with only a single temporary vector instead of two (vec1,vec2)


        v1_arr.fill(0.0, 0, 2 * (j + 1));
        v2_arr.fill(0.0, 0, 2 * (j + 1)); // http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/

        if (_math.math.abs(T_ij) >= _math.math.abs(T_ji)) {
          v1[_i3] = T_ij;
          v1[j] = _math.math.sub(λ1, T_ii);
          v2[_i3] = T_ij;
          v2[j] = _math.math.sub(λ2, T_ii);
        } else {
          v1[j] = T_ji;
          v1[_i3] = _math.math.sub(λ1, T_jj);
          v2[j] = T_ji;
          v2[_i3] = _math.math.sub(λ2, T_jj);
        }

        computeVec(λ1, v1, _i3);
        computeVec(λ2, v2, _i3);

        for (var _k2 = j + 1; _k2-- > 0;) {
          V_arr[2 * (V_off + N * _k2 + _i3) + 0] = v1_arr[2 * _k2 + 0];
          V_arr[2 * (V_off + N * _k2 + _i3) + 1] = v1_arr[2 * _k2 + 1];
          V_arr[2 * (V_off + N * _k2 + j) + 0] = v2_arr[2 * _k2 + 0];
          V_arr[2 * (V_off + N * _k2 + j) + 1] = v2_arr[2 * _k2 + 1];
        }

        Λ[Λ_off + _i3] = λ1;
        Λ[Λ_off + j] = λ2;
        --j;
      }
    } // COMPUTE COLUMN NORMS


    norm_sum.fill(0.0);
    norm_max.fill(0.0);

    for (var _i4 = 0; _i4 < N; _i4++) {
      for (var J = 0; J < N * 2; J++) {
        var _j = J >>> 1,
            V_ij = Math.abs(V_arr[2 * (V_off + N * _i4) + J]);

        if (V_ij > 0) {
          if (V_ij > norm_max[_j]) {
            var scale = norm_max[_j] / V_ij;
            norm_max[_j] = V_ij;
            norm_sum[_j] *= scale * scale;
          }

          var ratio = V_ij / norm_max[_j];
          norm_sum[_j] += ratio * ratio;
        }
      }
    }

    var norm = norm_sum;

    for (var _i5 = 0; _i5 < N; _i5++) {
      var max = norm_max[_i5];
      norm[_i5] = isFinite(max) ? Math.sqrt(norm_sum[_i5]) * max : max;
    } // NORMALIZE COLUMNS


    for (var _i6 = 0; _i6 < N; _i6++) {
      for (var _j2 = 0; _j2 < N; _j2++) {
        V[V_off + N * _i6 + _j2] = _math.math.div(V[V_off + N * _i6 + _j2], norm[_j2]);
      }
    }
  }

  return [new _nd_array.NDArray(Λ_shape, Λ), (0, _matmul.matmul2)(Q, new _nd_array.NDArray(V_shape, V))];
}

function schur_decomp(A) {
  A = (0, _nd_array.asarray)(A); // HESSENBERG DECOMPOSITION

  var N = A.shape[A.ndim - 1],
      _hessenberg_decomp = (0, _hessenberg.hessenberg_decomp)(A),
      _hessenberg_decomp2 = (0, _slicedToArray2["default"])(_hessenberg_decomp, 2),
      Q = _hessenberg_decomp2[0],
      H = _hessenberg_decomp2[1];

  A = undefined; // FRANCIS QR

  schur_qrfrancis_inplace(Q, H);
  return [Q, H];
}
/** Takes a Hessenberg Decomposition as input and performs a real Schur Decomposition IN-PLACE.
 *  Does not perform any scaling.
 *  Does not check Hessenberg property.
 */


function schur_qrfrancis_inplace(Q, H) {
  var N = Q.shape[Q.ndim - 1],
      DTypeArray = _dt.ARRAY_TYPES[Q.dtype],
      tmp = new DTypeArray(N);
  if (Q.shape[Q.ndim - 2] != N) throw new Error('Q is not square.');
  if (Q.ndim != H.ndim) throw new Error('Q.ndim != H.ndim.');

  for (var i = Q.ndim; i-- > 0;) {
    if (Q.shape[i] != H.shape[i]) throw new Error('Q.shape != H.shape.');
  }

  Q = Q.data;
  H = H.data;
  var TOL = Number.EPSILON;
  /** This function is recursively called to perform the Francis QR Algorithm, which is
   *  an implicit double-shift version of the QR Algorithm. This function only works on
   *  a subregion of H and an independent Q on each recursive call. That nested Q is
   *  applied to the remaining H and the outer Q once the nested schur decomposition
   *  is finished.
   *
   *  SEE: Prof. Dr. Peter Arbenz
   *       252-0504-00 G
   *       Numerical Methods for Solving Large Scale Eigenvalue Problems
   *       (Spring semester 2018)
   *       Chapter 3: The QR Algorithm
   *       http://people.inf.ethz.ch/arbenz/ewp/Lnotes/2010/chapter3.pdf 
   */

  var francis_qr = function francis_qr(Q, Q_off, Q_stride, H_off) {
    var h = function h(i, j) {
      return H[H_off + N * i + j];
    },
        is_zero = function is_zero(i) {
      if (Math.abs(h(i, i - 1)) > TOL * (Math.abs(h(i - 1, i - 1)) + Math.abs(h(i, i)))) // <- goes to else on NaN
        return false;else {
        H[H_off + N * i + i - 1] *= 0.0; // <- Handles NaN. If a value is that small, its digits are likely nonsense (due to cancellation error) so let's set it to zero.

        return true;
      }
    },

    /** Applies a two-sided given rotation.
     */
    giv = function giv(i, j, c, s) {
      if (j <= i) throw new Error('Assertion failed.'); // ROTATE ROWS IN H

      for (var k = Math.max(0, i - 1); k < Q_stride; k++) {
        var H_i = H[H_off + N * i + k],
            H_j = H[H_off + N * j + k];
        H[H_off + N * i + k] = s * H_j + c * H_i;
        H[H_off + N * j + k] = c * H_j - s * H_i;
      } // ROTATE COLUMNS IN H


      for (var _k3 = Math.min(Q_stride, j + 2); _k3-- > 0;) {
        var _H_i = H[H_off + N * _k3 + i],
            _H_j = H[H_off + N * _k3 + j];
        H[H_off + N * _k3 + i] = s * _H_j + c * _H_i;
        H[H_off + N * _k3 + j] = c * _H_j - s * _H_i;
      } // ROTATE ROWS IN Q


      for (var _k4 = Q_stride; _k4-- > 0;) {
        var Q_i = Q[Q_off + Q_stride * i + _k4],
            Q_j = Q[Q_off + Q_stride * j + _k4];
        Q[Q_off + Q_stride * i + _k4] = s * Q_j + c * Q_i;
        Q[Q_off + Q_stride * j + _k4] = c * Q_j - s * Q_i;
      }
    },

    /** Recursively performs the schur-decomposition of the sub-region [s,e).
     * 
     *  During iteration, the givens rotations are only applied to current
     *  subregion and are accumulated in a temporary matrix. Only after the
     *  Francis QR algorithm is done, the transformations are applied to 
     */
    recurse = function recurse(s, e) {
      if (e - s > Q_stride >>> 1) throw new Error('Assertion failed.'); // <- assert that memory is bounded by O(n)

      var n = e - s;
      if (n < 3) throw new Error('Assertion failed.');
      var q = new DTypeArray(n * n); // INIT q TO IDENTITY

      for (var _i7 = 0; _i7 < n; _i7++) {
        q[n * _i7 + _i7] = 1.0;
      } // RUN FRANCIS QR ON SUB-PROBLEM


      francis_qr(q, 0, n, H_off + N * s + s); // TRANSPOSE q

      for (var _i8 = 0; _i8 < n; _i8++) {
        for (var j = 0; j < _i8; j++) {
          var q_ij = q[n * _i8 + j];
          q[n * _i8 + j] = q[n * j + _i8];
          q[n * j + _i8] = q_ij;
        }
      } // APPLY q TO LEFT OF Q (Q' = q.T @ Q)


      for (var _i9 = 0; _i9 < Q_stride; _i9++) // <- each column in Q
      {
        tmp.fill(0.0, 0, n);

        for (var _j3 = 0; _j3 < n; _j3++) {
          for (var k = 0; k < n; k++) {
            tmp[k] += q[n * _j3 + k] * Q[Q_off + Q_stride * (s + _j3) + _i9];
          }
        }

        for (var _j4 = 0; _j4 < n; _j4++) {
          Q[Q_off + Q_stride * (s + _j4) + _i9] = tmp[_j4];
        }
      } // APPLY q TO LEFT OF H (H' = q.T @ H)


      for (var _i10 = e; _i10 < Q_stride; _i10++) // <- each column in H
      {
        tmp.fill(0.0, 0, n);

        for (var _j5 = 0; _j5 < n; _j5++) {
          for (var _k5 = 0; _k5 < n; _k5++) {
            tmp[_k5] += q[n * _j5 + _k5] * H[H_off + N * (s + _j5) + _i10];
          }
        }

        for (var _j6 = 0; _j6 < n; _j6++) {
          H[H_off + N * (s + _j6) + _i10] = tmp[_j6];
        }
      } // APPLY q TO RIGHT H (H" = H' @ q)


      for (var _i11 = 0; _i11 < s; _i11++) // <- each row in H
      {
        tmp.fill(0.0, 0, n);

        for (var _j7 = 0; _j7 < n; _j7++) {
          for (var _k6 = 0; _k6 < n; _k6++) {
            tmp[_k6] += q[n * _j7 + _k6] * H[H_off + N * _i11 + (s + _j7)];
          }
        }

        for (var _j8 = 0; _j8 < n; _j8++) {
          H[H_off + N * _i11 + (s + _j8)] = tmp[_j8];
        }
      }
    };

    var stuck_o_meter = 0,
        start = 0,
        end = Q_stride;

    while (true) {
      // DETECT ZEROS ON THE SUB-DIAGONAL AND SHRINK WORK SIZE ACCORDINGLY
      for (var done = false; !done;) {
        if (end - start < 3) return;
        if (end - start < Q_stride >>> 4) return recurse(start, end); // <- ZOOM IN IF SUBPROBLEM IS SMALL ENOUGH

        done = false;
        if (is_zero(start + 1)) start += 1;else if (is_zero(end - 1)) end -= 1;else if (is_zero(start + 2)) start += 2;else if (is_zero(end - 2)) end -= 2;else done = true;
        var mid = start + end >>> 1;

        for (var _i13 = start + 2; done && ++_i13 < end - 2;) {
          if (is_zero(_i13)) {
            done = false; // RUN NESTED FRANCIS QR ON THE SMALLER OF THE TWO 

            if (_i13 > mid) {
              recurse(_i13, end);
              end = _i13;
            } else {
              recurse(start, _i13);
              start = _i13;
            }
          }
        }

        if (!done) stuck_o_meter = 0;
      }

      stuck_o_meter += 1; // DETERMINE (DOUBLE) SHIFT FROM LOWER RIGHT 2x2 BLOCK

      var _i12 = end - 2,
          j = end - 1,
          tr = h(_i12, _i12) + h(j, j),
          det = h(_i12, _i12) * h(j, j) - h(_i12, j) * h(j, _i12); // FOR REAL EIGENVALUES LETS USE THE ONE THAT'S CLOSER TO THE LOWER RIGHT ENTRY


      if (tr * tr > 4 * det) {
        var sign = tr >= 0 ? +1.0 : -1.0,
            ev1 = 0.5 * (tr + sign * Math.sqrt(tr * tr - 4 * det)),
            ev2 = det * 2.0 / (tr + sign * Math.sqrt(tr * tr - 4 * det)); // <- Citardauq Formula
        // use the eigenvalue closer to A[j,j]
        // SEE: Bindel, Fall 2016, Matrix Computations (CS 6210)
        //      Notes for 2016-10-24
        //      https://www.cs.cornell.edu/~bindel/class/cs6210-f16/lec/2016-10-24.pdf

        if (Math.abs(h(j, j) - ev1) > Math.abs(h(j, j) - ev2)) ev1 = ev2;
        tr = ev1 * 2;
        det = ev1 * ev1;
      } // IF WE'RE STUCK, LET'S WIGGLE LIKE A FISH ON LAND (... well except maybe that fella: https://www.youtube.com/watch?v=fJLCSsnhLFc)
      // SEE: NUMERICAL RECIPES Webnote No. 16, Rev. 1
      //      Description of the QR Algorithm for Hessenberg Matrices
      //      http://numerical.recipes/webnotes/nr3web16.pdf


      if (stuck_o_meter % 16 == 0) {
        if (stuck_o_meter > 1e9) throw new Error('Too many iterations for a single eigenvalue.');
        tr = Math.abs(h(j, _i12)) + Math.abs(h(_i12, end - 3));
        det = tr * tr;
        tr *= Math.random() * 0.5 + 1.25;
      } // FIRST COLUMN OF DOUBLE SHIFTED MATRIX (H-sI)*(H-conj(s)I) = H² - 2*Re(s)*H + |s|²I


      _i12 = start + 0;
      j = start + 1;
      var k = start + 2,
          a1 = h(_i12, _i12) * h(_i12, _i12) + h(_i12, j) * h(j, _i12) - tr * h(_i12, _i12) + det,
          a2 = h(j, _i12) * (h(_i12, _i12) + h(j, j) - tr),
          a3 = h(j, _i12) * h(k, j); // APPLY "DOUBLE SHIFTED" GIVENS ROTATIONS

      for (var row = 2; row-- > 0; j = k, a2 = a3) {
        if (a2 != 0) {
          var norm = Math.hypot(a1, a2),
              c = a1 / norm,
              s = a2 / norm;
          giv(_i12, j, c, s); // if( Math.abs(c*a2 - s*a1) > 1e-8 ) throw new Error('Assertion failed.')

          a1 = c * a1 + s * a2;
        }
      } // REINSTATE HESSENBERG PROPERTY


      for (var col = start; col < end - 2; col++) {
        _i12 = col + 1;
        var J = Math.min(end, col + 4);

        for (j = col + 2; j < J; j++) {
          var H_i = h(_i12, col),
              H_j = h(j, col);
          if (H_j == 0) continue;

          var _norm = Math.hypot(H_i, H_j),
              _c = H_i / _norm,
              _s = H_j / _norm;

          giv(_i12, j, _c, _s); // if( Math.abs( h(j,col) ) > 1e-8 ) throw new Error('Assertion failed!')

          H[H_off + N * j + col] *= 0.0; // <- handles NaN
        }
      }
    }
  }; // BEGIN SCHUR DECOMPOSING MATRICES


  for (var off = 0; off < Q.length; off += N * N) {
    // TRANSPOSE Q
    for (var _i14 = 0; _i14 < N; _i14++) {
      for (var j = 0; j < _i14; j++) {
        var Q_ij = Q[off + N * _i14 + j];
        Q[off + N * _i14 + j] = Q[off + N * j + _i14];
        Q[off + N * j + _i14] = Q_ij;
      }
    } // RUN FRANCIS QR ALGORITHM


    francis_qr(Q, off, N, off); // BEGIN RESOLVE REAL-VALUED 2x2 BLOCKS

    for (var _j9 = 1; _j9 < N; _j9++) {
      var _i15 = _j9 - 1;

      if (H[off + N * _j9 + _i15] != 0) {
        var _ret = function () {
          // The goal is to find a givens rotation that Schur-decomposes a real-eigenvalue 2x2 matrix.
          // ┌                ┐ ┌            ┐ ┌                ┐   ┌      ┐
          // │ cos(α) -sin(α) │ │ H_ii  H_ij │ │ cos(α) -sin(α) │ ! │ λ₁ p │
          // │                │ │            │ │                │ = │      │
          // │ sin(α)  cos(α) │ │ H_ji  H_jj │ │ sin(α)  cos(α) │   │ 0  λ₂│ => 0 == (H_ji⋅cos(α) - H_ii⋅sin(α))⋅cos(α) + (H_jj⋅cos(α) - H_ij⋅sin(α))⋅sin(α)
          // └                ┘ └            ┘ └                ┘   └      ┘ => 0 == (H_jj-H_ii)⋅sin(2⋅α) + (H_ij+H_ji)⋅cos(2⋅α) + H_ji-H_ij =: f(α)
          //
          // A simple and very numerically accurate solution would be Binary Search. In order to do that, we have to bracket a solution. So let's determine the extrema of f(α).
          // 
          // f'(α_max) =!= 0 = 2*(H_jj-H_ii)⋅cos(2⋅α_max) - 2*(H_ij+H_ji)⋅sin(2⋅α_max) => α_max = atan2( H_jj-H_ii, H_ij+H_ji ) / 2 + n*π/2
          var H_ii = H[off + N * _i15 + _i15],
              H_ij = H[off + N * _i15 + _j9],
              H_ji = H[off + N * _j9 + _i15],
              H_jj = H[off + N * _j9 + _j9],
              tr = H_ii + H_jj,
              det = H_ii * H_jj - H_ij * H_ji;
          if (tr * tr < 4 * det) return "continue";
          var α_min = Math.atan2(H_jj - H_ii, H_ij + H_ji) / 2,
              α_max = α_min + Math.PI / 2 * (α_min <= 0 ? +1 : -1),
              α = (0, _root1d_bisect.root1d_bisect)(function (α) {
            return (H_ji * Math.cos(α) - H_ii * Math.sin(α)) * Math.cos(α) + (H_jj * Math.cos(α) - H_ij * Math.sin(α)) * Math.sin(α);
          }, α_min, α_max),
              c = Math.cos(α),
              s = Math.sin(α);

          for (var k = _i15; k < N; k++) // ROTATE ROWS IN H
          {
            var H_i = H[off + N * _i15 + k],
                H_j = H[off + N * _j9 + k];
            H[off + N * _i15 + k] = s * H_j + c * H_i;
            H[off + N * _j9 + k] = c * H_j - s * H_i;
          }

          for (var _k7 = _j9 + 1; _k7-- > 0;) // ROTATE COLUMNS IN H
          {
            var _H_i2 = H[off + N * _k7 + _i15],
                _H_j2 = H[off + N * _k7 + _j9];
            H[off + N * _k7 + _i15] = s * _H_j2 + c * _H_i2;
            H[off + N * _k7 + _j9] = c * _H_j2 - s * _H_i2;
          }

          for (var _k8 = N; _k8-- > 0;) // ROTATE ROWS IN Q
          {
            var Q_i = Q[off + N * _i15 + _k8],
                Q_j = Q[off + N * _j9 + _k8];
            Q[off + N * _i15 + _k8] = s * Q_j + c * Q_i;
            Q[off + N * _j9 + _k8] = c * Q_j - s * Q_i;
          }

          H[off + N * _j9 + _i15] *= 0.0; // <- Handles NaN. If a value is that small, its digits are likely nonsense (due to cancellation error) so let's set it to zero.
        }();

        if (_ret === "continue") continue;
      }
    } // END RESOLVE REAL-VALUED 2x2 BLOCKS
    // TRANSPOSE Q BACK


    for (var _i16 = 0; _i16 < N; _i16++) {
      for (var _j10 = 0; _j10 < _i16; _j10++) {
        var _Q_ij = Q[off + N * _i16 + _j10];
        Q[off + N * _i16 + _j10] = Q[off + N * _j10 + _i16];
        Q[off + N * _j10 + _i16] = _Q_ij;
      }
    }
  } // END SCHUR DECOMPOSING MATRICES

} // END schur_qrfrancis_inplace