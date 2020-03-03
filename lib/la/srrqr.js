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
exports._row_norm_update = _row_norm_update;
exports.srrqr_decomp_full = srrqr_decomp_full;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _dt = require("../dt");

var _math = _interopRequireDefault(require("../math"));

var _nd_array = require("../nd_array");

var _giv_rot = require("./_giv_rot");

var _norm5 = require("./norm");

var _rrqr = require("./rrqr");

var _transpose_inplace2 = require("./transpose_inplace");

// TODO: add economic srrqr
// FIXME: there seems to be an infinite-loop bug with close to zero matrices like:
//   [[0,0,0],
//    [0,0,0],
//    [Number.MIN_VALUE,0,0],
//    [0,Number.MIN_VALUE,0],
//    [0,0,Number.MIN_VALUE]]
function _row_norm_update(norm, AB, AB_off, M) {
  // https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Math/hypot#Polyfill
  AB_off |= 0;
  M <<= 1;

  for (var j = 0; j < M; j += 2) {
    var s = Math.abs(AB[AB_off + (j >>> 1)]);

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

function srrqr_decomp_full(X) {
  var opt = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : {};
  // Overview
  // --------
  //
  // [1] "EFFICIENT ALGORITHMS FOR COMPUTING A STRONG RANK-REVEALING QR FACTORIZATION"
  //      Ming Gu, Stanley C. Eisenstat,
  //      https://math.berkeley.edu/~mgu/MA273/Strong_RRQR.pdf
  //
  // "Weak" vs. "Strong" RRQR
  // ------------------------
  // While the RRQR works well in practice, there are constructed examples of matrices
  // for which the 
  //
  // L2R vs. BIN RRQR
  // ----------------
  // Instead of test every potential rank `k` from left to right, as suggested [1],
  // we use binary search to find the proper rank. This requires significantly less
  // "strong" swaps in most cases and few more swaps in very low-rank cases.
  X = (0, _nd_array.asarray)(X);
  if (X.ndim < 2) throw new Error('srrqr_decomp_full(A,opt): A must be at least 2D.');
  if (X.dtype.startsWith('complex')) throw new Error('srrqr_decomp_full(A,opt): Complex A not (yet) supported.'); //
  //

  var dtype = X.dtype === 'float32' ? 'float32' : 'float64',
      DTypeArray = _dt.ARRAY_TYPES[dtype],
      R_shape = X.shape,
      Q_shape = Int32Array.from(R_shape),
      P_shape = Q_shape.slice(0, -1),
      r_shape = R_shape.slice(0, -2),
      _R_shape$slice = R_shape.slice(-2),
      _R_shape$slice2 = (0, _slicedToArray2["default"])(_R_shape$slice, 2),
      M = _R_shape$slice2[0],
      N = _R_shape$slice2[1],
      L = Math.min(M, N); // <- M not M-1, because the last row still needs pivotization


  Q_shape[Q_shape.length - 1] = M;
  P_shape[P_shape.length - 1] = N; // At each step of the decomposition we have:
  //
  //   X[:,P] = Q @ R
  //
  // Where R is partially triangularized, P keeps track of the column swaps and
  // Q is an orthogonals matrix.
  //
  // As part of the SRRQR we test different ranks `k` that X might actually have.
  // For every `k`, R[k] consists of the following quadrants:
  //          ┏         ┓
  //          ┃A[k]┊B[k]┃
  //   R[k] = ┃┈┈┈┈┼┈┈┈┈┃
  //          ┃    ┊C[k]┃
  //          ┗         ┛
  //
  // Where A[k] is already triangularized. We try to find a column permutation P that
  // norm(C[k],'fro') becomes as small as possible. We achieve that by swapping the
  // columns that maximize det(A[k]). Since:
  //
  //   det(X) = det(R[k]) = det(A[k])*det(C[k])
  //
  // norm(C[k],'fro') becomes smaller and smaller using the so-called "strong" column
  // swaps. Such swaps are only performed until the factor of determinant increase is
  // above a certain threshold:
  //                          !
  //   det(A'[k]) / det(A[k]) > dtol
  //
  // Where `dtol` is some small-ish constant value above 1. Without this threshold, an
  // exponential amount of column swaps might be performed. With the threshold however,
  // at most log[f](sqrt(n)) column swaps are required.
  //
  // Gu and Eisenstat [1] have shown how to compute a matrix `W` of size (k,n-k) which,
  // for each possible column swap, contains the factors by which det(A[k]) would increase
  // or decrease, i.e.:
  //
  // W[i,j] = det(A[k] with columns i,j swapped) / det(A[k])
  //
  // In order to compute `W`, we need inv(A[k]) and A[k]\B[k]. To compute these efficiently,
  // we need to keep track of them in memory during iteration in a matrix `AB`.
  // Gu and Eisenstat [1] proposed "update" and "downdate" methods that can be used to
  // adjust AB during column swaps or changes of `k`.
  //
  // Once we have found a series of column permutations that minimized norm(C[k],'fro'),
  // we can check if C[k] can be considered zero:
  //             ?
  //  norm(C[k]) <= ztol
  //
  // If this is the case, we know that the rank is `k` or less. Using this test inside
  // of a binary search algorithm, we can find the correct rank by testing log2(n)
  // different values for `k`.
  //
  // At each step of the binary search we have a lower and upper inclusive limit `k0` and
  // `K` of the actual rank, while we are currently testing `k`. Depending on the result
  // of the aforementioned test the new binary search range is going to be either [k0+1,k]
  // or [k,K]. In the former case we need avoid downdating from a (nearly) rank-deficient A[k].
  // To achieve this we need to keep track of inv(A[k0]) and A[k0]\B[k0] in a matrix `AB0`.
  // This allows us to update from `AB0` instead of ever downdating `AB`.

  var R = DTypeArray.from(X.data);
  X = undefined;
  var r = new Int32Array(R.length / N / M),
      P = new Int32Array(R.length / M),
      // <- tracks column permutations
  Q = new DTypeArray(R.length / N * M),
      AB = new DTypeArray(M * N),
      // <- stores inv(A) and A\B in COLUMN MAJOR order
  AB0 = new DTypeArray(M * N),
      norm = new DTypeArray(N << 1),
      NORM = new _norm5.FrobeniusNorm(); // <─ underflow-safe representation of the column norm

  for (var i = r.length; i-- > 0;) {
    r[i] = Math.random() * 1024 - 1024;
  }

  var ztol = NaN,
      dtol = NaN; // If no "strong" column swaps have been performed yet, Q has a special
  // structure that allows faster givens rotations while eliminating columns
  // of R. As soon as a "strong" swap occours, Q is "dirty" from now on.

  var Q_dirty = false;
  var Q_off = 0,
      R_off = 0,
      P_off = 0; // `k0` and `K` are the lower and upper inclusive bounds of the binary
  // search range while `k` is the currently considered potential rank

  var k0 = 0,
      k = 0,
      K = 0; // PROCESS THE (TOLERANCE) OPTIONS

  var _ref = function () {
    var _opt = opt,
        _opt$dtol = _opt.dtol,
        dtol = _opt$dtol === void 0 ? 1.01 : _opt$dtol,
        _opt$ztol = _opt.ztol,
        ztol = _opt$ztol === void 0 ? undefined : _opt$ztol;

    if (dtol instanceof _nd_array.NDArray) {
      throw new Error("srrqr_decomp_full(A,opt): NDArray as opt.dtol not yet supported.");
    } else {
      if (!(dtol >= 1)) throw new Error("srrqr_decomp_full(A,opt): Invalid opt.dtol: ".concat(dtol, ". Must be >=1."));

      var _DTOL = Number(dtol);

      dtol = function dtol() {
        return _DTOL;
      };
    }

    if (ztol instanceof _nd_array.NDArray) {
      throw new Error("srrqr_decomp_full(A,opt): NDArray as opt.ztol not yet supported.");
    } else if (null == ztol) {
      var _eps = Math.sqrt((0, _dt.eps)(dtype));

      ztol = function ztol(i) {
        i *= M * N;
        var I = i + M * N;
        NORM.reset();

        for (; i < I; i++) {
          NORM.include(R[i]);
        }

        return _eps * NORM.result * Math.max(M, N); // <- should be a little stricter than rrqr_rank such they work well together
      };
    } else {
      if (!(ztol >= 0)) throw new Error("srrqr_decomp_full(A,opt): invalid opt.ztol: ".concat(ztol, ". Must be non-negative number."));

      var _ZTOL = Number(ztol);

      ztol = function ztol() {
        return _ZTOL;
      };
    }

    return [dtol, ztol];
  }(),
      _ref2 = (0, _slicedToArray2["default"])(_ref, 2),
      DTOL = _ref2[0],
      ZTOL = _ref2[1];

  opt = undefined;
  /* Updates inv(A[k]) and A[k]\B[k]. Used to update both `AB` and `AB0`.
   */

  var update = function update(AB, k) {
    //*DEBUG*/    // CHECK TRIANGULARITY OF R AND AB
    //*DEBUG*/    for( let j=0;   j < k; j++)
    //*DEBUG*/    for( let i=j; ++i < M;    )
    //*DEBUG*/    {
    //*DEBUG*/      if( !(R [R_off + N*i+j  ] === 0) ) throw new Error(`Assertion failed.`);
    //*DEBUG*/      if( !(AB[          i+j*M] === 0) ) throw new Error(`Assertion failed.`);
    //*DEBUG*/    }
    // UPDATE inv(A)
    {
      var R_kk = -R[R_off + N * k + k];
      AB[k + k * M] = -1 / R_kk;

      for (var _i = k; _i-- > 0;) {
        AB[_i + k * M] /= R_kk;
      }
    } // UPDATE A\B

    for (var j = k; ++j < N;) {
      for (var _i2 = 0; _i2 <= k; _i2++) {
        AB[_i2 + j * M] += AB[_i2 + k * M] * R[R_off + N * k + j];
      }
    }
  };
  /* Downdates inv(A) and A\B. Used to downdate AB and AB0 once
   * during each "strong" column swaps. AB0 is not downdate if
   * the column swap does not affect it.
   */


  var downdate = function downdate(AB, k) {
    //*DEBUG*/    // CHECK TRIANGULARITY OF R AND AB
    //*DEBUG*/    for( let j=0;   j < k; j++)
    //*DEBUG*/    for( let i=j; ++i < M;    )
    //*DEBUG*/    {
    //*DEBUG*/      if( !(R [R_off + N*i+j  ] === 0) ) throw new Error(`Assertion failed.`);
    //*DEBUG*/      if( !(AB[          i+j*M] === 0) ) throw new Error(`Assertion failed.`);
    //*DEBUG*/    }
    // DOWNDATE A\B
    for (var j = N; --j > k;) {
      AB[k + j * M] = 0;

      for (var _i3 = k; _i3-- > 0;) {
        AB[_i3 + j * M] -= AB[_i3 + k * M] * R[R_off + N * k + j];
      }
    } // DOWNDATE inv(A)


    AB[k + k * M] = 0;
    {
      var R_kk = -R[R_off + N * k + k];

      for (var _i4 = k; _i4-- > 0;) {
        AB[_i4 + k * M] *= R_kk;
      }
    }
  };
  /* Swaps column `p` and `k` and eliminates the (new) column `k` in R.
   * A\B and A0\B0 are updated accordingly. `p` must not be less than `k`.
   */


  var swap_elim = function swap_elim(p) {
    //*DEBUG*/    if(!(k <= p)) throw new Error('Assertion failed.');
    // swap columns in R
    for (var j = 0; j < M; j++) {
      var R_jk = R[R_off + N * j + k];
      R[R_off + N * j + k] = R[R_off + N * j + p];
      R[R_off + N * j + p] = R_jk;
    } // swap columns in A\B (which is column major)


    for (var _j = 0; _j < k; _j++) {
      var AB_jk = AB0[_j + k * M];
      AB0[_j + k * M] = AB0[_j + p * M];
      AB0[_j + p * M] = AB_jk;
    }

    for (var _j2 = 0; _j2 < k; _j2++) {
      var _AB_jk = AB[_j2 + k * M];
      AB[_j2 + k * M] = AB[_j2 + p * M];
      AB[_j2 + p * M] = _AB_jk;
    } // swap P


    var P_i = P[P_off + k];
    P[P_off + k] = P[P_off + p];
    P[P_off + p] = P_i; // RESET COLUMN NORM (INDEX k IS SET TO ZERO FOR THE NEXT RRQR)

    norm.fill(0.0); // ELIMINATE COLUMN k BELOW DIAGONAL (USING GIVEN ROTATIONS)

    var count = 0;

    for (var _j3 = k; ++_j3 < M;) {
      var kk = R_off + N * k + k,
          jk = R_off + N * _j3 + k,
          _R_jk = R[jk];

      if (0 !== _R_jk) {
        var R_kk = R[kk],
            _giv_rot_qr2 = (0, _giv_rot._giv_rot_qr)(R_kk, _R_jk),
            _giv_rot_qr3 = (0, _slicedToArray2["default"])(_giv_rot_qr2, 3),
            c = _giv_rot_qr3[0],
            s = _giv_rot_qr3[1],
            _norm2 = _giv_rot_qr3[2];

        R[jk] = 0;

        if (s !== 0) {
          R[kk] = _norm2;
          (0, _giv_rot._giv_rot_rows)(R, N - 1 - k, kk + 1, jk + 1, c, s);
          (0, _giv_rot._giv_rot_rows)(Q, Q_dirty ? M : _j3 + 1, // <- "strong" column swaps "contaminate" Q -> be fully rotate
          Q_off + M * k, Q_off + M * _j3, c, s);
        }
      }

      (0, _rrqr._norm_update)(norm, R, R_off + N * _j3, k + 1);
    }
  };
  /* Swaps column `k` and the column with the largest remaining norm.
   * Afterwards, this method eliminates the (new) column `k` in R.
   * A\B and A0\B0 are updated accordingly.
   */


  var piv_elim = function piv_elim() {
    //*DEBUG*/    check_col_norms();
    var p = -1,
        norm_p = -Infinity;

    for (var j = k; j < N; j++) {
      var norm_j = (0, _rrqr._norm)(norm, j);

      if (norm_p < norm_j) {
        p = j;
        norm_p = norm_j;
      }
    } //*DEBUG*/    if(!(norm_p >= 0))
    //*DEBUG*/      throw new Error('Assertion failed.');


    swap_elim(p);
  }; //*DEBUG*/  const check = (AB, k) =>
  //*DEBUG*/  {
  //*DEBUG*/    // check inv(A)
  //*DEBUG*/    for( let i=0; i < k; i++ )
  //*DEBUG*/    for( let j=0; j < k; j++ )
  //*DEBUG*/    {
  //*DEBUG*/      let sum = 0;
  //*DEBUG*/      for( let h=0; h < k; h++ )
  //*DEBUG*/        sum += AB[i+h*M] * R[R_off + N*h+j];
  //*DEBUG*/      if( ! (Math.abs(sum - (i===j)) <= 1e-4) )
  //*DEBUG*/        throw new Error(`${sum} != ${(i===j)*1}.`);
  //*DEBUG*/    }
  //*DEBUG*/    // check A\B
  //*DEBUG*/    for( let i=0; i < k; i++ )
  //*DEBUG*/    for( let j=k; j < N; j++ )
  //*DEBUG*/    {
  //*DEBUG*/      let sum = 0;
  //*DEBUG*/      for( let h=0; h < k; h++ )
  //*DEBUG*/        sum += AB[i+h*M] * R[R_off + N*h+j];
  //*DEBUG*/      if( ! (Math.abs(sum - AB[i+j*M]) <= 1e-4) )
  //*DEBUG*/        throw new Error(`${sum} != ${AB[i+j*M]}.`);
  //*DEBUG*/    }
  //*DEBUG*/  }
  //*DEBUG*/
  //*DEBUG*/
  //*DEBUG*/  const check_col_norms = () =>
  //*DEBUG*/  {
  //*DEBUG*/    for( let j=k; j < N; j++ )
  //*DEBUG*/    {
  //*DEBUG*/      const nj = _norm(norm, j),
  //*DEBUG*/            uj = Math.hypot(...function*(){
  //*DEBUG*/              for( let i=k; i < M; i++ )
  //*DEBUG*/                yield R[R_off + N*i+j];
  //*DEBUG*/            }());
  //*DEBUG*/      if( !(Math.abs(nj-uj) <= 1e-8) )
  //*DEBUG*/        throw new Error('Assertion failed.');
  //*DEBUG*/    }
  //*DEBUG*/  }

  /* Used to copy either AB to AB0 or AB0 to AB.
   */


  var copy = function copy(AB_src, AB_dst) {
    for (var j = 0; j < N; j++) {
      for (var _i5 = 0; _i5 < k; _i5++) {
        AB_dst[_i5 + j * M] = AB_src[_i5 + j * M];
      }
    }
  };
  /* Used to move column `p` of R to column `k` >= `p` using cyclic permutation.
   * This method is used to prepare a "strong" column swap. After the cycle
   * followed by a retriangulation, we can downdate AB (and sometimes AB0) such
   * the the column swap does no longer affect AB (and AB0).
   * 
   * The column swap in R looks roughly as follows:
   * 
   *      ┏                ┓       ┏                ┓
   *      ┃.   p x   x k   ┃       ┃.   x   x k p   ┃
   *      ┃ .  . .   . k   ┃       ┃ .  .   . k .   ┃
   *      ┃  . . .   . k   ┃       ┃  . .   . k .   ┃
   *      ┃   .. .   . k   ┃(cycle)┃   ..   . k .   ┃
   *   R: ┃    p x   x k   ┃  ==>  ┃    x   x k p   ┃
   *      ┃      x   x k   ┃       ┃    x   x k 0   ┃
   *      ┃        ⋱ x k   ┃       ┃      ⋱ x k .   ┃
   *      ┃          x k   ┃       ┃        x k .   ┃
   *      ┃            k   ┃       ┃          k .   ┃
   *      ┃              ⋱ ┃       ┃              ⋱ ┃
   *      ┗                ┛       ┗                ┛
   * 
   * The advantage of this cyclic permutation is that only (k-p) Givens rotations
   * are required to retriangulate R.
   *
   *           ┏              ┓       ┏              ┓
   *           ┃.             ┃       ┃.             ┃
   *           ┃ .  .       . ┃       ┃ .  .       . ┃
   *           ┃  . .       . ┃       ┃  . .       . ┃
   *           ┃   ..       . ┃(cycle)┃   ..       . ┃
   *   inv(A): ┃    p p … p p ┃  ==>  ┃    0 x   x x ┃
   *           ┃      x   x x ┃       ┃        ⋱ x x ┃
   *           ┃        ⋱ x x ┃       ┃          x x ┃
   *           ┃          x x ┃       ┃            k ┃
   *           ┃            k ┃       ┃    p p … p p ┃
   *           ┃              ┃       ┃              ┃
   *           ┗              ┛       ┗              ┛
   * 
   * Since the cycle is a column swap inside of triangulate region A of R, we
   * have to perform the same cycle as ROW permutation in inv(A).
   * 
   * 
   * 
   * Keep in mind that while retriangulating R, the same givens rotations in
   * the same order must be applied to the columns of inv(A) after which
   * inv(A) is triangular again as well.
   */


  var cycle = function cycle(AB, p, k) {
    for (var j = p; j < N; j++) {
      var AB_pj = AB[p + j * M];

      if (j < k) {
        for (var _i6 = p; _i6 < j; _i6++) {
          AB[_i6 + j * M] = AB[_i6 + 1 + j * M];
        }

        AB[j + j * M] = 0;
      } else for (var _i7 = p; _i7 < k; _i7++) {
        AB[_i7 + j * M] = AB[_i7 + 1 + j * M];
      }

      AB[k + j * M] = AB_pj;
    }
  };
  /* Computes the column norms of C
   */


  var update_col_norms = function update_col_norms() {
    norm.fill(0.0);

    for (var _i8 = k; _i8 < M; _i8++) {
      (0, _rrqr._norm_update)(norm, R, R_off + N * _i8, k);
    }
  };
  /* Computes Frobenius norm of C.
   */


  var norm_C = function norm_C() {
    //*DEBUG*/    check_col_norms();
    NORM.reset();

    for (var j = k; j < N; j++) {
      NORM.include((0, _rrqr._norm)(norm, j));
    }

    return NORM.result;
  };
  /* Adjusts the binary search range and moves `k` to the
   * middle of that new range. If `increase` is true the
   * binary search range is moved to the right of the
   * current `k`. 
   */


  var adjust_k = function adjust_k(increase) {
    //*DEBUG*/    if( 'boolean' !== typeof increase )
    //*DEBUG*/      throw new Error('Assertion failed.');
    //*DEBUG*/    if(!(k0 <= k     )) throw new Error('Assertion failed.');
    //*DEBUG*/    if(!(      k <= K)) throw new Error('Assertion failed.');
    //*DEBUG*/    check_col_norms();
    if (increase) {
      //*DEBUG*/      if(!(k < K)) throw new Error('Assertion failed.');
      piv_elim();
      update(AB, k++);
      copy(AB, AB0);
      k0 = k;
    } else {
      //*DEBUG*/      if(!(k0 < k)) throw new Error('Assertion failed.');
      //*DEBUG*/      if(!(K == k)) throw new Error('Assertion failed.');
      copy(AB0, AB);
      k = k0;
      update_col_norms();
    }

    var mid = k0 + K >>> 1;

    while (k < mid) {
      //*DEBUG*/      check_col_norms();
      if (norm_C() <= ztol) {
        // we have found a new upper bound for the rank
        K = k;

        if (k0 < k) {
          // if we can go back let's always go back binary search style
          copy(AB0, AB);
          k = k0;
          update_col_norms();
          mid = k0 + K >>> 1;
          increase = false;
          continue;
        } // if we can't go back we have found the correct rank, let's return


        break;
      }

      if (increase) piv_elim();
      update(AB, k++);
      if (!increase) update_col_norms(); //*DEBUG*/      check_col_norms();
    }
  }; //*DEBUG*/  const logdet_A = () =>
  //*DEBUG*/  {
  //*DEBUG*/    // CHECK TRIANGULARITY OF R AND AB
  //*DEBUG*/    for( let j=0;   j < k; j++)
  //*DEBUG*/    for( let i=j; ++i < M;    )
  //*DEBUG*/      if( !(R [R_off + N*i+j  ] === 0) ) throw new Error(`Assertion failed.`);
  //*DEBUG*/
  //*DEBUG*/    let logdet_A = 0;
  //*DEBUG*/    for( let i=0; i < k; i++ )
  //*DEBUG*/      logdet_A += Math.log2(Math.abs(R[R_off + N*i+i]));
  //*DEBUG*/    return logdet_A;
  //*DEBUG*/  };
  // FOR EACH MATRIX IN THE BATCH


  for (var r_off = 0; Q_off < Q.length; Q_off += M * M, R_off += M * N, P_off += N, r_off += 1) {
    // SCALE R // <- FIXME: rethink scaling regarding underflow and overflow
    NORM.reset();

    for (var _i9 = M * N; _i9-- > 0;) {
      NORM.include(R[R_off + _i9]);
    }

    var SCALE = NORM.max === 0 ? 1 : NORM.result;
    if (!isFinite(SCALE)) throw new Error('Assertion failed: ' + SCALE);
    if (!(0 < SCALE)) throw new Error('Assertion failed: ' + SCALE);
    if (1 !== SCALE) for (var _i10 = M * N; _i10-- > 0;) {
      R[R_off + _i10] /= SCALE;
    } // INIT P

    for (var _i11 = 0; _i11 < N; _i11++) {
      P[P_off + _i11] = _i11;
    } // INIT Q (TO IDENTITY)


    for (var _i12 = 0; _i12 < M; _i12++) {
      Q[Q_off + M * _i12 + _i12] = 1;
    } // INIT AB


    AB0.fill(0.0);
    AB.fill(0.0); // RETRIEVE TOLERANCES

    dtol = DTOL(r_off), ztol = ZTOL(r_off);
    if (!isFinite(dtol)) throw new Error("Assertion failed. Invalid dtol: ".concat(dtol, "."));
    if (!isFinite(ztol)) throw new Error("Assertion failed. Invalid ztol: ".concat(ztol, "."));
    if (!(1 <= dtol)) throw new Error("Assertion failed. Invalid dtol: ".concat(dtol, "."));
    if (!(0 <= ztol)) throw new Error("Assertion failed. Invalid ztol: ".concat(ztol, ".")); // INIT BINARY SEARCH BOUNDS

    k0 = k = 0;
    K = L;
    Q_dirty = false;
    update_col_norms();

    loop: for (;;) {
      //*DEBUG*/      check(AB, k );
      //*DEBUG*/      check(AB0,k0);
      //*DEBUG*/      if(!(k0 <= k)) throw new Error('Assertion failed.');
      //*DEBUG*/      if(!(   k<=N)) throw new Error('Assertion failed.');
      if (norm_C() <= ztol) {
        // WE HAVE FOUND A NEW UPPER BOUND FOR THE RANK
        K = k;

        if (k0 < k) {
          // IF WE CAN GO BACK, LET'S GO BACK BINARY SEARCH STYLE
          adjust_k(
          /*increase=*/
          false); //*DEBUG*/          check(AB, k );
          //*DEBUG*/          check(AB0,k0);
        } else if (k === N) break loop; // <- no more column swaps possible -> stop
        //*DEBUG*/        else if( k !== k0) throw new Error('Assertion failed.');
        // AT THIS POINT WE KNOW THE RANK BUT
        // LET'S STRONG SWAP AS MUCH AS POSSIBLE

      } // SEARCH BEST COLUMN SWAPS


      var p = -1,
          q = -1,
          F = -Infinity; // COMPUTE ROW NORMS OF inv(A)

      for (var j = 0; j < k; j++) {
        _row_norm_update(norm, AB, j * M, j + 1);
      } // FIND BEST "STRONG" COLUMN SWAP


      for (var _i13 = 0; _i13 < k; _i13++) {
        var _r = (0, _rrqr._norm)(norm, _i13);

        for (var _j4 = k; _j4 < N; _j4++) {
          var c = (0, _rrqr._norm)(norm, _j4);
          var f = Math.hypot(AB[_i13 + _j4 * M], _r * c);

          if (F < f) {
            F = f;
            p = _i13;
            q = _j4;
          }
        }
      } // IF NO GOOD COLUMN SWAP IS FOUND ANYMORE


      if (!(F > dtol)) {
        if (k0 >= K) {
          // WE HAVE FOUND THE EXACT RANK
          //*DEBUG*/          if( k0 !== K ) throw new Error('Assertion failed.');
          //*DEBUG*/          if( k0 !== k ) throw new Error('Assertion failed.');
          break loop;
        } else {
          // WE HAVE TO GO FURTHER (the current k is less than the rank)
          //*DEBUG*/          if(!(k < K))
          //*DEBUG*/            throw new Error('Assertion failed.');
          adjust_k(
          /*increase=*/
          true);
          continue loop;
        }
      } //*DEBUG*/      if( !(0 < k) )
      //*DEBUG*/        throw new Error('Assertion failed.');


      Q_dirty = true; //*DEBUG*/      const predict_logdet_A = logdet_A() + Math.log2(F);
      // IF inv(A0) IS AFFECTED BY COLUMN SWAP

      if (p < k0) {
        --k0; // MOVE COLUMN p TO k0 (VIA CYCLIC PERMUTATION) TODO use triangulary property to reduces Ops
        // CYCLE COLUMNS of R

        for (var _i14 = 0; _i14 <= k0; _i14++) {
          var R_ip = R[R_off + N * _i14 + p];

          for (var _j5 = Math.max(p, _i14 - 1); _j5 < k0; _j5++) {
            R[R_off + N * _i14 + _j5] = R[R_off + N * _i14 + (_j5 + 1)];
          }

          R[R_off + N * _i14 + k0] = R_ip;
        } // CYCLE ROWS OF inv(A)


        cycle(AB, p, k0);
        cycle(AB0, p, k0); // CYCLE P P

        {
          var P_p = P[P_off + p];

          for (var _j6 = p; _j6 < k0; _j6++) {
            P[P_off + _j6] = P[P_off + _j6 + 1];
          }

          P[P_off + k0] = P_p;
        } // RETRIANGULATE USING GIVENS ROTATIONS
        // (since cyclic permutation turned R from triangular to Hessenberg)

        for (var _i15 = p; _i15 < k0; _i15++) {
          var ii = R_off + N * _i15 + _i15,
              ji = R_off + N * (_i15 + 1) + _i15,
              R_ji = R[ji];

          if (0 !== R_ji) {
            var R_ii = R[ii],
                _giv_rot_qr4 = (0, _giv_rot._giv_rot_qr)(R_ii, R_ji),
                _giv_rot_qr5 = (0, _slicedToArray2["default"])(_giv_rot_qr4, 3),
                _c = _giv_rot_qr5[0],
                s = _giv_rot_qr5[1],
                _norm3 = _giv_rot_qr5[2];

            R[ji] = 0;

            if (s !== 0) {
              R[ii] = _norm3;
              (0, _giv_rot._giv_rot_rows)(R, N - 1 - _i15, ii + 1, ji + 1, _c, s);
              (0, _giv_rot._giv_rot_rows)(Q, M, Q_off + M * _i15, Q_off + M * (_i15 + 1), _c, s);
              (0, _giv_rot._giv_rot_rows)(AB, _i15 + 1, M * _i15, M * (_i15 + 1), _c, s);
              AB[k0 + (_i15 + 1) * M] = -s * AB[k0 + _i15 * M] + _c * AB[k0 + (_i15 + 1) * M];
              (0, _giv_rot._giv_rot_rows)(AB0, _i15 + 1, M * _i15, M * (_i15 + 1), _c, s);
              AB0[k0 + (_i15 + 1) * M] = -s * AB0[k0 + _i15 * M] + _c * AB0[k0 + (_i15 + 1) * M];
            }
          }

          AB[k0 + _i15 * M] = 0;
          AB0[k0 + _i15 * M] = 0;
        }

        downdate(AB0, k0);
        p = k0++;
      }

      --k; // <- go back one step since the swappend column needs retriangulation
      // MOVE COLUMN p TO k (VIA CYCLIC PERMUTATION) TODO use triangulary property to reduces Ops
      // CYCLE COLUMNS of R

      for (var _i16 = 0; _i16 <= k; _i16++) {
        var _R_ip = R[R_off + N * _i16 + p];

        for (var _j7 = Math.max(p, _i16 - 1); _j7 < k; _j7++) {
          R[R_off + N * _i16 + _j7] = R[R_off + N * _i16 + (_j7 + 1)];
        }

        R[R_off + N * _i16 + k] = _R_ip;
      } // CYCLE ROWS OF inv(A)


      cycle(AB, p, k); // CYCLE COLS OF inv(A0)

      for (var _i17 = 0; _i17 < k0; _i17++) {
        norm[_i17] = AB0[_i17 + p * M];
      }

      for (var _j8 = p; _j8 < k; _j8++) {
        for (var _i18 = 0; _i18 < k0; _i18++) {
          AB0[_i18 + _j8 * M] = AB0[_i18 + (_j8 + 1) * M];
        }
      }

      for (var _i19 = 0; _i19 < k0; _i19++) {
        AB0[_i19 + k * M] = norm[_i19];
      } // CYCLE P P


      {
        var _P_p = P[P_off + p];

        for (var _j9 = p; _j9 < k; _j9++) {
          P[P_off + _j9] = P[P_off + _j9 + 1];
        }

        P[P_off + k] = _P_p;
      } // RETRIANGULATE USING GIVENS ROTATIONS
      // (since cyclic permutation turned R from triangular to Hessenberg)

      for (var _i20 = p; _i20 < k; _i20++) {
        var _ii = R_off + N * _i20 + _i20,
            _ji = R_off + N * (_i20 + 1) + _i20,
            _R_ji = R[_ji];

        if (0 !== _R_ji) {
          var _R_ii = R[_ii],
              _giv_rot_qr6 = (0, _giv_rot._giv_rot_qr)(_R_ii, _R_ji),
              _giv_rot_qr7 = (0, _slicedToArray2["default"])(_giv_rot_qr6, 3),
              _c2 = _giv_rot_qr7[0],
              _s = _giv_rot_qr7[1],
              _norm4 = _giv_rot_qr7[2];

          R[_ji] = 0;

          if (_s !== 0) {
            R[_ii] = _norm4;
            (0, _giv_rot._giv_rot_rows)(R, N - 1 - _i20, _ii + 1, _ji + 1, _c2, _s);
            (0, _giv_rot._giv_rot_rows)(Q, M, Q_off + M * _i20, Q_off + M * (_i20 + 1), _c2, _s);
            (0, _giv_rot._giv_rot_rows)(AB, _i20 + 1, M * _i20, M * (_i20 + 1), _c2, _s); // rotate last row of inv(A)

            AB[k + (_i20 + 1) * M] = -_s * AB[k + _i20 * M] + _c2 * AB[k + (_i20 + 1) * M];
          }
        }

        AB[k + _i20 * M] = 0; // <- aside from rounding errors, inv(A) is going to be triangular again
      }

      downdate(AB, k);
      swap_elim(q);

      if (p < k0) {
        //*DEBUG*/        if(!(p === k0-1)) throw new Error('Assertion failed.');
        update(AB0, p);
      }

      update(AB, k++); //*DEBUG*/      check(AB, k );
      //*DEBUG*/      check(AB0,k0);
      //*DEBUG*/      // CHECK DET PREDICTION
      //*DEBUG*/      if( !(Math.abs(logdet_A() - predict_logdet_A) <= 1e-6) )
      //*DEBUG*/        throw new Error(`${logdet_A()} != ${predict_logdet_A}`); 
    } // UNSCALE R


    if (1 !== SCALE) for (var _i21 = M * N; _i21-- > 0;) {
      R[R_off + _i21] *= SCALE;
    } // WRITE RANK

    r[r_off] = k; // TRANSPOSE Q (was transposed for cache alignment)

    (0, _transpose_inplace2._transpose_inplace)(M, Q, Q_off);
  }

  return [new _nd_array.NDArray(Q_shape, Q), new _nd_array.NDArray(R_shape, R), new _nd_array.NDArray(P_shape, P), new _nd_array.NDArray(r_shape, r)];
}