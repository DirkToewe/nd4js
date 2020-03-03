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
exports.svd_jac_classic = svd_jac_classic;

var _regenerator = _interopRequireDefault(require("@babel/runtime/regenerator"));

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _dt = require("../dt");

var _nd_array = require("../nd_array");

var _matmul = require("./matmul");

var _qr = require("./qr");

var _transpose_inplace = require("./transpose_inplace");

var _svd_jac_utils = require("./_svd_jac_utils");

var _giv_rot = require("./_giv_rot");

function svd_jac_classic(A) {
  A = (0, _nd_array.asarray)(A);
  if (A.dtype.startsWith('complex')) throw new Error('svd_jac_1sided(A): A.dtype must be float.');
  var shape = A.shape,
      N = shape[shape.length - 2]; // QR DECOMPOSE RECTANGULAR MATRICES AND WORK ON (QUADRATIC) R

  {
    var M = shape[shape.length - 1]; // if A is not square use QR Decomposition

    if (N > M) {
      var _qr_decomp = (0, _qr.qr_decomp)(A),
          _qr_decomp2 = (0, _slicedToArray2["default"])(_qr_decomp, 2),
          Q = _qr_decomp2[0],
          R = _qr_decomp2[1],
          _svd_jac_classic = svd_jac_classic(R),
          _svd_jac_classic2 = (0, _slicedToArray2["default"])(_svd_jac_classic, 3),
          _U = _svd_jac_classic2[0],
          _sv = _svd_jac_classic2[1],
          _V = _svd_jac_classic2[2];

      return [(0, _matmul.matmul2)(Q, _U), _sv, _V];
    }

    if (N < M) {
      var _qr_decomp3 = (0, _qr.qr_decomp)(A.T),
          _qr_decomp4 = (0, _slicedToArray2["default"])(_qr_decomp3, 2),
          _Q = _qr_decomp4[0],
          _R = _qr_decomp4[1],
          _svd_jac_classic3 = svd_jac_classic(_R),
          _svd_jac_classic4 = (0, _slicedToArray2["default"])(_svd_jac_classic3, 3),
          _U2 = _svd_jac_classic4[0],
          _sv2 = _svd_jac_classic4[1],
          _V2 = _svd_jac_classic4[2];

      (0, _transpose_inplace.transpose_inplace)(_V2);
      return [_V2, _sv2, (0, _matmul.matmul2)(_Q, _U2).T];
    }
  } // ALLOCATE RESULT DATA

  var DType = A.dtype === 'float32' ? 'float32' : 'float64',
      DTypeArray = _dt.ARRAY_TYPES[DType],
      TOL = (0, _dt.eps)(DType) * N,
      U = DTypeArray.from(A.data);
  A = undefined; // <- potentially allow GC

  var S = new DTypeArray(N * N),
      // <- tempory storage for decomposition
  V = new DTypeArray(U.length),
      sv = new DTypeArray(U.length / N),
      diag = new DTypeArray(N),
      ord = Int32Array.from({
    length: N
  }, function (_, i) {
    return i;
  });
  if (1 > N) throw new Error('Assertion failed.');

  if (1 === N) {
    for (var i = U.length; i-- > 0;) {
      if (U[i] < +0.0) {
        U[i] *= -1.0;
        sv[i] = -1.0;
      } else sv[i] = +1.0;
    }

    return [new _nd_array.NDArray(shape, sv), new _nd_array.NDArray(shape.slice(0, -1), U), new _nd_array.NDArray(shape, V.fill(1))];
  } //
  // BUILD TRIANGLE TREE
  //
  // the size of each level of the tree


  var treeSizes = Int32Array.from( /*#__PURE__*/_regenerator["default"].mark(function _callee() {
    var n;
    return _regenerator["default"].wrap(function _callee$(_context) {
      while (1) {
        switch (_context.prev = _context.next) {
          case 0:
            n = N;

          case 1:
            if (!(n > 1)) {
              _context.next = 7;
              break;
            }

            n = n + 1 >>> 1;
            _context.next = 5;
            return n;

          case 5:
            _context.next = 1;
            break;

          case 7:
          case "end":
            return _context.stop();
        }
      }
    }, _callee);
  })());
  var treeData = new DTypeArray(treeSizes.reduce(function (len, N) {
    return len + (N * N + N >>> 1);
  }, 0)); //*DEBUG*/  const piv = (i,j) => {
  //*DEBUG*/    const S_ij = S[N*i+j],
  //*DEBUG*/          S_ji = S[N*j+i];
  //*DEBUG*/    return S_ij*S_ij + S_ji*S_ji;
  //*DEBUG*/  };

  var find_pivot = function find_pivot() {
    var k = 0,
        l = 0,
        i = 0,
        j = 0;

    var _loop = function _loop(_off, h) {
      var val = function val(k, l) {
        k += i;
        l += j;
        off = _off;
        return treeData[_off + (k * k + k >>> 1) + l];
      };

      var n = treeSizes[h];
      _off -= n * n + n >>> 1;
      var max = val(0, 0);
      k = i;
      l = j;

      if (i + 1 < n) {
        var v10 = val(1, 0);

        if (v10 > max) {
          max = v10;
          k = i + 1;
          l = j;
        }

        var v11 = val(1, 1);

        if (v11 > max) {
          max = v11;
          k = i + 1;
          l = j + 1;
        }
      }

      if (i != j) {
        var v01 = val(0, 1);

        if (v01 > max) {
          max = v01;
          k = i;
          l = j + 1;
        }
      }

      i = 2 * k;
      j = 2 * l;
      off = _off;
    };

    for (var off = treeData.length, h = treeSizes.length; h-- > 0;) {
      _loop(off, h);
    }

    var max = -Infinity;

    var hyp = function hyp(s, t) {
      s += i;
      t += j;
      var S_st = S[N * s + t],
          S_ts = S[N * t + s];
      return S_st * S_st + S_ts * S_ts;
    };

    if (i + 1 < N) {
      var h10 = hyp(1, 0);

      if (h10 > max) {
        max = h10;
        k = i + 1;
        l = j;
      }
    }

    if (j + 1 < N) {
      var h01 = hyp(0, 1);

      if (h01 > max) {
        max = h01;
        k = i;
        l = j + 1;
      }
    }

    if (i != j) {
      var h00 = hyp(0, 0);

      if (h00 > max) {
        max = h00;
        k = i;
        l = j;
      }

      if (i + 1 < N) {
        var h11 = hyp(1, 1);

        if (h11 > max) {
          k = i + 1;
          l = j + 1;
        }
      }
    }

    return [l, k];
  };
  /** updates the specified row in the triangle tree.
   */


  var update_row = function update_row(row) {
    row = row >>> 1 << 1; // <- round down to pow2
    // build bottom tree level

    for (var _i = row; _i < row + 2 && _i < N; _i++) {
      var r = _i >>> 1,
          _off2 = r * r + r >>> 1;

      for (var j = 0; j < row + 1; j++) {
        var x = S[N * _i + j],
            y = S[N * j + _i],
            S_ij = _i === j ? -Infinity : x * x + y * y,
            k = _off2 + (j >>> 1);
        treeData[k] = _i % 2 || j % 2 ? Math.max(treeData[k], S_ij) : S_ij;
      }
    } // build remaining tree levels


    for (var _off3 = 0, h = 1; h < treeSizes.length; h++) {
      var _N = treeSizes[h - 1],
          OFF = _off3;
      _off3 += _N * _N + _N >>> 1;
      row = row >>> 2 << 1;

      for (var _R2 = row; _R2 < row + 2 && _R2 < _N; _R2++) {
        var ROFF = OFF + (_R2 * _R2 + _R2 >>> 1),
            _r = _R2 >>> 1,
            roff = _off3 + (_r * _r + _r >>> 1);

        for (var C = 0; C <= _R2; C++) {
          var _k = roff + (C >>> 1),
              T = treeData[ROFF + C];

          treeData[_k] = _R2 % 2 || C % 2 ? Math.max(treeData[_k], T) : T;
        }
      }
    }
  };
  /** updates the specified column in the triangle tree.
   */


  var update_col = function update_col(col) {
    col = col >>> 1 << 1; // <- round down to pow2
    // build bottom tree level

    var J = Math.min(col + 2, N);

    for (var _i2 = col; _i2 < N; _i2++) {
      var r = _i2 >>> 1,
          _off4 = r * r + r >>> 1;

      for (var j = col; j < J; j++) {
        var x = S[N * _i2 + j],
            y = S[N * j + _i2],
            S_ij = _i2 === j ? -Infinity : x * x + y * y,
            k = _off4 + (j >>> 1);
        treeData[k] = _i2 % 2 || j % 2 ? Math.max(treeData[k], S_ij) : S_ij;
      }
    } // build remaining tree levels


    for (var _off5 = 0, h = 1; h < treeSizes.length; h++) {
      col = col >>> 2 << 1;

      var _N2 = treeSizes[h - 1],
          _J = Math.min(col + 2, _N2),
          OFF = _off5;

      _off5 += _N2 * _N2 + _N2 >>> 1;

      for (var _R3 = col; _R3 < _N2; _R3++) {
        var ROFF = OFF + (_R3 * _R3 + _R3 >>> 1),
            _r2 = _R3 >>> 1,
            roff = _off5 + (_r2 * _r2 + _r2 >>> 1);

        for (var C = col; C < _J && C <= _R3; C++) {
          var _k2 = roff + (C >>> 1),
              T = treeData[ROFF + C];

          treeData[_k2] = _R3 % 2 || C % 2 ? Math.max(treeData[_k2], T) : T;
        }
      }
    }
  };

  for (var UV_off = 0, sv_off = 0; sv_off < sv.length; UV_off += N * N, sv_off += N) {
    // MOVE FROM U TO S
    for (var _i3 = 0; _i3 < N; _i3++) {
      for (var j = 0; j < N; j++) {
        S[N * _i3 + j] = U[UV_off + N * _i3 + j];
        U[UV_off + N * _i3 + j] = +(_i3 === j);
      }
    }

    ; // INIT V TO IDENTITY

    for (var _i4 = 0; _i4 < N; _i4++) {
      for (var _j = 0; _j < N; _j++) {
        V[UV_off + N * _i4 + _j] = +(_i4 === _j);
      }
    } // INIT TRIANGLE TREE


    for (var _i5 = 0; _i5 < N; _i5 += 2) {
      update_row(_i5);
    } //
    // (CLASSICAL) JACOBI SVD ITERATIONS
    //


    for (var s = -1, t = -1;;) {
      // FIND THE OFF DIAGONAL PAIR WITH THE LARGEST HYPOTHENUSE
      var _find_pivot = find_pivot(),
          _find_pivot2 = (0, _slicedToArray2["default"])(_find_pivot, 2),
          l = _find_pivot2[0],
          k = _find_pivot2[1];

      if (l >= k) throw new Error('Assertion failed.');
      if (s == k && t == l) break; // <- DRY PRINCIPLE

      s = k;
      t = l; //*DEBUG*/      // CHECK THAT THIS IS TRULY THE MAXIMUM (TODO: COMMENT OUT)
      //*DEBUG*/      for( let i=1; i < N; i++ )
      //*DEBUG*/      for( let j=0; j < i; j++ )
      //*DEBUG*/        if( ! (piv(i,j) <= piv(k,l)) ) throw new Error(`Assertion failed: ${i}, ${j}.`);

      var S_kk = S[N * k + k],
          S_kl = S[N * k + l],
          S_lk = S[N * l + k],
          S_ll = S[N * l + l]; // stopping criterion inspiredy by:
      //  "Jacobi's Method is More Accurate than QR"
      //   by James Demmel
      //   SIAM J. Matrix Anal. Appl, vol. 13, pp. 1204-1245, 1992

      if (!(Math.max(Math.abs(S_kl), Math.abs(S_lk)) > Math.sqrt(Math.abs(S_kk * S_ll)) * TOL)) break;

      var _svd_jac_angles2 = (0, _svd_jac_utils._svd_jac_angles)(S_ll, S_lk, S_kl, S_kk),
          _svd_jac_angles3 = (0, _slicedToArray2["default"])(_svd_jac_angles2, 4),
          cα = _svd_jac_angles3[0],
          sα = _svd_jac_angles3[1],
          cβ = _svd_jac_angles3[2],
          sβ = _svd_jac_angles3[3]; // ROTATE S


      (0, _giv_rot._giv_rot_rows)(S, N, N * l, N * k, cα, sα);
      (0, _giv_rot._giv_rot_cols)(S, N, l, k, cβ, sβ); // ENTRIES (k,l) AND (l,k) ARE REMAINDERS (CANCELLATION ERROR) FROM ELIMINATION => SHOULD BE SAFELY ZEROABLE

      S[N * k + l] = 0.0;
      S[N * l + k] = 0.0; // UPDATE TRIANGLE TREE ROWS AND COLUMNS

      update_row(k);
      update_row(l);
      update_col(k);
      update_col(l); // ROTATE U & V

      (0, _giv_rot._giv_rot_rows)(U, N, UV_off + N * l, UV_off + N * k, cα, sα);
      (0, _giv_rot._giv_rot_rows)(V, N, UV_off + N * l, UV_off + N * k, cβ, -sβ);
    }

    (0, _svd_jac_utils._svd_jac_post)(N, U, S, V, UV_off, sv, sv_off, ord);
  }

  return [new _nd_array.NDArray(shape, U), new _nd_array.NDArray(shape.slice(0, -1), sv), new _nd_array.NDArray(shape, V)];
}