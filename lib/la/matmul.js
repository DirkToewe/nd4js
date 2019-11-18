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
exports.matmul2 = matmul2;
exports.matmul = matmul;

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _regenerator = _interopRequireDefault(require("@babel/runtime/regenerator"));

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _nd_array = require("../nd_array");

var _dt = require("../dt");

var _math = _interopRequireDefault(require("../math"));

var _ref =
/*#__PURE__*/
_regenerator["default"].mark(function _callee() {
  var mk_matmul2;
  return _regenerator["default"].wrap(function _callee$(_context) {
    while (1) {
      switch (_context.prev = _context.next) {
        case 0:
          mk_matmul2 = function mk_matmul2(addProduct) {
            return new Function('A_shape', 'A', 'B_shape', 'B', 'C_shape', 'C', "\n      const I = A_shape[A_shape.length-2],\n            K = A_shape[A_shape.length-1],\n            J = B_shape[B_shape.length-1]; \n      let\n        a = 0, aStride = 1,\n        b = 0, bStride = 1,\n        c = 0;\n\n      const loops = new Int32Array(C_shape.length-2);\n      for( let d=0; d >= 0; )\n        if( d === C_shape.length-2 ) {\n          aStride = I*K;\n          bStride = K*J;\n          for( const cEnd = c + I*J; c < cEnd; c += J, b -= bStride ) {\n          for( const aEnd = a + K  ; a < aEnd; c -= J, a++ ) {\n          for( const bEnd = b +   J; b < bEnd; c += 1, b++ ) {\n            ".concat(addProduct, "\n          }}}\n          b += bStride;\n          d -= 1;\n        }\n        else\n        {\n          if( loops[d]++ > 0 ) {\n            if( loops[d] > C_shape[d] ) {\n              aStride *= A_shape[ d - C_shape.length + A_shape.length ] || 1;\n              bStride *= B_shape[ d - C_shape.length + B_shape.length ] || 1;\n              loops[d--] = 0;\n              continue;\n            }\n            if( ! (A_shape[ d - C_shape.length + A_shape.length ] > 1) ) a -= aStride;\n            if( ! (B_shape[ d - C_shape.length + B_shape.length ] > 1) ) b -= bStride;\n          }\n          ++d;\n        }\n    "));
          };

          _context.next = 3;
          return mk_matmul2('C[c] += A[a]*B[b];');

        case 3:
          _context.next = 5;
          return mk_matmul2("\n    C[2*c+0] += A[a] * B[2*b+0];\n    C[2*c+1] += A[a] * B[2*b+1];\n  ");

        case 5:
          _context.next = 7;
          return mk_matmul2("\n    C[2*c+0] += B[b] * A[2*a+0];\n    C[2*c+1] += B[b] * A[2*a+1];\n  ");

        case 7:
          _context.next = 9;
          return mk_matmul2("\n    C[2*c+0] += B[2*b+0] * A[2*a+0]  -  B[2*b+1] * A[2*a+1];\n    C[2*c+1] += B[2*b+0] * A[2*a+1]  +  B[2*b+1] * A[2*a+0];\n  ");

        case 9:
          _context.next = 11;
          return mk_matmul2('C[c] = this.add(C[c], this.mul(A[a],B[b]));').bind(_math["default"]);

        case 11:
        case "end":
          return _context.stop();
      }
    }
  }, _callee);
})(),
    _ref2 = (0, _slicedToArray2["default"])(_ref, 5),
    matmul2_RR = _ref2[0],
    matmul2_RC = _ref2[1],
    matmul2_CR = _ref2[2],
    matmul2_CC = _ref2[3],
    matmul2_ELSE = _ref2[4];

function matmul2(a, b) {
  a = (0, _nd_array.asarray)(a);
  b = (0, _nd_array.asarray)(b);
  if (a.ndim < 2) throw new Error('A must be at least 2D.');
  if (b.ndim < 2) throw new Error('B must be at least 2D.');

  var _a$shape$slice = a.shape.slice(-2),
      _a$shape$slice2 = (0, _slicedToArray2["default"])(_a$shape$slice, 2),
      I = _a$shape$slice2[0],
      K = _a$shape$slice2[1],
      J = b.shape[b.ndim - 1];

  if (b.shape[b.ndim - 2] != K) throw new Error('The last dimension of A and the 2nd to last dimension of B do not match.');
  var ndim = Math.max(a.ndim, b.ndim),
      shape = Int32Array.from({
    length: ndim
  }, function () {
    return 1;
  });
  shape[ndim - 2] = I;
  shape[ndim - 1] = J; // FIND COMMON (BROADCASTED) SHAPE

  for (var _i = 0, _arr = [a, b]; _i < _arr.length; _i++) {
    var arr = _arr[_i];

    for (var i = ndim - 2, j = arr.ndim - 2; i-- > 0 && j-- > 0;) {
      if (1 === shape[i]) shape[i] = arr.shape[j];else if (shape[i] != arr.shape[j] && arr.shape[j] != 1) throw new Error('Shapes are not broadcast-compatible.');
    }
  } // GENERATE RESULT DATA


  var DTypeArray = _dt.ARRAY_TYPES[(0, _dt.super_dtype)(a.dtype, b.dtype)]; // INTEGER TYPE ID (0: scalar, 1: complex, 2: object)


  var c = new DTypeArray(shape.reduce(function (m, n) {
    return m * n;
  }, 1));
  if (c instanceof Array) c.fill(0);
  var pairing = 0;

  for (var _i2 = 0, _arr2 = [a, b]; _i2 < _arr2.length; _i2++) {
    var ab = _arr2[_i2];
    pairing *= 3;

    switch (ab.dtype) {
      default:
        pairing += 1;

      case 'complex128':
        pairing += 1;

      case 'int32':
      case 'float32':
      case 'float64':
    }
  }

  switch (pairing) {
    case 0:
      matmul2_RR(a.shape, a.data, b.shape, b.data, shape, c);
      break;

    case 1:
      matmul2_RC(a.shape, a.data, b.shape, b.data._array, shape, c._array);
      break;

    case 3:
      matmul2_CR(a.shape, a.data._array, b.shape, b.data, shape, c._array);
      break;

    case 4:
      matmul2_CC(a.shape, a.data._array, b.shape, b.data._array, shape, c._array);
      break;

    default:
      matmul2_ELSE(a.shape, a.data, b.shape, b.data, shape, c);
      break;
  }

  return new _nd_array.NDArray(shape, c);
}

function matmul() {
  for (var _len = arguments.length, matrices = new Array(_len), _key = 0; _key < _len; _key++) {
    matrices[_key] = arguments[_key];
  }

  matrices = matrices.map(_nd_array.asarray);
  if (matrices.length == 1) return matrices[0];
  if (matrices.length == 2) return matmul2.apply(void 0, (0, _toConsumableArray2["default"])(matrices));
  /** Returns the number of floating point operations necessary to
   *  matrix multiply two arrays of the given shapes.
   */

  function nOps(shapeA, shapeB) {
    var _shapeA$slice = shapeA.slice(-2),
        _shapeA$slice2 = (0, _slicedToArray2["default"])(_shapeA$slice, 2),
        I = _shapeA$slice2[0],
        K = _shapeA$slice2[1],
        J = shapeB[shapeB.length - 1];

    if (shapeB[shapeB.length - 2] != K) throw new Error('Shape mismatch.');
    var ndim = Math.max(shapeA.length, shapeB.length),
        shape = Int32Array.from({
      length: ndim
    }, function () {
      return 1;
    });
    shape[ndim - 2] = I;
    shape[ndim - 1] = J; // FIND COMMON (BROADCASTED) SHAPE

    for (var _i3 = 0, _arr3 = [shapeA, shapeB]; _i3 < _arr3.length; _i3++) {
      var shp = _arr3[_i3];

      for (var i = ndim - 2, j = shp.length - 2; i-- > 0 && j-- > 0;) {
        if (1 === shape[i]) shape[i] = shp[j];else if (shape[i] != shp[j] && shp[j] != 1) throw new Error('Shapes are not broadcast-compatible.');
      }
    }

    return [shape.reduce(function (a, b) {
      return a * b;
    }, 1) * K, shape];
  }

  ; // https://en.wikipedia.org/wiki/Matrix_chain_multiplication
  // https://www.geeksforgeeks.org/matrix-chain-multiplication-dp-8/
  // opMap[i][j] (j >= i) caches the optimal number of FLOPs required to multiply matrices[i] up to (including) matrices[j]

  var opMap = Array.from({
    length: matrices.length
  }, function () {
    return [];
  }); // initialized opMap

  for (var i = 0; i < matrices.length; i++) {
    opMap[i][i] = [0, matrices[i].shape];
  } // compute remaining opMap


  for (var len = 2; len <= matrices.length; len++) {
    for (var _i4 = 0; _i4 <= matrices.length - len; _i4++) {
      var minFlops = Infinity,
          minShape = void 0;

      for (var j = 1; j < len; j++) {
        var _opMap$_i = (0, _slicedToArray2["default"])(opMap[_i4][_i4 + j - 1], 2),
            lFlops = _opMap$_i[0],
            lShape = _opMap$_i[1];

        var _opMap = (0, _slicedToArray2["default"])(opMap[_i4 + j][_i4 + len - 1], 2),
            rFlops = _opMap[0],
            rShape = _opMap[1];

        var _nOps = nOps(lShape, rShape),
            _nOps2 = (0, _slicedToArray2["default"])(_nOps, 2),
            flops = _nOps2[0],
            shape = _nOps2[1];

        flops += lFlops + rFlops;

        if (flops < minFlops) {
          minFlops = flops;
          minShape = shape; // <- the shape should always be the same so this is not strictly necessary
        }
      }

      if (minShape === undefined) throw new Error('Integer overflow (too many FLOPs).');
      opMap[_i4][_i4 + len - 1] = [minFlops, minShape];
    }
  } // compute the result using the minimal number of FLOPs using opMap


  function product(from, to) {
    if (from == to) return matrices[from];
    var minFlops = Infinity,
        minIdx;

    for (var _i5 = from; _i5 < to; _i5++) {
      var _opMap$from$_i = (0, _slicedToArray2["default"])(opMap[from][_i5], 2),
          _lFlops = _opMap$from$_i[0],
          _lShape = _opMap$from$_i[1];

      var _opMap$to = (0, _slicedToArray2["default"])(opMap[_i5 + 1][to], 2),
          _rFlops = _opMap$to[0],
          _rShape = _opMap$to[1];

      var _nOps3 = nOps(_lShape, _rShape),
          _nOps4 = (0, _slicedToArray2["default"])(_nOps3, 2),
          _flops = _nOps4[0],
          _ = _nOps4[1];

      _flops += _lFlops + _rFlops;

      if (_flops < minFlops) {
        minFlops = _flops;
        minIdx = _i5;
      }
    }

    return matmul2(product(from, minIdx), product(minIdx + 1, to));
  }

  return product(0, matrices.length - 1);
}