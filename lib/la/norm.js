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
exports.norm = norm;
exports.FrobeniusNorm = void 0;

var _classCallCheck2 = _interopRequireDefault(require("@babel/runtime/helpers/classCallCheck"));

var _createClass2 = _interopRequireDefault(require("@babel/runtime/helpers/createClass"));

var _nd_array = require("../nd_array");

function _createForOfIteratorHelper(o) { if (typeof Symbol === "undefined" || o[Symbol.iterator] == null) { if (Array.isArray(o) || (o = _unsupportedIterableToArray(o))) { var i = 0; var F = function F() {}; return { s: F, n: function n() { if (i >= o.length) return { done: true }; return { done: false, value: o[i++] }; }, e: function e(_e) { throw _e; }, f: F }; } throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method."); } var it, normalCompletion = true, didErr = false, err; return { s: function s() { it = o[Symbol.iterator](); }, n: function n() { var step = it.next(); normalCompletion = step.done; return step; }, e: function e(_e2) { didErr = true; err = _e2; }, f: function f() { try { if (!normalCompletion && it["return"] != null) it["return"](); } finally { if (didErr) throw err; } } }; }

function _unsupportedIterableToArray(o, minLen) { if (!o) return; if (typeof o === "string") return _arrayLikeToArray(o, minLen); var n = Object.prototype.toString.call(o).slice(8, -1); if (n === "Object" && o.constructor) n = o.constructor.name; if (n === "Map" || n === "Set") return Array.from(o); if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n)) return _arrayLikeToArray(o, minLen); }

function _arrayLikeToArray(arr, len) { if (len == null || len > arr.length) len = arr.length; for (var i = 0, arr2 = new Array(len); i < len; i++) { arr2[i] = arr[i]; } return arr2; }

var FrobeniusNorm = /*#__PURE__*/function () {
  function FrobeniusNorm() {
    (0, _classCallCheck2["default"])(this, FrobeniusNorm);
    this.sum = 0.0;
    this.max = 0.0;
    Object.seal(this);
  }

  (0, _createClass2["default"])(FrobeniusNorm, [{
    key: "reset",
    value: function reset() {
      this.sum = this.max = 0;
    }
  }, {
    key: "include",
    value: function include(x) {
      x = Math.abs(x);

      if (x !== 0) {
        if (this.max < x) {
          var s = this.max / x;
          this.sum *= s * s;
          this.max = x;
        }

        x /= this.max;
        this.sum += x * x;
      }
    }
  }, {
    key: "resultIncl",
    value: function resultIncl(x) {
      x = Math.abs(x);
      var sum = this.sum,
          max = this.max;

      if (x !== 0) {
        if (max < x) {
          var s = max / x;
          sum *= s * s;
          max = x;
        }

        x /= max;
        sum += x * x;
      }

      return isFinite(max) ? Math.sqrt(sum) * max : max;
    }
  }, {
    key: "result",
    get: function get() {
      return isFinite(this.max) ? Math.sqrt(this.sum) * this.max : this.max;
    }
  }]);
  return FrobeniusNorm;
}();

exports.FrobeniusNorm = FrobeniusNorm;

function norm(A) {
  var ord = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 'fro';
  var axis = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : undefined;
  A = (0, _nd_array.asarray)(A);
  var norm;

  switch (ord) {
    case 'fro':
      norm = new FrobeniusNorm();
      break;

    default:
      throw new Error("norm(A,ord,axis): Unsupported ord: ".concat(ord, "."));
  }

  if (null == axis) {
    var _iterator = _createForOfIteratorHelper(A.data),
        _step;

    try {
      for (_iterator.s(); !(_step = _iterator.n()).done;) {
        var x = _step.value;
        norm.include(x);
      }
    } catch (err) {
      _iterator.e(err);
    } finally {
      _iterator.f();
    }

    return norm.result;
  }

  throw new Error('norm(A,ord,axis): axis argument not yet supported.');
} //    if( null == reducer )
//    {
//      if( null    ==  dtype      ) return this.data.reduce(axes)
//      if('string' === typeof axes) return this.data.reduce(dtype)
//      reducer = dtype; dtype = undefined
//    }
//    if( null == dtype ) dtype = 'object'
//
//    if( ! is_subdtype(this.dtype, dtype) )
//      throw new Error('New dtype must be a super-dtype.')
//
//    const oldNDim = this.ndim
//
//    if( 'number' === typeof axes ) axes = [axes]
//    if( axes instanceof NDArray ) {
//      if( ! is_subdtype(axes.dtype,'int32') ) throw new Error(`Invalid dtype ${axes.dtype} for axes.`) 
//      if( axes.ndim === 1 )
//        axes = axes.data
//      else
//        throw new Error('Only 1D nd.Array allowed for axes.')
//    }
//    axes = new Set( function*() {
//      for( let ax of axes ) {
//        if( 0 > ax )  ax += oldNDim
//        if( 0 > ax || ax >= oldNDim ) throw new Error('Reduction axis '+ax+' out of bounds.')
//        yield +ax
//      }
//    }() )
//
//    const
//      oldShape= this.shape,
//      newShape= oldShape.filter( (size,i) => ! axes.has(i) ),
//      newData = new ARRAY_TYPES[dtype]( newShape.reduce((a,b) => a*b, 1) ),
//      oldData = this.data
//    let
//      newIdx = 0,
//      oldIdx = 0
//
//    function fill(d, reduce)
//    {
//      if( oldShape.length === d )
//      {
//        if(reduce) newData[newIdx]=reducer(oldData[oldIdx++], newData[newIdx])
//        else       newData[newIdx]=        oldData[oldIdx++]
//        ++newIdx
//      }
//      else if( axes.has(d) ) {
//        const idx = newIdx
//        for( let j=oldShape[d]; j-- > 0; reduce=true ) // <- copy over very first value and only then start reduction
//        {
//          newIdx = idx
//          fill(d+1, reduce)
//        }
//      }
//      else for( let j=oldShape[d]; j-- > 0; )
//        fill(d+1, reduce)
//    }
//    fill(0,false)
//
//    if( newIdx !=   newData.length ) throw new Error([newIdx,   newData.length])
//    if( oldIdx != this.data.length ) throw new Error([oldIdx, this.data.length])
//
//    return new NDArray(newShape,newData)