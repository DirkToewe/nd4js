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
exports.array = array;
exports.asarray = asarray;
exports.NDArray = void 0;

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _classCallCheck2 = _interopRequireDefault(require("@babel/runtime/helpers/classCallCheck"));

var _possibleConstructorReturn2 = _interopRequireDefault(require("@babel/runtime/helpers/possibleConstructorReturn"));

var _createClass2 = _interopRequireDefault(require("@babel/runtime/helpers/createClass"));

var _inherits2 = _interopRequireDefault(require("@babel/runtime/helpers/inherits"));

var _wrapNativeSuper2 = _interopRequireDefault(require("@babel/runtime/helpers/wrapNativeSuper"));

var _regenerator = _interopRequireDefault(require("@babel/runtime/regenerator"));

var _typeof2 = _interopRequireDefault(require("@babel/runtime/helpers/typeof"));

var _zip_elems = require("./zip_elems");

var _dt = require("./dt");

function array(dtype, content) {
  var _marked = /*#__PURE__*/_regenerator["default"].mark(shape);

  if (null == content) {
    content = dtype;
    dtype = undefined;
  }

  if (('object' == (0, _typeof2["default"])(content) || 'function' == typeof content) && 'shape' in content && 'data' in content) {
    var _data = content.data;
    if (null != dtype && !(_data instanceof _dt.ARRAY_TYPES[dtype])) _data = _dt.ARRAY_TYPES[dtype].from(_data);else _data = _data.slice();
    return new NDArray(content.shape, _data);
  }

  function shape(content) {
    var len;
    return _regenerator["default"].wrap(function shape$(_context) {
      while (1) {
        switch (_context.prev = _context.next) {
          case 0:
            if (!('object' === (0, _typeof2["default"])(content) && 'length' in content)) {
              _context.next = 7;
              break;
            }

            len = content.length;

            if (!(typeof len !== 'number' || !(len % 1 === 0) || len === 0)) {
              _context.next = 4;
              break;
            }

            throw 'Illegal argument(s).';

          case 4:
            _context.next = 6;
            return len;

          case 6:
            return _context.delegateYield(shape(content[0]), "t0", 7);

          case 7:
          case "end":
            return _context.stop();
        }
      }
    }, _marked);
  }

  shape = Int32Array.from(shape(content));
  if (null == dtype) dtype = function dtype(d, content) {
    if (d === shape.length) return (0, _dt.dtypeof)(content);else {
      var j = shape[d] - 1,
          dt = dtype(d + 1, content[j]);

      for (; j >= 0; j--) {
        dt = (0, _dt.super_dtype)(dt, dtype(d + 1, content[j]));
      }

      return dt;
    }
  }(0, content);
  (0, _dt._check_dtype)(dtype);
  var data = new _dt.ARRAY_TYPES[dtype](shape.reduce(function (a, b) {
    return a * b;
  }, 1));
  var idx = 0;

  function fill(d, content) {
    if (d === shape.length) data[idx++] = content;else {
      if (content.length !== shape[d]) throw new Error('Nested array not axis-aligned.');

      for (var j = 0; j < shape[d]; j++) {
        fill(d + 1, content[j]);
      }
    }
  }

  fill(0, content);
  return new NDArray(shape, data);
}

function asarray(dtype, arrayLike) {
  if (null == arrayLike) {
    arrayLike = dtype;
    dtype = undefined;
  }

  if (dtype != null && dtype !== 'float' && !(dtype in _dt.ARRAY_TYPES)) throw new Error("asarray(dtype, arrayLike): Invalid dtype: '".concat(dtype, "'."));

  if (arrayLike instanceof NDArray) {
    if (dtype == null || dtype === arrayLike.dtype || dtype === 'float' && arrayLike.dtype.startsWith('float')) return arrayLike;
    return new NDArray(arrayLike.shape, _dt.ARRAY_TYPES[dtype].from(arrayLike.data));
  }

  if (dtype === 'float') dtype = 'float64';
  return array(dtype, arrayLike);
}

var NDArray = /*#__PURE__*/function (_Function) {
  (0, _inherits2["default"])(NDArray, _Function);
  (0, _createClass2["default"])(NDArray, null, [{
    key: "name",
    get: function get() {
      return 'nd.Array';
    } //
    // FACTORY METHODS
    //

  }]);

  function NDArray(shape, data) {
    var _this;

    (0, _classCallCheck2["default"])(this, NDArray);
    if (!(shape instanceof Int32Array)) throw new Error('Shape must be Int32Array.');
    if (shape.some(function (s) {
      return s < 1;
    })) throw new Error("Invalid shape: ".concat(shape, "."));
    if (data.length != shape.reduce(function (x, y) {
      return x * y;
    }, 1)) throw new Error("Shape [".concat(shape, "] does not match array length of ").concat(data.length, "."));

    var self = function self() {
      for (var _len = arguments.length, indices = new Array(_len), _key = 0; _key < _len; _key++) {
        indices[_key] = arguments[_key];
      }

      return self.data[self._flat_idx(indices)];
    };

    Object.setPrototypeOf(self, NDArray.prototype);
    Object.freeze(shape.buffer);
    self.shape = shape;
    self.data = data;
    Object.seal(self);
    return (0, _possibleConstructorReturn2["default"])(_this, self);
  } //
  // GENERAL
  //


  (0, _createClass2["default"])(NDArray, [{
    key: "set",
    //    TODO implements this in combination with a PROXY
    //    get length() {
    //      return this.shape[0];
    //    }
    value: function set(indices, value) {
      this.data[this._flat_idx(indices)] = value;
    }
  }, {
    key: "modify",
    value: function modify(indices, modifier) {
      var i = this._flat_idx(indices);

      this.data[i] = modifier.apply(void 0, [this.data[i]].concat((0, _toConsumableArray2["default"])(indices)));
    }
  }, {
    key: "_flat_idx",
    value: function _flat_idx(indices) {
      var shape = this.shape;
      if (indices.length != shape.length) throw new Error("Multi-index [".concat(indices, "] does not have expected length of ").concat(shape.length, "."));
      var flat_idx = 0,
          stride = 1;

      for (var i = shape.length; i-- > 0; stride *= shape[i]) {
        var idx = indices[i];
        if (idx % 1 != 0) throw new Error("Multi-index [".concat(indices, "] contains non-integer entries."));
        if (idx < 0) idx += shape[i];
        if (idx < 0 || idx >= shape[i]) throw new Error("Multi-index [".concat(indices, "] out of bounds [").concat(shape, "]."));
        flat_idx += idx * stride;
      }

      return flat_idx;
    }
  }, {
    key: Symbol["for"]('nodejs.util.inspect.custom'),
    value: function value(depth, options) {
      return this.toString();
    }
  }, {
    key: "toString",
    value: function toString() {
      var max_len_fn = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : function (d, shape) {
        return d >= shape.length - 2 ? 6 : 6;
      };
      return function () {
        var _marked2 = /*#__PURE__*/_regenerator["default"].mark(entries),
            _marked3 = /*#__PURE__*/_regenerator["default"].mark(str);

        if (!(max_len_fn instanceof Function)) {
          var max_len = max_len_fn * 1;

          max_len_fn = function max_len_fn() {
            return max_len;
          };
        }

        var shape = this.shape,
            data = this.data,
            strides = new Int32Array(shape.length);
        strides[strides.length - 1] = 1;

        for (var i = strides.length; --i > 0;) {
          strides[i - 1] = strides[i] * shape[i];
        }
        /** Collects the String representations of all (displayed) entries.
         */


        function entries(d, idx) {
          var max_len, last, first, data_idx, j, _j, _j2;

          return _regenerator["default"].wrap(function entries$(_context2) {
            while (1) {
              switch (_context2.prev = _context2.next) {
                case 0:
                  max_len = max_len_fn(d, shape), last = max_len >>> 1, first = max_len - last;

                  if (!(shape.length === d)) {
                    _context2.next = 16;
                    break;
                  }

                  data_idx = data[idx];
                  _context2.t0 = data_idx;
                  _context2.next = _context2.t0 === null ? 6 : _context2.t0 === undefined ? 9 : 12;
                  break;

                case 6:
                  _context2.next = 8;
                  return 'null';

                case 8:
                  return _context2.abrupt("break", 14);

                case 9:
                  _context2.next = 11;
                  return 'undefined';

                case 11:
                  return _context2.abrupt("break", 14);

                case 12:
                  _context2.next = 14;
                  return data_idx.toString();

                case 14:
                  _context2.next = 37;
                  break;

                case 16:
                  if (!(shape[d] > max_len + 1)) {
                    _context2.next = 31;
                    break;
                  }

                  j = 0;

                case 18:
                  if (!(j < first)) {
                    _context2.next = 23;
                    break;
                  }

                  return _context2.delegateYield(entries(d + 1, idx + j * strides[d]), "t1", 20);

                case 20:
                  j++;
                  _context2.next = 18;
                  break;

                case 23:
                  _j = shape[d] - last;

                case 24:
                  if (!(_j < shape[d])) {
                    _context2.next = 29;
                    break;
                  }

                  return _context2.delegateYield(entries(d + 1, idx + _j * strides[d]), "t2", 26);

                case 26:
                  _j++;
                  _context2.next = 24;
                  break;

                case 29:
                  _context2.next = 37;
                  break;

                case 31:
                  _j2 = 0;

                case 32:
                  if (!(_j2 < shape[d])) {
                    _context2.next = 37;
                    break;
                  }

                  return _context2.delegateYield(entries(d + 1, idx + _j2 * strides[d]), "t3", 34);

                case 34:
                  _j2++;
                  _context2.next = 32;
                  break;

                case 37:
                case "end":
                  return _context2.stop();
              }
            }
          }, _marked2);
        }

        entries = (0, _toConsumableArray2["default"])(entries(0, 0)); // pad all entries to the same string length

        if (this.ndim > 1) {
          var padLen = entries.reduce(function (a, b) {
            return Math.max(a, b.length);
          }, 0);
          entries = entries.map(function (e) {
            return e.padStart(padLen);
          });
        }

        var iEntries = 0;

        function str(indent, d, idx) {
          var prefix, infix, suffix, max_len, last, first, j, _j3, _j4;

          return _regenerator["default"].wrap(function str$(_context3) {
            while (1) {
              switch (_context3.prev = _context3.next) {
                case 0:
                  if (!(shape.length === d)) {
                    _context3.next = 4;
                    break;
                  }

                  _context3.next = 3;
                  return entries[iEntries++];

                case 3:
                  return _context3.abrupt("return");

                case 4:
                  indent += ' ';
                  prefix = '[ ', infix = ', ', suffix = ' ]';

                  if (d < shape.length - 1) {
                    infix = ',\n' + indent;
                    prefix = '[';
                    suffix = ']';
                  }

                  _context3.next = 9;
                  return prefix;

                case 9:
                  max_len = max_len_fn(d, shape), last = max_len >>> 1, first = max_len - last;

                  if (!(shape[d] > max_len + 1)) {
                    _context3.next = 31;
                    break;
                  }

                  j = 0;

                case 12:
                  if (!(j < first)) {
                    _context3.next = 19;
                    break;
                  }

                  return _context3.delegateYield(str(indent, d + 1, idx + j * strides[d]), "t0", 14);

                case 14:
                  _context3.next = 16;
                  return infix;

                case 16:
                  j++;
                  _context3.next = 12;
                  break;

                case 19:
                  _context3.next = 21;
                  return "...".concat(shape[d] - max_len, " more...");

                case 21:
                  _j3 = shape[d] - last;

                case 22:
                  if (!(_j3 < shape[d])) {
                    _context3.next = 29;
                    break;
                  }

                  _context3.next = 25;
                  return infix;

                case 25:
                  return _context3.delegateYield(str(indent, d + 1, idx + _j3 * strides[d]), "t1", 26);

                case 26:
                  _j3++;
                  _context3.next = 22;
                  break;

                case 29:
                  _context3.next = 40;
                  break;

                case 31:
                  _j4 = 0;

                case 32:
                  if (!(_j4 < shape[d])) {
                    _context3.next = 40;
                    break;
                  }

                  if (!(_j4 > 0)) {
                    _context3.next = 36;
                    break;
                  }

                  _context3.next = 36;
                  return infix;

                case 36:
                  return _context3.delegateYield(str(indent, d + 1, idx + _j4 * strides[d]), "t2", 37);

                case 37:
                  _j4++;
                  _context3.next = 32;
                  break;

                case 40:
                  _context3.next = 42;
                  return suffix;

                case 42:
                case "end":
                  return _context3.stop();
              }
            }
          }, _marked3);
        }

        return (0, _toConsumableArray2["default"])(str('', 0, 0)).join('');
      }.apply(this);
    }
  }, {
    key: "toNestedArray",
    value: function toNestedArray() {
      var data = this.data,
          shape = this.shape;
      var flat_idx = 0;

      var toNested = function toNested(d) {
        return d === shape.length ? data[flat_idx++] : Array.from({
          length: shape[d]
        }, function () {
          return toNested(d + 1);
        });
      };

      return toNested(0);
    } //
    // ITERATION
    //

  }, {
    key: Symbol.iterator,
    value: /*#__PURE__*/_regenerator["default"].mark(function value() {
      var shape, stride, i;
      return _regenerator["default"].wrap(function value$(_context4) {
        while (1) {
          switch (_context4.prev = _context4.next) {
            case 0:
              shape = this.shape.slice(1), stride = shape.reduce(function (a, b) {
                return a * b;
              }, 1);
              i = 0;

            case 2:
              if (!(i < this.data.length)) {
                _context4.next = 7;
                break;
              }

              _context4.next = 5;
              return new NDArray(shape, this.data.slice(i, i += stride));

            case 5:
              _context4.next = 2;
              break;

            case 7:
            case "end":
              return _context4.stop();
          }
        }
      }, value, this);
    })
  }, {
    key: "forEach",
    value: function forEach(consumer) {
      var len = this.shape[0],
          shape = this.shape.slice(1),
          stride = shape.reduce(function (a, b) {
        return a * b;
      }, 1);

      for (var i = 0, idx = 0; i < len; i++, idx += stride) {
        consumer(new NDArray(shape, this.data.slice(idx, idx + stride)), i);
      }
    }
  }, {
    key: "elems",
    value: /*#__PURE__*/_regenerator["default"].mark(function elems() {
      var _marked4, shape, data, multi_idx, flat_idx, elems;

      return _regenerator["default"].wrap(function elems$(_context6) {
        while (1) {
          switch (_context6.prev = _context6.next) {
            case 0:
              elems = function _ref(d) {
                return _regenerator["default"].wrap(function elems$(_context5) {
                  while (1) {
                    switch (_context5.prev = _context5.next) {
                      case 0:
                        if (!(d === shape.length)) {
                          _context5.next = 5;
                          break;
                        }

                        _context5.next = 3;
                        return [multi_idx.slice(), data[flat_idx++]];

                      case 3:
                        _context5.next = 11;
                        break;

                      case 5:
                        multi_idx[d] = 0;

                      case 6:
                        if (!(multi_idx[d] < shape[d])) {
                          _context5.next = 11;
                          break;
                        }

                        return _context5.delegateYield(elems(d + 1), "t0", 8);

                      case 8:
                        multi_idx[d]++;
                        _context5.next = 6;
                        break;

                      case 11:
                      case "end":
                        return _context5.stop();
                    }
                  }
                }, _marked4);
              };

              _marked4 = /*#__PURE__*/_regenerator["default"].mark(elems);
              shape = this.shape, data = this.data, multi_idx = new Int32Array(shape.length); // <- index in result

              flat_idx = 0;
              return _context6.delegateYield(elems(0), "t0", 5);

            case 5:
            case "end":
              return _context6.stop();
          }
        }
      }, elems, this);
    })
  }, {
    key: "forElems",
    value: function forElems(consumer) {
      var shape = this.shape,
          data = this.data,
          multi_idx = new Int32Array(shape.length); // <- index in result

      var flat_idx = 0;

      function forEntries(d) {
        if (d === shape.length) consumer.apply(void 0, [data[flat_idx++]].concat((0, _toConsumableArray2["default"])(multi_idx)));else for (multi_idx[d] = 0; multi_idx[d] < shape[d]; multi_idx[d]++) {
          forEntries(d + 1);
        }
      }

      forEntries(0);
    } //
    // TRANSFORMATION
    //

  }, {
    key: "valueOf",
    value: function valueOf() {
      if (this.shape.length === 0) return this.data[0];
      return this;
    }
  }, {
    key: "mapElems",
    value: function mapElems(dtype, mapper) {
      if (null == mapper) {
        if (dtype == null) return new NDArray(this.shape, this.data.slice());

        if (dtype instanceof Function) {
          mapper = dtype;
          dtype = undefined;
        }
      }

      if (null == mapper) mapper = function mapper(x) {
        return x;
      };
      return (0, _zip_elems.zip_elems)([this], dtype, mapper);
    }
  }, {
    key: "transpose",
    value: function transpose() // <- TODO allow for ellipsis '...' as input
    {
      for (var _len2 = arguments.length, axes = new Array(_len2), _key2 = 0; _key2 < _len2; _key2++) {
        axes[_key2] = arguments[_key2];
      }

      axes = (0, _toConsumableArray2["default"])(axes);
      var newShape = Int32Array.from(this.shape),
          ndim = newShape.length;
      var strides = new Int32Array(ndim);
      if (ndim > 0) strides[ndim - 1] = 1;

      for (var i = ndim; --i > 0;) {
        strides[i - 1] = newShape[i] * strides[i];
      } // BY DEFAULT THE LAST 2 AXES ARE SWAPPED


      if (axes.length == 0) {
        if (ndim <= 1) return this.sliceElems();
        newShape[ndim - 2] = this.shape[ndim - 1];
        newShape[ndim - 1] = this.shape[ndim - 2];
        var tmp = strides[ndim - 2];
        strides[ndim - 2] = strides[ndim - 1];
        strides[ndim - 1] = tmp;
      } // 
      else {
          var set = new Set(axes),
              _strides = strides;
          strides = new Int32Array(ndim);
          if (set.size != axes.length) throw new Error('Duplicate axes are not allowed.');

          for (var _i = axes.length; _i-- > 0;) {
            var j = axes[_i];
            if (0 > j || j >= ndim) throw new Error('Axis out of bounds.');
            newShape[_i] = this.shape[j];
            strides[_i] = _strides[j];
          } // COMPLETE WITH REMAINING/MISSING/IMPLIED INDICES


          for (var _i2 = 0, _j5 = set.size; _i2 < ndim; _i2++) {
            if (!set.has(_i2)) {
              newShape[_j5] = this.shape[_i2];
              strides[_j5++] = _strides[_i2];
            }
          }
        }

      var oldData = this.data,
          newData = new oldData.__proto__.constructor(oldData.length);
      var newI = oldData.length;

      function copy(d, oldI) {
        if (ndim > d) for (var _i3 = newShape[d]; _i3-- > 0; oldI -= strides[d]) {
          copy(d + 1, oldI);
        } else newData[--newI] = oldData[oldI];
      }

      copy(0, newI - 1);
      return new NDArray(newShape, newData);
    }
  }, {
    key: "reshape",
    value: function reshape() {
      for (var _len3 = arguments.length, shape = new Array(_len3), _key3 = 0; _key3 < _len3; _key3++) {
        shape[_key3] = arguments[_key3];
      }

      shape = Int32Array.from(shape);
      var len = this.data.length;
      var rest = 1,
          infer = -1;

      for (var i = shape.length; --i >= 0;) {
        if (shape[i] < 0) {
          if (shape[i] !== -1) throw new Error('Invalid dimension: ' + shape[i]);
          if (infer !== -1) throw new Error('At most on dimension may be -1.');
          infer = i;
        } else rest *= shape[i];
      }

      if (infer !== -1) {
        if (0 !== len % rest) throw new Error('Shape cannot be inferred.');
        shape[infer] = len / rest;
      }

      return new NDArray(shape, this.data);
    }
  }, {
    key: "reduceElems",
    value: function reduceElems(axes, dtype, reducer) {
      if (null == reducer) {
        if (null == dtype) return this.data.reduce(axes);
        if ('string' === typeof axes) return this.data.reduce(dtype);
        reducer = dtype;
        dtype = undefined;
      }

      if (null == dtype) dtype = 'object';
      if (!(0, _dt.is_subdtype)(this.dtype, dtype)) throw new Error('New dtype must be a super-dtype.');
      var oldNDim = this.ndim;
      if ('number' === typeof axes) axes = [axes];

      if (axes instanceof NDArray) {
        if (!(0, _dt.is_subdtype)(axes.dtype, 'int32')) throw new Error("Invalid dtype ".concat(axes.dtype, " for axes."));
        if (axes.ndim === 1) axes = axes.data;else throw new Error('Only 1D nd.Array allowed for axes.');
      }

      axes = new Set( /*#__PURE__*/_regenerator["default"].mark(function _callee() {
        var _iteratorNormalCompletion, _didIteratorError, _iteratorError, _iterator, _step, ax;

        return _regenerator["default"].wrap(function _callee$(_context7) {
          while (1) {
            switch (_context7.prev = _context7.next) {
              case 0:
                _iteratorNormalCompletion = true;
                _didIteratorError = false;
                _iteratorError = undefined;
                _context7.prev = 3;
                _iterator = axes[Symbol.iterator]();

              case 5:
                if (_iteratorNormalCompletion = (_step = _iterator.next()).done) {
                  _context7.next = 15;
                  break;
                }

                ax = _step.value;
                if (0 > ax) ax += oldNDim;

                if (!(0 > ax || ax >= oldNDim)) {
                  _context7.next = 10;
                  break;
                }

                throw new Error('Reduction axis ' + ax + ' out of bounds.');

              case 10:
                _context7.next = 12;
                return +ax;

              case 12:
                _iteratorNormalCompletion = true;
                _context7.next = 5;
                break;

              case 15:
                _context7.next = 21;
                break;

              case 17:
                _context7.prev = 17;
                _context7.t0 = _context7["catch"](3);
                _didIteratorError = true;
                _iteratorError = _context7.t0;

              case 21:
                _context7.prev = 21;
                _context7.prev = 22;

                if (!_iteratorNormalCompletion && _iterator["return"] != null) {
                  _iterator["return"]();
                }

              case 24:
                _context7.prev = 24;

                if (!_didIteratorError) {
                  _context7.next = 27;
                  break;
                }

                throw _iteratorError;

              case 27:
                return _context7.finish(24);

              case 28:
                return _context7.finish(21);

              case 29:
              case "end":
                return _context7.stop();
            }
          }
        }, _callee, null, [[3, 17, 21, 29], [22,, 24, 28]]);
      })());
      var oldShape = this.shape,
          newShape = oldShape.filter(function (size, i) {
        return !axes.has(i);
      }),
          newData = new _dt.ARRAY_TYPES[dtype](newShape.reduce(function (a, b) {
        return a * b;
      }, 1)),
          oldData = this.data;
      var newIdx = 0,
          oldIdx = 0;

      function fill(d, reduce) {
        if (oldShape.length === d) {
          if (reduce) newData[newIdx] = reducer(oldData[oldIdx++], newData[newIdx]);else newData[newIdx] = oldData[oldIdx++];
          ++newIdx;
        } else if (axes.has(d)) {
          var idx = newIdx;

          for (var j = oldShape[d]; j-- > 0; reduce = true) // <- copy over very first value and only then start reduction
          {
            newIdx = idx;
            fill(d + 1, reduce);
          }
        } else for (var _j6 = oldShape[d]; _j6-- > 0;) {
          fill(d + 1, reduce);
        }
      }

      fill(0, false);
      if (newIdx != newData.length) throw new Error([newIdx, newData.length]);
      if (oldIdx != this.data.length) throw new Error([oldIdx, this.data.length]);
      return new NDArray(newShape, newData);
    }
  }, {
    key: "sliceElems",
    value: function sliceElems() {
      for (var _len4 = arguments.length, slices = new Array(_len4), _key4 = 0; _key4 < _len4; _key4++) {
        slices[_key4] = arguments[_key4];
      }

      var oldNdim = this.ndim,
          oldShape = this.shape;
      var nEllipses = 0,
          nNew = 0,
          nReduced = 0,
          // <- dimensions of this matrix that to not appear in the sliced one
      nSliced = 0; // some basic acCOUNTing

      for (var _i4 = 0, _slices = slices; _i4 < _slices.length; _i4++) {
        var slc = _slices[_i4];
        if (slc === '...') ++nEllipses;else if (slc === 'new') ++nNew;else if (typeof slc === 'number') {
          if (!(slc % 1 == 0)) throw new Error('Only integral indices allowed.');else ++nReduced;
        } else {
          var len = slc.length;
          if (typeof len !== 'number' || !(len % 1 == 0) || len > 3 || len == 1) throw new Error('Illegal argument(s).');
          ++nSliced;
        }
      }

      if (nEllipses > 1) throw new Error('Highlander! There can be only one ... (ellipsis).');
      if (nEllipses == 0) slices.push('...');
      if (nReduced + nSliced > oldNdim) throw new Error('Cannot sliced more dimensions than available.'); // fill in omitted slice arguments and replace ellipsis

      var off = 0;
      var strides = [],
          newShape = Int32Array.from( /*#__PURE__*/_regenerator["default"].mark(function _callee2() {
        var d, stride, i, slc, _i5, _slc, _slc2, _slc2$, start, _slc2$2, end, _slc2$3, step;

        return _regenerator["default"].wrap(function _callee2$(_context8) {
          while (1) {
            switch (_context8.prev = _context8.next) {
              case 0:
                d = oldNdim, stride = 1;
                i = slices.length;

              case 2:
                if (!(--i >= 0)) {
                  _context8.next = 54;
                  break;
                }

                slc = slices[i];
                _context8.t0 = slc;
                _context8.next = _context8.t0 === 'new' ? 7 : _context8.t0 === '...' ? 11 : 21;
                break;

              case 7:
                strides.push(0);
                _context8.next = 10;
                return 1;

              case 10:
                return _context8.abrupt("break", 52);

              case 11:
                _i5 = nReduced + nSliced;

              case 12:
                if (!(_i5 < oldNdim)) {
                  _context8.next = 20;
                  break;
                }

                strides.push(stride);
                stride *= oldShape[--d];
                _context8.next = 17;
                return oldShape[d];

              case 17:
                _i5++;
                _context8.next = 12;
                break;

              case 20:
                return _context8.abrupt("break", 52);

              case 21:
                --d;

                if (!('number' === typeof slc)) {
                  _context8.next = 29;
                  break;
                }

                if (0 > slc) slc += oldShape[d];

                if (!(0 > slc || slc >= oldShape[d])) {
                  _context8.next = 26;
                  break;
                }

                throw new Error('Index out of bounds.');

              case 26:
                off += stride * slc;
                _context8.next = 51;
                break;

              case 29:
                _slc = slc, _slc2 = (0, _slicedToArray2["default"])(_slc, 3), _slc2$ = _slc2[0], start = _slc2$ === void 0 ? undefined : _slc2$, _slc2$2 = _slc2[1], end = _slc2$2 === void 0 ? undefined : _slc2$2, _slc2$3 = _slc2[2], step = _slc2$3 === void 0 ? 1 : _slc2$3;
                if (end < 0) end += oldShape[d];
                if (start < 0) start += oldShape[d];

                if (!(0 > start || start >= oldShape[d])) {
                  _context8.next = 34;
                  break;
                }

                throw new Error('Slice start out of bounds.');

              case 34:
                if (!(step == 0)) {
                  _context8.next = 36;
                  break;
                }

                throw new Error('Stride/step cannot be zero.');

              case 36:
                if (!(step < 0)) {
                  _context8.next = 43;
                  break;
                }

                if (start == null) start = oldShape[d] - 1;
                if (end == null) end = -1;

                if (!(-1 > end || end >= oldShape[d])) {
                  _context8.next = 41;
                  break;
                }

                throw new Error('Slice end out of bounds.');

              case 41:
                _context8.next = 47;
                break;

              case 43:
                if (start == null) start = 0;
                if (end == null) end = oldShape[d];

                if (!(0 > end || end > oldShape[d])) {
                  _context8.next = 47;
                  break;
                }

                throw new Error('Slice end out of bounds.');

              case 47:
                strides.push(stride * step);
                off += stride * start;
                _context8.next = 51;
                return 1 + Math.trunc((-Math.sign(step) + end - start) / step);

              case 51:
                stride *= oldShape[d];

              case 52:
                _context8.next = 2;
                break;

              case 54:
              case "end":
                return _context8.stop();
            }
          }
        }, _callee2);
      })());
      var oldData = this.data,
          newData = new _dt.ARRAY_TYPES[this.dtype](newShape.reduce(function (a, b) {
        return a * b;
      }, 1));
      var newIdx = 0;

      function fill(d, oldIdx) {
        if (0 === d) newData[newIdx++] = oldData[oldIdx];else for (var i = newShape[--d]; i-- > 0; oldIdx += strides[d]) {
          fill(d, oldIdx);
        }
      }

      fill(newShape.length, off);
      newShape.reverse();
      return new NDArray(newShape, newData);
    } //    convolve(options, conv_op) // <- IDEA: img.convolve({shape=[4,4,3], paddings=['same','same','none']}, (values, ...indices) => values.sum() )

  }, {
    key: "ndim",
    get: function get() {
      return this.shape.length;
    }
  }, {
    key: "dtype",
    get: function get() {
      for (var _i6 = 0, _Object$entries = Object.entries(_dt.ARRAY_TYPES); _i6 < _Object$entries.length; _i6++) {
        var _Object$entries$_i = (0, _slicedToArray2["default"])(_Object$entries[_i6], 2),
            dtype = _Object$entries$_i[0],
            ArrayType = _Object$entries$_i[1];

        if (this.data instanceof ArrayType) return dtype;
      }

      throw new Error('Data type could not determined.');
    }
  }, {
    key: "T",
    get: function get() {
      if (this.ndim == 0) return this.sliceElems(); // <- TODO maybe [[]]

      if (this.ndim == 1) return this.sliceElems([], 'new');
      return this.transpose();
    }
  }, {
    key: "H",
    get: function get() {
      var result = this.T;
      if (this.dtype != 'object' && this.dtype != 'complex128') return result;
      return result.mapElems(math.conj);
    }
  }]);
  return NDArray;
}( /*#__PURE__*/(0, _wrapNativeSuper2["default"])(Function));

exports.NDArray = NDArray;