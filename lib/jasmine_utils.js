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
exports.forEachItemIn = exports.CUSTOM_MATCHERS = void 0;

var _toConsumableArray2 = _interopRequireDefault(require("@babel/runtime/helpers/toConsumableArray"));

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _nd_array = require("./nd_array");

var _math = _interopRequireDefault(require("./math"));

function unravel(flat, shape) {
  var result = Int32Array.from(shape);

  for (var i = result.length; i-- > 0;) {
    result[i] = flat % shape[i];
    flat = Math.trunc(flat / shape[i]);
  }

  return result;
}

var toBeBand = function toBeBand(label, lower, upper) {
  return function (util, customEq) {
    return {
      compare: function compare(act) {
        var _ref = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : {},
            _ref$tol = _ref.tol,
            tol = _ref$tol === void 0 ? 0 : _ref$tol;

        var LEN = act.data.length,
            _act$shape$slice = act.shape.slice(-2),
            _act$shape$slice2 = (0, _slicedToArray2["default"])(_act$shape$slice, 2),
            M = _act$shape$slice2[0],
            N = _act$shape$slice2[1];

        for (var flat = 0; flat < LEN; flat++) {
          var i = Math.trunc(flat / N) % M,
              j = flat % N;
          if ((i - j > lower || j - i > upper) && !(_math["default"].abs(act.data[flat]) <= tol)) return {
            pass: false,
            message: "Expected\n".concat(act, "\nto be ").concat(label, ", but actual(").concat(unravel(flat, act.shape), ") = ").concat(act.data[flat], " \u2249 0.")
          };
        }

        return {
          pass: true
        };
      }
    };
  };
};

var CUSTOM_MATCHERS = {
  toBeAllCloseTo: function toBeAllCloseTo(util, customEq) {
    return {
      compare: function compare(act, exp) {
        var _ref2 = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : {},
            _ref2$rtol = _ref2.rtol,
            rtol = _ref2$rtol === void 0 ? 1e-5 : _ref2$rtol,
            _ref2$atol = _ref2.atol,
            atol = _ref2$atol === void 0 ? 1e-8 : _ref2$atol;

        act = (0, _nd_array.asarray)(act);
        exp = (0, _nd_array.asarray)(exp);
        if (Object.is(act, exp)) console.warn('Actual and expected are identical.');
        var ndim = Math.max(act.ndim, exp.ndim);
        var shape = Array.from({
          length: ndim
        }, function () {
          return 1;
        });

        for (var a = act.ndim, e = exp.ndim, i = ndim; i-- > 0;) {
          var as = act.shape[--a] || 1,
              es = exp.shape[--e] || 1;
          if (as === 1) shape[i] = es;else if (es !== 1 && es !== as) return {
            pass: false,
            message: "Expected shape [".concat(act.shape, "] to be broadcast-compatible to [").concat(exp.shape, "].")
          };else shape[i] = as;
        }

        var flatAct = 0,
            strideAct,
            flatExp = 0,
            strideExp;

        var is_close = function is_close(x, y) {
          var tol = atol + rtol * _math["default"].max(_math["default"].abs(x), _math["default"].abs(y));

          return _math["default"].abs(_math["default"].sub(x, y)) <= tol;
        };

        function visit(d) {
          if (ndim === d) {
            if (!is_close(act.data[flatAct], exp.data[flatExp])) throw "actual:\n".concat(act, "\nexpected:\n").concat(exp, "\nactual(").concat(unravel(flatAct, act.shape), ") = ").concat(act.data[flatAct], " \u2249 ").concat(exp.data[flatExp], " = expected(").concat(unravel(flatExp, exp.shape), ")");
            flatAct += strideAct = 1;
            flatExp += strideExp = 1;
            return;
          }

          for (var _i = 0;;) {
            visit(d + 1);
            if (++_i === shape[d]) break;
            if (!(act.shape[d - ndim + act.ndim] > 1)) flatAct -= strideAct;
            if (!(exp.shape[d - ndim + exp.ndim] > 1)) flatExp -= strideExp;
          }

          strideAct *= act.shape[d - ndim + act.ndim] || 1;
          strideExp *= exp.shape[d - ndim + exp.ndim] || 1;
        }

        try {
          visit(0);
          return {
            pass: true
          };
        } catch (message) {
          return {
            pass: false,
            message: message
          };
        }
      }
    };
  },
  toBeCloseTo: function toBeCloseTo(util, customEq) {
    return {
      compare: function compare(act, exp) {
        var _ref3 = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : {},
            _ref3$rtol = _ref3.rtol,
            rtol = _ref3$rtol === void 0 ? 1e-5 : _ref3$rtol,
            _ref3$atol = _ref3.atol,
            atol = _ref3$atol === void 0 ? 1e-8 : _ref3$atol;

        var tol = atol + rtol * _math["default"].max(_math["default"].abs(act), _math["default"].abs(exp));

        var result = {
          pass: _math["default"].abs(_math["default"].sub(act, exp)) <= tol
        };
        if (!result.pass) result.message = "Expected actual=".concat(act, " to be close to expected=").concat(exp, ".");
        return result;
      }
    };
  },
  toBeDiagonal: toBeBand('diagonal', 0, 0),
  toBeLowerBidiagonal: toBeBand('lower bidiagonal', +1, 0),
  toBeUpperBidiagonal: toBeBand('upper bidiagonal', 0, +1),
  toBeTridiagonal: toBeBand('tridiagonal', +1, +1),
  toBeLowerTriangular: toBeBand('lower triangular', Infinity, 0),
  toBeUpperTriangular: toBeBand('upper triangular', 0, Infinity),
  toBeLowerHessenberg: toBeBand('lower hessenberg', Infinity, +1),
  toBeUpperHessenberg: toBeBand('upper hessenberg', +1, Infinity)
};
exports.CUSTOM_MATCHERS = CUSTOM_MATCHERS;

function toStr(item) {
  if (ARRAY_TYPES.some(function (Type) {
    return item instanceof Type;
  })) {
    var lf = '',
        infix = ', ',
        result;

    if (item.length > 9) {
      var head = Array.from(item.slice(0, 4), toStr),
          tail = Array.from(item.slice(-4), toStr);

      if ([].concat((0, _toConsumableArray2["default"])(head), (0, _toConsumableArray2["default"])(tail)).some(function (s) {
        return /\n/.test(s);
      })) {
        infix = '\n,\n';
        lf = '\n';
      }

      result = "".concat(head.join(infix), ",").concat(lf || ' ', "... ").concat(item.length - 8, " more...").concat(lf || ' ', ", ").concat(tail.join(infix));
    } else {
      result = Array.from(item, toStr);

      if (result.some(function (s) {
        return /\n/.test(s);
      })) {
        infix = '\n,\n';
        lf = '\n';
      }

      result = result.join(infix);
    }

    return item instanceof Array ? "[".concat(lf).concat(result).concat(lf, "]") : "".concat(item.constructor.name, ".of(").concat(lf).concat(result).concat(lf, ")");
  }

  return 'string' === typeof item ? "\"".concat(item, "\"") : "".concat(item);
}

var ARRAY_TYPES = [Array, Int8Array, Uint8Array, Uint8ClampedArray, Float32Array, Int16Array, Uint16Array, Float64Array, Int32Array, Uint32Array];

var forEachItemIn = function forEachItemIn(items) {
  return {
    it: function (_it) {
      function it(_x, _x2, _x3) {
        return _it.apply(this, arguments);
      }

      it.toString = function () {
        return _it.toString();
      };

      return it;
    }(function (description, specFn, timeout) {
      if ('string' !== typeof description) throw new Error('description must be string.');
      var item, i;
      var spec = it(description, function () {
        return new Promise(function (resolve, reject) {
          var iter = items[Symbol.iterator]();
          i = -1;
          var intervalID = setInterval(function () {
            try {
              var t0 = performance.now();

              for (;;) {
                var _iter$next = iter.next(),
                    value = _iter$next.value,
                    done = _iter$next.done;

                ++i;

                if (done) {
                  if (0 === i) throw new Error('forEachItemIn: items may not be empty.'); //              console.log(`${i} tests ↴`)

                  clearInterval(intervalID);
                  return resolve();
                }

                specFn(item = value);
                if (250 < performance.now() - t0) return;
              }
            } catch (err) {
              clearInterval(intervalID);
              return reject(err);
            }
          }, 0);
        });
      }, timeout);

      spec.addExpectationResult = function (passed, data, isError) {
        if (passed) return;

        try {
          if (!passed) {
            var msg = function msg(_msg) {
              var str = toStr(item);
              if (/\n/.test(str)) str = '\n' + str;
              return "For item[".concat(i, "] = ").concat(str, ":\n").concat(_msg);
            };

            if ('message' in data) data.message = msg(data.message);else data.error.message = msg(data.error.message);
          }
        } catch (err) {
          console.error(err);
        } finally {
          Object.getPrototypeOf(spec).addExpectationResult.call(spec, passed, data, isError);
        }
      };

      return spec;
    })
  };
};

exports.forEachItemIn = forEachItemIn;