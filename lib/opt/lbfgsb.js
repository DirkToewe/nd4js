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
exports.min_lbfgsb_gen = min_lbfgsb_gen;

var _regenerator = _interopRequireDefault(require("@babel/runtime/regenerator"));

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _nd_array = require("../nd_array");

var _float64_utils = require("../dt/float64_utils");

var _line_search_error = require("./line_search/line_search_error");

var _more_thuente_u = require("./line_search/more_thuente_u123");

var _lbfgsb_solver = require("./_lbfgsb_solver");

var _marked = /*#__PURE__*/_regenerator["default"].mark(min_lbfgsb_gen);

function _createForOfIteratorHelper(o, allowArrayLike) { var it; if (typeof Symbol === "undefined" || o[Symbol.iterator] == null) { if (Array.isArray(o) || (it = _unsupportedIterableToArray(o)) || allowArrayLike && o && typeof o.length === "number") { if (it) o = it; var i = 0; var F = function F() {}; return { s: F, n: function n() { if (i >= o.length) return { done: true }; return { done: false, value: o[i++] }; }, e: function e(_e) { throw _e; }, f: F }; } throw new TypeError("Invalid attempt to iterate non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method."); } var normalCompletion = true, didErr = false, err; return { s: function s() { it = o[Symbol.iterator](); }, n: function n() { var step = it.next(); normalCompletion = step.done; return step; }, e: function e(_e2) { didErr = true; err = _e2; }, f: function f() { try { if (!normalCompletion && it["return"] != null) it["return"](); } finally { if (didErr) throw err; } } }; }

function _unsupportedIterableToArray(o, minLen) { if (!o) return; if (typeof o === "string") return _arrayLikeToArray(o, minLen); var n = Object.prototype.toString.call(o).slice(8, -1); if (n === "Object" && o.constructor) n = o.constructor.name; if (n === "Map" || n === "Set") return Array.from(o); if (n === "Arguments" || /^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(n)) return _arrayLikeToArray(o, minLen); }

function _arrayLikeToArray(arr, len) { if (len == null || len > arr.length) len = arr.length; for (var i = 0, arr2 = new Array(len); i < len; i++) { arr2[i] = arr[i]; } return arr2; }

// REFERENCES
// ----------
// .. [1] "On the Limited Memory BFGS Method for Large Scale Optimization."
//         Dong C. Liu and Jorge Nocedal
function min_lbfgsb_gen(fg, x0, bounds) {
  var _ref,
      _ref$historySize,
      historySize,
      _ref$lineSearch,
      lineSearch,
      _ref$updateTol,
      updateTol,
      _ref$scaleInit,
      scaleInit,
      _x0$shape,
      L,
      shape,
      i,
      x,
      _fg,
      _fg2,
      f,
      g,
      solver,
      indices,
      dx,
      dg,
      step,
      _loop,
      _args2 = arguments;

  return _regenerator["default"].wrap(function min_lbfgsb_gen$(_context2) {
    while (1) {
      switch (_context2.prev = _context2.next) {
        case 0:
          _ref = _args2.length > 3 && _args2[3] !== undefined ? _args2[3] : {}, _ref$historySize = _ref.historySize, historySize = _ref$historySize === void 0 ? 8 : _ref$historySize, _ref$lineSearch = _ref.lineSearch, lineSearch = _ref$lineSearch === void 0 ? (0, _more_thuente_u.more_thuente_u123)({
            fRed: 1e-3,
            gRed: 0.8
          }) : _ref$lineSearch, _ref$updateTol = _ref.updateTol, updateTol = _ref$updateTol === void 0 ? 1e-14 : _ref$updateTol, _ref$scaleInit = _ref.scaleInit, scaleInit = _ref$scaleInit === void 0 ? 1e-3 : _ref$scaleInit;

          if (!(historySize % 1 !== 0)) {
            _context2.next = 3;
            break;
          }

          throw new Error('min_lbfgsb_gen(fg, x0, opt): opt.historySize must be an integer.');

        case 3:
          if (historySize > 0) {
            _context2.next = 5;
            break;
          }

          throw new Error('min_lbfgsb_gen(fg, x0, opt): opt.historySize must be positive.');

        case 5:
          if (fg instanceof Function) {
            _context2.next = 7;
            break;
          }

          throw new Error('min_lbfgsb_gen(fg,x0,bounds): fg must be a function.');

        case 7:
          x0 = (0, _nd_array.array)('float64', x0);
          bounds = (0, _nd_array.asarray)(bounds);
          _x0$shape = (0, _slicedToArray2["default"])(x0.shape, 1), L = _x0$shape[0], shape = x0.shape;

          if (!(x0.ndim !== 1)) {
            _context2.next = 12;
            break;
          }

          throw new Error('min_lbfgsb_gen(fg,x0,bounds): x0.ndim must be 1.');

        case 12:
          if (!(bounds.ndim !== 2)) {
            _context2.next = 14;
            break;
          }

          throw new Error('min_lbfgsb_gen(fg,x0,bounds): bounds.ndim must be 2.');

        case 14:
          if (!(bounds.shape[0] !== L)) {
            _context2.next = 16;
            break;
          }

          throw new Error('min_lbfgsb_gen(fg,x0,bounds): bounds.shape[0] must be x0.shape[0].');

        case 16:
          if (!(bounds.shape[1] !== 2)) {
            _context2.next = 18;
            break;
          }

          throw new Error('min_lbfgsb_gen(fg,x0,bounds): bounds.shape[1] must be 2.');

        case 18:
          updateTol *= 1;
          scaleInit *= 1;

          if (0 <= updateTol) {
            _context2.next = 22;
            break;
          }

          throw new Error('min_lbfgsb_gen(fg,x0,bounds,opt): opt.updateTol most be non-negative.');

        case 22:
          if (0 < scaleInit) {
            _context2.next = 24;
            break;
          }

          throw new Error('min_lbfgsb_gen(fg,x0,bounds,opt): opt.scaleInit most be greater than 0.');

        case 24:
          bounds = bounds.data;
          updateTol *= 1;
          i = bounds.length;

        case 27:
          if (!((i -= 2) >= 0)) {
            _context2.next = 32;
            break;
          }

          if (bounds[i + 0] <= bounds[i + 1]) {
            _context2.next = 30;
            break;
          }

          throw new Error('Assertion failed.');

        case 30:
          _context2.next = 27;
          break;

        case 32:
          x = x0;
          x0 = undefined;
          _fg = fg(x), _fg2 = (0, _slicedToArray2["default"])(_fg, 2), f = _fg2[0], g = _fg2[1];
          f *= 1;
          g = (0, _nd_array.array)('float64', g);

          if (isFinite(f)) {
            _context2.next = 39;
            break;
          }

          throw new Error('min_lbfgsb_gen(fg, x0, bounds): fg must return [float,float[]].');

        case 39:
          if (g.dtype.startsWith('float')) {
            _context2.next = 41;
            break;
          }

          throw new Error('min_lbfgsb_gen(fg, x0, bounds): fg must return [float,float[]].');

        case 41:
          lineSearch = lineSearch(fg);
          solver = new _lbfgsb_solver.LBFGSB_Solver(historySize, L);
          solver.scale = 1 / scaleInit; // TODO: Implement scaling as described in [1]

          indices = Int32Array.from({
            length: L
          }, function (_, i) {
            return i;
          }), dx = new Float64Array(L), dg = new Float64Array(L);

          step = function step() {
            var g_arg = g.data.slice(),
                x_arg = x.data.slice(),
                n_dir = new Float64Array(L);

            var _solver$compute_cauch = solver.compute_cauchyGeneralized(g_arg, x_arg, bounds, indices,
            /*complete=*/
            false),
                _solver$compute_cauch2 = (0, _slicedToArray2["default"])(_solver$compute_cauch, 2),
                n_fix = _solver$compute_cauch2[0],
                df = _solver$compute_cauch2[1];

            solver.compute_subspace_Hv(g_arg, n_dir, indices, n_fix); //*DEBUG*/    if( 0 !== n_fix ) throw new Error('Not yet tested.');                        
            //*DEBUG*/    if( 0 !== df ) throw new Error('Not yet tested.');
            //*DEBUG*/    if( ! g_arg.every( (gi,i) => gi === g.data[i] ) ) throw new Error('Not yet tested.');

            /*DEBUG*/

            if (!indices.subarray(0, n_fix).every(function (i) {
              return n_dir[i] === 0;
            }))
              /*DEBUG*/
              throw new Error('Assertion failed.');

            var αMax = function () {
              var αMax = Infinity;

              var _iterator = _createForOfIteratorHelper(indices.subarray(n_fix)),
                  _step;

              try {
                for (_iterator.s(); !(_step = _iterator.n()).done;) {
                  var _i = _step.value;

                  if (0 !== n_dir[_i]) {
                    var end = bounds[2 * _i + (n_dir[_i] < 0)];
                    var travel = (x_arg[_i] - end) / n_dir[_i];
                    if (!(0 <= travel)) throw new Error('Assertion failed: ' + travel);
                    var n = 0;

                    while (Math.sign(n_dir[_i]) * (x_arg[_i] - travel * n_dir[_i] - end) < 0) {
                      if (++n > 2) throw new Error('Assertion failed.');
                      travel = (0, _float64_utils.nextDown)(travel);
                    }

                    if (n_dir[_i] < 0 && x_arg[_i] - travel * n_dir[_i] > end) throw new Error('Assertion failed.');
                    if (n_dir[_i] > 0 && x_arg[_i] - travel * n_dir[_i] < end) throw new Error('Assertion failed.');
                    αMax = Math.min(αMax, travel);
                  }
                }
              } catch (err) {
                _iterator.e(err);
              } finally {
                _iterator.f();
              }

              return αMax;
            }(); // TODO: g_arg is only the quasi newton approximation of the gradient at x_arg.
            //       Alternatively we could compute the gradient at x_arg insteal which might
            //       improve robustness.


            var X, F, G;

            try {
              var _lineSearch = lineSearch(new _nd_array.NDArray(shape, x_arg.slice()), f + df, new _nd_array.NDArray(shape, g_arg), new _nd_array.NDArray(shape, n_dir),
              /*αMin=*/
              undefined,
              /*α0  =*/
              undefined,
              /*αMax=*/
              αMax);

              var _lineSearch2 = (0, _slicedToArray2["default"])(_lineSearch, 3);

              X = _lineSearch2[0];
              F = _lineSearch2[1];
              G = _lineSearch2[2];
            } catch (err) {
              if (err instanceof _line_search_error.LineSearchBoundReachedError || err instanceof _line_search_error.LineSearchBisectionError) {
                X = err.x;
                F = err.f;
                G = err.g;
              } else if (err instanceof _line_search_error.LineSearchNoProgressError && df < 0) {
                X = new _nd_array.NDArray(shape, x_arg);

                var _fg3 = fg(new _nd_array.NDArray(shape, x_arg.slice()));

                var _fg4 = (0, _slicedToArray2["default"])(_fg3, 2);

                F = _fg4[0];
                G = _fg4[1];
              } else throw err;
            }

            F *= 1;
            X = (0, _nd_array.array)('float64', X);
            G = (0, _nd_array.array)('float64', G);

            for (var _i2 = L; _i2-- > 0;) {
              dx[_i2] = X.data[_i2] - x.data[_i2];
            }

            for (var _i3 = L; _i3-- > 0;) {
              dg[_i3] = G.data[_i3] - g.data[_i3];
            }

            if (dx.every(function (x) {
              return x === 0;
            })) throw new _line_search_error.LineSearchNoProgressError();
            var mx = dg.reduce(function (mx, x) {
              return Math.max(mx, Math.abs(x));
            }, 0); // <- avoid underflow

            var gg = 0,
                gx = 0;

            for (var _i4 = dg.length; _i4-- > 0;) {
              var xi = dx[_i4] / mx,
                  gi = dg[_i4] / mx;
              gg += gi * gi;
              gx += xi * gi;
            }

            if (updateTol * gg < gx) {
              solver.update(dx, dg); // solver.scale = Math.max(solver.scale, gg/gx);

              solver.scale = gg / gx;
            }

            x = X;
            f = F;
            g = G;
          };

          _loop = /*#__PURE__*/_regenerator["default"].mark(function _loop() {
            var γ;
            return _regenerator["default"].wrap(function _loop$(_context) {
              while (1) {
                switch (_context.prev = _context.next) {
                  case 0:
                    if (x instanceof _nd_array.NDArray) {
                      _context.next = 2;
                      break;
                    }

                    throw new Error('Assertion failed.');

                  case 2:
                    if (g instanceof _nd_array.NDArray) {
                      _context.next = 4;
                      break;
                    }

                    throw new Error('Assertion failed.');

                  case 4:
                    if (!(x.shape[0] !== L)) {
                      _context.next = 6;
                      break;
                    }

                    throw new Error('Assertion#7 failed.');

                  case 6:
                    if (!(g.shape[0] !== L)) {
                      _context.next = 8;
                      break;
                    }

                    throw new Error('Assertion#8 failed.');

                  case 8:
                    if (!(x.ndim !== 1)) {
                      _context.next = 10;
                      break;
                    }

                    throw new Error('Assertion#4 failed.');

                  case 10:
                    if (!(g.ndim !== 1)) {
                      _context.next = 12;
                      break;
                    }

                    throw new Error('Assertion#6 failed.');

                  case 12:
                    if (!isNaN(f)) {
                      _context.next = 14;
                      break;
                    }

                    throw new Error('Assertion#5 failed.');

                  case 14:
                    γ = g.data.slice();
                    x.data.forEach(function (xi, i) {
                      var _bounds$subarray = bounds.subarray(2 * i, 2 + 2 * i),
                          _bounds$subarray2 = (0, _slicedToArray2["default"])(_bounds$subarray, 2),
                          x_min = _bounds$subarray2[0],
                          x_max = _bounds$subarray2[1],
                          gi = γ[i];

                      if (0 < gi && xi === x_min) γ[i] = 0;else if (0 > gi && xi === x_max) γ[i] = 0;
                    });
                    _context.next = 18;
                    return [new _nd_array.NDArray(shape, x.data.slice()), f, new _nd_array.NDArray(shape, γ), new _nd_array.NDArray(shape, g.data.slice())];

                  case 18:
                    _context.prev = 18;
                    step();
                    return _context.abrupt("break", 32);

                  case 23:
                    _context.prev = 23;
                    _context.t0 = _context["catch"](18);

                    if (!(_context.t0 instanceof _line_search_error.LineSearchError && solver.m > 0)) {
                      _context.next = 29;
                      break;
                    }

                    solver.forget(solver.m + 1 >>> 1);
                    _context.next = 30;
                    break;

                  case 29:
                    throw _context.t0;

                  case 30:
                    _context.next = 18;
                    break;

                  case 32:
                  case "end":
                    return _context.stop();
                }
              }
            }, _loop, null, [[18, 23]]);
          });

        case 47:
          return _context2.delegateYield(_loop(), "t0", 48);

        case 48:
          _context2.next = 47;
          break;

        case 50:
        case "end":
          return _context2.stop();
      }
    }
  }, _marked);
}