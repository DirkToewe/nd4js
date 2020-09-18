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
exports.min_lbfgs_gen = min_lbfgs_gen;
exports.lsq_lbfgs_gen = lsq_lbfgs_gen;
exports.fit_lbfgs_gen = fit_lbfgs_gen;

var _regenerator = _interopRequireDefault(require("@babel/runtime/regenerator"));

var _slicedToArray2 = _interopRequireDefault(require("@babel/runtime/helpers/slicedToArray"));

var _nd_array = require("../nd_array");

var _line_search_error = require("./line_search/line_search_error");

var _lbfgs_solver = require("./_lbfgs_solver");

var _more_thuente_abc = require("./line_search/more_thuente_abc");

var _marked = /*#__PURE__*/_regenerator["default"].mark(min_lbfgs_gen),
    _marked2 = /*#__PURE__*/_regenerator["default"].mark(lsq_lbfgs_gen),
    _marked3 = /*#__PURE__*/_regenerator["default"].mark(fit_lbfgs_gen);

// REFERENCES
// ----------
// .. [1] "On the Limited Memory BFGS Method for Large Scale Optimization."
//         Dong C. Liu and Jorge Nocedal
// .. [2] "Updating Quasi-Newton Matrices with Limited Storage"
//         Jorge Nocedal 
//         MATHEMATICS OF COMPUTATION, VOL. 35, NO. 151, JULY 1980, PAGES 773-78
//         https://courses.engr.illinois.edu/ece544na/fa2014/nocedal80.pdf
// TODO: Implement scaling as described in [1].
function min_lbfgs_gen(fg, x0) {
  var _ref,
      _ref$historySize,
      historySize,
      _ref$lineSearch,
      lineSearch,
      _ref$updateTol,
      updateTol,
      _ref$negDir,
      negDir0,
      _ref$scaling,
      scaling,
      x,
      _x$shape,
      L,
      shape,
      _fg,
      _fg2,
      f,
      g,
      solver,
      dx,
      dg,
      step,
      _args = arguments;

  return _regenerator["default"].wrap(function min_lbfgs_gen$(_context) {
    while (1) {
      switch (_context.prev = _context.next) {
        case 0:
          _ref = _args.length > 2 && _args[2] !== undefined ? _args[2] : {}, _ref$historySize = _ref.historySize, historySize = _ref$historySize === void 0 ? 8 : _ref$historySize, _ref$lineSearch = _ref.lineSearch, lineSearch = _ref$lineSearch === void 0 ? (0, _more_thuente_abc.more_thuente_abc)() : _ref$lineSearch, _ref$updateTol = _ref.updateTol, updateTol = _ref$updateTol === void 0 ? 1e-14 : _ref$updateTol, _ref$negDir = _ref.negDir0, negDir0 = _ref$negDir === void 0 ? function (g) {
            return g;
          } : _ref$negDir, _ref$scaling = _ref.scaling, scaling = _ref$scaling === void 0 ? function () {
            var scale = 1024;
            return {
              update: function update(dx, dg) {
                if (dx.ndim !== 1) throw new Error('Assertion failed.');
                if (dg.ndim !== 1) throw new Error('Assertion failed.');
                dx = dx.data;
                dg = dg.data;
                var len = dx.length;
                if (len !== dg.length) throw new Error('Assertion failed.');
                var mx = dg.reduce(function (mx, x) {
                  return Math.max(mx, Math.abs(x));
                }, 0); // <- avoid underflow

                var gg = 0,
                    gx = 0;

                for (var i = dg.length; i-- > 0;) {
                  var gi = dg[i] / mx,
                      xi = dx[i] / mx;
                  gg += gi * gi;
                  gx += xi * gi;
                }

                scale = gg / gx; // scale = Math.max(scale, gg/gx);

                if (!isFinite(scale)) throw new Error('Assertion failed.');
              },
              apply: function apply(negDir) {
                if (negDir.ndim !== 1) throw new Error('Assertion failed.');
                var d = negDir.data;

                for (var i = d.length; i-- > 0;) {
                  d[i] /= scale;
                }

                return negDir;
              }
            };
          }() : _ref$scaling;

          if (!(historySize % 1 !== 0)) {
            _context.next = 3;
            break;
          }

          throw new Error('min_lbfgs_gen(fg, x0, opt): opt.historySize must be an integer.');

        case 3:
          if (historySize > 0) {
            _context.next = 5;
            break;
          }

          throw new Error('min_lbfgs_gen(fg, x0, opt): opt.historySize must be positive.');

        case 5:
          x = (0, _nd_array.array)('float64', x0);
          x0 = undefined;

          if (!(x.ndim !== 1)) {
            _context.next = 9;
            break;
          }

          throw new Error('min_lbfgs_gen(fg, x0, opt): x0.ndim must be 1.');

        case 9:
          _x$shape = (0, _slicedToArray2["default"])(x.shape, 1), L = _x$shape[0], shape = x.shape;
          _fg = fg(x), _fg2 = (0, _slicedToArray2["default"])(_fg, 2), f = _fg2[0], g = _fg2[1];
          g = (0, _nd_array.array)('float64', g);
          f *= 1;

          if (!isNaN(f)) {
            _context.next = 15;
            break;
          }

          throw new Error('min_lbfgs_gen(fg, x0, opt): fg must return [float,float[]].');

        case 15:
          if (!(g.ndim !== 1)) {
            _context.next = 17;
            break;
          }

          throw new Error('min_lbfgs_gen(fg, x0, opt): fg must return [float,float[]].');

        case 17:
          if (!(g.shape[0] !== L)) {
            _context.next = 19;
            break;
          }

          throw new Error('min_lbfgs_gen(fg, x0, opt): fg must return [float,float[]].');

        case 19:
          lineSearch = lineSearch(fg);
          solver = new _lbfgs_solver.LBFGS_Solver(historySize, L);
          dx = new Float64Array(L), dg = new Float64Array(L);

          step = function step() {
            var negDir = Float64Array.from(g.data);
            solver.compute_Hv_phase1(negDir);
            negDir = negDir0(new _nd_array.NDArray(shape, negDir));
            negDir = (0, _nd_array.array)('float64', negDir);
            negDir = scaling.apply(negDir);
            negDir = (0, _nd_array.array)('float64', negDir);
            solver.compute_Hv_phase2(negDir.data);
            var X, F, G;
            var update = true;

            try {
              var _lineSearch = lineSearch(new _nd_array.NDArray(shape, x.data.slice()), f, new _nd_array.NDArray(shape, g.data.slice()), negDir);

              var _lineSearch2 = (0, _slicedToArray2["default"])(_lineSearch, 3);

              X = _lineSearch2[0];
              F = _lineSearch2[1];
              G = _lineSearch2[2];
            } catch (err) {
              if (!(err instanceof _line_search_error.LineSearchBisectionError)) throw err;
              X = err.x;
              F = err.f;
              G = err.g;
              update = false;
            }

            X = (0, _nd_array.array)('float64', X);
            G = (0, _nd_array.array)('float64', G);
            F *= 1;
            if (X.ndim !== 1) throw new Error('Assertion failed.');
            if (G.ndim !== 1) throw new Error('Assertion failed.');
            if (X.shape[0] !== L) throw new Error('Assertion failed.');
            if (G.shape[0] !== L) throw new Error('Assertion failed.');
            if (isNaN(F)) throw new Error('Assertion failed.');

            for (var i = L; i-- > 0;) {
              dg[i] = G.data[i] - g.data[i];
            }

            for (var _i = L; _i-- > 0;) {
              dx[_i] = X.data[_i] - x.data[_i];
            }

            if (dx.every(function (x) {
              return x === 0;
            })) throw new _line_search_error.LineSearchNoProgressError();

            update = update || function () {
              var mx = dg.reduce(function (mx, x) {
                return Math.max(mx, Math.abs(x));
              }, 0); // <- avoid underflow

              var gg = 0,
                  gx = 0;

              for (var _i2 = dg.length; _i2-- > 0;) {
                var gi = dg[_i2] / mx,
                    xi = dx[_i2] / mx;
                gg += gi * gi;
                gx += xi * gi;
              }

              return updateTol * gg < gx;
            }();

            if (update) {
              solver.update(dx, dg);
              scaling.update(new _nd_array.NDArray(shape, dx.slice()), new _nd_array.NDArray(shape, dg.slice()));
            }

            x = X;
            f = F;
            g = G;
          };

        case 23:
          if (x instanceof _nd_array.NDArray) {
            _context.next = 25;
            break;
          }

          throw new Error('Assertion failed.');

        case 25:
          if (g instanceof _nd_array.NDArray) {
            _context.next = 27;
            break;
          }

          throw new Error('Assertion failed.');

        case 27:
          if (!isNaN(f)) {
            _context.next = 29;
            break;
          }

          throw new Error('Assertion#5 failed.');

        case 29:
          if (!(x.ndim !== 1)) {
            _context.next = 31;
            break;
          }

          throw new Error('Assertion#4 failed.');

        case 31:
          if (!(g.ndim !== 1)) {
            _context.next = 33;
            break;
          }

          throw new Error('Assertion#6 failed.');

        case 33:
          if (!(x.shape[0] !== L)) {
            _context.next = 35;
            break;
          }

          throw new Error('Assertion#7 failed.');

        case 35:
          if (!(g.shape[0] !== L)) {
            _context.next = 37;
            break;
          }

          throw new Error('Assertion#8 failed.');

        case 37:
          _context.next = 39;
          return [new _nd_array.NDArray(shape, x.data.slice()), f, new _nd_array.NDArray(shape, g.data.slice())];

        case 39:
          _context.prev = 39;
          step();
          return _context.abrupt("break", 53);

        case 44:
          _context.prev = 44;
          _context.t0 = _context["catch"](39);

          if (!(_context.t0 instanceof _line_search_error.LineSearchError && solver.m > 0)) {
            _context.next = 50;
            break;
          }

          solver.forget(solver.m + 1 >>> 1);
          _context.next = 51;
          break;

        case 50:
          throw _context.t0;

        case 51:
          _context.next = 39;
          break;

        case 53:
          _context.next = 23;
          break;

        case 55:
        case "end":
          return _context.stop();
      }
    }
  }, _marked, null, [[39, 44]]);
}

function lsq_lbfgs_gen(fJ, //: (x0: float[nInputs]) => [f: float[nSamples], j: float[nSamples,nInputs]] - THE OPTIMIZATION FUNCTION AND ITS JACOBIAN
x0) {
  var _ref2,
      _ref2$historySize,
      historySize,
      _ref2$lineSearch,
      lineSearch,
      _ref2$updateTol,
      updateTol,
      x,
      _fJ,
      _fJ2,
      f,
      j,
      _j$shape,
      K,
      L,
      shape_x,
      shape_f,
      shape_j,
      shape_e,
      eg,
      last_x,
      last_f,
      last_j,
      negDir0,
      dx,
      dg,
      solver,
      _eg,
      _eg2,
      e,
      g,
      step,
      _args2 = arguments;

  return _regenerator["default"].wrap(function lsq_lbfgs_gen$(_context2) {
    while (1) {
      switch (_context2.prev = _context2.next) {
        case 0:
          _ref2 = _args2.length > 2 && _args2[2] !== undefined ? _args2[2] : {}, _ref2$historySize = _ref2.historySize, historySize = _ref2$historySize === void 0 ? 8 : _ref2$historySize, _ref2$lineSearch = _ref2.lineSearch, lineSearch = _ref2$lineSearch === void 0 ? (0, _more_thuente_abc.more_thuente_abc)() : _ref2$lineSearch, _ref2$updateTol = _ref2.updateTol, updateTol = _ref2$updateTol === void 0 ? 1e-13 : _ref2$updateTol;

          if (!(historySize % 1 !== 0)) {
            _context2.next = 3;
            break;
          }

          throw new Error('min_lbfgs_gen(fg, x0, opt): opt.historySize must be an integer.');

        case 3:
          if (historySize > 0) {
            _context2.next = 5;
            break;
          }

          throw new Error('min_lbfgs_gen(fg, x0, opt): opt.historySize must be positive.');

        case 5:
          // param | meaning
          // ------|-------------------------
          //  x    | function input
          //  f    | function output f(x)
          //  j    | jacobian f'(x)
          //  e    | mean squared error (MSE)
          //  g    | gradient of MSE
          x = (0, _nd_array.array)('float64', x0);
          x0 = undefined;

          if (!(x.ndim !== 1)) {
            _context2.next = 9;
            break;
          }

          throw new Error('min_lbfgs_gen(fg, x0, opt): x0.ndim must be 1.');

        case 9:
          _fJ = fJ(x), _fJ2 = (0, _slicedToArray2["default"])(_fJ, 2), f = _fJ2[0], j = _fJ2[1];
          f = (0, _nd_array.array)('float64', f);
          j = (0, _nd_array.array)('float64', j);

          if (!(f.ndim !== 1)) {
            _context2.next = 14;
            break;
          }

          throw new Error('lsq_lbfgs_gen(fg, x0: float[L], opt): fg must return [float[K],float[K,L]].');

        case 14:
          if (!(j.ndim !== 2)) {
            _context2.next = 16;
            break;
          }

          throw new Error('lsq_lbfgs_gen(fg, x0: float[L], opt): fg must return [float[K],float[K,L]].');

        case 16:
          _j$shape = (0, _slicedToArray2["default"])(j.shape, 2), K = _j$shape[0], L = _j$shape[1], shape_x = x.shape, shape_f = f.shape, shape_j = j.shape, shape_e = Int32Array.of();

          if (!(f.shape[0] !== K)) {
            _context2.next = 19;
            break;
          }

          throw new Error('lsq_lbfgs_gen(fg, x0: float[L], opt): fg must return [float[K],float[K,L]].');

        case 19:
          if (!(x.shape[0] !== L)) {
            _context2.next = 21;
            break;
          }

          throw new Error('lsq_lbfgs_gen(fg, x0: float[L], opt): fg must return [float[K],float[K,L]].');

        case 21:
          fJ = function () {
            var FJ = fJ;
            return function (x) {
              x = new _nd_array.NDArray(shape_x, x.data);

              var _FJ = FJ(x),
                  _FJ2 = (0, _slicedToArray2["default"])(_FJ, 2),
                  f = _FJ2[0],
                  j = _FJ2[1];

              f = (0, _nd_array.array)('float64', f);
              j = (0, _nd_array.array)('float64', j);
              if (f.data.some(isNaN)) throw new Error('Assertion failed.');
              if (j.data.some(isNaN)) throw new Error('Assertion failed.');
              if (f.ndim !== 1) throw new Error('lsq_lbfgs_gen(fg, x0: float[L], opt): fJ(x: float[L]) must return [float[K],float[K,L]].');
              if (j.ndim !== 2) throw new Error('lsq_lbfgs_gen(fg, x0: float[L], opt): fJ(x: float[L]) must return [float[K],float[K,L]].');
              if (f.shape[0] !== K) throw new Error('lsq_lbfgs_gen(fg, x0: float[L], opt): fJ(x: float[L]) must return [float[K],float[K,L]].');
              if (j.shape[0] !== K) throw new Error('lsq_lbfgs_gen(fg, x0: float[L], opt): fJ(x: float[L]) must return [float[K],float[K,L]].');
              if (j.shape[1] !== L) throw new Error('lsq_lbfgs_gen(fg, x0: float[L], opt): fJ(x: float[L]) must return [float[K],float[K,L]].');
              return [f, j];
            };
          }();

          eg = function eg(f, j) {
            var F = f.data,
                J = j.data;
            var e = 0,
                G = new Float64Array(L);

            for (var i = K; i-- > 0;) {
              e += F[i] * F[i] / K;
            }

            for (var _i3 = K; _i3-- > 0;) {
              for (var _j = L; _j-- > 0;) {
                G[_j] += 2 * F[_i3] * J[L * _i3 + _j] / K;
              }
            }

            return [e, new _nd_array.NDArray(shape_x, G)];
          };

          lineSearch = lineSearch(function (x) {
            last_x = (0, _nd_array.array)('float64', x);
            x = undefined;
            if (last_x.ndim !== 1) throw new Error('Assertion failed.');
            if (last_x.shape[0] !== L) throw new Error('Assertion failed.');

            var _fJ3 = fJ(last_x);

            var _fJ4 = (0, _slicedToArray2["default"])(_fJ3, 2);

            last_f = _fJ4[0];
            last_j = _fJ4[1];
            return eg(last_f, last_j);
          });

          negDir0 = function negDir0(negDir) {
            // OPTION 1: USE CAUCHY POINT OF THE GAUSS NEWTON APPROXIMATION
            var G = g.data,
                J = j.data; // compute the minimum along negDir (similar to cauchy point)

            var a = 0;

            for (var i = L; i-- > 0;) {
              a += G[i] * negDir[i];
            }

            if (0 < a) {
              var b = 0;

              for (var _i4 = K; _i4-- > 0;) {
                var Jd = 0;

                for (var _j2 = L; _j2-- > 0;) {
                  Jd += J[L * _i4 + _j2] * negDir[_j2];
                }

                b += Jd * Jd * 2 / K;
              }

              if (!(0 < b)) throw new Error('Assertion failed.');
              a /= b;

              for (var _i5 = L; _i5-- > 0;) {
                negDir[_i5] *= a;
              }
            } // // OPTION 1: USE GAUSS NEWTON APPROXIMATION
            // const nd = lstsq(
            //   matmul2(j.T, j),
            //   new NDArray(Int32Array.of(L,1), negDir)
            // );
            // negDir.set(nd.data);
            // for( let i=L; i-- > 0; )
            //   negDir[i] *= K/2;

          };

          dx = new Float64Array(L), dg = new Float64Array(L), solver = new _lbfgs_solver.LBFGS_Solver(historySize, L);
          _eg = eg(f, j), _eg2 = (0, _slicedToArray2["default"])(_eg, 2), e = _eg2[0], g = _eg2[1];

          step = function step() {
            var negDir = Float64Array.from(g.data);
            solver.compute_Hv_phase1(negDir);
            negDir0(negDir);
            solver.compute_Hv_phase2(negDir);
            var X,
                E,
                G,
                update = true;

            try {
              var _lineSearch3 = lineSearch(new _nd_array.NDArray(shape_x, x.data.slice()), e, new _nd_array.NDArray(shape_x, g.data.slice()), new _nd_array.NDArray(shape_x, negDir));

              var _lineSearch4 = (0, _slicedToArray2["default"])(_lineSearch3, 3);

              X = _lineSearch4[0];
              E = _lineSearch4[1];
              G = _lineSearch4[2];
            } catch (err) {
              if (!(err instanceof _line_search_error.LineSearchBisectionError)) throw err;
              X = err.x;
              E = err.f;
              G = err.g;
              update = false;
            }

            x = (0, _nd_array.array)('float64', x);
            e *= 1;
            g = (0, _nd_array.array)('float64', g);
            if (x.ndim !== 1) throw new Error('Assertion failed.');
            if (g.ndim !== 1) throw new Error('Assertion failed.');
            if (x.shape[0] !== L) throw new Error('Assertion failed.');
            if (g.shape[0] !== L) throw new Error('Assertion failed.');
            if (isNaN(e)) throw new Error('Assertion failed.');
            var F = last_f,
                J = last_j; // only recompute f,j is necessary

            if (last_x.data.some(function (xi, i) {
              return xi !== X.data[i];
            })) {
              var _fJ5 = fJ(X);

              var _fJ6 = (0, _slicedToArray2["default"])(_fJ5, 2);

              F = _fJ6[0];
              J = _fJ6[1];
            }

            for (var i = L; i-- > 0;) {
              dx[i] = X.data[i] - x.data[i];
            }

            for (var _i6 = L; _i6-- > 0;) {
              dg[_i6] = G.data[_i6] - g.data[_i6];
            }

            if (dx.every(function (x) {
              return x === 0;
            })) throw new _line_search_error.LineSearchNoProgressError();

            update = update || function () {
              var mx = dg.reduce(function (mx, x) {
                return Math.max(mx, Math.abs(x));
              }, 0); // <- avoid underflow

              var gg = 0,
                  gx = 0;

              for (var _i7 = dg.length; _i7-- > 0;) {
                var gi = dg[_i7] / mx,
                    xi = dx[_i7] / mx;
                gg += gi * gi;
                gx += xi * gi;
              }

              return updateTol * gg < gx;
            }();

            if (update) solver.update(dx, dg);
            x = X;
            e = E;
            g = G;
            f = F;
            j = J;
          };

        case 28:
          _context2.next = 30;
          return [new _nd_array.NDArray(shape_x, x.data.slice()), e, new _nd_array.NDArray(shape_x, g.data.slice()), new _nd_array.NDArray(shape_f, f.data.slice()), new _nd_array.NDArray(shape_j, j.data.slice())];

        case 30:
          _context2.prev = 30;
          step();
          return _context2.abrupt("break", 44);

        case 35:
          _context2.prev = 35;
          _context2.t0 = _context2["catch"](30);

          if (!(_context2.t0 instanceof _line_search_error.LineSearchError && solver.m > 0)) {
            _context2.next = 41;
            break;
          }

          solver.forget(solver.m + 1 >>> 1);
          _context2.next = 42;
          break;

        case 41:
          throw _context2.t0;

        case 42:
          _context2.next = 30;
          break;

        case 44:
          _context2.next = 28;
          break;

        case 46:
        case "end":
          return _context2.stop();
      }
    }
  }, _marked2, null, [[30, 35]]);
}

function fit_lbfgs_gen(x, //: float[nSamples,nDim]
y, //: float[nSamples]
fg, //: (params: float[nParam]) => (x: float[nDim]) => [f: float, g: float[nParam]],
p0, //: float[nParam]
opt) {
  var FloatArray, _p0$shape, P, _x$shape2, M, N, x_shape, R, J, RJ, fJ;

  return _regenerator["default"].wrap(function fit_lbfgs_gen$(_context3) {
    while (1) {
      switch (_context3.prev = _context3.next) {
        case 0:
          if (fg instanceof Function) {
            _context3.next = 2;
            break;
          }

          throw new Error('fit_lbfgs_gen(x,y, fg, p0): fg must be a function.');

        case 2:
          x = (0, _nd_array.array)(x); // <- TODO allow non-ndarray x ?

          y = (0, _nd_array.array)('float64', y);
          p0 = (0, _nd_array.array)('float64', p0);
          FloatArray = Float64Array;

          if (!(x.ndim !== 2)) {
            _context3.next = 8;
            break;
          }

          throw new Error('fit_lbfgs_gen(x,y, fg, p0): x.ndim must be 2.');

        case 8:
          if (!(y.ndim !== 1)) {
            _context3.next = 10;
            break;
          }

          throw new Error('fit_lbfgs_gen(x,y, fg, p0): y.ndim must be 1.');

        case 10:
          if (!(p0.ndim !== 1)) {
            _context3.next = 12;
            break;
          }

          throw new Error('fit_lbfgs_gen(x,y, fg, p0): p0.ndim must be 1.');

        case 12:
          _p0$shape = (0, _slicedToArray2["default"])(p0.shape, 1), P = _p0$shape[0], _x$shape2 = (0, _slicedToArray2["default"])(x.shape, 2), M = _x$shape2[0], N = _x$shape2[1], x_shape = Int32Array.of(N);

          if (!(M != y.shape[0])) {
            _context3.next = 15;
            break;
          }

          throw new Error('fit_lbfgs_gen(x,y, fg, p0): x.shape[0] must be equal to y.shape[0].');

        case 15:
          x = x.data.slice(); // <- TODO: TypedArray.subarray could be used here

          y = y.data; // if x is frozen, no protection copies are necessary while passing it to fg

          Object.freeze(x.buffer);
          x = Array.from({
            length: M
          }, function (_, i) {
            return Object.freeze(new _nd_array.NDArray(x_shape, x.subarray(N * i, N * i + 1)));
          });
          Object.freeze(x);
          R = new FloatArray(M), J = new FloatArray(M * P), RJ = [new _nd_array.NDArray(Int32Array.of(M), R), new _nd_array.NDArray(Int32Array.of(M, P), J)];

          fJ = function fJ(p) {
            var fgp = fg(p);
            if (!(fgp instanceof Function)) throw new Error();

            for (var i = M; i-- > 0;) {
              var _fgp = fgp(x[i]),
                  _fgp2 = (0, _slicedToArray2["default"])(_fgp, 2),
                  _f = _fgp2[0],
                  _g = _fgp2[1];

              _f = (0, _nd_array.asarray)(_f);
              _g = (0, _nd_array.asarray)(_g);
              if (_f.ndim !== 0) throw new Error('fit_lbfgs_gen(x,y, fg, p0): fg must have signature float[nParams] => float[nDim] => [float, float[nParams]].');
              if (_g.ndim !== 1) throw new Error('fit_lbfgs_gen(x,y, fg, p0): fg must have signature float[nParams] => float[nDim] => [float, float[nParams]].');
              if (_g.shape[0] !== P) throw new Error('fit_lbfgs_gen(x,y, fg, p0): fg must have signature float[nParams] => float[nDim] => [float, float[nParams]].');
              _f = _f.data[0];
              _g = _g.data;
              if (isNaN(_f)) throw new Error('fit_lbfgs_gen(x,y, fg, p0): NaN encountered.');
              if (_g.some(isNaN)) throw new Error('fit_lbfgs_gen(x,y, fg, p0): NaN encountered.');
              R[i] = _f - y[i];

              for (var _j3 = P; _j3-- > 0;) {
                J[P * i + _j3] = _g[_j3];
              }
            }

            return RJ;
          };

          return _context3.delegateYield(lsq_lbfgs_gen(fJ, p0, opt), "t0", 23);

        case 23:
        case "end":
          return _context3.stop();
      }
    }
  }, _marked3);
}