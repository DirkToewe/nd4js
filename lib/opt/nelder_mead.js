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
exports.min_nelder_mead_gen = min_nelder_mead_gen;

var _regenerator = _interopRequireDefault(require("@babel/runtime/regenerator"));

var _nd_array = require("../nd_array");

var _simplex = require("../geom/simplex");

var _matmul2 = require("../la/matmul");

var _norm = require("../la/norm");

var _min_max = require("../iter/min_max");

var _alea_rng = require("../rand/alea_rng");

var _optimization_error = require("./optimization_error");

var _marked = /*#__PURE__*/_regenerator["default"].mark(min_nelder_mead_gen);

// REFERENCES
// ----------
// .. [1] "A  simplex  method  for  function  minimization."
//         J. A. Nelder and R. Mead
var RNG = new _alea_rng.AleaRNG('opt/nelder_mead.js');

function min_nelder_mead_gen(f, x0) {
  var _ref,
      _ref$scale,
      scale,
      _ref$reflect,
      reflect,
      _ref$expand,
      expand,
      _ref$worstRatio,
      worstRatio,
      _ref$contractOuter,
      contractOuter,
      _ref$contractInner,
      contractInner,
      _ref$shrink,
      shrink,
      _ref$rng,
      rng,
      N,
      shape_x,
      verts_f,
      centroid,
      x,
      distance,
      _matmul,
      verts_x,
      i,
      j,
      _i,
      b,
      _i2,
      _j,
      x_ij,
      f0,
      _b,
      f_best,
      f_worst,
      _i3,
      _j2,
      _j3,
      dist_min,
      dist_max,
      _i4,
      _dist,
      dist,
      _i5,
      f_reflect,
      _i6,
      _i7,
      f_expand,
      _i8,
      _i9,
      f_contract,
      _i10,
      _i11,
      _f_contract,
      _i12,
      stuck,
      _i13,
      _j4,
      _x_ij,
      c,
      _args = arguments;

  return _regenerator["default"].wrap(function min_nelder_mead_gen$(_context) {
    while (1) {
      switch (_context.prev = _context.next) {
        case 0:
          distance = function _distance(x, y, y_off) {
            if (!(y_off >= 0)) throw new Error('Assertion failed.');
            if (!(y_off <= y.length - N)) throw new Error('Assertion failed.');
            if (y_off % N !== 0) throw new Error('Assertion failed.');
            var norm = new _norm.FrobeniusNorm();

            for (var i = N; i-- > 0;) {
              norm.include(x[i] - y[i + y_off]);
            }

            return norm.result;
          };

          _ref = _args.length > 2 && _args[2] !== undefined ? _args[2] : {}, _ref$scale = _ref.scale, scale = _ref$scale === void 0 ? 0.001 : _ref$scale, _ref$reflect = _ref.reflect, reflect = _ref$reflect === void 0 ? 1 : _ref$reflect, _ref$expand = _ref.expand, expand = _ref$expand === void 0 ? Math.E / 8 : _ref$expand, _ref$worstRatio = _ref.worstRatio, worstRatio = _ref$worstRatio === void 0 ? 1e-10 : _ref$worstRatio, _ref$contractOuter = _ref.contractOuter, contractOuter = _ref$contractOuter === void 0 ? 2 / (1 + Math.sqrt(5)) : _ref$contractOuter, _ref$contractInner = _ref.contractInner, contractInner = _ref$contractInner === void 0 ? 2 / (1 + Math.sqrt(5)) : _ref$contractInner, _ref$shrink = _ref.shrink, shrink = _ref$shrink === void 0 ? 2.5 / Math.PI : _ref$shrink, _ref$rng = _ref.rng, rng = _ref$rng === void 0 ? RNG : _ref$rng;

          if (f instanceof Function) {
            _context.next = 4;
            break;
          }

          throw new Error('min_nelder_mead_gen(f,x0): f has to be a Function.');

        case 4:
          if (0 < scale) {
            _context.next = 6;
            break;
          }

          throw new Error('min_nelder_mead_gen(f,x0,opt={}): opt.scale must be greater than 0.');

        case 6:
          if (0 < reflect) {
            _context.next = 8;
            break;
          }

          throw new Error('min_nelder_mead_gen(f,x0,opt={}): opt.reflect must be greater than 0.');

        case 8:
          if (0 <= expand) {
            _context.next = 10;
            break;
          }

          throw new Error('min_nelder_mead_gen(f,x0,opt={}): opt.expand must be non-negative.');

        case 10:
          if (!(0 < expand)) console.warn('min_nelder_mead_gen(f,x0,opt={}): opt.expand should be greater than 0.');

          if (0 < contractOuter) {
            _context.next = 13;
            break;
          }

          throw new Error('min_nelder_mead_gen(f,x0,opt={}): opt.contractOuter must be greater than 0.');

        case 13:
          if (1 >= contractOuter) {
            _context.next = 15;
            break;
          }

          throw new Error('min_nelder_mead_gen(f,x0,opt={}): opt.contractOuter must not be greater than 1.');

        case 15:
          if (!(1 > contractOuter)) console.warn('min_nelder_mead_gen(f,x0,opt={}): opt.contractOuter must should be less than 1.');

          if (0 < contractInner) {
            _context.next = 18;
            break;
          }

          throw new Error('min_nelder_mead_gen(f,x0,opt={}): opt.contractInner must be greater than 0.');

        case 18:
          if (1 >= contractInner) {
            _context.next = 20;
            break;
          }

          throw new Error('min_nelder_mead_gen(f,x0,opt={}): opt.contractInner must not be greater than 1.');

        case 20:
          if (!(1 > contractInner)) console.warn('min_nelder_mead_gen(f,x0,opt={}): opt.contractInner must should be less than 1.');

          if (0 < shrink) {
            _context.next = 23;
            break;
          }

          throw new Error('min_nelder_mead_gen(f,x0,opt={}): opt.shrink must be greater than 0.');

        case 23:
          if (1 >= shrink) {
            _context.next = 25;
            break;
          }

          throw new Error('min_nelder_mead_gen(f,x0,opt={}): opt.shrink must not be greater than 1.');

        case 25:
          if (!(1 > shrink)) console.warn('min_nelder_mead_gen(f,x0,opt={}): opt.shrink must should be less than 1.');
          x0 = (0, _nd_array.array)('float64', x0);

          if (!(x0.ndim !== 1)) {
            _context.next = 29;
            break;
          }

          throw new Error('min_nelder_mead_gen(f,x0): x0.ndim has to be 1.');

        case 29:
          N = x0.shape[0], shape_x = x0.shape;
          verts_f = new Float64Array(N + 1), centroid = new Float64Array(N), x = x0.data;
          x0 = undefined;

          f = function () {
            var F = f;
            return function (x) {
              var off = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0;
              if (!(off >= 0)) throw new Error('Assertion failed.');
              if (!(off <= x.length - N)) throw new Error('Assertion failed.');
              if (off % N !== 0) throw new Error('Assertion failed.');
              var f = 1 * F(new _nd_array.NDArray(shape_x, x.slice(off, off + N)));
              if (isNaN(f)) throw new Error('Assertion failed.');
              return f;
            };
          }();

          // INIT SIMPLEX
          _matmul = (0, _matmul2.matmul2)((0, _simplex.regular_simplex)(N), rng.ortho(N)), verts_x = _matmul.data;

          for (i = N + 1; i-- > 0;) {
            for (j = N; j-- > 0;) {
              verts_x[N * i + j] = verts_x[N * i + j] * scale + x[j];
            }
          }

          for (_i = N + 1; _i-- > 0;) {
            verts_f[_i] = f(verts_x, N * _i);
          }

          ;
          b = (0, _min_max.argmin)(verts_f); // YIELD STARTING POINT (and simplex)

          _context.next = 40;
          return [new _nd_array.NDArray(shape_x, verts_x.slice(N * b, N * (b + 1))), verts_f[b]];

        case 40:
          // FIND WORST VERTEX (and move it to row 0)
          ;
          _i2 = (0, _min_max.argmax)(verts_f);

          for (_j = N; _j-- > 0;) {
            x_ij = verts_x[N * _i2 + _j];
            verts_x[N * _i2 + _j] = verts_x[_j];
            verts_x[_j] = x_ij;
          }

          f0 = verts_f[0];
          verts_f[0] = verts_f[_i2];
          verts_f[_i2] = f0;
          _b = (0, _min_max.argmin)(verts_f.subarray(1)) + 1, f_best = verts_f[_b], f_worst = verts_f[0]; // COMPUTE CENTROID OF REMAINING VERTICES (write to x2)

          centroid.fill(0.0);

          for (_i3 = N + 1; _i3-- > 1;) {
            for (_j2 = N; _j2-- > 0;) {
              centroid[_j2] += verts_x[N * _i3 + _j2];
            }
          }

          for (_j3 = N; _j3-- > 0;) {
            centroid[_j3] /= N;
          }

          dist_min = +Infinity, dist_max = -Infinity;

          for (_i4 = N + 1; _i4-- > 1;) {
            _dist = distance(centroid, verts_x, N * _i4);
            dist_min = Math.min(_dist, dist_min);
            dist_max = Math.max(_dist, dist_max);
          }

          dist = 0.5 * distance(centroid, verts_x, 0); // if( !(0 < dist) )
          //   throw new OptimizationNoProgressError();

          if (!(0 < dist)) {
            _context.next = 64;
            break;
          }

          // REFLECT
          for (_i5 = N; _i5-- > 0;) {
            x[_i5] = centroid[_i5] + reflect * (centroid[_i5] - verts_x[_i5]);
          }

          f_reflect = f(x);

          if (!(f_reflect < f_worst)) {
            _context.next = 64;
            break;
          }

          for (_i6 = N; _i6-- > 0;) {
            verts_x[_i6] = x[_i6];
          }

          verts_f[0] = f_reflect;

          if (!(f_reflect < f_best)) {
            _context.next = 63;
            break;
          }

          // EXPAND (if reflection new best)
          if (dist * worstRatio <= dist_min) {
            for (_i7 = N; _i7-- > 0;) {
              x[_i7] = centroid[_i7] + expand * (x[_i7] - centroid[_i7]);
            }

            f_expand = f(x);

            if (f_expand < f_reflect) {
              for (_i8 = N; _i8-- > 0;) {
                verts_x[_i8] = x[_i8];
              }

              verts_f[0] = f_expand;
            }
          } // YIELD NEW BEST POINT (and simplex)


          _context.next = 63;
          return [new _nd_array.NDArray(shape_x, verts_x.slice(0, N)), verts_f[0]];

        case 63:
          return _context.abrupt("continue", 107);

        case 64:
          if (!(dist >= dist_max * worstRatio)) {
            _context.next = 85;
            break;
          }

          // OUTER CONTRACTION
          ;

          for (_i9 = N; _i9-- > 0;) {
            x[_i9] = centroid[_i9] + contractOuter * (centroid[_i9] - verts_x[_i9]);
          }

          f_contract = f(x);

          if (!(f_contract < f_worst)) {
            _context.next = 75;
            break;
          }

          for (_i10 = N; _i10-- > 0;) {
            verts_x[_i10] = x[_i10];
          }

          verts_f[0] = f_contract;

          if (!(f_contract < f_best)) {
            _context.next = 74;
            break;
          }

          _context.next = 74;
          return [new _nd_array.NDArray(shape_x, verts_x.slice(0, N)), verts_f[0]];

        case 74:
          return _context.abrupt("continue", 107);

        case 75:
          // INNER CONTRACTION
          ;

          for (_i11 = N; _i11-- > 0;) {
            x[_i11] = centroid[_i11] - contractInner * (centroid[_i11] - verts_x[_i11]);
          }

          _f_contract = f(x);

          if (!(_f_contract < f_worst)) {
            _context.next = 85;
            break;
          }

          for (_i12 = N; _i12-- > 0;) {
            verts_x[_i12] = x[_i12];
          }

          verts_f[0] = _f_contract;

          if (!(_f_contract < f_best)) {
            _context.next = 84;
            break;
          }

          _context.next = 84;
          return [new _nd_array.NDArray(shape_x, verts_x.slice(0, N)), verts_f[0]];

        case 84:
          return _context.abrupt("continue", 107);

        case 85:
          stuck = true; // SHRINK

          _i13 = N + 1;

        case 87:
          if (!(_i13-- > 0)) {
            _context.next = 101;
            break;
          }

          if (!(_i13 !== _b)) {
            _context.next = 99;
            break;
          }

          _j4 = N;

        case 90:
          if (!(_j4-- > 0)) {
            _context.next = 98;
            break;
          }

          _x_ij = verts_x[N * _i13 + _j4] + (1 - shrink) * (verts_x[N * _b + _j4] - verts_x[N * _i13 + _j4]);
          stuck = stuck && verts_x[N * _i13 + _j4] === _x_ij;
          verts_x[N * _i13 + _j4] = _x_ij;

          if (isFinite(_x_ij)) {
            _context.next = 96;
            break;
          }

          throw new Error('Assertion failed.');

        case 96:
          _context.next = 90;
          break;

        case 98:
          verts_f[_i13] = f(verts_x, N * _i13);

        case 99:
          _context.next = 87;
          break;

        case 101:
          if (!stuck) {
            _context.next = 103;
            break;
          }

          throw new _optimization_error.OptimizationNoProgressError();

        case 103:
          c = (0, _min_max.argmin)(verts_f);

          if (!(verts_f[c] < f_best)) {
            _context.next = 107;
            break;
          }

          _context.next = 107;
          return [new _nd_array.NDArray(shape_x, verts_x.slice(N * c, N * (c + 1))), verts_f[c]];

        case 107:
          _context.next = 40;
          break;

        case 109:
        case "end":
          return _context.stop();
      }
    }
  }, _marked);
}