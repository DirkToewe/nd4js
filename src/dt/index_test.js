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

import {Complex} from './complex'
import {is_subdtype, super_dtype, dtypeof, eps} from '.'

describe('dt', () => {
  it("is_subdtype adheres to int32 < float32 < float64 < complex128 < object", () => {
    expect( is_subdtype(     'int32','object') ).toBe(true);
    expect( is_subdtype(   'float32','object') ).toBe(true);
    expect( is_subdtype(   'float64','object') ).toBe(true);
    expect( is_subdtype('complex128','object') ).toBe(true);
    expect( is_subdtype(    'object','object') ).toBe(true);

    expect( is_subdtype(     'int32','complex128') ).toBe(true);
    expect( is_subdtype(   'float32','complex128') ).toBe(true);
    expect( is_subdtype(   'float64','complex128') ).toBe(true);
    expect( is_subdtype('complex128','complex128') ).toBe(true);
    expect( is_subdtype(    'object','complex128') ).toBe(false);

    expect( is_subdtype(     'int32','float64') ).toBe(true);
    expect( is_subdtype(   'float32','float64') ).toBe(true);
    expect( is_subdtype(   'float64','float64') ).toBe(true);
    expect( is_subdtype('complex128','float64') ).toBe(false);
    expect( is_subdtype(    'object','float64') ).toBe(false);

    expect( is_subdtype(     'int32','float32') ).toBe(true);
    expect( is_subdtype(   'float32','float32') ).toBe(true);
    expect( is_subdtype(   'float64','float32') ).toBe(false);
    expect( is_subdtype('complex128','float32') ).toBe(false);
    expect( is_subdtype(    'object','float32') ).toBe(false);

    expect( is_subdtype(     'int32','int32') ).toBe(true);
    expect( is_subdtype(   'float32','int32') ).toBe(false);
    expect( is_subdtype(   'float64','int32') ).toBe(false);
    expect( is_subdtype('complex128','int32') ).toBe(false);
    expect( is_subdtype(    'object','int32') ).toBe(false);
  });

  it("super_dtype adheres to int32 < float32 < float64 < complex128 < object", () => {
    expect( super_dtype('int32', 'int32') ).toBe('int32');
  
    expect( super_dtype('int32',  'float32') ).toBe('float32');
    expect( super_dtype('float32','float32') ).toBe('float32');
    expect( super_dtype('float32',  'int32') ).toBe('float32');
  
    expect( super_dtype('int32',  'float64') ).toBe('float64');
    expect( super_dtype('float32','float64') ).toBe('float64');
    expect( super_dtype('float64','float64') ).toBe('float64');
    expect( super_dtype('float64',  'int32') ).toBe('float64');
    expect( super_dtype('float64','float32') ).toBe('float64');
  
    expect( super_dtype(     'int32','complex128') ).toBe('complex128');
    expect( super_dtype(   'float32','complex128') ).toBe('complex128');
    expect( super_dtype(   'float64','complex128') ).toBe('complex128');
    expect( super_dtype('complex128','complex128') ).toBe('complex128');
    expect( super_dtype('complex128',     'int32') ).toBe('complex128');
    expect( super_dtype('complex128',   'float32') ).toBe('complex128');
    expect( super_dtype('complex128',   'float64') ).toBe('complex128');
  
    expect( super_dtype(     'int32','object') ).toBe('object');
    expect( super_dtype(   'float32','object') ).toBe('object');
    expect( super_dtype(   'float64','object') ).toBe('object');
    expect( super_dtype('complex128','object') ).toBe('object');
    expect( super_dtype('object',    'object') ).toBe('object');
    expect( super_dtype('object',     'int32') ).toBe('object');
    expect( super_dtype('object',   'float32') ).toBe('object');
    expect( super_dtype('object',   'float64') ).toBe('object');
    expect( super_dtype('object','complex128') ).toBe('object');
  });

  it("dtypeof returns the proper dtype", () => {
    expect( dtypeof( 1337)              ).toBe(     'int32');
    expect( dtypeof(new Complex(1337))  ).toBe(     'int32');
    expect( dtypeof(1.337)              ).toBe(   'float64');
    expect( dtypeof(new Complex(1.337)) ).toBe(   'float64');
    expect( dtypeof(new Complex(1,2))   ).toBe('complex128');
    expect( dtypeof('str')              ).toBe(    'object');
    expect( dtypeof({x:2})              ).toBe(    'object');
    expect( dtypeof([1,2])              ).toBe(    'object');
  });

  it('eps works correctly', () => {
    expect(             eps('float64')   + 1  ).toBeGreaterThan(1)
    expect(             eps('float64')/2 + 1  ).toBe(1)
    expect( Math.fround(eps('float32')   + 1) ).toBeGreaterThan(1)
    expect( Math.fround(eps('float32')/2 + 1) ).toBe(1)
    expect(          eps('complex128')   + 1  ).toBeGreaterThan(1)
    expect(          eps('complex128')/2 + 1  ).toBe(1)
  })
});
