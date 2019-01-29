'''
This file is part of ND.JS.

ND.JS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ND.JS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ND.JS. If not, see <http://www.gnu.org/licenses/>.
'''

import base64, io, numpy as np


def main():

  def shapes():
    yield tuple()
    for l in range(1,4):
      yield l,
      for m in range(1,4):
        yield l,m
        for n in range(1,4):
          yield l,m,n

  def test_cases():

    for shape in shapes():
      for dtype in ['<i4','>i4','<f4','>f4','<f8','>f8','<c16','>c16']:
        jstype = {
          'i4' :      'int32',
          'f4' :    'float32',
          'f8' :    'float64',
          'c16': 'complex128'
        }[dtype[1:]]
    
        for f_order in [False, True]:
          a = np.random.uniform(-1e3,+1e3,shape)
          if dtype[1] == 'c':
            a = a + 1j * np.random.uniform(-1e3,+1e3,shape)
    
          a = a.astype(dtype)
          if f_order and shape != ():
            a = np.asfortranarray(a)
    
          with io.BytesIO() as bytes:
            np.save(bytes, a)
            bytes.seek(0)
            bytes = bytes.read()
            bytes = base64.encodebytes(bytes).decode('ASCII').replace('\n','')
    
          toStr = str
          if dtype[1] == 'c':
            toStr = lambda x: 'new Complex(%s,%s)' % (x.real,x.imag)
    
          string = '''
            new NDArray( Int32Array.of({SHAPE}), ARRAY_TYPES["{DTYPE}"].of({DATA}) )
          '''.format(
            SHAPE = ','.join(   str(s) for s in shape ),
            DTYPE = jstype,
            DATA  = ','.join( toStr(d) for d in a.flatten() )
          ).strip()
    
          yield bytes,string

  with open('./npy_test_data.js', mode='w') as out:
    out.write('''
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

import {NDArray} from '../nd_array'
import {ARRAY_TYPES, Complex} from '../dt'


export function* npy_test_data()
{
''')

    for bytes,string in test_cases():
      out.write('  yield ["%s",\n         %s]\n' % (bytes,string))

    out.write('}\n')


if '__main__' == __name__:
  main()
