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

import numpy as np, base64, io, scipy as sci, scipy.linalg as la


def to_b64( A ):
  with io.BytesIO() as bytes:
    assert A.dtype != '<i8'
    np.save(bytes, A, allow_pickle=False)
    bytes.seek(0)
    bytes = bytes.read()
    return base64.encodebytes(bytes).decode('ASCII').replace('\n','')


def main():
  print("import {b64_decode_gen}  from '../io/b64'")
  print("import {npy_deserialize} from '../io/npy'")
  print()
  print()
  print('export function* pldlp_decomp_test_data()')
  print('{')

  for n in range(1,17):
    for repeat in range(16):
      A = np.random.rand(n,n)*4 - 2
      A = (A + A.T) / 2

      L,D,P = la.ldl(A)
      L = L[P,:]
      assert np.all( np.tril(L) == L ), L
      P = P.astype( np.dtype('<i4', align=False) ) 

      if 0 < repeat:
        print()

      print('  yield [')
      for a in [A,L,D,P]:
        print("    '%s'," % (to_b64(a),))
      print('  ].map(b64_decode_gen).map(npy_deserialize);')

  print('}')


if '__main__' == __name__:
  main()

