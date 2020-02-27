'''
This file is part of ND4JS.

ND4JS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ND4JS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ND4JS. If not, see <http://www.gnu.org/licenses/>.
'''

import numpy as np
import tensorflow as tf
import tensorflow.linalg as tfla


def solve_test_case( J, f ):
  assert J.ndim == 2
  assert f.ndim == 2

  assert J.shape[0] == f.shape[0]
  assert         1  == f.shape[1]

  with tf.Graph().as_default():
    J = tf.constant(name='J', value=J)
#     JT= tf.transpose(J, name='JT')
    f = tf.constant(name='f', value=f)
    d = tf.norm(J, name='d', axis=0)
    D = tf.diag(d, name='D')

    JD = J / d
  
    λ = tf.placeholder(name='lambda', dtype=tf.float64, shape=[])

    Dx = tfla.lstsq(JD, -f, l2_regularizer=λ)
  
    r = tfla.norm(Dx)
    dr, = tf.gradients(r, λ)
  
    with tf.Session() as sess:
      inputs = {λ: 0.0}
      return sess.run([r,dr], feed_dict=inputs)


def generate_test_data():  
  np.random.seed(1337)

  overdet = []
  underdet = []

  for m in range(1,17):
    for n in range(1,17):
      J = np.random.rand(m,n)*4-2
      f = np.random.rand(m,1)*4-2

      r,dr = solve_test_case(J,f)

      if m >= n:
        overdet.append([J,f, r,dr])
  
      if m <= n:
        underdet.append([J,f, r,dr])

  print('\n'*7)

  print('''export function* computeMinGlobal_overdet_gen()
{''')
  for [J,f, r,dr] in overdet:
    print('''  yield [
    new NDArray( Int32Array.of({m}),   Float64Array.of({}) ),
    new NDArray( Int32Array.of({m},{n}), Float64Array.of({}) ),
    {},
    {}
  ];'''.format(
      ','.join( str(x) for x in f.flatten() ),
      ','.join( str(x) for x in J.flatten() ),
      r, dr,
      m=J.shape[0],
      n=J.shape[1]
    ))

  print('''}

export function* computeMinGlobal_underdet_gen()
{''')
  for [J,f, r,dr] in underdet:
    print('''  yield [
    new NDArray( Int32Array.of({m}),   Float64Array.of({}) ),
    new NDArray( Int32Array.of({m},{n}), Float64Array.of({}) ),
    {},
    {}
  ];'''.format(
      ','.join( str(x) for x in f.flatten() ),
      ','.join( str(x) for x in J.flatten() ),
      r, dr,
      m=J.shape[0],
      n=J.shape[1]
    ))
  print('}')

if __name__ == '__main__':
  generate_test_data()
