'use strict'

const nd = require('./nd')

const A = nd.array(
  [[1,2],
   [3,4]]
)
console.log('ND:\n'+A.toString())

console.log('3+4:', nd.math.add(3,4))
nd.math.add = (x,y) => x*y
console.log('3*4:', nd.math.add(3,4))

console.log('nd.dt:', nd.dt)
console.log('nd.la:', nd.la)
