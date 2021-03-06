const nd = require('nd4js@1.3.0') // <- 1st line is auto-generated

const a = nd.array(
  [[10,   0,  0],
   [ 0,   0,100],
   [ 0,1000,  0]]
)
console.log('A =\n' + a)

const a2 = a.mapElems(x => x*2)
console.log('2a = \n' + a2)

const v = nd.tabulate([3,1], (i,j) => i+1)
console.log('v =\n' + v)

const apv = nd.zip_elems([a,v], (x,y) => x+y)
console.log('A + v = \n' + apv)

const av = nd.la.matmul(a,v)
console.log('A v =\n' + av)
