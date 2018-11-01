const tf = require('@tensorflow/tfjs')
const fs = require('fs')

const { performance } = require('perf_hooks');

const bandPartWithBackend = ( a, numLower, numUpper ) =>
{
  if( numLower%1 !== 0 ) throw new Error(`bandPart(): numLower=${numLower} is no integer.`);
  if( numUpper%1 !== 0 ) throw new Error(`bandPart(): numUpper=${numUpper} is no integer.`);
  if( a.rank < 2 ) throw new Error(`bandPart(): a.rank = ${a.rank} < 2.`)

  if( numLower < 0 ) numLower = a.shape[a.rank-2];
  if( numUpper < 0 ) numUpper = a.shape[a.rank-1];

  const B_shape = Array.from(a.shape),
       [M,N]   = a.shape.slice(-2),
        dtype  = a.dtype,
        wordsPerElem = dtype.startsWith('complex') ? 2 : 1,
        B = a.dataSync().slice(); a = undefined;

  Object.freeze(B_shape);

  for( let off=0; off < B.length; off += M*N )
  for( let i=0; i < M; i++ )
  for( let j=0; j < N; j++ )
  if( j-i > numUpper ||
      i-j > numLower )
    for( let k=0; k < wordsPerElem; k++ )
      B[wordsPerElem*(off + N*i+j) + k] = 0;

  return tf.Tensor.make(B_shape,{values: B},dtype);
}

const bandPartWithoutBackend = ( a, numLower, numUpper ) =>
{
  if( numLower%1 !== 0 ) throw new Error(`bandPart(): numLower=${numLower} is no integer.`);
  if( numUpper%1 !== 0 ) throw new Error(`bandPart(): numUpper=${numUpper} is no integer.`);
  if( a.rank < 2 ) throw new Error(`bandPart(): a.rank = ${a.rank} < 2.`)

  const [M,N] = a.shape.slice(-2);

  if( numLower < 0 ) numLower = M;
  if( numUpper < 0 ) numUpper = N;

  return tf.tidy( () => {
    numLower = tf.scalar(numLower, a.dtype);
    numUpper = tf.scalar(numUpper, a.dtype);
  
    const i = tf.range(0, M, 1, a.dtype).reshape([-1,1]),
          j = tf.range(0, N, 1, a.dtype);
  
    const inBand = tf.logicalAnd(
      tf.sub(i,j).lessEqual(numLower),
      tf.sub(j,i).lessEqual(numUpper)
    );
  
    return tf.mul(a, inBand.cast(a.dtype) );
  });
}

/** Quick benchmark of TFJS 0.13.3.
  */
async function main()
{
  tf.setBackend('cpu', true);

  for( const [name,bandPart] of [
    ['withBackend',   bandPartWithBackend   ],
    ['withoutBackend',bandPartWithoutBackend]
  ])
  {
    console.log(`Testing '${name}'...`)

    const Ns = [],
          Ts = [];

    for( let N=0; N++ < 160; )
    {
      const a = tf.randomUniform([N,N,N],-1,+1);

      const t0 = performance.now();

      tf.tidy( () => {
        const b = bandPart(a,-1,0);
        b.dataSync();
      });

      const dt = performance.now() - t0;

      if( N%16 == 0 )
        console.log(`[${name}] N =${N.toString().padStart(4)}: ${(dt/1000).toFixed(3)}sec`);

      Ns.push(N);
      Ts.push(dt/1000);
    }

    const results = JSON.stringify({
      n : Ns,
      dt: Ts
    });

    await new Promise( f => fs.writeFile('benchmarks/tf_bandPart_'+name+'.json', results, 'utf-8', f ) );
  }

  console.log('DONE!');
}
main().catch( err => console.error(err) );
