const tf = require('@tensorflow/tfjs')
const fs = require('fs')

const { performance } = require('perf_hooks');

/** Quick benchmark of TFJS 0.13.3.
  */
function main()
{
  tf.setBackend('cpu', false);

  const Ns = [],
        Ts = [];
  for( let N=1; N < 512; N++ )
    tf.tidy( () => {
      const a = tf.randomUniform([N,N],-1,+1);
  
      const t0 = performance.now();
  
      tf.tidy( () => {
        const [q,r] = tf.linalg.qr(a);
        q.dataSync();
        r.dataSync();
      });
  
      const dt = performance.now() - t0;
  
      console.log(`N = ${N}: ${(dt/1000).toFixed(3)}sec`);
  
      Ns.push(N);
      Ts.push(dt);
    });

  const results = JSON.stringify({
    n : Ns,
    dt: Ts
  });
  fs.writeFile('benchmarks/tf_qr_results.json', results, 'utf-8', () => console.log('DONE!') );
}
main()
