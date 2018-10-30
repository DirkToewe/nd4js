const nd = require('../nd.js');
const fs = require('fs');

const { performance } = require('perf_hooks');

/** 
  * 
  * @returns
  */
function main()
{
  const Ns = [],
        Ts = [];
  for( let N=1; N < 512; N++ )
  {
    const a = nd.tabulate([N,N], 'float64', () => Math.random()*2 - 1 );

    const t0 = performance.now();

    const [q,r] = nd.la.qr_decomp(a);

    const dt = performance.now() - t0;

    console.log(`N = ${N}: ${(dt/1000).toFixed(3)}sec`);

    Ns.push(N);
    Ts.push(dt);
  }

  const results = JSON.stringify({
    n : Ns,
    dt: Ts
  });
  fs.writeFile('benchmarks/nd_qr_results.json', results, 'utf-8', () => console.log('DONE!') );
}
main()
