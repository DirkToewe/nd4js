'strict mode';

const {performance} = require('perf_hooks');

const nd = require('./lib');

const {TrustRegionSolverODR   } = require('./lib/opt/_trust_region_solver_odr'   ),
      {TrustRegionSolverODR_NY} = require('./lib/opt/_trust_region_solver_odr_ny');
const { solve } = require('./lib/la/solve');

console.log({
  TrustRegionSolverODR,
  TrustRegionSolverODR_NY
});

console.log('  MX  |  NX  |  NP  |    ODR    |    ONY    |  RATIO');
console.log('------|------|------|-----------|-----------|---------');

// for( const MX of [...function*(){ for( let MX=1; MX < 8192; MX = Math.max(MX+1, MX*1.25 | 0) ) yield MX }()].reverse() )
for( let MX=1; MX < 8192          ; MX = Math.max(MX+1, MX*1.25 | 0) )
for( let NP=1; NP < 64 && NP <= MX; NP = Math.max(NP+1, NP*1.25 | 0) )
for( let NX=1; NX < 64            ; NX = Math.max(NX+1, NX*1.25 | 0) )
{
  const fgg = (p,x) => {
    return [
      nd.tabulate([MX],    'float64', () => 0),
      nd.tabulate([MX,NP], 'float64', () => 0),
      nd.tabulate([MX,NX], 'float64', () => 0)
    ];
  };

  const [odr,ony] = function(){
    const p0 = nd.tabulate(   [NP], 'float64', () => 0),
         dx0 = nd.tabulate([MX,NX], 'float64', () => 0);

    return [
      new TrustRegionSolverODR   (fgg, p0,dx0),
      new TrustRegionSolverODR_NY(fgg, p0,dx0)
    ]
  }();

  const {M,N} = odr;

  let nTest = 128,
      tOdr  = 0.0,
      tOny  = 0.0;

  for( let repeat=0; repeat++ < nTest; )
  {
    for( let i=MX*NX; i-- > 0; ) odr.J11[i] = ony.J11[i] = Math.random()*1.8 + 0.1;
    for( let i=MX*NX; i-- > 0; ) odr.J21[i] = ony.J21[i] = Math.random()*4 - 2;
    for( let i=MX*NP; i-- > 0; ) odr.J22[i] = ony.J22[i] = Math.random()*4 - 2;
    for( let i=M    ; i-- > 0; ) odr. F0[i] = ony. F0[i] = Math.random()*4 - 2;

    odr.rank = -1;
    ony.rank = -1;
    ony.prepared = false;

    let t0 = performance.now();
     odr.computeNewton();
    tOdr += performance.now() - t0;

       t0 = performance.now();
     ony.computeNewton();
    tOny += performance.now() - t0;

    const TOL = nd.la.norm(odr.newton_dX) * Number.EPSILON;

    for( let i=N; i-- > 0; )
    {
      const delta = Math.abs(
        odr.newton_dX[i] -
        ony.newton_dX[i]
      );

      if( ! (delta <= TOL) )
        throw new Error('Assertion failed.');
    }
  }

  console.log( MX.toString().padStart(5)
      + ' |' + NX.toString().padStart(5)
      + ' |' + NP.toString().padStart(5)
      + ' |' + (tOdr/nTest).toFixed(3).padStart(8)+'ms'
      + ' |' + (tOny/nTest).toFixed(3).padStart(8)+'ms'
      + ' |' + (tOny/tOdr ).toFixed(3).padStart(8)
  );
}