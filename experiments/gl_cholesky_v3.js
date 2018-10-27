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

/** This script contains a proof-of-concept implementation of the
  * Cholesky Decomposition in WebGL2.
  */

const gl = function(){
  const canvas = document.createElement('canvas');
  canvas.width  = 1;
  canvas.height = 1;
  const gl = canvas.getContext("webgl2",{
    antialias: false,
    stencil: false,
    alpha: false
  });

  console.log( 'GLSL Version:', gl.getParameter(gl.SHADING_LANGUAGE_VERSION) );

  if( ! gl.getExtension('EXT_color_buffer_float') )
    throw new Error('HDR rendering not supported.');

  return gl;
}();

const gl_LLT = function(){
  const VERT_SHADER = `\
    #version 300 es
    precision highp float;

    // INPUTS
    in vec2 pos;

    // OUTPUTS

    // SHADER
    void main() {
      gl_Position = vec4( pos, 0.0, 1.0 );
    }
  `;

  const FRAG_SHADER = `\
    #version 300 es
    precision highp float;

    // INPUTS
    uniform highp int size;
    uniform highp sampler2D matrix;

    // OUTPUTS
    out highp float newVal;

    float m( int i, int j )
    {
      return texelFetch( matrix, ivec2(j,i), /*lod=*/0 ).x;
    }
    
    void main()
    {
      int i = int(gl_FragCoord.y),
          j = int(gl_FragCoord.x);
      float sum = 0.0,
            rem = 0.0;

      for( int k=0; k <= min(i,j); k++ )
      { // https://en.wikipedia.org/wiki/Kahan_summation_algorithm#The_algorithm
        float r = m(i,k)*m(j,k) - rem,
              s = sum  + r;
        rem =(s - sum) - r;
        sum = s;
      }

      newVal = sum;
    }
  `;

   //
  // INIT
 //
  const mkShader = (code,type) => {
    const shader = gl.createShader(type);
    gl.shaderSource(shader, code);
    gl.compileShader(shader);

    const log = gl.getShaderInfoLog(shader);
    console.log(log);

    const compiled = gl.getShaderParameter(shader, gl.COMPILE_STATUS);

    if(compiled) return shader;

    gl.deleteShader(shader);
    throw new Error('Could not compile shader.');
  }

  const vertShader = mkShader(VERT_SHADER, gl.  VERTEX_SHADER);
  const fragShader = mkShader(FRAG_SHADER, gl.FRAGMENT_SHADER);

  const program = gl.createProgram();
  gl.attachShader(program, vertShader);
  gl.attachShader(program, fragShader);
  gl.linkProgram(program);
  gl. useProgram(program);

  const posLoc = gl. getAttribLocation(program, 'pos' ),
       sizeLoc = gl.getUniformLocation(program, 'size');

  const posBuf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
  gl.bufferData(gl.ARRAY_BUFFER, Float32Array.of(
   -1,+1,
   +1,+1,
   -1,-1,
   +1,-1
  ), gl.STATIC_DRAW);

  let inTex = gl.createTexture(),
     outTex = gl.createTexture();
  for( const tex of [inTex,outTex] )
  {
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  }

  const frameBuf = gl.createFramebuffer();
  gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

  return S => {
    if( S.ndim != 2 ) throw new Error('S is not a matrix.');
    const [M,N] = S.shape.slice(-2);
    if( M != N ) throw new Error('S must be quadratic.');
    if( S.dtype != 'float32' ) throw new Error("S.dtype must be 'float32'.");

    gl.useProgram(program);

    gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
    gl.enableVertexAttribArray(posLoc);
    gl.vertexAttribPointer(posLoc, 2, gl.FLOAT, false, 0, 0);

    gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

     //
    // UPLOAD
   //
    console.time('UPLOAD');

    gl.uniform1i(sizeLoc,N);

    gl.bindTexture(gl.TEXTURE_2D, inTex);
    gl.texImage2D(
      /*target=*/gl.TEXTURE_2D,
      /*levelOfDetail=*/0,
      /*internalFormat=*/gl.R32F,
      /*width,height=*/N,N,
      /*border=*/0,
      /*format=*/gl.RED,
      /*type=*/gl.FLOAT,
      /*srcData=*/S.data
    );
    gl.bindTexture(gl.TEXTURE_2D, outTex);
    gl.texImage2D(
      /*target=*/gl.TEXTURE_2D,
      /*levelOfDetail=*/0,
      /*internalFormat=*/gl.R32F,
      /*width,height=*/N,N,
      /*border=*/0,
      /*format=*/gl.RED,
      /*type=*/gl.FLOAT,
      /*srcData=*/null
    );
    gl.activeTexture(gl.TEXTURE0 + 0);
    gl.viewport(0,0,N,N);
    gl.drawBuffers([gl.COLOR_ATTACHMENT0]);

    gl.finish();
    console.timeEnd('UPLOAD');

     //
    // COMPUTE
   //
    console.time('COMPUTATION');

    // SET WRITE TO MATRIX 2
    gl.framebufferTexture2D(
      /*target=*/gl.DRAW_FRAMEBUFFER,
      /*attachment=*/gl.COLOR_ATTACHMENT0,
      /*texTarget=*/gl.TEXTURE_2D,
      /*texture=*/outTex,
      /*levelOfDetail=*/0
    );

    // SET READ FROM MATRIX 1
    gl.bindTexture(gl.TEXTURE_2D, inTex);

    // COMPUTE
    gl.drawArrays(gl.TRIANGLE_STRIP, 0,4);

    gl.finish();
    console.timeEnd('COMPUTATION');

     //
    // DOWNLOAD
   //
    let outArr = new Float32Array(4*N*N);
    console.time('DOWNLOAD')
    gl.readBuffer(gl.COLOR_ATTACHMENT0);
    gl.readPixels(
      /*x,y=*/0,0,
      /*w,h=*/N,N,
      /*format=*/gl.RGBA,
      /*type=*/gl.FLOAT,
      /*writeTo=*/outArr
    );
    console.timeEnd('DOWNLOAD')
    outArr = outArr.filter( (_,i) => i%4 == 0 );

    return new nd.Array(S.shape,outArr)
  };
}();

const gl_cholesky/*: nd.Array => nd.Array*/ = function() {
  const VERT_SHADER = `\
    #version 300 es
    precision highp float;

    // INPUTS
    in float cornerIndex;
    uniform sampler2D matrix;

    uniform float size,  // <- number of rows/columns in the matrix
//                  NaN,
                  offset;// <- index before column that's currently finalized

    // OUTPUTS
    flat out highp float m_kk;

    // SHADER
    void main() {
      m_kk = texelFetch( matrix, ivec2(offset+1.5), /*lod=*/0 ).x;
      m_kk = sqrt(m_kk);

      // 1
      //  |\
      //  |_\
      // 2   3
      vec2 index;
      if( 1.0 == cornerIndex ) { index = vec2(offset,   offset - 0.5); }
      if( 2.0 == cornerIndex ) { index = vec2(offset,     size + 1.0); }
      if( 3.0 == cornerIndex ) { index = vec2(size + 1.5, size + 1.0); }
      gl_Position = vec4( index/size*2.0 - 1.0, 0.0, 1.0 );
    }
  `;

  const FRAG_SHADER = `\
    #version 300 es
    precision highp float;

    // INPUTS
    flat in highp float m_kk;
    uniform float offset;
    uniform sampler2D matrix;

    // OUTPUTS
    out float newVal;

    float m( int i, int j )
    {
      return texelFetch( matrix, ivec2(j,i), /*lod=*/0 ).x;
    }
    
    void main()
    {
      int i = int(gl_FragCoord.y),
          j = int(gl_FragCoord.x),
          k = int(offset + 1.5);

      float m_ik = m(i,k) / m_kk,
            m_jk = m(j,k) / m_kk,
            m_ij = m(i,j);

           if( j <  k )  newVal = m_ij; // <- copy column left to the currently finalized column
      else if( j == k ) {
           if( i == k )  newVal = m_kk; // <- diagonal element
           else          newVal = m_ik; // <- elements below diagonal element
      }    else          newVal = m_ij - m_ik*m_jk; // <- remaining elements
    }
  `;

   //
  // INIT
 //
  const mkShader = (code,type) => {
    const shader = gl.createShader(type);
    gl.shaderSource(shader, code);
    gl.compileShader(shader);

    const log = gl.getShaderInfoLog(shader);
    console.log(log);

    const compiled = gl.getShaderParameter(shader, gl.COMPILE_STATUS);

    if(compiled) return shader;

    gl.deleteShader(shader);
    throw new Error('Could not compile shader.');
  }

  const vertShader = mkShader(VERT_SHADER, gl.  VERTEX_SHADER);
  const fragShader = mkShader(FRAG_SHADER, gl.FRAGMENT_SHADER);

  const program = gl.createProgram();
  gl.attachShader(program, vertShader);
  gl.attachShader(program, fragShader);
  gl.linkProgram(program);
  gl.useProgram(program);

  const cornerLoc = gl. getAttribLocation(program, 'cornerIndex'),
        offsetLoc = gl.getUniformLocation(program, 'offset'),
          sizeLoc = gl.getUniformLocation(program, 'size'  ),
           nanLoc = gl.getUniformLocation(program, 'NaN'  );

  const cornerBuf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, cornerBuf);
  gl.bufferData(gl.ARRAY_BUFFER, Float32Array.of(1,2,3), gl.STATIC_DRAW);

  let matrix1Tex = gl.createTexture(),
      matrix2Tex = gl.createTexture();
  for( const tex of [matrix1Tex,matrix2Tex] )
  {
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  }

  const frameBuf = gl.createFramebuffer();
  gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);
  gl.uniform1f(nanLoc, NaN);

  return S => {
    if( S.ndim != 2 ) throw new Error('S is not a matrix.');
    const [M,N] = S.shape.slice(-2);
    if( M != N ) throw new Error('S must be quadratic.');
    if( S.dtype != 'float32' ) throw new Error("S.dtype must be 'float32'.");

    gl.useProgram(program);

    gl.bindBuffer(gl.ARRAY_BUFFER, cornerBuf);
    gl.enableVertexAttribArray(cornerBuf);
    gl.vertexAttribPointer(cornerBuf, 1, gl.FLOAT, false, 0, 0);

    gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

     //
    // UPLOAD
   //
    console.time('upload');

    gl.uniform1f(sizeLoc,N);

    gl.bindTexture(gl.TEXTURE_2D, matrix1Tex);
    gl.texImage2D(
      /*target=*/gl.TEXTURE_2D,
      /*levelOfDetail=*/0,
      /*internalFormat=*/gl.R32F,
      /*width,height=*/N,N,
      /*border=*/0,
      /*format=*/gl.RED,
      /*type=*/gl.FLOAT,
      /*srcData=*/S.data
    );
    gl.bindTexture(gl.TEXTURE_2D, matrix2Tex);
    gl.texImage2D(
      /*target=*/gl.TEXTURE_2D,
      /*levelOfDetail=*/0,
      /*internalFormat=*/gl.R32F,
      /*width,height=*/N,N,
      /*border=*/0,
      /*format=*/gl.RED,
      /*type=*/gl.FLOAT,
      /*srcData=*/null
    );
    gl.activeTexture(gl.TEXTURE0 + 0);
    gl.viewport(0,0,N,N);
    gl.drawBuffers([gl.COLOR_ATTACHMENT0]);

    gl.finish();
    console.timeEnd('upload');

     //
    // COMPUTE
   //
    console.time('computation');

    for( let off = -1; off < N-1; off++ )
    {
      // SET OFFSET
      gl.uniform1f(offsetLoc,off);

      // SET WRITE TO MATRIX 2
      gl.framebufferTexture2D(
        /*target=*/gl.DRAW_FRAMEBUFFER,
        /*attachment=*/gl.COLOR_ATTACHMENT0,
        /*texTarget=*/gl.TEXTURE_2D,
        /*texture=*/matrix2Tex,
        /*levelOfDetail=*/0
      );

      // SET READ FROM MATRIX 1
      gl.bindTexture(gl.TEXTURE_2D, matrix1Tex);

      // COMPUTE
      gl.drawArrays(gl.TRIANGLES, 0,3);

      [matrix1Tex,matrix2Tex] = [matrix2Tex,matrix1Tex];
    }

    gl.finish();
    console.timeEnd('computation');

     //
    // DOWNLOAD
   //
    let matrixOut = new Float32Array(4*N*N);
    console.time('download')
    gl.readBuffer(gl.COLOR_ATTACHMENT0);
    gl.readPixels(
      /*x,y=*/0,0,
      /*w,h=*/N,N,
      /*format=*/gl.RGBA,
      /*type=*/gl.FLOAT,
      /*writeTo=*/matrixOut
    );
    console.timeEnd('download')
    matrixOut = matrixOut.filter( (_,i) => i%4 == 0 );

    return new nd.Array(S.shape,matrixOut)
  };
}();

function main()
{
  const sleep = dt => new Promise( f => setTimeout(f,dt) );

  const is_close = (x,y) => {
    const atol = 1e-2,
          rtol = 1e-2,
           tol = atol + rtol * Math.max(
            Math.abs(x),
            Math.abs(y)
          );
    return Math.abs(x-y) <= tol;
  };

  async function test_LLT()
  {
    for( let run=0; ++run < 1024; )
    {
      await sleep();

      console.log('\nRUN:',run);
      let
        randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from,
        shape = 128;
      shape = [shape,shape]; 
      const L = nd.tabulate(shape,'float32',(...indices) => {
        const [i,j] = indices.slice(-2);
        if( i==j ) return Math.random()*1.0 + 0.5;
        if( i< j ) return 0;
        return Math.random()*0.2 - 0.1;
      });
      let LLT = nd.la.matmul2(L,L.T);
      LLT = nd.Array.from(LLT,'float32')

      const llT = gl_LLT(L);

      nd.Array.from([LLT,llT], 'float32', (x,y,...indices) => {
        if( ! is_close(x,y) ) {
          let msg = '{\n'+LLT+'\n} expected but {\n'+llT+'\n} encountered.\n'+x+' != '+y+' at index ['+indices+']';
          throw new Error(msg);
        }
      });

    }
  }
//  test_LLT();


  async function test_cholesky()
  {
    for( let run=0; ++run < 1024; )
    {
      await sleep();

      console.log('\nRUN:',run);
      let
        randInt = (from,until) => Math.floor(Math.random()*(until-from)) + from,
        shape = 2*1024;
      shape = [shape,shape]; 
      const L = nd.tabulate(shape,'float32',(...indices) => {
        const [i,j] = indices.slice(-2);
        if( i==j ) return Math.random()*1.0 + 0.5;
        if( i< j ) return 0;
        return Math.random()*0.2 - 0.1;
      });
      let LLT = gl_LLT(L);

      const label = `N = ${shape[0].toString().padStart(4)} `;
      console.time(label);
      let l = gl_cholesky(LLT);
      console.timeEnd(label); 

      l = nd.la.tril(l);

      nd.Array.from([L,l], 'float32', (x,y,...indices) => {
        if( ! is_close(x,y) ) {
          let msg = '{\n'+L+'\n} expected but {\n'+l+'\n} encountered.\n'+x+' != '+y+' at index ['+indices+']';
          throw new Error(msg);
        }
      });
    }
  }
  test_cholesky();
}

main();
