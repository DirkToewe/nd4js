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

  gl.disable(gl.RASTERIZER_DISCARD);
  gl.disable(gl.DEPTH_TEST);
  gl.disable(gl.CULL_FACE);
  gl.disable(gl.DITHER);
  gl.disable(gl.BLEND);

  console.log( 'GLSL Version:', gl.getParameter(gl.SHADING_LANGUAGE_VERSION) );

  if( ! gl.getExtension('EXT_color_buffer_float') )
    throw new Error('HDR rendering not supported.');

  return gl;
}();


const mkProgram = (...shaderCode) => {

  if( shaderCode.length != 2 )
    throw new Error('mkProgram(): exactly two shader code strings expected.');
    
  const program = gl.createProgram(),
        shaders = [];

  try {
    for( let i=0; i < shaderCode.length; i++ )
    {
      const shaderType = i == 0 ? gl.  VERTEX_SHADER
                                : gl.FRAGMENT_SHADER;

      const shader = gl.createShader(shaderType);
      shaders.push(shader);
      gl.shaderSource(shader, shaderCode[i]);
      gl.compileShader(shader);
      gl.attachShader(program, shader);

      const log = gl.getShaderInfoLog(shader);
      console.log(log);

      const compiled = gl.getShaderParameter(shader, gl.COMPILE_STATUS);

      if(compiled) continue;

      const msg = gl.getShaderInfoLog(shader);

      const lines = shaderCode[i].split(/\n/);

      for( let i=0; i < lines.length; i++ ) {
        lines[i] = (i+1).toString().padStart(4) + ':  ' + lines[i];
      }

      throw new Error('Could not compile shader:\n' + msg + '\nSource:\n' + lines.join('\n') );
    }

    gl.linkProgram(program);

    const linked = gl.getProgramParameter(program, gl.LINK_STATUS);

    if(linked) return program;

    var msg = gl.getProgramInfoLog(program);

    throw new Error('Could not link program:\n' + msg);
  }
  catch(err) {
    for( const shader of shaders )
      gl.deleteShader(shader);
    gl.deleteProgram(program);
    throw err;
  }
};


const gl_matmul = function(){

  const program = mkProgram(
     //
    // VERTEX SHADER
   //
    `#version 300 es
      precision highp float;
      precision highp int;

      // INPUTS
      in vec2 pos;

      // OUTPUTS

      // SHADER
      void main() {
        gl_Position = vec4( pos, 0.0, 1.0 );
      }
    `,
     //
    // MATMUL SHADER
   //
    `#version 300 es
      precision highp float;
      precision highp int;

      // INPUTS
      uniform highp int innerSize;
      uniform highp sampler2D matrix_A;
      uniform highp sampler2D matrix_B;

      // OUTPUTS
      out highp vec4 rgba;

      highp float A( int i, int j ) { return texelFetch( matrix_A, ivec2(j,i), /*lod=*/0 ).x; }
      highp float B( int i, int j ) { return texelFetch( matrix_B, ivec2(j,i), /*lod=*/0 ).x; }

      highp vec4 float_to_rgba( highp float val )
      {
        // ZERO
        if( 0.0 == val )
          return vec4( vec3(0), (1.0/val < 0.0) ? 128.0/255.0 : 0.0 );

        // NaN
        if( isnan(val) )
          return vec4( vec3(255), 127 ) / 255.0;

        highp int sign = val < 0.0 ? 128 : 0;

        // +- INFINITY
        if( isinf(val) )
          return vec4(0, 0, 128, float(sign | 127)) / 255.0;

        val = abs(val);

        highp float ex = floor(log2(val));

        highp int exBits = int(ex) + 127;

        val *= exp2( 23.0 - ex );

        highp int mantissa = int(val) - (1<<23);

        highp ivec4 bytes = ivec4(
          mantissa>> 0 & 0xFF,
          mantissa>> 8 & 0xFF,
          mantissa>>16 | 0xFF & exBits<<7,
                  sign |        exBits>>1
        );

        return vec4(bytes) / 255.0;
      }
      
      void main()
      {
        highp int i = int(gl_FragCoord.y),
                  j = int(gl_FragCoord.x);

        // https://en.wikipedia.org/wiki/Kahan_summation_algorithm#The_algorithm
        highp float C_ij = 0.0,
                    c_ij = 0.0;

        for( int k=0; k < innerSize; k++ )
        {
          highp float s = A(i,k)*B(k,j) - c_ij,
                      S = C_ij + s;
          c_ij = (S - C_ij) - s;
          C_ij =  S;
        }

        rgba = float_to_rgba(C_ij);
      }
    `
  );

   //
  // INIT
 //
  gl.useProgram(program);

  const   pos_loc = gl. getAttribLocation(program, 'pos'),
    innerSize_loc = gl.getUniformLocation(program, 'innerSize'),
     matrix_A_loc = gl.getUniformLocation(program, 'matrix_A'),
     matrix_B_loc = gl.getUniformLocation(program, 'matrix_B');

  const posBuf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
  gl.bufferData(gl.ARRAY_BUFFER, Float32Array.of(
   -1,+1,
   +1,+1,
   -1,-1,
   +1,-1
  ), gl.STATIC_DRAW);

  const [matrix_A_tex,
         matrix_B_tex,
         matrix_C_tex] = function*(){
    for( let i=3; i-- > 0; )
    {
      const tex = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, tex);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      yield tex;
    }
  }();

  const frameBuf = gl.createFramebuffer();
  gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

  return (A,B, timeCallback) => {
    if( A.ndim != 2 ) throw new Error('A is not a matrix.');
    if( B.ndim != 2 ) throw new Error('B is not a matrix.');
    if( A.dtype != 'float32' ) throw new Error("A.dtype must be 'float32'.");
    if( B.dtype != 'float32' ) throw new Error("B.dtype must be 'float32'.");

    const [I,K] = A.shape.slice(-2),
          [L,J] = B.shape.slice(-2);
    if( K != L ) throw new Error('A.shape[-1] != B.shape[-2]');

    gl.useProgram(program);

    gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
    gl.enableVertexAttribArray(pos_loc);
    gl.vertexAttribPointer(pos_loc, 2, gl.FLOAT, false, 0, 0);

    gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

     //
    // UPLOAD
   //
//    console.time('UPLOAD');

    gl.bindTexture(gl.TEXTURE_2D, matrix_A_tex);
    gl.texImage2D(
      /*target=*/gl.TEXTURE_2D,
      /*levelOfDetail=*/0,
      /*internalFormat=*/gl.R32F,
      /*width,height=*/K,I,
      /*border=*/0,
      /*format=*/gl.RED,
      /*type=*/gl.FLOAT,
      /*srcData=*/A.data
    );
    gl.bindTexture(gl.TEXTURE_2D, matrix_B_tex);
    gl.texImage2D(
      /*target=*/gl.TEXTURE_2D,
      /*levelOfDetail=*/0,
      /*internalFormat=*/gl.R32F,
      /*width,height=*/J,K,
      /*border=*/0,
      /*format=*/gl.RED,
      /*type=*/gl.FLOAT,
      /*srcData=*/B.data
    );

    gl.finish();
    const t0 = performance.now();

    gl.bindTexture(gl.TEXTURE_2D, matrix_C_tex);
    gl.texImage2D(
      /*target=*/gl.TEXTURE_2D,
      /*levelOfDetail=*/0,
      /*internalFormat=*/gl.RGBA,
      /*width,height=*/J,I,
      /*border=*/0,
      /*format=*/gl.RGBA,
      /*type=*/gl.UNSIGNED_BYTE,
      /*srcData=*/null
    );
    gl.viewport(0,0,J,I);
    gl.drawBuffers([gl.COLOR_ATTACHMENT0]);

//    gl.finish();
//    console.timeEnd('UPLOAD');

     //
    // COMPUTE
   //
//    console.time('COMPUTATION');

    gl.uniform1i(innerSize_loc, K);
    gl.uniform1i(matrix_A_loc, 0); // texture unit 0
    gl.uniform1i(matrix_B_loc, 1); // texture unit 1

    // SET WRITE TO MATRIX C
    gl.framebufferTexture2D(
      /*target=*/gl.DRAW_FRAMEBUFFER,
      /*attachment=*/gl.COLOR_ATTACHMENT0,
      /*texTarget=*/gl.TEXTURE_2D,
      /*texture=*/matrix_C_tex,
      /*levelOfDetail=*/0
    );

    // READ FROM MATRICES A,B
    gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, matrix_A_tex);
    gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, matrix_B_tex);

    // COMPUTE
    gl.drawArrays(gl.TRIANGLE_STRIP, 0,4);

    if( null != timeCallback ) {
      gl.finish();
      const dt = performance.now() - t0;
      timeCallback(dt);
//    console.timeEnd('COMPUTATION');
    }

     //
    // DOWNLOAD
   //
    let outArr = new Float32Array(I*J);
//    console.time('DOWNLOAD')
    gl.readBuffer(gl.COLOR_ATTACHMENT0);
    gl.readPixels(
      /*x,y=*/0,0,
      /*w,h=*/J,I,
      /*format=*/gl.RGBA,
      /*type=*/gl.UNSIGNED_BYTE,
      /*writeTo=*/new Uint8Array(outArr.buffer)
    );
//    console.timeEnd('DOWNLOAD')

    return new nd.Array(Int32Array.of(I,J), outArr)
  };
}();


const gl_matmul_block1x4 = function(){

  const program = mkProgram(
     //
    // VERTEX SHADER
   //
    `#version 300 es
      precision highp float;
      precision highp int;

      // INPUTS
      in vec2 pos;

      // OUTPUTS

      // SHADER
      void main() {
        gl_Position = vec4( pos, 0.0, 1.0 );
      }
    `,
     //
    // MATMUL SHADER
   //
    `#version 300 es
      precision highp float;
      precision highp int;

      // INPUTS
      uniform highp int innerSize;
      uniform highp sampler2D matrix_A;
      uniform highp sampler2D matrix_B;

      // OUTPUTS
      out highp vec4 C_ij;

      highp vec4 A( int i, int j ) { return texelFetch( matrix_A, ivec2(j,i), /*lod=*/0 ); }
      highp mat4 B( int i, int j ) {
        return mat4(
          texelFetch( matrix_B, ivec2(j,i*4+0), /*lod=*/0 ),
          texelFetch( matrix_B, ivec2(j,i*4+1), /*lod=*/0 ),
          texelFetch( matrix_B, ivec2(j,i*4+2), /*lod=*/0 ),
          texelFetch( matrix_B, ivec2(j,i*4+3), /*lod=*/0 )
        );
      }
      
      void main()
      {
        highp int i = int(gl_FragCoord.y),
                  j = int(gl_FragCoord.x);

        // https://en.wikipedia.org/wiki/Kahan_summation_algorithm#The_algorithm
                   C_ij = vec4(0);
        highp vec4 c_ij = vec4(0);

        for( int k=0; k < innerSize; k++ )
        {
          highp vec4 s = B(k,j)*A(i,k) - c_ij,
                     S = C_ij + s;
          c_ij = (S - C_ij) - s;
          C_ij =  S;
        }
      }
    `
  );

   //
  // INIT
 //
  const   pos_loc = gl. getAttribLocation(program, 'pos'),
    innerSize_loc = gl.getUniformLocation(program, 'innerSize'),
     matrix_A_loc = gl.getUniformLocation(program, 'matrix_A'),
     matrix_B_loc = gl.getUniformLocation(program, 'matrix_B');

  const posBuf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
  gl.bufferData(gl.ARRAY_BUFFER, Float32Array.of(
   -1,+1,
   +1,+1,
   -1,-1,
   +1,-1
  ), gl.STATIC_DRAW);

  const [matrix_A_tex,
         matrix_B_tex,
         matrix_C_tex] = function*(){
    for( let i=3; i-- > 0; )
    {
      const tex = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, tex);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      yield tex;
    }
  }();

  const frameBuf = gl.createFramebuffer();

  return (A,B, timeCallback) => {
    if( A.ndim != 2 ) throw new Error('A is not a matrix.');
    if( B.ndim != 2 ) throw new Error('B is not a matrix.');
    if( A.dtype != 'float32' ) throw new Error("A.dtype must be 'float32'.");
    if( B.dtype != 'float32' ) throw new Error("B.dtype must be 'float32'.");

    const [I,K] = A.shape.slice(-2),
          [L,J] = B.shape.slice(-2);
    if( K != L ) throw new Error('A.shape[-1] != B.shape[-2]');

    if( K%4 != 0 ) throw new Error();
    if( J%4 != 0 ) throw new Error();

    gl.useProgram(program);

    gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
    gl.enableVertexAttribArray(pos_loc);
    gl.vertexAttribPointer(pos_loc, 2, gl.FLOAT, false, 0, 0);

    gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

     //
    // UPLOAD
   //
    gl.bindTexture(gl.TEXTURE_2D, matrix_A_tex);
    gl.texImage2D(
      /*target=*/gl.TEXTURE_2D,
      /*levelOfDetail=*/0,
      /*internalFormat=*/gl.RGBA32F,
      /*width,height=*/K>>>2, I,
      /*border=*/0,
      /*format=*/gl.RGBA,
      /*type=*/gl.FLOAT,
      /*srcData=*/A.data
    );
    gl.bindTexture(gl.TEXTURE_2D, matrix_B_tex);
    gl.texImage2D(
      /*target=*/gl.TEXTURE_2D,
      /*levelOfDetail=*/0,
      /*internalFormat=*/gl.RGBA32F,
      /*width,height=*/J>>>2, K,
      /*border=*/0,
      /*format=*/gl.RGBA,
      /*type=*/gl.FLOAT,
      /*srcData=*/B.data
    );

    gl.finish();
    const t0 = performance.now();

    gl.bindTexture(gl.TEXTURE_2D, matrix_C_tex);
    gl.texImage2D(
      /*target=*/gl.TEXTURE_2D,
      /*levelOfDetail=*/0,
      /*internalFormat=*/gl.RGBA32F,
      /*width,height=*/J>>>2,I,
      /*border=*/0,
      /*format=*/gl.RGBA,
      /*type=*/gl.FLOAT,
      /*srcData=*/null
    );
    gl.viewport(0,0,J>>>2,I);
    gl.drawBuffers([gl.COLOR_ATTACHMENT0]);

//    gl.finish();
//    console.timeEnd('UPLOAD');

     //
    // COMPUTE
   //
//    console.time('UPLOAD');

    gl.uniform1i(innerSize_loc, K>>>2);
    gl.uniform1i(matrix_A_loc, 0); // texture unit 0
    gl.uniform1i(matrix_B_loc, 1); // texture unit 1
//    console.time('COMPUTATION');

    // SET WRITE TO MATRIX C
    gl.framebufferTexture2D(
      /*target=*/gl.DRAW_FRAMEBUFFER,
      /*attachment=*/gl.COLOR_ATTACHMENT0,
      /*texTarget=*/gl.TEXTURE_2D,
      /*texture=*/matrix_C_tex,
      /*levelOfDetail=*/0
    );

    // READ FROM MATRICES A,B
    gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, matrix_A_tex);
    gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, matrix_B_tex);

    // COMPUTE
    gl.drawArrays(gl.TRIANGLE_STRIP, 0,4);

    if( null != timeCallback ) {
      gl.finish();
      const dt = performance.now() - t0;
      timeCallback(dt);
//    console.timeEnd('COMPUTATION');
    }

     //
    // DOWNLOAD
   //
    let outArr = new Float32Array(I*J);
//    console.time('DOWNLOAD')
    gl.readBuffer(gl.COLOR_ATTACHMENT0);
    gl.readPixels(
      /*x,y=*/0,0,
      /*w,h=*/J>>>2,I,
      /*format=*/gl.RGBA,
      /*type=*/gl.FLOAT,
      /*writeTo=*/outArr
    );
//    console.timeEnd('DOWNLOAD')

    return new nd.Array(Int32Array.of(I,J), outArr)
  };
}();


const gl_matmul_block4x4_v1 = function(){

  const vertShader = `\
    #version 300 es
    precision highp float;
    precision highp int;

    // INPUTS
    in vec2 pos;

    // OUTPUTS

    // SHADER
    void main() {
      gl_Position = vec4( pos, 0.0, 1.0 );
    }
  `;

  const program_mul = mkProgram(
    vertShader,
     //
    // MATMUL SHADER
   //
    `#version 300 es
      precision highp float;
      precision highp int;

      // INPUTS
      uniform highp int innerSize;
      uniform highp sampler2D matrix_A;
      uniform highp sampler2D matrix_B;

      // OUTPUTS
      out highp vec4 block4x4[4];

      highp mat4 read( highp sampler2D matrix, int i, int j ) {
         return mat4(
          texelFetch( matrix, ivec2(j,i*4+0), /*lod=*/0 ),
          texelFetch( matrix, ivec2(j,i*4+1), /*lod=*/0 ),
          texelFetch( matrix, ivec2(j,i*4+2), /*lod=*/0 ),
          texelFetch( matrix, ivec2(j,i*4+3), /*lod=*/0 )
        );
      }

      highp mat4 A( int i, int j ) { return read(matrix_A, i,j); }
      highp mat4 B( int i, int j ) { return read(matrix_B, i,j); }
      
      void main()
      {
        highp int i = int(gl_FragCoord.y),
                  j = int(gl_FragCoord.x);

        // https://en.wikipedia.org/wiki/Kahan_summation_algorithm#The_algorithm
        highp mat4 C_ij = mat4(0),
                   c_ij = mat4(0);

        for( int k=0; k < innerSize; k++ )
        {
          highp mat4 s = B(k,j)*A(i,k) - c_ij,
                     S = C_ij + s;
          c_ij = (S - C_ij) - s;
          C_ij =  S;
        }

        block4x4[0] = C_ij[0];
        block4x4[1] = C_ij[1];
        block4x4[2] = C_ij[2];
        block4x4[3] = C_ij[3];
      }
    `
  );

  const program_post = mkProgram(
    vertShader,
     //
    // MATMUL SHADER
   //
    `#version 300 es
      precision highp float;
      precision highp int;

      // INPUTS
      uniform highp sampler2D C0;
      uniform highp sampler2D C1;
      uniform highp sampler2D C2;
      uniform highp sampler2D C3;

      // OUTPUTS
      out highp vec4 C_ij;

      highp vec4 read( highp sampler2D matrix, int i, int j ) {
        return texelFetch( matrix, ivec2(j,i), /*lod=*/0 );
      }
      
      void main()
      {
        highp int i = int(gl_FragCoord.y),
                  j = int(gl_FragCoord.x),
                  k = i%4;
        i /= 4;
        switch(k) {
          case 0: C_ij = read(C0,i,j); break;
          case 1: C_ij = read(C1,i,j); break;
          case 2: C_ij = read(C2,i,j); break;
          case 3: C_ij = read(C3,i,j); break;
        }
      }
    `
  );

   //
  // INIT
 //
  const   pos_mul = gl. getAttribLocation(program_mul, 'pos'),
    innerSize_mul = gl.getUniformLocation(program_mul, 'innerSize'),
     matrix_A_mul = gl.getUniformLocation(program_mul, 'matrix_A'),
     matrix_B_mul = gl.getUniformLocation(program_mul, 'matrix_B');

  const   pos_post = gl. getAttribLocation(program_post, 'pos'),
    matrix_C0_post = gl.getUniformLocation(program_post, 'C0'),
    matrix_C1_post = gl.getUniformLocation(program_post, 'C1'),
    matrix_C2_post = gl.getUniformLocation(program_post, 'C2'),
    matrix_C3_post = gl.getUniformLocation(program_post, 'C3');

  const posBuf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
  gl.bufferData(gl.ARRAY_BUFFER, Float32Array.of(
   -1,+1,
   +1,+1,
   -1,-1,
   +1,-1
  ), gl.STATIC_DRAW);

  const [matrix_A_tex,
         matrix_B_tex,
         matrix_C0_tex,
         matrix_C1_tex,
         matrix_C2_tex,
         matrix_C3_tex,
         matrix_C_tex] = function*(){
    for( let i=7; i-- > 0; )
    {
      const tex = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, tex);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      yield tex;
    }
  }();

  const frameBuf = gl.createFramebuffer();
  gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

  return (A,B, timeCallback) => {
    if( A.ndim != 2 ) throw new Error('A is not a matrix.');
    if( B.ndim != 2 ) throw new Error('B is not a matrix.');
    if( A.dtype != 'float32' ) throw new Error("A.dtype must be 'float32'.");
    if( B.dtype != 'float32' ) throw new Error("B.dtype must be 'float32'.");

    const [I,K] = A.shape.slice(-2),
          [L,J] = B.shape.slice(-2);
    if( K != L ) throw new Error('A.shape[-1] != B.shape[-2]');

    if( I%4 != 0 ) throw new Error();
    if( K%4 != 0 ) throw new Error();
    if( J%4 != 0 ) throw new Error();

     //
    // UPLOAD
   //
//    console.time('UPLOAD');

    gl.bindTexture(gl.TEXTURE_2D, matrix_A_tex);
    gl.texImage2D(
      /*target=*/gl.TEXTURE_2D,
      /*levelOfDetail=*/0,
      /*internalFormat=*/gl.RGBA32F,
      /*width,height=*/K>>>2, I,
      /*border=*/0,
      /*format=*/gl.RGBA,
      /*type=*/gl.FLOAT,
      /*srcData=*/A.data
    );
    gl.bindTexture(gl.TEXTURE_2D, matrix_B_tex);
    gl.texImage2D(
      /*target=*/gl.TEXTURE_2D,
      /*levelOfDetail=*/0,
      /*internalFormat=*/gl.RGBA32F,
      /*width,height=*/J>>>2, K,
      /*border=*/0,
      /*format=*/gl.RGBA,
      /*type=*/gl.FLOAT,
      /*srcData=*/B.data
    );

    gl.finish();
    const t0 = performance.now();

    for( const matrix of [
      matrix_C0_tex,
      matrix_C1_tex,
      matrix_C2_tex,
      matrix_C3_tex
    ] )
    {
      gl.bindTexture(gl.TEXTURE_2D, matrix);
      gl.texImage2D(
        /*target=*/gl.TEXTURE_2D,
        /*levelOfDetail=*/0,
        /*internalFormat=*/gl.RGBA32F,
        /*width,height=*/J>>>2,I>>>2,
        /*border=*/0,
        /*format=*/gl.RGBA,
        /*type=*/gl.FLOAT,
        /*srcData=*/null
      );
    }
    gl.bindTexture(gl.TEXTURE_2D, matrix_C_tex);
    gl.texImage2D(
      /*target=*/gl.TEXTURE_2D,
      /*levelOfDetail=*/0,
      /*internalFormat=*/gl.RGBA32F,
      /*width,height=*/J>>>2,I,
      /*border=*/0,
      /*format=*/gl.RGBA,
      /*type=*/gl.FLOAT,
      /*srcData=*/null
    );

//    gl.finish();
//    console.timeEnd('UPLOAD');

     //
    // COMPUTE
   //
//    console.time('COMPUTATION');


    gl.useProgram(program_mul);

    gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

    gl.viewport(0,0,J>>>2,I>>>2);
    gl.drawBuffers([
      gl.COLOR_ATTACHMENT0,
      gl.COLOR_ATTACHMENT1,
      gl.COLOR_ATTACHMENT2,
      gl.COLOR_ATTACHMENT3
    ]);
    gl.uniform1i(innerSize_mul, K>>>2);

    gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
    gl.enableVertexAttribArray(pos_mul);
    gl.vertexAttribPointer(pos_mul, 2, gl.FLOAT, false, 0, 0);

    // SET WRITE TO MATRIX C
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, matrix_C0_tex, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT1, gl.TEXTURE_2D, matrix_C1_tex, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT2, gl.TEXTURE_2D, matrix_C2_tex, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT3, gl.TEXTURE_2D, matrix_C3_tex, /*lod=*/0);

    // READ FROM MATRICES A,B
    gl.uniform1i(matrix_A_mul, 0); // texture unit 0
    gl.uniform1i(matrix_B_mul, 1); // texture unit 1
    gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, matrix_A_tex);
    gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, matrix_B_tex);

    // COMPUTE
    gl.drawArrays(gl.TRIANGLE_STRIP, 0,4);

//    gl.finish();
//    console.timeEnd('COMPUTATION');

     //
    // POSTPROCESS
   //
//    console.time('POSTPROCESS');

    gl.useProgram(program_post);

    gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

    gl.viewport(0,0,J>>>2,I);
    gl.drawBuffers([
      gl.COLOR_ATTACHMENT0
    ]);

    gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
    gl.enableVertexAttribArray(pos_post);
    gl.vertexAttribPointer(pos_post, 2, gl.FLOAT, false, 0, 0);

    // SET WRITE TO MATRIX C
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, matrix_C_tex, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT1, gl.TEXTURE_2D, null, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT2, gl.TEXTURE_2D, null, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT3, gl.TEXTURE_2D, null, /*lod=*/0);

    // READ FROM MATRICES C0,C1,C2,C3
    gl.uniform1i(matrix_C0_post, 0);
    gl.uniform1i(matrix_C1_post, 1);
    gl.uniform1i(matrix_C2_post, 2);
    gl.uniform1i(matrix_C3_post, 3);
    gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, matrix_C0_tex);
    gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, matrix_C1_tex);
    gl.activeTexture(gl.TEXTURE2); gl.bindTexture(gl.TEXTURE_2D, matrix_C2_tex);
    gl.activeTexture(gl.TEXTURE3); gl.bindTexture(gl.TEXTURE_2D, matrix_C3_tex);

    // COMPUTE
    gl.drawArrays(gl.TRIANGLE_STRIP, 0,4);

    if( null != timeCallback ) {
      gl.finish();
      const dt = performance.now() - t0;
      timeCallback(dt);
//    console.timeEnd('POSTPROCESS');
    }


     //
    // DOWNLOAD
   //
    let outArr = new Float32Array(I*J);
//    console.time('DOWNLOAD')
    gl.readBuffer(gl.COLOR_ATTACHMENT0);
    gl.readPixels(
      /*x,y=*/0,0,
      /*w,h=*/J>>>2,I,
      /*format=*/gl.RGBA,
      /*type=*/gl.FLOAT,
      /*writeTo=*/outArr
    );
//    console.timeEnd('DOWNLOAD')

    return new nd.Array(Int32Array.of(I,J), outArr)
  };
}();


const gl_matmul_block2x2_v1 = function(){

  const vertShader = `\
    #version 300 es
    precision highp float;
    precision highp int;

    // INPUTS
    in vec2 pos;

    // OUTPUTS

    // SHADER
    void main() {
      gl_Position = vec4( pos, 0.0, 1.0 );
    }
  `;

  const program_mul = mkProgram(
    vertShader,
     //
    // MATMUL SHADER
   //
    `#version 300 es
      precision highp float;
      precision highp int;

      // INPUTS
      uniform highp int innerSize;
      uniform highp sampler2D matrix_A;
      uniform highp sampler2D matrix_B;

      // OUTPUTS
      out highp vec4 block2x2;

      highp mat2 read( highp sampler2D matrix, int i, int j ) {
         return mat2(
          texelFetch( matrix, ivec2(j,i*2+0), /*lod=*/0 ).xy,
          texelFetch( matrix, ivec2(j,i*2+1), /*lod=*/0 ).xy
        );
      }

      highp mat2 A( int i, int j ) { return read(matrix_A, i,j); }
      highp mat2 B( int i, int j ) { return read(matrix_B, i,j); }
      
      void main()
      {
        highp int i = int(gl_FragCoord.y),
                  j = int(gl_FragCoord.x);

        // https://en.wikipedia.org/wiki/Kahan_summation_algorithm#The_algorithm
        highp mat2 C_ij = mat2(0),
                   c_ij = mat2(0);

        for( int k=0; k < innerSize; k++ )
        {
          highp mat2 s = B(k,j)*A(i,k) - c_ij,
                     S = C_ij + s;
          c_ij = (S - C_ij) - s;
          C_ij =  S;
        }

        block2x2 = vec4( C_ij[0], C_ij[1] );
      }
    `
  );

  const program_post = mkProgram(
    vertShader,
     //
    // MATMUL SHADER
   //
    `#version 300 es
      precision highp float;
      precision highp int;

      // INPUTS
      uniform highp sampler2D matrix_C;

      // OUTPUTS
      out highp vec4 C_ij;

      highp vec4 C( int i, int j ) {

        highp int k = i % 2 * 2;
        i /= 2;

        highp vec4 c1 = texelFetch( matrix_C, ivec2(0+2*j,i), /*lod=*/0 ),
                   c2 = texelFetch( matrix_C, ivec2(1+2*j,i), /*lod=*/0 );

        return vec4(
          c1[k], c1[k+1],
          c2[k], c2[k+1]
        );
      }
      
      void main()
      {
        highp int i = int(gl_FragCoord.y),
                  j = int(gl_FragCoord.x);
        C_ij = C(i,j);
      }
    `
  );

   //
  // INIT
 //
  const   pos_mul = gl. getAttribLocation(program_mul, 'pos'),
    innerSize_mul = gl.getUniformLocation(program_mul, 'innerSize'),
     matrix_A_mul = gl.getUniformLocation(program_mul, 'matrix_A'),
     matrix_B_mul = gl.getUniformLocation(program_mul, 'matrix_B');

  const  pos_post = gl. getAttribLocation(program_post, 'pos'),
    matrix_C_post = gl.getUniformLocation(program_post, 'matrix_C');

  const posBuf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
  gl.bufferData(gl.ARRAY_BUFFER, Float32Array.of(
   -1,+1,
   +1,+1,
   -1,-1,
   +1,-1
  ), gl.STATIC_DRAW);

  const [A_tex,
         B_tex,
         C_tex,
         D_tex] = function*(){
    for( let i=4; i-- > 0; )
    {
      const tex = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, tex);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      yield tex;
    }
  }();

  const frameBuf = gl.createFramebuffer();
  gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

  return (A,B, timeCallback) => {
    if( A.ndim != 2 ) throw new Error('A is not a matrix.');
    if( B.ndim != 2 ) throw new Error('B is not a matrix.');
    if( A.dtype != 'float32' ) throw new Error("A.dtype must be 'float32'.");
    if( B.dtype != 'float32' ) throw new Error("B.dtype must be 'float32'.");

    const [I,K] = A.shape.slice(-2),
          [L,J] = B.shape.slice(-2);
    if( K != L ) throw new Error('A.shape[-1] != B.shape[-2]');

    if( I%2 != 0 ) throw new Error();
    if( K%2 != 0 ) throw new Error();
    if( J%4 != 0 ) throw new Error();

     //
    // UPLOAD
   //
//    console.time('UPLOAD');

    gl.bindTexture(gl.TEXTURE_2D, A_tex); gl.texImage2D(gl.TEXTURE_2D, /*lod=*/0, /*internal=*/gl.RG32F,   /*w,h=*/K>>>1, I,     /*border=*/0, /*format=*/gl.RG,   /*typ=*/gl.FLOAT, /*srcDat=*/A.data);
    gl.bindTexture(gl.TEXTURE_2D, B_tex); gl.texImage2D(gl.TEXTURE_2D, /*lod=*/0, /*internal=*/gl.RG32F,   /*w,h=*/J>>>1, K,     /*border=*/0, /*format=*/gl.RG,   /*typ=*/gl.FLOAT, /*srcDat=*/B.data);

    gl.finish();
    const t0 = performance.now();

    gl.bindTexture(gl.TEXTURE_2D, C_tex); gl.texImage2D(gl.TEXTURE_2D, /*lod=*/0, /*internal=*/gl.RGBA32F, /*w,h=*/J>>>1, I>>>1, /*border=*/0, /*format=*/gl.RGBA, /*typ=*/gl.FLOAT, /*srcDat=*/null  );
    gl.bindTexture(gl.TEXTURE_2D, D_tex); gl.texImage2D(gl.TEXTURE_2D, /*lod=*/0, /*internal=*/gl.RGBA32F, /*w,h=*/J>>>2, I,     /*border=*/0, /*format=*/gl.RGBA, /*typ=*/gl.FLOAT, /*srcDat=*/null  );

//    gl.finish();
//    console.timeEnd('UPLOAD');

     //
    // COMPUTE
   //
//    console.time('COMPUTATION');

    gl.useProgram(program_mul);

    gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

    gl.viewport(0,0,J>>>1,I>>>1);
    gl.drawBuffers([
      gl.COLOR_ATTACHMENT0
    ]);
    gl.uniform1i(innerSize_mul, K>>>1);

    gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
    gl.enableVertexAttribArray(pos_mul);
    gl.vertexAttribPointer(pos_mul, 2, gl.FLOAT, false, 0, 0);

    // SET WRITE TO MATRIX C
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, C_tex, /*lod=*/0);

    // READ FROM MATRICES A,B
    gl.uniform1i(matrix_A_mul, 0);
    gl.uniform1i(matrix_B_mul, 1);
    gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, A_tex);
    gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, B_tex);

    // COMPUTE
    gl.drawArrays(gl.TRIANGLE_STRIP, 0,4);

//    gl.finish();
//    console.timeEnd('COMPUTATION');

     //
    // POSTPROCESS
   //
//    console.time('POSTPROCESS');

    gl.useProgram(program_post);

    gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

    gl.viewport(0,0,J>>>2,I);
    gl.drawBuffers([
      gl.COLOR_ATTACHMENT0
    ]);

    gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
    gl.enableVertexAttribArray(pos_post);
    gl.vertexAttribPointer(pos_post, 2, gl.FLOAT, false, 0, 0);

    // SET WRITE TO MATRIX D
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, D_tex, /*lod=*/0);

    // READ FROM MATRIX C
    gl.uniform1i(matrix_C_post, 0);
    gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, C_tex);

    // COMPUTE
    gl.drawArrays(gl.TRIANGLE_STRIP, 0,4);

    if( null != timeCallback ) {
      gl.finish();
      const dt = performance.now() - t0;
      timeCallback(dt);
//    console.timeEnd('POSTPROCESS');
    }

     //
    // DOWNLOAD
   //
    let outArr = new Float32Array(I*J);
//    console.time('DOWNLOAD')
    gl.readBuffer(gl.COLOR_ATTACHMENT0);
    gl.readPixels(
      /*x,y=*/0,0,
      /*w,h=*/J>>>2,I,
      /*format=*/gl.RGBA,
      /*type=*/gl.FLOAT,
      /*writeTo=*/outArr
    );
//    console.timeEnd('DOWNLOAD')

    return new nd.Array(Int32Array.of(I,J), outArr)
  };
}();


const gl_matmul_block2x2_v2 = function(){

  const vertShader = `\
    #version 300 es
    precision highp float;
    precision highp int;

    // INPUTS
    in vec2 pos;

    // OUTPUTS

    // SHADER
    void main() {
      gl_Position = vec4( pos, 0.0, 1.0 );
    }
  `;

  const program_pre = mkProgram(
    vertShader,
     //
    // MATMUL SHADER
   //
    `#version 300 es
      precision highp float;
      precision highp int;

      // INPUTS
      uniform highp sampler2D matrix_AB;

      // OUTPUTS
      out highp vec4 block2x2;

      highp vec4 AB( int i, int j ) {
         return vec4(
          texelFetch( matrix_AB, ivec2(j,i*2+0), /*lod=*/0 ).xy,
          texelFetch( matrix_AB, ivec2(j,i*2+1), /*lod=*/0 ).xy
        );
      }
      
      void main()
      {
        highp int i = int(gl_FragCoord.y),
                  j = int(gl_FragCoord.x);
        block2x2 = AB(i,j);
      }
    `
  );

  const program_mul = mkProgram(
    vertShader,
     //
    // MATMUL SHADER
   //
    `#version 300 es
      precision highp float;
      precision highp int;

      // INPUTS
      uniform highp int innerSize;
      uniform highp sampler2D matrix_A;
      uniform highp sampler2D matrix_B;

      // OUTPUTS
      out highp vec4 block2x2;

      highp mat2 read( highp sampler2D matrix, int i, int j ) {
        highp vec4 m_ij = texelFetch( matrix, ivec2(j,i), /*lod=*/0 );
        return mat2(
          m_ij.xy,
          m_ij.zw
        );
      }

      highp mat2 A( int i, int j ) { return read(matrix_A, i,j); }
      highp mat2 B( int i, int j ) { return read(matrix_B, i,j); }
      
      void main()
      {
        highp int i = int(gl_FragCoord.y),
                  j = int(gl_FragCoord.x);

        // https://en.wikipedia.org/wiki/Kahan_summation_algorithm#The_algorithm
        highp mat2 C_ij = mat2(0),
                   c_ij = mat2(0);

        for( int k=0; k < innerSize; k++ )
        {
          highp mat2 s = B(k,j)*A(i,k) - c_ij,
                     S = C_ij + s;
          c_ij = (S - C_ij) - s;
          C_ij =  S;
        }

        block2x2 = vec4( C_ij[0], C_ij[1] );
      }
    `
  );

  const program_post = mkProgram(
    vertShader,
     //
    // MATMUL SHADER
   //
    `#version 300 es
      precision highp float;
      precision highp int;

      // INPUTS
      uniform highp sampler2D matrix_C;

      // OUTPUTS
      out highp vec4 C_ij;

      highp vec4 C( int i, int j ) {

        highp int k = i % 2 * 2;
        i /= 2;

        highp vec4 c1 = texelFetch( matrix_C, ivec2(0+2*j,i), /*lod=*/0 ),
                   c2 = texelFetch( matrix_C, ivec2(1+2*j,i), /*lod=*/0 );

        return vec4(
          c1[k], c1[k+1],
          c2[k], c2[k+1]
        );
      }
      
      void main()
      {
        highp int i = int(gl_FragCoord.y),
                  j = int(gl_FragCoord.x);
        C_ij = C(i,j);
      }
    `
  );

   //
  // INIT
 //
  const    pos_pre = gl. getAttribLocation(program_mul, 'pos'),
     matrix_AB_pre = gl.getUniformLocation(program_mul, 'matrix_AB');

  const   pos_mul = gl. getAttribLocation(program_mul, 'pos'),
    innerSize_mul = gl.getUniformLocation(program_mul, 'innerSize'),
     matrix_A_mul = gl.getUniformLocation(program_mul, 'matrix_A'),
     matrix_B_mul = gl.getUniformLocation(program_mul, 'matrix_B');

  const  pos_post = gl. getAttribLocation(program_post, 'pos'),
    matrix_C_post = gl.getUniformLocation(program_post, 'matrix_C');

  const posBuf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
  gl.bufferData(gl.ARRAY_BUFFER, Float32Array.of(
   -1,+1,
   +1,+1,
   -1,-1,
   +1,-1
  ), gl.STATIC_DRAW);

  const [a_tex,
         b_tex,
         A_tex,
         B_tex,
         C_tex,
         D_tex] = function*(){
    for( let i=6; i-- > 0; )
    {
      const tex = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, tex);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      yield tex;
    }
  }();

  const frameBuf = gl.createFramebuffer();
  gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

  return (A,B, timeCallback) => {
    if( A.ndim != 2 ) throw new Error('A is not a matrix.');
    if( B.ndim != 2 ) throw new Error('B is not a matrix.');
    if( A.dtype != 'float32' ) throw new Error("A.dtype must be 'float32'.");
    if( B.dtype != 'float32' ) throw new Error("B.dtype must be 'float32'.");

    const [I,K] = A.shape.slice(-2),
          [L,J] = B.shape.slice(-2);
    if( K != L ) throw new Error('A.shape[-1] != B.shape[-2]');

    if( I%2 != 0 ) throw new Error();
    if( K%2 != 0 ) throw new Error();
    if( J%4 != 0 ) throw new Error();

     //
    // UPLOAD
   //
//    console.time('UPLOAD');

    gl.bindTexture(gl.TEXTURE_2D, a_tex); gl.texImage2D(gl.TEXTURE_2D, /*lod=*/0, /*internal=*/gl.RG32F,   /*w,h=*/K>>>1, I,     /*border=*/0, /*read=*/gl.RG,   /*typ=*/gl.FLOAT, /*srcDat=*/A.data);
    gl.bindTexture(gl.TEXTURE_2D, b_tex); gl.texImage2D(gl.TEXTURE_2D, /*lod=*/0, /*internal=*/gl.RG32F,   /*w,h=*/J>>>1, K,     /*border=*/0, /*read=*/gl.RG,   /*typ=*/gl.FLOAT, /*srcDat=*/B.data);

    gl.finish();
    const t0 = performance.now();

    gl.bindTexture(gl.TEXTURE_2D, A_tex); gl.texImage2D(gl.TEXTURE_2D, /*lod=*/0, /*internal=*/gl.RGBA32F, /*w,h=*/K>>>1, I>>>1, /*border=*/0, /*read=*/gl.RGBA, /*typ=*/gl.FLOAT, /*srcDat=*/null  );
    gl.bindTexture(gl.TEXTURE_2D, B_tex); gl.texImage2D(gl.TEXTURE_2D, /*lod=*/0, /*internal=*/gl.RGBA32F, /*w,h=*/J>>>1, K>>>1, /*border=*/0, /*read=*/gl.RGBA, /*typ=*/gl.FLOAT, /*srcDat=*/null  );
    gl.bindTexture(gl.TEXTURE_2D, C_tex); gl.texImage2D(gl.TEXTURE_2D, /*lod=*/0, /*internal=*/gl.RGBA32F, /*w,h=*/J>>>1, I>>>1, /*border=*/0, /*read=*/gl.RGBA, /*typ=*/gl.FLOAT, /*srcDat=*/null  );
    gl.bindTexture(gl.TEXTURE_2D, D_tex); gl.texImage2D(gl.TEXTURE_2D, /*lod=*/0, /*internal=*/gl.RGBA32F, /*w,h=*/J>>>2, I,     /*border=*/0, /*read=*/gl.RGBA, /*typ=*/gl.FLOAT, /*srcDat=*/null  );

//    gl.finish();
//    console.timeEnd('UPLOAD');

     //
    // PREPROCESS A
   //
//    console.time('PREPROCESS_A');

    gl.useProgram(program_pre);

    gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

    gl.viewport(0,0,K>>>1,I>>>1);
    gl.drawBuffers([
      gl.COLOR_ATTACHMENT0
    ]);

    gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
    gl.enableVertexAttribArray(pos_pre);
    gl.vertexAttribPointer(pos_pre, 2, gl.FLOAT, false, 0, 0);

    // SET WRITE TO MATRIX A
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, A_tex, /*lod=*/0);

    // READ FROM MATRIX a
    gl.uniform1i(matrix_AB_pre, 0);
    gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, a_tex);

    // COMPUTE
    gl.drawArrays(gl.TRIANGLE_STRIP, 0,4);

//    gl.finish();
//    console.timeEnd('PREPROCESS_A');

     //
    // PREPROCESS B
   //
//    console.time('PREPROCESS_B');

    gl.useProgram(program_pre);

    gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

    gl.viewport(0,0,J>>>1,K>>>1);
    gl.drawBuffers([
      gl.COLOR_ATTACHMENT0
    ]);

    gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
    gl.enableVertexAttribArray(pos_pre);
    gl.vertexAttribPointer(pos_pre, 2, gl.FLOAT, false, 0, 0);

    // SET WRITE TO MATRIX A
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, B_tex, /*lod=*/0);

    // READ FROM MATRIX a
    gl.uniform1i(matrix_AB_pre, 0);
    gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, b_tex);

    // COMPUTE
    gl.drawArrays(gl.TRIANGLE_STRIP, 0,4);

//    gl.finish();
//    console.timeEnd('PREPROCESS_B');

     //
    // COMPUTE
   //
//    console.time('COMPUTATION');

    gl.useProgram(program_mul);

    gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

    gl.viewport(0,0,J>>>1,I>>>1);
    gl.drawBuffers([
      gl.COLOR_ATTACHMENT0
    ]);
    gl.uniform1i(innerSize_mul, K>>>1);

    gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
    gl.enableVertexAttribArray(pos_mul);
    gl.vertexAttribPointer(pos_mul, 2, gl.FLOAT, false, 0, 0);

    // SET WRITE TO MATRIX C
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, C_tex, /*lod=*/0);

    // READ FROM MATRICES A,B
    gl.uniform1i(matrix_A_mul, 0);
    gl.uniform1i(matrix_B_mul, 1);
    gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, A_tex);
    gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, B_tex);

    // COMPUTE
    gl.drawArrays(gl.TRIANGLE_STRIP, 0,4);

//    gl.finish();
//    console.timeEnd('COMPUTATION');

     //
    // POSTPROCESS
   //
//    console.time('POSTPROCESS');

    gl.useProgram(program_post);

    gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

    gl.viewport(0,0,J>>>2,I);
    gl.drawBuffers([
      gl.COLOR_ATTACHMENT0
    ]);

    gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
    gl.enableVertexAttribArray(pos_post);
    gl.vertexAttribPointer(pos_post, 2, gl.FLOAT, false, 0, 0);

    // SET WRITE TO MATRIX D
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, D_tex, /*lod=*/0);

    // READ FROM MATRIX C
    gl.uniform1i(matrix_C_post, 0);
    gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, C_tex);

    // COMPUTE
    gl.drawArrays(gl.TRIANGLE_STRIP, 0,4);

    if( null != timeCallback ) {
      gl.finish();
      const dt = performance.now() - t0;
      timeCallback(dt);
//    console.timeEnd('POSTPROCESS');
    }

     //
    // DOWNLOAD
   //
    let outArr = new Float32Array(I*J);
//    console.time('DOWNLOAD')
    gl.readBuffer(gl.COLOR_ATTACHMENT0);
    gl.readPixels(
      /*x,y=*/0,0,
      /*w,h=*/J>>>2,I,
      /*format=*/gl.RGBA,
      /*type=*/gl.FLOAT,
      /*writeTo=*/outArr
    );
//    console.timeEnd('DOWNLOAD')

    return new nd.Array(Int32Array.of(I,J), outArr)
  };
}();


const gl_matmul_block4x4_v2 = function(){

  const vertShader = `\
    #version 300 es
    precision highp float;
    precision highp int;

    // INPUTS
    in vec2 pos;

    // OUTPUTS

    // SHADER
    void main() {
      gl_Position = vec4( pos, 0.0, 1.0 );
    }
  `;

  const program_mul = mkProgram(
    vertShader,
     //
    // MATMUL SHADER
   //
    `#version 300 es
      precision highp float;
      precision highp int;

      // INPUTS
      uniform highp int innerSize;
      uniform highp sampler2D matrix_A;
      uniform highp sampler2D matrix_B;

      // OUTPUTS
      out highp vec4 block4x4[4];

      highp mat4 read( highp sampler2D matrix, int i, int j ) {
         return mat4(
          texelFetch( matrix, ivec2(j,i*4+0), /*lod=*/0 ),
          texelFetch( matrix, ivec2(j,i*4+1), /*lod=*/0 ),
          texelFetch( matrix, ivec2(j,i*4+2), /*lod=*/0 ),
          texelFetch( matrix, ivec2(j,i*4+3), /*lod=*/0 )
        );
      }

      highp mat4 A( int i, int j ) { return read(matrix_A, i,j); }
      highp mat4 B( int i, int j ) { return read(matrix_B, i,j); }
      
      void main()
      {
        highp int i = int(gl_FragCoord.y),
                  j = int(gl_FragCoord.x);

        // https://en.wikipedia.org/wiki/Kahan_summation_algorithm#The_algorithm
        highp mat4 C_ij = mat4(0),
                   c_ij = mat4(0);

        for( int k=0; k < innerSize; k++ )
        {
          highp mat4 s = B(k,j)*A(i,k) - c_ij,
                     S = C_ij + s;
          c_ij = (S - C_ij) - s;
          C_ij =  S;
        }

        block4x4[0] = C_ij[0];
        block4x4[1] = C_ij[1];
        block4x4[2] = C_ij[2];
        block4x4[3] = C_ij[3];
      }
    `
  );

  const program_post = mkProgram(
    vertShader,
     //
    // MATMUL SHADER
   //
    `#version 300 es
      precision highp float;
      precision highp int;

      // INPUTS
      uniform highp sampler2D C[4];

      // OUTPUTS
      out highp vec4 C_ij;

      highp vec4 read( highp sampler2D matrix, int i, int j ) {
        return texelFetch( matrix, ivec2(j,i), /*lod=*/0 );
      }
      
      void main()
      {
        highp int i = int(gl_FragCoord.y),
                  j = int(gl_FragCoord.x),
                  k = i%4;
        i /= 4;
        switch(k) {
          case 0: C_ij = read(C[0], i,j); break;
          case 1: C_ij = read(C[1], i,j); break;
          case 2: C_ij = read(C[2], i,j); break;
          case 3: C_ij = read(C[3], i,j); break;
        }
      }
    `
  );

   //
  // INIT
 //
  const   pos_mul = gl. getAttribLocation(program_mul, 'pos'),
    innerSize_mul = gl.getUniformLocation(program_mul, 'innerSize'),
     matrix_A_mul = gl.getUniformLocation(program_mul, 'matrix_A'),
     matrix_B_mul = gl.getUniformLocation(program_mul, 'matrix_B');

  const   pos_post = gl. getAttribLocation(program_post, 'pos'),
    matrix_C0_post = gl.getUniformLocation(program_post, 'C[0]'),
    matrix_C1_post = gl.getUniformLocation(program_post, 'C[1]'),
    matrix_C2_post = gl.getUniformLocation(program_post, 'C[2]'),
    matrix_C3_post = gl.getUniformLocation(program_post, 'C[3]');

  const posBuf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
  gl.bufferData(gl.ARRAY_BUFFER, Float32Array.of(
   -1,+1,
   +1,+1,
   -1,-1,
   +1,-1
  ), gl.STATIC_DRAW);

  const [matrix_A_tex,
         matrix_B_tex,
         matrix_C0_tex,
         matrix_C1_tex,
         matrix_C2_tex,
         matrix_C3_tex,
         matrix_C_tex] = function*(){
    for( let i=7; i-- > 0; )
    {
      const tex = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, tex);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      yield tex;
    }
  }();

  const frameBuf = gl.createFramebuffer();
  gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

  return (A,B, timeCallback) => {
    if( A.ndim != 2 ) throw new Error('A is not a matrix.');
    if( B.ndim != 2 ) throw new Error('B is not a matrix.');
    if( A.dtype != 'float32' ) throw new Error("A.dtype must be 'float32'.");
    if( B.dtype != 'float32' ) throw new Error("B.dtype must be 'float32'.");

    const [I,K] = A.shape.slice(-2),
          [L,J] = B.shape.slice(-2);
    if( K != L ) throw new Error('A.shape[-1] != B.shape[-2]');

    if( I%4 != 0 ) throw new Error();
    if( K%4 != 0 ) throw new Error();
    if( J%4 != 0 ) throw new Error();

     //
    // UPLOAD
   //
//    console.time('UPLOAD');

    gl.bindTexture(gl.TEXTURE_2D, matrix_A_tex);
    gl.texImage2D(
      /*target=*/gl.TEXTURE_2D,
      /*levelOfDetail=*/0,
      /*internalFormat=*/gl.RGBA32F,
      /*width,height=*/K>>>2, I,
      /*border=*/0,
      /*format=*/gl.RGBA,
      /*type=*/gl.FLOAT,
      /*srcData=*/A.data
    );
    gl.bindTexture(gl.TEXTURE_2D, matrix_B_tex);
    gl.texImage2D(
      /*target=*/gl.TEXTURE_2D,
      /*levelOfDetail=*/0,
      /*internalFormat=*/gl.RGBA32F,
      /*width,height=*/J>>>2, K,
      /*border=*/0,
      /*format=*/gl.RGBA,
      /*type=*/gl.FLOAT,
      /*srcData=*/B.data
    );

    gl.finish();
    const t0 = performance.now();

    for( const matrix of [
      matrix_C0_tex,
      matrix_C1_tex,
      matrix_C2_tex,
      matrix_C3_tex
    ] )
    {
      gl.bindTexture(gl.TEXTURE_2D, matrix);
      gl.texImage2D(
        /*target=*/gl.TEXTURE_2D,
        /*levelOfDetail=*/0,
        /*internalFormat=*/gl.RGBA32F,
        /*width,height=*/J>>>2,I>>>2,
        /*border=*/0,
        /*format=*/gl.RGBA,
        /*type=*/gl.FLOAT,
        /*srcData=*/null
      );
    }
    gl.bindTexture(gl.TEXTURE_2D, matrix_C_tex);
    gl.texImage2D(
      /*target=*/gl.TEXTURE_2D,
      /*levelOfDetail=*/0,
      /*internalFormat=*/gl.RGBA32F,
      /*width,height=*/J>>>2,I,
      /*border=*/0,
      /*format=*/gl.RGBA,
      /*type=*/gl.FLOAT,
      /*srcData=*/null
    );

//    gl.finish();
//    console.timeEnd('UPLOAD');

     //
    // COMPUTE
   //
//    console.time('COMPUTATION');


    gl.useProgram(program_mul);

    gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

    gl.viewport(0,0,J>>>2,I>>>2);
    gl.drawBuffers([
      gl.COLOR_ATTACHMENT0,
      gl.COLOR_ATTACHMENT1,
      gl.COLOR_ATTACHMENT2,
      gl.COLOR_ATTACHMENT3
    ]);
    gl.uniform1i(innerSize_mul, K>>>2);

    gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
    gl.enableVertexAttribArray(pos_mul);
    gl.vertexAttribPointer(pos_mul, 2, gl.FLOAT, false, 0, 0);

    // SET WRITE TO MATRIX C
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, matrix_C0_tex, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT1, gl.TEXTURE_2D, matrix_C1_tex, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT2, gl.TEXTURE_2D, matrix_C2_tex, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT3, gl.TEXTURE_2D, matrix_C3_tex, /*lod=*/0);

    // READ FROM MATRICES A,B
    gl.uniform1i(matrix_A_mul, 0); // texture unit 0
    gl.uniform1i(matrix_B_mul, 1); // texture unit 1
    gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, matrix_A_tex);
    gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, matrix_B_tex);

    // COMPUTE
    gl.drawArrays(gl.TRIANGLE_STRIP, 0,4);

//    gl.finish();
//    console.timeEnd('COMPUTATION');

     //
    // POSTPROCESS
   //
//    console.time('POSTPROCESS');

    gl.useProgram(program_post);

    gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

    gl.viewport(0,0,J>>>2,I);
    gl.drawBuffers([
      gl.COLOR_ATTACHMENT0
    ]);

    gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
    gl.enableVertexAttribArray(pos_post);
    gl.vertexAttribPointer(pos_post, 2, gl.FLOAT, false, 0, 0);

    // SET WRITE TO MATRIX C
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, matrix_C_tex, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT1, gl.TEXTURE_2D, null, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT2, gl.TEXTURE_2D, null, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT3, gl.TEXTURE_2D, null, /*lod=*/0);

    // READ FROM MATRICES C0,C1,C2,C3
    gl.uniform1i(matrix_C0_post, 0);
    gl.uniform1i(matrix_C1_post, 1);
    gl.uniform1i(matrix_C2_post, 2);
    gl.uniform1i(matrix_C3_post, 3);
    gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, matrix_C0_tex);
    gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, matrix_C1_tex);
    gl.activeTexture(gl.TEXTURE2); gl.bindTexture(gl.TEXTURE_2D, matrix_C2_tex);
    gl.activeTexture(gl.TEXTURE3); gl.bindTexture(gl.TEXTURE_2D, matrix_C3_tex);

    // COMPUTE
    gl.drawArrays(gl.TRIANGLE_STRIP, 0,4);

    if( null != timeCallback ) {
      gl.finish();
      const dt = performance.now() - t0;
      timeCallback(dt);
//    console.timeEnd('POSTPROCESS');
    }


     //
    // DOWNLOAD
   //
    let outArr = new Float32Array(I*J);
//    console.time('DOWNLOAD')
    gl.readBuffer(gl.COLOR_ATTACHMENT0);
    gl.readPixels(
      /*x,y=*/0,0,
      /*w,h=*/J>>>2,I,
      /*format=*/gl.RGBA,
      /*type=*/gl.FLOAT,
      /*writeTo=*/outArr
    );
//    console.timeEnd('DOWNLOAD')

    return new nd.Array(Int32Array.of(I,J), outArr)
  };
}();


const gl_matmul_block4x16 = function(){

  const vertShader = `\
    #version 300 es
    precision highp float;
    precision highp int;

    // INPUTS
    in vec2 pos;

    // OUTPUTS

    // SHADER
    void main() {
      gl_Position = vec4( pos, 0.0, 1.0 );
    }
  `;

  const program_mul = mkProgram(
    vertShader,
     //
    // MATMUL SHADER
   //
    `#version 300 es
      precision highp float;
      precision highp int;

      // INPUTS
      uniform highp int innerSize;
      uniform highp sampler2D matrix_A;
      uniform highp sampler2D matrix_B;

      // OUTPUTS
      out highp vec4 block4x4[8];

      highp mat4 A( int i, int j ) {
        return mat4(
          texelFetch( matrix_A, ivec2(j,i*4+0), /*lod=*/0 ),
          texelFetch( matrix_A, ivec2(j,i*4+1), /*lod=*/0 ),
          texelFetch( matrix_A, ivec2(j,i*4+2), /*lod=*/0 ),
          texelFetch( matrix_A, ivec2(j,i*4+3), /*lod=*/0 )
        );
      }
      
      void main()
      {
        highp int i = int(gl_FragCoord.y),
                  j = int(gl_FragCoord.x);
        j *= 2;

        // https://en.wikipedia.org/wiki/Kahan_summation_algorithm#The_algorithm
        highp mat4 C0_ij = mat4(0),
                   c0_ij = mat4(0),
                   C1_ij = mat4(0),
                   c1_ij = mat4(0);

        for( int k=0; k < innerSize; k++ )
        {
          highp mat4 A_ik = A(i,k),
                    B0_kj,
                    B1_kj;

          // READ IN B AS CONSEQUTIVELY AS POSSIBLE
          B0_kj[0] = texelFetch(matrix_B, ivec2(j+0,k*4+0), 0); B1_kj[0] = texelFetch(matrix_B, ivec2(j+1,k*4+0), 0);
          B0_kj[1] = texelFetch(matrix_B, ivec2(j+0,k*4+1), 0); B1_kj[1] = texelFetch(matrix_B, ivec2(j+1,k*4+1), 0);
          B0_kj[2] = texelFetch(matrix_B, ivec2(j+0,k*4+2), 0); B1_kj[2] = texelFetch(matrix_B, ivec2(j+1,k*4+2), 0);
          B0_kj[3] = texelFetch(matrix_B, ivec2(j+0,k*4+3), 0); B1_kj[3] = texelFetch(matrix_B, ivec2(j+1,k*4+3), 0);

          highp mat4 s,S;

          s = B0_kj*A_ik - c0_ij,
          S = C0_ij + s;
          c0_ij = (S - C0_ij) - s;
          C0_ij =  S;

          s = B1_kj*A_ik - c1_ij,
          S = C1_ij + s;
          c1_ij = (S - C1_ij) - s;
          C1_ij =  S;
        }

        block4x4[0] = C0_ij[0];
        block4x4[1] = C0_ij[1];
        block4x4[2] = C0_ij[2];
        block4x4[3] = C0_ij[3];
        block4x4[4] = C1_ij[0];
        block4x4[5] = C1_ij[1];
        block4x4[6] = C1_ij[2];
        block4x4[7] = C1_ij[3];
      }
    `
  );

  const program_post = mkProgram(
    vertShader,
     //
    // MATMUL SHADER
   //
    `#version 300 es
      precision highp float;
      precision highp int;

      // INPUTS
      uniform highp sampler2D C[8];

      // OUTPUTS
      out highp vec4 C_ij;

      highp vec4 read( highp sampler2D matrix, int i, int j ) {
        return texelFetch( matrix, ivec2(j,i), /*lod=*/0 );
      }
      
      void main()
      {
        highp int i = int(gl_FragCoord.y),
                  j = int(gl_FragCoord.x),
                  k = i%4 + j%2*4;
        i /= 4;
        j /= 2;
        switch(k) {
          case 0: C_ij = read(C[0], i,j); break;
          case 1: C_ij = read(C[1], i,j); break;
          case 2: C_ij = read(C[2], i,j); break;
          case 3: C_ij = read(C[3], i,j); break;
          case 4: C_ij = read(C[4], i,j); break;
          case 5: C_ij = read(C[5], i,j); break;
          case 6: C_ij = read(C[6], i,j); break;
          case 7: C_ij = read(C[7], i,j); break;
        }
      }
    `
  );

   //
  // INIT
 //
  const   pos_mul = gl. getAttribLocation(program_mul, 'pos'),
    innerSize_mul = gl.getUniformLocation(program_mul, 'innerSize'),
     matrix_A_mul = gl.getUniformLocation(program_mul, 'matrix_A'),
     matrix_B_mul = gl.getUniformLocation(program_mul, 'matrix_B');

  const   pos_post = gl. getAttribLocation(program_post, 'pos'),
    matrix_C0_post = gl.getUniformLocation(program_post, 'C[0]'),
    matrix_C1_post = gl.getUniformLocation(program_post, 'C[1]'),
    matrix_C2_post = gl.getUniformLocation(program_post, 'C[2]'),
    matrix_C3_post = gl.getUniformLocation(program_post, 'C[3]'),
    matrix_C4_post = gl.getUniformLocation(program_post, 'C[4]'),
    matrix_C5_post = gl.getUniformLocation(program_post, 'C[5]'),
    matrix_C6_post = gl.getUniformLocation(program_post, 'C[6]'),
    matrix_C7_post = gl.getUniformLocation(program_post, 'C[7]');

  const posBuf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
  gl.bufferData(gl.ARRAY_BUFFER, Float32Array.of(
   -1,+1,
   +1,+1,
   -1,-1,
   +1,-1
  ), gl.STATIC_DRAW);

  const [matrix_A_tex,
         matrix_B_tex,
         matrix_C0_tex,
         matrix_C1_tex,
         matrix_C2_tex,
         matrix_C3_tex,
         matrix_C4_tex,
         matrix_C5_tex,
         matrix_C6_tex,
         matrix_C7_tex,
         matrix_C_tex] = function*(){
    for( let i=11; i-- > 0; )
    {
      const tex = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, tex);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      yield tex;
    }
  }();

  const frameBuf = gl.createFramebuffer();
  gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

  return (A,B, timeCallback) => {
    if( A.ndim != 2 ) throw new Error('A is not a matrix.');
    if( B.ndim != 2 ) throw new Error('B is not a matrix.');
    if( A.dtype != 'float32' ) throw new Error("A.dtype must be 'float32'.");
    if( B.dtype != 'float32' ) throw new Error("B.dtype must be 'float32'.");

    const [I,K] = A.shape.slice(-2),
          [L,J] = B.shape.slice(-2);
    if( K != L ) throw new Error('A.shape[-1] != B.shape[-2]');

    if( I%4 != 0 ) throw new Error();
    if( K%4 != 0 ) throw new Error();
    if( J%8 != 0 ) throw new Error();

     //
    // UPLOAD
   //
//    console.time('UPLOAD');

    gl.bindTexture(gl.TEXTURE_2D, matrix_A_tex);
    gl.texImage2D(
      /*target=*/gl.TEXTURE_2D,
      /*levelOfDetail=*/0,
      /*internalFormat=*/gl.RGBA32F,
      /*width,height=*/K>>>2, I,
      /*border=*/0,
      /*format=*/gl.RGBA,
      /*type=*/gl.FLOAT,
      /*srcData=*/A.data
    );
    gl.bindTexture(gl.TEXTURE_2D, matrix_B_tex);
    gl.texImage2D(
      /*target=*/gl.TEXTURE_2D,
      /*levelOfDetail=*/0,
      /*internalFormat=*/gl.RGBA32F,
      /*width,height=*/J>>>2, K,
      /*border=*/0,
      /*format=*/gl.RGBA,
      /*type=*/gl.FLOAT,
      /*srcData=*/B.data
    );

    gl.finish();
    const t0 = performance.now();

    for( const matrix of [
      matrix_C0_tex,
      matrix_C1_tex,
      matrix_C2_tex,
      matrix_C3_tex,
      matrix_C4_tex,
      matrix_C5_tex,
      matrix_C6_tex,
      matrix_C7_tex,
    ] )
    {
      gl.bindTexture(gl.TEXTURE_2D, matrix);
      gl.texImage2D(
        /*target=*/gl.TEXTURE_2D,
        /*levelOfDetail=*/0,
        /*internalFormat=*/gl.RGBA32F,
        /*width,height=*/J>>>3,I>>>2,
        /*border=*/0,
        /*format=*/gl.RGBA,
        /*type=*/gl.FLOAT,
        /*srcData=*/null
      );
    }
    gl.bindTexture(gl.TEXTURE_2D, matrix_C_tex);
    gl.texImage2D(
      /*target=*/gl.TEXTURE_2D,
      /*levelOfDetail=*/0,
      /*internalFormat=*/gl.RGBA32F,
      /*width,height=*/J>>>2,I,
      /*border=*/0,
      /*format=*/gl.RGBA,
      /*type=*/gl.FLOAT,
      /*srcData=*/null
    );

//    gl.finish();
//    console.timeEnd('UPLOAD');

     //
    // COMPUTE
   //
//    console.time('COMPUTATION');


    gl.useProgram(program_mul);

    gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

    gl.viewport(0,0,J>>>3,I>>>2);
    gl.drawBuffers([
      gl.COLOR_ATTACHMENT0,
      gl.COLOR_ATTACHMENT1,
      gl.COLOR_ATTACHMENT2,
      gl.COLOR_ATTACHMENT3,
      gl.COLOR_ATTACHMENT4,
      gl.COLOR_ATTACHMENT5,
      gl.COLOR_ATTACHMENT6,
      gl.COLOR_ATTACHMENT7
    ]);
    gl.uniform1i(innerSize_mul, K>>>2);

    gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
    gl.enableVertexAttribArray(pos_mul);
    gl.vertexAttribPointer(pos_mul, 2, gl.FLOAT, false, 0, 0);

    // SET WRITE TO MATRIX C
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, matrix_C0_tex, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT1, gl.TEXTURE_2D, matrix_C1_tex, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT2, gl.TEXTURE_2D, matrix_C2_tex, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT3, gl.TEXTURE_2D, matrix_C3_tex, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT4, gl.TEXTURE_2D, matrix_C4_tex, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT5, gl.TEXTURE_2D, matrix_C5_tex, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT6, gl.TEXTURE_2D, matrix_C6_tex, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT7, gl.TEXTURE_2D, matrix_C7_tex, /*lod=*/0);

    // READ FROM MATRICES A,B
    gl.uniform1i(matrix_A_mul, 0); // texture unit 0
    gl.uniform1i(matrix_B_mul, 1); // texture unit 1
    gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, matrix_A_tex);
    gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, matrix_B_tex);

    // COMPUTE
    gl.drawArrays(gl.TRIANGLE_STRIP, 0,4);

//    gl.finish();
//    console.timeEnd('COMPUTATION');

     //
    // POSTPROCESS
   //
//    console.time('POSTPROCESS');

    gl.useProgram(program_post);

    gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

    gl.viewport(0,0,J>>>2,I);
    gl.drawBuffers([
      gl.COLOR_ATTACHMENT0
    ]);

    gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
    gl.enableVertexAttribArray(pos_post);
    gl.vertexAttribPointer(pos_post, 2, gl.FLOAT, false, 0, 0);

    // SET WRITE TO MATRIX C
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, matrix_C_tex, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT1, gl.TEXTURE_2D, null, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT2, gl.TEXTURE_2D, null, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT3, gl.TEXTURE_2D, null, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT4, gl.TEXTURE_2D, null, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT5, gl.TEXTURE_2D, null, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT6, gl.TEXTURE_2D, null, /*lod=*/0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT7, gl.TEXTURE_2D, null, /*lod=*/0);

    // READ FROM MATRICES C0,C1,C2,C3
    gl.uniform1i(matrix_C0_post, 0);
    gl.uniform1i(matrix_C1_post, 1);
    gl.uniform1i(matrix_C2_post, 2);
    gl.uniform1i(matrix_C3_post, 3);
    gl.uniform1i(matrix_C4_post, 4);
    gl.uniform1i(matrix_C5_post, 5);
    gl.uniform1i(matrix_C6_post, 6);
    gl.uniform1i(matrix_C7_post, 7);
    gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, matrix_C0_tex);
    gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, matrix_C1_tex);
    gl.activeTexture(gl.TEXTURE2); gl.bindTexture(gl.TEXTURE_2D, matrix_C2_tex);
    gl.activeTexture(gl.TEXTURE3); gl.bindTexture(gl.TEXTURE_2D, matrix_C3_tex);
    gl.activeTexture(gl.TEXTURE4); gl.bindTexture(gl.TEXTURE_2D, matrix_C4_tex);
    gl.activeTexture(gl.TEXTURE5); gl.bindTexture(gl.TEXTURE_2D, matrix_C5_tex);
    gl.activeTexture(gl.TEXTURE6); gl.bindTexture(gl.TEXTURE_2D, matrix_C6_tex);
    gl.activeTexture(gl.TEXTURE7); gl.bindTexture(gl.TEXTURE_2D, matrix_C7_tex);

    // COMPUTE
    gl.drawArrays(gl.TRIANGLE_STRIP, 0,4);

    if( null != timeCallback ) {
      gl.finish();
      const dt = performance.now() - t0;
      timeCallback(dt);
//    console.timeEnd('POSTPROCESS');
    }


     //
    // DOWNLOAD
   //
    let outArr = new Float32Array(I*J);
//    console.time('DOWNLOAD')
    gl.readBuffer(gl.COLOR_ATTACHMENT0);
    gl.readPixels(
      /*x,y=*/0,0,
      /*w,h=*/J>>>2,I,
      /*format=*/gl.RGBA,
      /*type=*/gl.FLOAT,
      /*writeTo=*/outArr
    );
//    console.timeEnd('DOWNLOAD')

    return new nd.Array(Int32Array.of(I,J), outArr)
  };
}();


function main()
{
  const sleep = dt => new Promise( f => setTimeout(f,dt) );

  const timeit = job => {
    const t0 = performance.now();
    job();
    return performance.now() - t0;
  };

  const is_close = (x,y) => {
    const atol = 1e-3,
          rtol = 1e-2,
           tol = atol + rtol * Math.max(
            Math.abs(x),
            Math.abs(y)
          );
    return Math.abs(x-y) <= tol;
  };

  async function test( method_dict )
  {
    for( let run=0; ++run <= 8; )
    {
      console.log(`Run${run.toString().padStart(4)}`);

      const N = 128;
      const I = 4*( Math.floor(Math.random()*N) + 1 );
      const K = 4*( Math.floor(Math.random()*N) + 1 );
      const J = 8*( Math.floor(Math.random()*N) + 1 );

//      console.log('  I,K,J:', [I,K,J]);

      const a = nd.tabulate([I,K], 'float32', (i,j) => Math.random()*8-4 ),
            b = nd.tabulate([K,J], 'float32', (i,j) => Math.random()*8-4 ),
            C = nd.la.matmul2(a,b);

      for( const [name,gl_matmul] of Object.entries(method_dict) )
      {
        console.log(' ',name);
        await sleep();

 //       console.log(' ', name);
        const c = gl_matmul(a,b);

  //      console.log('A:'); console.log(a.toString());
  //      console.log('B:'); console.log(b.toString());
  //      console.log('C:'); console.log(c.toString());

        nd.Array.from([C,c], 'float32', (x,y,...indices) => {
          if( ! is_close(x,y) ) {
            let msg = '{\n'+C+'\n} expected but {\n'+c+'\n} encountered.\n'+x+' != '+y+' at index ['+indices+']';
            throw new Error(msg);
          }
        });
      }
    }
  }
/*
  test({
 //   gl_matmul,
    gl_matmul_block1x4,
    gl_matmul_block2x2_v1,
    gl_matmul_block2x2_v2,
    gl_matmul_block4x4_v1,
    gl_matmul_block4x4_v2,
    gl_matmul_block4x16
  });
*/
  async function benchmark( method_dict )
  {
    const plot = document.createElement('div');
    document.body.appendChild(plot);

    const sizes = [],
          times = {};

    for( const name of Object.keys(method_dict) )
      times[name] = [];

    // CREATE PLOT
    {
      const data = [];

      for( const [name,time] of Object.entries(times) )
        data.push({
          type: 'scattergl',
          mode: 'lines',
          name,
          x: sizes,
          y: time
        });

      const layout = {
        title: '<b>WebGL Matrix Multiplication Benchmark</b><br><i>Warning:</i> Benchmark quite demanding in the end.',
        height: 900,
        xaxis: { title: 'size' },
        yaxis: { title: 'time [msec]' }
      };

      Plotly.plot(plot,{data,layout});
    }

    for( let N=8; N <= 2*1024; N += 8 )
    {
      const [I,K,J] = [N,N,N];

      sizes.push(N);

      const A = nd.tabulate([I,K], 'float32', (i,j) => Math.random()*2-1 ),
            B = nd.tabulate([K,J], 'float32', (i,j) => Math.random()*2-1 );
      let   C;

      for( const [name,gl_matmul] of Object.entries(method_dict) )
      {
        await sleep(16); // <- maybe sleeping a little allows the GPU to flush/cool
        C = gl_matmul(A,B, dt => times[name].push(dt) );
      }

      // UPDATE PLOT
      Plotly.restyle(plot, {
        y: Object.values(times).map( x => x.slice() )
      });
    }
  }

  benchmark({
    gl_matmul,
    gl_matmul_block1x4,
    gl_matmul_block2x2_v1,
    gl_matmul_block2x2_v2,
    gl_matmul_block4x4_v1,
    gl_matmul_block4x4_v2,
    gl_matmul_block4x16
  });
}

main();
