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

/** A veritable WebGL Rube-Goldberg-Machine: converts js_float32->gl_rgba->gl_float32->gl_rgba->js_float32... yay!
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

const gl_rgba_to_float/*: nd.Array => nd.Array*/ = function() {
  const VERT_SHADER = `\
    precision highp float;

    // INPUTS
    attribute highp vec2 pos2d;
    uniform highp vec2 size;

    // OUTPUTS
    varying highp vec2 index;

    // SHADER
    void main() {
      index = (pos2d+1.0) / 2.0;
      gl_Position = vec4( pos2d, 0.0, 1.0 );
    }
  `;

  const FRAG_SHADER = `\
    precision highp float;
    precision highp int;

    // INPUTS
    varying highp vec2 index;
    uniform sampler2D matrix;
    uniform float NaN;

    // OUTPUTS

    highp float rgba_to_float( highp vec4 rgba )
    { // https://en.wikipedia.org/wiki/IEEE_754-1985#Representation_of_numbers
      //
      // S: sign bit (1 means negative)
      // E: exponent bits
      // M: mantissa bits
      //
      // |bytes[3] |bytes[2] |bytes[1] |bytes[0] |
      // |SEEE EEEE|EMMM MMMM|MMMM MMMM|MMMM MMMM|

      highp ivec4 bytes = ivec4(rgba*255.0);

      highp int sign = 1 - (bytes[3]/128*2),
            mantissa =(bytes[2]/128 == 0 ? 128*256*256 : 0)
                     + bytes[2] * 256*256
                     + bytes[1] * 256
                     + bytes[0],
            exponent = bytes[3]*2
                     - bytes[3]/128*256
                     + bytes[2]/128;

      if( exponent == 255 && mantissa != 128*256*256 )
        return NaN;

      return float(sign*mantissa) * exp2( float(exponent-150) );
    }

    bool isNaN( highp float x ) {
      return ! ( x >= -0.0 || x <= +0.0 );
    }

    bool isInf( highp float x ) {
      return 0.0 == 1.0/x;
    }

    highp vec4 float_to_rgba( highp float val )
    {
      // ZERO
      if( 0.0 == val ) return vec4( vec3(0), (1.0/val < 0.0) ? 128.0/255.0 : 0.0 );

      // NaN
      if( isNaN(val) ) return vec4( vec3(255), 127) / 255.0;

      highp int sign = val < 0.0 ? 128 : 0;

      // +- INFINITY
      if( isInf(val) ) return vec4(0, 0, 128, sign + 127) / 255.0;

      val = abs(val);

      // USE BINARY SEARCH TO FIND THE EXPONENT
      highp int lo = -127,
                hi = +128;
      for( int i = 0; i < 8; i++ ) {
        highp int mid = (lo+hi)/2;
        highp float midVal = val*exp2( float(-mid) );
        if( midVal >= 1.0 ) lo = mid;
        if( midVal <  2.0 ) hi = mid;
      };
      highp int exponent = lo + 127;

      val *= exp2( float(23-lo) );

      highp int mantissa = int(val) - 128*256*256;

      highp vec4 bytes = vec4(
        mantissa           - mantissa/ 256         *256,
        mantissa/ 256      - mantissa/(256*256)    *256,
        mantissa/(256*256) - mantissa/(256*256*128)*128 + exponent*128 - exponent/2*256,
                                                   sign + exponent/2
      );
      bytes /= 255.0;

      return bytes;
    }

    void main()
    {
      gl_FragColor = float_to_rgba( 2.0*rgba_to_float(
        texture2D(matrix, index)
      ) );
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

  const pos2dLoc = gl. getAttribLocation(program, 'pos2dLoc'),
         sizeLoc = gl.getUniformLocation(program, 'size'),
          NaNLoc = gl.getUniformLocation(program, 'NaN');

  const pos2dBuf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, pos2dBuf);
  gl.bufferData(gl.ARRAY_BUFFER, Float32Array.of(
    -1,+1,
    +1,+1,
    -1,-1,
    +1,-1
  ), gl.STATIC_DRAW);

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

  return S => {
    if( S.ndim != 2 ) throw new Error('S is not a matrix.');
    const [M,N] = S.shape.slice(-2);
    if( M != N ) throw new Error('S must be quadratic.');
    if( S.dtype != 'float32' ) throw new Error("S.dtype must be 'float32'.");

    gl.useProgram(program);

    gl.bindBuffer(gl.ARRAY_BUFFER, pos2dBuf);
    gl.enableVertexAttribArray(pos2dBuf);
    gl.vertexAttribPointer(pos2dBuf, 2, gl.FLOAT, false, 0, 0);

    gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuf);

     //
    // UPLOAD
   //
    console.time('upload');

    gl.bindTexture(gl.TEXTURE_2D, matrix1Tex);
    gl.texImage2D(
      /*target=*/gl.TEXTURE_2D,
      /*levelOfDetail=*/0,
      /*internalFormat=*/gl.RGBA,
      /*width,height=*/N,N,
      /*border=*/0,
      /*format=*/gl.RGBA,
      /*type=*/gl.UNSIGNED_BYTE,
      /*srcData=*/new Uint8Array(S.data.buffer)
    );
    gl.bindTexture(gl.TEXTURE_2D, matrix2Tex);
    gl.texImage2D(
      /*target=*/gl.TEXTURE_2D,
      /*levelOfDetail=*/0,
      /*internalFormat=*/gl.RGBA,
      /*width,height=*/N,N,
      /*border=*/0,
      /*format=*/gl.RGBA,
      /*type=*/gl.UNSIGNED_BYTE,
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

    gl.uniform1f(NaNLoc, NaN);
    gl.uniform2f(sizeLoc, N,N);

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
    gl.drawArrays(gl.TRIANGLE_STRIP, 0,4);

    [matrix1Tex,matrix2Tex] = [matrix2Tex,matrix1Tex];

    gl.finish();
    console.timeEnd('computation');

     //
    // DOWNLOAD
   //
    let matrixOut = new Float32Array(N*N);
    console.time('download')
    gl.readBuffer(gl.COLOR_ATTACHMENT0);
    gl.readPixels(
      /*x,y=*/0,0,
      /*w,h=*/N,N,
      /*format=*/gl.RGBA,
      /*type=*/gl.UNSIGNED_BYTE,
      /*writeTo=*/new Uint8Array(matrixOut.buffer)
    );
    console.timeEnd('download')

    return new nd.Array(S.shape,matrixOut)
  };
}();

async function main()
{
  const sleep = dt => new Promise( f => setTimeout(f,dt) );

  const [gl_NaN] = gl_rgba_to_float( nd.array('float32',[[NaN]]) ).data;
  if( ! isNaN(gl_NaN) ) throw new Error(`Assertion error: ${gl_NaN} != NaN`);

  const gl_infPos = gl_rgba_to_float( nd.array('float32',[[+Infinity]]) ).data;
  const gl_infNeg = gl_rgba_to_float( nd.array('float32',[[-Infinity]]) ).data;

  if( ! ( ! isNaN(gl_infPos) && ! isFinite(gl_infPos) && gl_infPos > 0 ) ) throw new Error(`${gl_infPos} != +Inf`);
  if( ! ( ! isNaN(gl_infNeg) && ! isFinite(gl_infNeg) && gl_infNeg < 0 ) ) throw new Error(`${gl_infNeg} != -Inf`);

  const a = nd.tabulate([5,5], 'float32', (i,j) => 10*i+j+11 );
  const b = gl_rgba_to_float(a);
  console.log( a.toString() );
  console.log( b.toString() );

  for( let run=0; ++run < 128; )
  {
    await sleep(); // <- keeps the browser responsive

    console.log(`Run${run.toString().padStart(4)}`);
    const N = Math.trunc( Math.random()*1024 ) + 1;
    const a = nd.tabulate([N,N], 'float32', (i,j) => Math.random()*2e37 - 1e37 );
    const b = gl_rgba_to_float(a);
    nd.Array.from([a,b], 'float32', (a,b) => {
      if( 2.0*a !== b ) throw new Error(`Assertion error: ${a} != ${b}`);
    });
  }
}

main();
