const     path = require('path'),
       webpack = require('webpack'),
  WorkerPlugin = require('worker-plugin');
 
//   ClosurePlugin = require('closure-webpack-plugin');

const cfg_module = {
  rules: [{
    test: /\.js$/,
    include: path.resolve(__dirname, 'src'),
    // exclude: /(.*_test\.js)|(.*\.py)|(.*_test_data\.js)/,
    use: {
      loader: 'babel-loader',
      options: {
        presets: ['@babel/preset-env'],
        plugins: [
          '@babel/transform-runtime'
        ]
      }
    }
  }]
}

module.exports = [
  {
    mode: 'development',
    entry: './src/help.js',
    output: {
      path: path.resolve(__dirname, 'dist'),
      filename: 'nd.js',
      library: 'nd',
      libraryTarget: 'umd',
      globalObject: "this"
    },
    plugins: [
      new WorkerPlugin({
        globalObject: 'this'
      })
    ],
    module: cfg_module,
  },
  {
    mode: 'production',
    entry: './src/index.js',
    output: {
      path: path.resolve(__dirname, 'dist'),
      filename: 'nd.min.js',
      library: 'nd',
      libraryTarget: 'umd',
      globalObject: "this"
    },
    plugins: [
      new WorkerPlugin({
        globalObject: 'this'
      })
    ],
    module: cfg_module
  }
]
