'use strict';

const glob = require('glob');

module.exports = api => {
  api.cache(true);

  return {
    presets: ['@babel/preset-env'],
    plugins: ['@babel/transform-runtime'],
    only: glob.sync("src/**/*.js").filter(
      path => !path.endsWith('_test.js') &&
              !path.endsWith('_test_data.js')
    ),
    sourceType: 'module',
    sourceMaps: false,

// TODO: get source maps to work and enable the following configs:
//    sourceMaps: 'inline',
//    minified: true,
//    comments: false,
  };
};
