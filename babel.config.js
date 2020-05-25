'use strict';

const glob = require('glob');

module.exports = api => {
  api.cache(true);

  return {
    presets: ['@babel/preset-env'],
    plugins: ['@babel/transform-runtime'],
    only: glob.sync("src/**/*.js").filter( path => {
      path = path.split(/[\\/]/).pop();

      return !path.  endsWith('_test.js')
          && !path.  endsWith('_test_data.js')
          && !path.  endsWith('_test_utils.js')
          && !path.startsWith('_generic_test')
          &&  path !== 'jasmine_utils.js';
    }),
    sourceType: 'module',
    sourceMaps: false,

// TODO: get source maps to work and enable the following configs:
//    sourceMaps: 'inline',
//    minified: true,
//    comments: false,
  };
};
