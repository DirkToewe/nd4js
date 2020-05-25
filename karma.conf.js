const os = require('os');

const HOUR = 60*60*1000

const browsers = [
  'FirefoxHeadless',
  'ChromeHeadless'
];

const tested_files = [
  'src/**/*_test.js',
];

const reporters = [
//  'spec'
  'progress'
];


const nExecutors = Math.min( browsers.some(b => b.startsWith('Firefox')) ? 12 : 24, // <- Firefox can't seem handle more than 12
                   Math.max( 1, Math.floor(os.cpus().length / browsers.length) ));

module.exports = config => {
  config.set({
    basePath: '',

    frameworks: [
      'parallel',
      'jasmine'
    ],

//    browserDisconnectTolerance: 1e6,
    browserDisconnectTimeout: 1*HOUR,
    browserNoActivityTimeout: 1*HOUR,
        browserSocketTimeout: 1*HOUR,
          processKillTimeout: 1*HOUR,
//              captureTimeout: 1*HOUR,

    plugins: [
      'karma-parallel',
      'karma-chrome-launcher',
      'karma-firefox-launcher',
      'karma-jasmine',
      'karma-webpack',
      'karma-spec-reporter',
      'karma-sourcemap-loader',
    ],

    parallelOptions: {
      executors: nExecutors,
//      shardStrategy: 'round-robin'
    },

    client: {
      jasmine: {
        random: false,
        failFast: true,
        oneFailurePerSpec: true,
        stopSpecOnExpectationFailure: true,
                 timeoutInterval: 8*HOUR,
          defaultTimeoutInterval: 8*HOUR,
        DEFAULT_TIMEOUT_INTERVAL: 8*HOUR
      }
    },

    webpack: {
      mode: 'development',
      devtool: 'inline-source-map'
    },

    files: tested_files,

    exclude: [
    ],

    preprocessors: {
      'src/**/*.js': ['webpack', 'sourcemap']
    },

    reporters,

    specReporter: {
      showSpecTiming: true
    },

    browsers,

    port: 9876,
    logLevel: config.LOG_INFO,
    autoWatch: true,
    singleRun: true,
    concurrency: Infinity,
    colors: true
  })
}
