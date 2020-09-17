const       os = require('os'),
  WorkerPlugin = require('worker-plugin');

const SEC = 1000,
      MIN = 60*SEC,
     HOUR = 60*MIN;

const browsers = [
  'FirefoxHeadless',
  'ChromeHeadless'
];

const tested_files = [
  'src/**/*_test.js',
];

const reporters = [
//  'spec'
//  'progress'
  'progress-bar'
];


const nExecutors = Math.min( browsers.some(b => b.includes('Firefox')) ? 12 : 24, // <- Firefox can't seem handle more than 12
                   Math.max( 1, Math.floor(os.cpus().length / browsers.length) - 1 ));

module.exports = config => {
  config.set({
    basePath: '',

    frameworks: [
      'parallel',
      'jasmine'
    ],

    // browserDisconnectTolerance: 10*SEC,
    browserDisconnectTimeout: 5*MIN,
    browserNoActivityTimeout: 5*MIN,
        browserSocketTimeout: 5*MIN,
          processKillTimeout: 5*MIN,
//              captureTimeout: 1*HOUR,

    autoWatch: false,
//    reportSlowerThan: 5*SEC,

    plugins: [
      'karma-parallel',
      'karma-chrome-launcher',
      'karma-firefox-launcher',
      'karma-jasmine',
      'karma-webpack',
      'karma-spec-reporter',
      'karma-sourcemap-loader',
      require('./progress_bar_reporter')
    ],

    parallelOptions: {
      executors: nExecutors,
//      shardStrategy: 'round-robin'
    },

    client: {
      jasmine: {
        random: false,
//        random: true,
//        failFast: true,
        oneFailurePerSpec: true,
        stopSpecOnExpectationFailure: true,
                 timeoutInterval: 8*HOUR,
          defaultTimeoutInterval: 8*HOUR,
        DEFAULT_TIMEOUT_INTERVAL: 8*HOUR
      }
    },

    webpack: {
      mode: 'development',
      devtool: 'inline-source-map',
      plugins: [
        new WorkerPlugin({
          globalObject: 'this'
        })
      ]
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

    customLaunchers: {
      MyChromeHeadless: {
        base: "ChromeHeadless",
        flags: [
          "--js-flags=\"--max_old_space_size=2048 --max_semi_space_size=1024\""
        ]
      }
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
