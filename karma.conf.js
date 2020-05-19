const HOUR = 60*60*1000

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
      executors: 11,
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

    files: [
//      'src/opt/line_search/*_test.js',
//      'src/**/lbfgs_test.js',
//      'src/**/lbfgsb_test.js',
//      'src/opt/**/*_test.js',
      'src/**/*_test.js',
    ],

    exclude: [
    ],

    preprocessors: {
      'src/**/*.js': ['webpack', 'sourcemap']
    },

    reporters: [
//      'spec'
      'progress'
    ],

    specReporter: {
      showSpecTiming: true 
    },

    browsers: [
      'FirefoxHeadless',
      'ChromeHeadless',
    ],

    port: 9876,
    logLevel: config.LOG_INFO,
    autoWatch: true,
    singleRun: true,
    concurrency: Infinity,
    colors: true
  })
}
