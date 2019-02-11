// Karma configuration
// Generated on Tue Dec 04 2018 16:07:10 GMT+0100 (Central European Standard Time)

const HOUR = 60*60*1000

module.exports = function(config) {
  config.set({
    basePath: '',

    frameworks: ['jasmine'],

//    browserDisconnectTolerance: 1e6,
    browserDisconnectTimeout: 24*HOUR,
    browserNoActivityTimeout: 24*HOUR,
        browserSocketTimeout: 24*HOUR,
//              captureTimeout: 8*HOUR,

    plugins: [
      'karma-chrome-launcher',
      'karma-firefox-launcher',
      'karma-jasmine',
      'karma-webpack',
      'karma-spec-reporter',
      'karma-sourcemap-loader',
    ],

    client: {
      jasmine: {
        random: false,
        failFast: true,
        oneFailurePerSpec: false,
               timeoutInterval: 8*HOUR,
        defaultTimeoutInterval: 8*HOUR
      }
    },

    webpack: {
      mode: 'development',
      devtool: 'inline-source-map'
    },

    files: [
      'src/**/*_test.js'
    ],

    exclude: [
    ],

    preprocessors: {
      'src/**/*.js': ['webpack', 'sourcemap']
    },

    reporters: [
      'spec'
//      'progress'
    ],

    specReporter: {
      showSpecTiming: true 
    },

    browsers: [
      'Chrome',
      'Firefox',
    ],

    port: 9876,
    logLevel: config.LOG_INFO,
    autoWatch: true,
    singleRun: true,
    concurrency: Infinity,
    colors: true
  })
}
