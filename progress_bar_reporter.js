'use strict';

/* This file is part of ND4JS.
 *
 * ND4JS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ND4JS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ND4JS. If not, see <http://www.gnu.org/licenses/>.
 */

require('colors');


const SUBSTEPS_UTF8 = [
  ' ',
  '▏',
  '▎',
  '▍',
  '▌',
  '▋',
  '▊',
  '▉',
  '█'
].join('');
const SUBSTEPS_ASCII = [
  ' ',
  '-',
  '=',
  '#'
].join('');


const formatTime = milisecs => {
  const h = Math.floor(milisecs / (60*60*1000)),
        m = Math.floor(milisecs / (   60*1000)) - h*60,
        s =            milisecs / (      1000)  - h*60*60 - m*60;
//  return `${h}:${m.toString().padStart(2,'0')}:${s.toFixed(3).padStart(6,'0')}`;
  return `${h}:${m.toString().padStart(2,'0')}:${Math.round(s).toString().padStart(2,'0')}`;
}


function ProgressBarReporter(baseReporterDecorator, config)
{
  baseReporterDecorator(this);

  const {
    barLength = 32,
     useAscii = false
  } = config.progressBarReporter || {};

  const SUBSTEPS = useAscii
    ? SUBSTEPS_ASCII
    : SUBSTEPS_UTF8;

      this.USE_COLORS = !!config.colors;
  if( this.USE_COLORS )
    Object.assign(this,{
      LOG_SINGLE_BROWSER    :    '%s: ' + '%s'.cyan + '\n',
      LOG_MULTI_BROWSER     : '%s %s: ' + '%s'.cyan + '\n',
      SPEC_FAILURE          : '%s %s FAILED'.red + '\n',
      SPEC_SLOW             : '%s SLOW %s: %s'.yellow + '\n',
      ERROR                 : '%s ERROR'.red + '\n',
      FINISHED_CRASH        : ' CRASH'.red,
      FINISHED_DISCONNECTED : ' DISCONNECTED'.red,
      FINISHED_ERROR        : ' ERROR'.red,
      FINISHED_SUCCESS      : ' SUCCESS'.green
    });
  else
    Object.assign(this,{
      LOG_SINGLE_BROWSER    :    '%s: ' + '%s' + '\n',
      LOG_MULTI_BROWSER     : '%s %s: ' + '%s' + '\n',
      SPEC_FAILURE          : '%s %s FAILED' + '\n',
      SPEC_SLOW             : '%s SLOW %s: %s' + '\n',
      ERROR                 : '%s ERROR' + '\n',
      FINISHED_CRASH        : ' CRASH',
      FINISHED_DISCONNECTED : ' DISCONNECTED',
      FINISHED_ERROR        : ' ERROR',
      FINISHED_SUCCESS      : ' SUCCESS'
    });

  this.renderBrowser = function(browser) {
    const     results = browser.lastResult,
      nDone = results.success + results.failed,
      nTotal= results.total   - results.skipped;

    const crashed = results.totalTime && nDone !== nTotal,
       browserLen = this._browsers.reduce((len,b) => Math.max(len,b.toString().length), 0);

    // Add Browser Info
    let msg = `${browser.toString().padEnd(browserLen)}`;
    if( this.USE_COLORS )
      msg = msg.green;

    // Add progress bar
    ;{
      const                             nSub = SUBSTEPS.length-1,
        p = Math.floor( nDone*barLength*nSub / nTotal );
      let bar = SUBSTEPS[nSub].repeat(p/nSub | 0);
      if( 0 !== p%nSub )
        bar += SUBSTEPS[p%nSub];
      bar = bar.padEnd(barLength,' ');

      if( this.USE_COLORS ) {                               bar = bar.bgBlack;
             if(results.failed || results.error || crashed) bar = bar.red;
        else if( nDone === nTotal                         ) bar = bar.green;
        else                                                bar = bar.blue;
      }

      let lBound = useAscii ? '[' : '▕',
          rBound = useAscii ? ']' : '▏';
      if( this.USE_COLORS && useAscii ) {
        lBound = lBound.bgBlack;
        rBound = rBound.bgBlack;
      }
      msg += lBound;
      msg += bar;
      msg += rBound;
    }

    // Add test progress numbers
    msg += `${nDone.toString().padStart(4)}/${results.total.toString().padEnd(4)}`;

    let skips = '',
        fails = '';
    if(results.failed ) fails = ` ${results.failed } FAIL${results.failed  > 1 ? 'S' : ''}`;
    if(results.skipped) skips = ` ${results.skipped} SKIP${results.skipped > 1 ? 'S' : ''}`;
    if( this.USE_COLORS ) {
      fails = fails.red;
      skips = skips.yellow;
    }

    if( browser.isConnected ) {
      msg += `(${formatTime(results.totalTime || results.netTime)})`;

           if( results.disconnected               ) msg += this.FINISHED_DISCONNECTED;
      else if( results.error  && !results.failed  ) msg += this.FINISHED_ERROR;
      else if( crashed                            ) msg += this.FINISHED_CRASH;
      else if(!results.failed && nDone === nTotal ) msg += this.FINISHED_SUCCESS;
    }

    msg += fails;
    msg += skips;

    return msg
  }

  this.writeCommonMsg = function (msg) {
    this.write(this._remove() + msg + this._render());
  }

  this.specSuccess = function () {
    this.writeCommonMsg('');
  }

  this.onBrowserStart = function (browser) {
    this._browsers.push(browser)
    this.writeCommonMsg('');
  }

  this.onBrowserComplete = function () {
    this.writeCommonMsg('');
  }

  this.onRunStart = function () {
    this._browsers = []
    this._nRendered = 0;

    const stdout_write = process.stdout.write;
    this._stdout_write = Object.hasOwnProperty(process.stdout, 'write') ? stdout_write : null;
    const self = this;
    process.stdout.write = function(msg, ...args) {
      return stdout_write.apply(this, [self._remove() +  msg + self._render(), ...args]);
    };
  }

  this.onRunComplete = function () {
    console.assert('_stdout_write' in this);
    if( null != this._stdout_write ) process.stdout = this._stdout_write; else delete process.stdout.write;
    delete this._stdout_write;
  }

  this._remove = function () {
    let cmd = '';
    for( let i=this._nRendered; i-- > 0; )
      cmd += '\x1B[1A' + '\x1B[2K';
    return cmd;
  }

  this._render = function () {
    this._nRendered =  this._browsers.length;
    const sorted = [...this._browsers].sort( (a,b) => {
      const x = a.toString(),
            y = b.toString(),   result = (x>y) - (x<y);
      if( 0 !== result ) return result;
      a = a.lastResult.total;
      b = b.lastResult.total;
      return b-a;
    });
    return sorted.map( b => this.renderBrowser(b) + '\n' ).join('');
  }
};

ProgressBarReporter.$inject = ['baseReporterDecorator', 'config'];

module.exports = {
  'reporter:progress-bar': ['type', ProgressBarReporter]
};
