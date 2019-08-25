'use strict'

const fs = require('fs')

const version = require('./package.json').version
const license = `\
/* ND.JS v${version}
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
 `

for( const path of ['./dist/nd.js', './dist/nd.min.js'] )
{
  console.log(`Adding license to '${path}'...`);

  fs.readFile(path, 'utf-8', (err, contents) => {
    if(err)
      throw err
    const out = fs.createWriteStream(path)
    out.write(license,  'utf-8')
    out.write(contents, 'utf-8')
    out.end()
  });
}
