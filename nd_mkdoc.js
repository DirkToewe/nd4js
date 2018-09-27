{
  'use strict';

  const
    fs = require('fs'),
    nd = require('./nd.js'),
    html = nd.createDocHTML();
  fs.writeFile('./docs/doc.html', html.serialize(), err => { if(err) return console.log(err) });
}
