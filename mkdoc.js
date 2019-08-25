'use strict';
{
  console.log('Generating docs...');

  const nd = require('./dist/nd.js')
  const version = require('./package.json').version

  function create_doc_jsdom()
  {
    const
      jsdom= require('jsdom'),
      dom  = new jsdom.JSDOM('<!DOCTYPE html>'),
      doc  = dom.window.document,
      meta = doc.createElement('meta'),
      title= doc.createElement('title'),
      ul   = doc.createElement('ul'),
      h1   = doc.createElement('h1');
    h1   .innerHTML = `ND.JS v${version} Documentation`
    title.innerHTML = `ND.JS v${version} Documentation`
    meta.setAttribute('charset', 'utf-8')
    doc.lang = 'en'
    doc.head.appendChild(meta)
    doc.head.appendChild(title)
    doc.body.style = 'font-family:monospace'
    doc.body.appendChild(h1)
    doc.body.appendChild(ul)
    
    const ND = {...nd};
    delete ND.Array
    
    function addSection( key, val )
    {
      const
        h2 = doc.createElement('h2'),
        pre= doc.createElement('pre'),
        li = doc.createElement('li'),
        a  = doc.createElement('a');
    
      h2.id = key.replace(' ','_').replace('[','').replace(']','')
      a.innerHTML = key
      a.setAttribute('href', '#'+h2.id)
      ul.appendChild(li)
      li.appendChild(a)
    
      h2 .innerHTML = key
      pre.innerHTML = nd.help_str(val)
      doc.body.appendChild(h2)
      doc.body.appendChild(pre)
    }
    
    for( const [key,val] of Object.entries(ND) )
      addSection('nd.'+key,val)

    for( const [key,val] of Object.entries(ND.la) )
      if( key != '__doc__' )
        addSection('nd.la.'+key,val)

    for( const [key,val] of Object.entries(ND.opt) )
      if( key != '__doc__' )
        addSection('nd.opt.'+key,val)

    for( const [key,val] of Object.entries(ND.opt.line_search) )
      if( key != '__doc__' )
        addSection('nd.opt.line_search.'+key,val)

    for( const [key,val] of Object.entries(ND.io) )
      if( key != '__doc__' )
        addSection('nd.io.'+key,val)

    for( const [key,val] of Object.entries(ND.dt) )
      if( key != '__doc__' )
        addSection('nd.dt.'+key,val)

    const classes = [
      'nd.NDArray',
      'nd.dt.Complex'
    ]
    // ADD CLASS DOCS
    for( const name of classes )
    {
      const clazz = name.split('.').slice(1).reduce((obj,name) => obj[name], nd)
      
      for( const key of Object.getOwnPropertyNames(clazz) )
        if( ! key.startsWith('__') )
          try {
            addSection(name+'.'+key+' [STATIC]', clazz[key] )
          }
          catch(err) {}
    
      for( const key of Object.getOwnPropertyNames(clazz.prototype) )
        if( ! key.startsWith('__') )
          try {
            addSection(name+'.'+key, clazz.prototype[key] )
          }
          catch(err) {}
    }
    
    return dom
  }
  
  const
    fs = require('fs'),
    html = create_doc_jsdom();
  fs.writeFile('./docs/doc.html', html.serialize(), err => { if(err) return console.log(err) });
}
