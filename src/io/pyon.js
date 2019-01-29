'use strict';

/* This file is part of ND.JS.
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


export function pyon_parse( char_seq )
{
  // next: returns the next character
  // skip: returns the next non-whitespace character
  const [next,skip,err] = function(){
    let end = false,
        counter = 0;
    const iter = char_seq[Symbol.iterator](),
          next = () => {
            const {value, done} = iter.next();
            ++counter;
            if(done) {
              if(end) throw new Error('pyon_parse: Character sequence ended unexpectedly.');
              end = true;
              return null;
            }
            return value;
          },
          skip = () => {
            for( let char;; )
              switch( char = next() ) {
                default: return char;
                case '\f':
                case '\n':
                case '\r':
                case '\t':
                case '\v':
                case  ' ':
              }
          },
          err = (encountered, expected) => {
            const prefix = expected == null ? 'Invalid character' : `Expected ${expected}, but `;
            throw new Error(`pyon_parse: ${prefix} '${encountered}' encountered as ${counter}-th character.`);
          };
    return [next, skip, err];
  }()

  function any() {
    switch(char) {
      default: err(char);
      case '{': return dict();
      case '"':
      case "'": return str();
      case '(':
      case '[': return list();
      case 'T': return True();
      case 'F': return False();
      case 'N': return None();
      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7': case '.':
      case '8': case '+':
      case '9': case '-': return num();
    }
  }

  function num()
  {
    let num = '';

    loop: for( ; ; char=next() )
      switch(char) {
        case '\f':
        case '\n':
        case '\r':
        case '\t':
        case '\v':
        case  ' ': char = skip();
        default  : break loop;
        case '0':
        case '1':
        case '2':
        case '3': case 'a': case 'A':
        case '4': case 'b': case 'B':
        case '5': case 'c': case 'C':
        case '6': case 'd': case 'D':
        case '7': case 'e': case 'E': case '.':
        case '8': case 'f': case 'F': case '+':
        case '9': case 'e': case 'E': case '-': num += char;
      }

    num = num.trim()*1;
    if( isNaN(num) ) err(num);
    return num;
  }

  function dict()
  {
    const result = {};
    if( '{' !== char ) throw new Error('Assertion failed.');
    for(;;) {
      if( '}' === (char = skip()) ) {
        char = skip();
        return result;
      }

      const key = any();

      if( ':' !== char )
        err(char, end);
      char = skip();

      result[key] = any();

      switch(char) {
        default: err(char, ",' or '}");
        case ',': continue;
        case '}':
          char = skip();
          return result;
      }
    }
  }

  function list()
  {
    const end = function(){
      switch(char) {
        case '(': return ')';
        case '[': return ']';
        default: throw new Error('Assertion failed.');
      }
    }();

    const result = [];

    for(;;) {
      switch( char = skip() ) {
        case end:
          char = skip();
          return result;
        default:
      }

      result.push( any() );

      switch(char) {
        default: err(char, `,' or '${end}`);
        case ',': continue;
        case end:
          char = skip();
          return result;
      }
    }
  }

  function True()
  {
    if( 'T' !== char ) throw new Error('Assertion failed.');
    for( const c of 'rue' ) {
      char = next();
      if( c !== char )
        err(char, c);
    }
    char = skip();
    return true;
  }

  function False()
  {
    if( 'F' !== char ) throw new Error('Assertion failed.');
    for( const c of 'alse' ) {
      char = next();
      if( c !== char )
        err(char, c);
    }
    char = skip();
    return false;
  }

  function str()
  {
    const END = char;
    switch(END) {
      default: throw new Error('Assertion failed.');
      case "'":
      case '"':
    }
    let result = '';
    for(;;)
      switch( char = next() )
      {
        case END:
          char = skip();
          return result;
        case '\\':
          switch( char = next() )
          {
            case 'u':
              let codePoint = 0;
              for( let i=0; i < 4; i++ )
              {
                codePoint <<= 8;
                switch( char = next() ) {
                  default: err(char, 'Hexadecimal digit');
                  case 'f': case 'F': ++codePoint;
                  case 'e': case 'E': ++codePoint;
                  case 'd': case 'D': ++codePoint;
                  case 'c': case 'C': ++codePoint;
                  case 'b': case 'B': ++codePoint;
                  case 'a': case 'A': ++codePoint;
                  case '9':           ++codePoint;
                  case '8':           ++codePoint;
                  case '7':           ++codePoint;
                  case '6':           ++codePoint;
                  case '5':           ++codePoint;
                  case '4':           ++codePoint;
                  case '3':           ++codePoint;
                  case '2':           ++codePoint;
                  case '1':           ++codePoint;
                  case '0':           ++codePoint;
                  
                }
              }
              result += String.fromCodePoint(codePoint);
              continue;
            case 'b': result += '\b'; continue;
            case 'f': result += '\f'; continue;
            case 'n': result += '\n'; continue;
            case 'r': result += '\r'; continue;
            case 't': result += '\t'; continue;
            case '"':
            case "'":
            case "\\": // FALL THROUGH
          }
        default: result += char;
      }
  }

  function None()
  {
    if( 'N' !== char ) throw new Error('Assertion failed.');
    for( const c of 'one' ) {
      char = next();
      if( c !== char )
        err(char, c);
    }
    char = skip();
    return null;
  }

  let char = skip();
  return any();
}
