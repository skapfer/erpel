
#include "tinyconf.h"
#include <map>
#include <fstream>
#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <stdexcept>
#include <stdlib.h>
#include <memory>
#include <limits.h>
using namespace std;

typedef map <string, string> submap_type;
typedef map <string, submap_type> map_type;

namespace {
    inline void nop () {
    }

    void syntax_error (const char *fmt, ...) {
        char msg[201];
        va_list va;
        va_start (va, fmt);
        vsnprintf (msg, 200, fmt, va);
        msg[200] = '\0';
        throw std::runtime_error (msg);
    }

    // your standard trim utility function
    string trim (const string &s) {
        // trim at beginning
        int i;
        for (i = 0; i != (int)s.size (); ++i) {
            if (s[i] == ' ' || s[i] == '\n' || s[i] == '\t')
                continue;
            else
                break;
        }
        // trim at end
        int j;
        for (j = (int)s.size (); j != 0; --j) {
            if (s[j-1] == ' ' || s[j-1] == '\n' || s[j-1] == '\t')
                continue;
            else
                break;
        }
        return string (s, i, j); // handles overlaps
    }

    // discard current line of input
    void skip_line (std::istream &is) {
        string dummy;
        getline (is, dummy);
    }

    bool read_key (submap_type *m, std::istream &is) {
        // check for end-of-section and end-of-file
        is >> ws;
        switch (is.peek ()) {
        case EOF:
        case '[':
            return false;
        case '#':
            skip_line (is);
            return read_key (m, is);
        }
        // read one line, split at the =
        string line;
        getline (is, line);
        string::size_type seppos = line.find ('=', 0);
        if (seppos == string::npos)
            syntax_error ("expected = in key");
        // insert into map
        string key = trim (string (line, 0, seppos));
        string value = trim (string (line, seppos+1, string::npos));
        (*m)[key] = value;
        return true;
    }

    bool read_section (map_type *m, std::istream &is) {
        // section header
        is >> ws;
        switch (is.peek ()) {
        case EOF:
            return false;
        case '#':
            skip_line (is);
            return read_section (m, is);
        case '[':
            is.get ();
            break;
        default:
            syntax_error ("expected section header [, got %c", is.peek ());
        }
        std::string section;
        getline (is, section, ']');
        if (section.find ('\n') != std::string::npos)
            syntax_error ("expected ], got \\n");
        if (!is)
            syntax_error ("expected ], got EOF");
        // read keys
        while (read_key (&(*m)[section], is))
            nop ();
        return true;
    }
}

void Configuration::init (const std::string &conffilename) {
    auto_ptr <map_type> m (new map_type);
    // read configuration
    ifstream is (conffilename.c_str ());
    if (!is)
        syntax_error ("Unable to open file \"%s\"", conffilename.c_str ());
    while (read_section (m.get (), is))
        nop ();

    // conf. read completely, no more exceptions possible.
    term ();
    mem = (void *)(m.release ());
}

void Configuration::init (const Configuration &other) {
    map_type *m = new map_type ( *(map_type *)other.mem );
    term ();
    mem = (void *)m;
}

void Configuration::term () {
    if (mem) {
        delete (map_type *)mem;
        mem = 0;
    }
}

bool Configuration::has_section (const string_t &section) const {
    map_type *m = (map_type *)mem;
    map_type::const_iterator iter = m->find (section);
    return iter != m->end ();
}

bool Configuration::has_key (const string_t &section, const string_t &key) const {
    map_type *m = (map_type *)mem;
    map_type::const_iterator iter = m->find (section);
    if (iter == m->end ())
        return false;
    submap_type::const_iterator iter2 = iter->second.find (key);
    return iter2 != iter->second.end ();
}

string Configuration::value (const string_t &section, const string_t &key) const {
    map_type *m = (map_type *)mem;
    map_type::const_iterator iter = m->find (section);
    if (iter == m->end ())
        throw std::runtime_error ("No such section \"" + section + "\"");
    submap_type::const_iterator iter2 = iter->second.find (key);
    if (iter2 == iter->second.end ())
        throw std::runtime_error ("No such key \"" + key
                                  + "\" in section \"" + section + "\"");
    return iter2->second;
}

string Configuration::string (const string_t &section, const string_t &key) const {
    return value (section, key);
}

string Configuration::string (const string_t &section, const string_t &key, const string_t &def) const {
    if (has_key (section, key))
        return value (section, key);
    else
        return def;
}

bool Configuration::boolean (const string_t &section, const string_t &key) const {
    string_t v = value (section, key);
    if (v == "false" || v == "FALSE" || v == "0")
        return false;
    else if (v == "true" || v == "TRUE" || v == "1")
        return true;
    throw std::runtime_error ("Invalid boolean \"" + v + "\"");
}

int Configuration::integer (const string_t &section, const string_t &key) const {
    string_t v = value (section, key);
    errno = 0;
    long x = strtol (v.c_str (), NULL, 0);
    if (errno != 0 || x > INT_MAX || x < INT_MIN)
        throw std::runtime_error ("Invalid integer \"" + v + "\"");
    return x;
}

double Configuration::floating (const string_t &section, const string_t &key) const {
    string_t v = value (section, key);
    errno = 0;
    double x = strtod (v.c_str (), NULL);
    if (errno != 0)
        throw std::runtime_error ("Invalid floating-point number \"" + v + "\"");
    else
        return x;
}

#ifndef NDEBUG
void Configuration::dump (ostream &os) const {
    os << "Configuration " << (void *)this << "\n";
    map_type *m = (map_type *)mem;
    map_type::const_iterator iter;
    for (iter = m->begin (); iter != m->end (); ++iter) {
        os << "[" << iter->first << "]\n";
        submap_type::const_iterator iter2 = iter->second.begin ();
        submap_type::const_iterator end2  = iter->second.end ();
        for (; iter2 != end2; ++iter2)
            os << iter2->first << " = " << iter2->second << "\n";
    }
}
#endif

