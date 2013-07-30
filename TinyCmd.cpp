
#include "config.hpp"
#include "TinyCmd.hpp"
#include <stdarg.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <errno.h>
#include <stdio.h>

typedef std::string string_t;

namespace {
    string_t cpp_sprintf (const string_t &msg, ...) {
        va_list va;
        va_start (va, msg);
        int bufsz = vsnprintf (0, 0, msg.c_str (), va);
        va_end (va);
        va_start (va, msg);
        string_t ret = string_t (bufsz+1, ' ');
#ifndef NDEBUG
        int bufsz2 =
#endif
        vsnprintf (&ret[0], bufsz+1, msg.c_str (), va);
        assert (bufsz == bufsz2);
        va_end (va);
        assert (ret[ret.size ()-1] == '\0');
        ret.erase (ret.size ()-1);
        return ret;
    }

    bool split_opt (string_t *name, string_t *value, const string_t &in) {
        string_t::size_type equals = in.find ("=");
        if (equals == in.npos) {
            return false;
        } else {
            name->assign (in, 0, equals);
            value->assign (in, equals+1, in.npos);
            return true;
        }
    }

    string_t strip (const string_t &str) {
        string_t ret = str;
        while (!ret.empty () && isspace (int (ret[0])))
            ret.assign (ret, 1, ret.npos);
        while (!str.empty () && isspace (int (ret[ret.size () - 1])))
            ret.assign (ret, 0, ret.size () - 1);
        return ret;
    }
}

ConfigSource::Exception::Exception (const string_t &msg)
    : std::runtime_error (msg) { }

CommandLine::CommandLine (int, char **argv) {
    ++argv;

    while (*argv) {
        string_t arg = *argv;

        if (arg.substr (0, 2) != "--") {
            throw Exception (cpp_sprintf ("non-option on commandline (NYI): %s", *argv));
        }

        arg.erase (0, 2);
        ++argv;

        string_t optname, optval;

        if (!split_opt (&optname, &optval, arg)) {
            optname = arg;
            if (!argv) {
                throw Exception (cpp_sprintf ("option without value: --%s", optname.c_str ()));
            }

            optval = *argv++;
        }

        if (was_given (optname)) {
            throw Exception (cpp_sprintf ("duplicated option: --%s", optname.c_str ()));
        }

        data_[optname] = OptionData_ (optval);
    }
}

ConfigSource::OptionData_::OptionData_ (const string_t &value_)
    : value (value_), used (false) { }

bool ConfigSource::was_given (const string_t &name) const {
    return data_.find (name) != data_.end ();
}

string_t ConfigSource::string (const string_t &name) {
    if (!was_given (name)) {
        throw Exception (cpp_sprintf ("missing option --%s", name.c_str ()));
    } else {
        data_[name].used = true;
        return data_[name].value;
    }
}

string_t ConfigSource::string (const string_t &name, const string_t &default_value) {
    if (was_given (name))
        return string (name);
    else
        return default_value;
}

double ConfigSource::floating (const string_t &name, double default_value) {
    if (was_given (name))
        return floating (name);
    else
        return default_value;
}

bool ConfigSource::boolean (const string_t &name, bool default_value) {
    if (was_given (name))
        return boolean (name);
    else
        return default_value;
}

int ConfigSource::integer (const string_t &name, int default_value) {
    if (was_given (name))
        return integer (name);
    else
        return default_value;
}

// slightly paranoid string-to-float
double ConfigSource::floating (const string_t &name) {
    string_t value = strip (string (name));
    for (;;) {
        if (value == "inf")
            return INFINITY;
        if (value == "-inf")
            return -INFINITY;
        if (value == "nan")
            return NAN;
        if (value == "")
            break;
        const char *beginptr = value.c_str ();
        char *endptr = 0;
        errno = 0;
        double ret = strtod (beginptr, &endptr);
        if (endptr != beginptr + value.size ())
            break;
        if (errno == ERANGE)
            throw Exception (cpp_sprintf ("value %s cannot be represented in a double-precision float", beginptr));
        if (errno != 0)
            break;
        return ret;
    }
    throw Exception (cpp_sprintf ("invalid floating point literal \"%s\" (option %s)", value.c_str (), name.c_str ()));
}

// slightly paranoid string-to-int
int ConfigSource::integer (const string_t &name) {
    string_t value = strip (string (name));
    for (;;) {
        if (value == "")
            break;
        const char *beginptr = value.c_str ();
        char *endptr = 0;
        errno = 0;
        // 10 forces base-10 numbers, since octal numbers surprise new users sometimes
        long ret = strtol (beginptr, &endptr, 10);
        if (endptr != beginptr + value.size ())
            break;
        if (errno == ERANGE)
            throw Exception (cpp_sprintf ("value %s cannot be represented in a long integer", beginptr));
        if (errno != 0)
            break;
        if (ret > INT_MAX || ret < INT_MIN)
            throw Exception (cpp_sprintf ("value %s cannot be represented in a integer", beginptr));
        return int (ret);
    }
    throw Exception (cpp_sprintf ("invalid floating point literal \"%s\" (option %s)", value.c_str (), name.c_str ()));
}

bool ConfigSource::boolean (const string_t &name) {
    string_t value = strip (string (name));
    if (value == "true" || value == "TRUE" || value == "True" || value == "1")
        return true;
    if (value == "false" || value == "FALSE" || value == "False" || value == "0")
        return false;
    throw Exception (cpp_sprintf ("invalid boolean literal \"%s\" (option %s)", value.c_str (), name.c_str ()));
}


int ConfigSource::num_unused_options () const {
    int ret = 0;
    
    std::map <string_t, OptionData_>::const_iterator it;
    for (it = data_.begin (); it != data_.end (); ++it) {
        if (it->second.used == false)
            ++ret;
    }

    return ret;
}

void ConfigSource::check_no_unused_options () const {
    std::map <string_t, OptionData_>::const_iterator it;
    for (it = data_.begin (); it != data_.end (); ++it) {
        if (it->second.used == false) {
            throw Exception (cpp_sprintf ("unused option --%s=%s", 
                             it->first.c_str (),
                             it->second.value.c_str ()));
        }
    }
}
