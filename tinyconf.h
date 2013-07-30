#ifndef TINYCONF_H_INCLUDED 
#define TINYCONF_H_INCLUDED 

#include "config.hpp"
#include <string>
#ifndef NDEBUG
#include <ostream>
#endif

class Configuration {
public:
    typedef std::string string_t;

    Configuration ();
    Configuration (const Configuration &);
    Configuration (const string_t &conffilename);
    Configuration &operator= (const Configuration &);
    ~Configuration ();
    void init (const string_t &conffilename);

    bool has_section (const string_t &section) const;
    bool has_key     (const string_t &section, const string_t &key) const;

    string_t string (const string_t &section, const string_t &key) const;
    string_t string (const string_t &section, const string_t &key, const string_t &default_) const;
    bool    boolean (const string_t &section, const string_t &key) const;
    int     integer (const string_t &section, const string_t &key) const;
    double floating (const string_t &section, const string_t &key) const;

#ifndef NDEBUG
    void dump (std::ostream &) const;
#endif
private:
    void *mem;
    string_t value (const string_t &section, const string_t &key) const;
    void init (const Configuration &other);
    void term ();
};

//
// inline implementation
//
inline Configuration::Configuration () {
    mem = 0;
}

inline Configuration::Configuration (const std::string &conffilename) {
    mem = 0;
    init (conffilename);
}

inline Configuration::Configuration (const Configuration &other) {
    mem = 0;
    init (other);
}

inline Configuration &Configuration::operator= (const Configuration &other) {
    init (other);
    return *this;
}

inline Configuration::~Configuration () {
    term ();
}

#endif /* TINYCONF_H_INCLUDED */
