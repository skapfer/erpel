#ifndef TINYCMD_HPP_INCLUDED
#define TINYCMD_HPP_INCLUDED

#include <string>
#include <stdexcept>
#include <map>

class ConfigSource {
protected:
    typedef std::string string_t;
public:
    virtual ~ConfigSource () {}

    bool was_given  (const string_t &name) const;

    // extract option values
    string_t string (const string_t &name);
    string_t string (const string_t &name, const string_t &default_value);
    // do the boring conversion work
    double floating (const string_t &name);
    double floating (const string_t &name, double default_value);
    int     integer (const string_t &name);
    int     integer (const string_t &name, int default_value);
    bool    boolean (const string_t &name);
    bool    boolean (const string_t &name, bool default_value);

    int num_unused_options () const;

    // check that every argument was used at least once (catch typos)
    void check_no_unused_options () const;

    class Exception : public std::runtime_error {
    public:
        Exception (const string_t &);
    };

protected:
    struct OptionData_ {
        OptionData_ () { }
        OptionData_ (const string_t &);
        string_t value;
        bool used;
    };
    std::map <string_t, OptionData_> data_;
    typedef std::map <string_t, OptionData_>::iterator iter_t_;
};

class CommandLine : public ConfigSource {
public:
    CommandLine (int argc, char **argv);
    virtual ~CommandLine () {}
};


#endif // TINYCMD_HPP_INCLUDED
