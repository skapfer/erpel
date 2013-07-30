// weird functions copied from old code
#ifndef SHARED_HPP_INCLUDED
#define SHARED_HPP_INCLUDED

#include "config.hpp"
#include <utility>
#include <string>
#include <stdint.h>
#include <stdio.h>
#include "ublas.hpp"
#include <boost/multi_array.hpp>

double string_to_double (std::string);
std::pair <double,double> split_double_pair (std::string);
bool ends_with (std::string x, std::string end);

extern std::string outfile_prefix;
std::string cpp_sprintf (const std::string &msg, ...);
std::string make_filename (std::string name, ...);
void die (std::string name, ...) does_not_return;
void set_nice ();

void parse_XxYxZ (int *, std::string);

double inner_prod (const std::vector <double> &, const std::vector <double> &);

double lame_from_bulk_shear (double bulk, double shear);
double bulk_from_lame_shear (double lame, double shear);
double poisson_from_lame_shear (double lame, double shear);

// perf timer is not thread-safe!
class perf_timer_class {
public:
    perf_timer_class (const char *name);
    ~perf_timer_class ();

    double average_usec () const;
    double total_usec () const;
    const std::string &name () const;

    static uint64_t hptime ();
    void dump (FILE *) const;
    static void dump_all (FILE *);

    typedef std::vector <perf_timer_class *>::const_iterator iterator;
    static iterator begin (), end ();
private:
    uint64_t acc_;
    std::string name_;
    int count_;
    friend class perf_timer;
};

class perf_timer {
public:
    perf_timer (perf_timer_class &x) : my_ (x) { x.acc_ -= perf_timer_class::hptime (); }
    ~perf_timer () { my_.count_++; my_.acc_ += perf_timer_class::hptime (); }
private:
    perf_timer_class &my_;
};

template <typename T>
static size_t mem_use (const boost::multi_array <T, 3> &ary) {
    return sizeof (T) * ary.shape ()[0] * ary.shape ()[1] * ary.shape ()[2];
}

template <typename T>
static size_t mem_use (const boost::numeric::ublas::vector <T> &ary) {
    return sizeof (T) * ary.size ();
}

template <typename T>
static size_t mem_use (const boost::numeric::ublas::matrix <T> &ary) {
    return sizeof (T) * ary.size1 () * ary.size2 ();
}

template <typename T>
static size_t mem_use (const std::vector <std::vector <T> > &ary) {
    size_t ret = size_t ();
    for (size_t i = size_t (); i != ary.size (); ++i)
        ret += ary[i].size () * sizeof (T);
    return ret + sizeof (std::vector <T>) * ary.size ();
}

template <typename T>
static size_t mem_use (const std::vector <T> &ary) {
    return sizeof (T) * ary.size ();
}

inline
int max3 (int a, int b, int c) {
    return std::max (a, std::max (b, c));
}

FILE *bzopen (std::string filename, std::string mode);

#endif // SHARED_HPP_INCLUDED
