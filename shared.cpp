#include "shared.hpp"
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdexcept>
#include <errno.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <limits.h>

typedef std::string string_t;

// why not export it?
//namespace {
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
//}

double string_to_double (string_t value) {
    const char *beginptr = value.c_str ();
    char *endptr = 0;
    errno = 0;
    double ret = strtod (beginptr, &endptr);
    if (endptr != beginptr + value.size ())
        throw std::runtime_error (cpp_sprintf ("value %s is not a valid floating-point number", beginptr));
    if (errno == ERANGE)
        throw std::runtime_error (cpp_sprintf ("value %s cannot be represented in a double-precision float", beginptr));
    if (errno != 0)
        throw std::runtime_error (cpp_sprintf ("value %s is not a valid floating-point number", beginptr));
    return ret;
}

std::pair <double,double> split_double_pair (string_t dp) {
    string_t::size_type colon = dp.find (':');
    if (colon == dp.npos)
        throw std::runtime_error ("invalid moduli \"" + dp + "\"");
    double first = string_to_double (string_t (dp, 0, colon));
    double secnd = string_to_double (string_t (dp, colon+1, dp.npos));
    return std::make_pair (first, secnd);
}

bool ends_with (string_t x, string_t end) {
    string_t::size_type pos = x.rfind (end);
    return end.size () + pos == x.size ();
}

std::string outfile_prefix;

std::string make_filename (std::string msg, ...) {
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
    return outfile_prefix + ret;
}

void die (std::string msg, ...) {
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
    fprintf (stderr, "%s\n", ret.c_str ());
    exit (-1);
}

static void ignore (int) {  // avoid stupid warning about ignored result
}

void set_nice () {
    ignore (nice (10));
}

static std::string chop_x (std::string &x) {
    std::string::size_type end = x.find ('x');
    if (end == x.npos) {
        die ("[parse_XxYxZ] invalid string, contains too few x");
    } else {
        std::string ret;
        ret.assign (x, 0u, end);
        x.assign (x, end+1, x.npos);
        return ret;
    }
}

static unsigned long bulletproof_a2ui (std::string value) {
    const char *beginptr = value.c_str ();
    char *endptr = 0;
    errno = 0;
    unsigned long ret = strtoul (beginptr, &endptr, 10);
    if (endptr != beginptr + value.size ())
        throw std::runtime_error (cpp_sprintf ("value %s is not a valid unsigned integer", beginptr));
    if (errno == ERANGE)
        throw std::runtime_error (cpp_sprintf ("value %s cannot be represented in a unsigned long", beginptr));
    if (errno != 0)
        throw std::runtime_error (cpp_sprintf ("value %s is not a valid unsigned integer", beginptr));
    return ret;
}

void parse_XxYxZ (int *out, std::string in) {
    std::string eins = chop_x (in);
    std::string zwei = chop_x (in);
    std::string drei = in;

    unsigned long x1 = bulletproof_a2ui (eins);
    unsigned long x2 = bulletproof_a2ui (zwei);
    unsigned long x3 = bulletproof_a2ui (drei);

    if (x1 < 1 || x1 > INT_MAX)
        die ("[parse_XxYxZ] %s is either too small or too large", eins.c_str ());
    if (x2 < 1 || x2 > INT_MAX)
        die ("[parse_XxYxZ] %s is either too small or too large", zwei.c_str ());
    if (x3 < 1 || x3 > INT_MAX)
        die ("[parse_XxYxZ] %s is either too small or too large", drei.c_str ());

    *out++ = int (x1);
    *out++ = int (x2);
    *out++ = int (x3);
}

double inner_prod (const std::vector <double> &lhs, const std::vector <double> &rhs) {
    if (lhs.size () != rhs.size ()) {
        throw std::runtime_error ("inner_prod: sizes of vectors don't match");
    }

    double ret = 0.;
    for (int i = 0; i != (int)lhs.size (); ++i)
        ret += lhs[i]*rhs[i];
    return ret;
}

double lame_from_bulk_shear (double bulk, double shear) {
    return bulk - (2./3) * shear;
}

double bulk_from_lame_shear (double lame, double shear) {
    return lame + (2./3) * shear;
}

double poisson_from_lame_shear (double lame, double shear) {
    return lame * .5 / (lame+shear);
}

static time_t hptime_off = 0;

uint64_t perf_timer_class::hptime () {
    struct timeval x;
    if (gettimeofday (&x, 0) != 0) {
        perror ("gettimeofday");
        abort ();
    }

    if (hptime_off == 0)
        hptime_off = x.tv_sec;

    uint64_t ret = x.tv_sec - ::hptime_off;
    ret *= 1000000;
    ret += x.tv_usec;
    return ret;
}

void perf_timer_class::dump (FILE *fp) const {
    fprintf (fp, "timer class %s\n   average: %f msec,\n    total: %f sec\n",
             name_.c_str (), average_usec ()/1e3, total_usec ()/1e6);
}

void perf_timer_class::dump_all (FILE *fp) {
    for (iterator it = begin (); it != end (); ++it)
        (*it)->dump (fp);
}

double perf_timer_class::average_usec () const {
    return double (acc_)/count_;
}

double perf_timer_class::total_usec () const {
    return double (acc_);
}

const std::string &perf_timer_class::name () const {
    return name_;
}

static
std::vector <perf_timer_class *> &
perf_timer_list () {
    static std::vector <perf_timer_class *> bla;
    return bla;
}

perf_timer_class::perf_timer_class (const char *name)
    : name_(name) {

    acc_ = 0;
    count_ = 0;
    perf_timer_list ().push_back (this);
}

perf_timer_class::~perf_timer_class () {
    std::vector <perf_timer_class *>::iterator it;
    it = std::find (perf_timer_list ().begin (), perf_timer_list ().end (), this);
    perf_timer_list ().erase (it);
}

perf_timer_class::iterator perf_timer_class::begin () {
    return perf_timer_list ().begin ();
}

perf_timer_class::iterator perf_timer_class::end () {
    return perf_timer_list ().end ();
}
