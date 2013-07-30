// a few hacks to make things happen fast.
#ifndef PERF_HPP_INC
#define PERF_HPP_INC

#include <string.h>


namespace std {
    template <>
    inline
    void fill <double  *, double> (double *st, double *nd, const double &val) {
        if (val == 0.)
            memset (st, 0, (nd-st) * sizeof (double));
        else
            while (st != nd)
                *st++ = val;
    }
}

namespace from_blas {
    // just provide a prototype
    extern "C"
    void daxpy (int n, double alpha, const double *x, int incx, double *y, int incy);
    extern "C"
    double ddot (int n, const double *x, int incx, const double *y, int incy);
    extern "C"
    int idamax (int n, const double *x, int incx);
}

void add_to_vector (boost::numeric::ublas::vector <double> *dst,
                    double pref,
                    const boost::numeric::ublas::vector <double> &incr) {
#ifdef HAVE_BLAS
    from_blas::daxpy (dst->size (), pref, &incr(0), 1, &(*dst)(0), 1);
#else
    noalias (*dst) += pref * incr;
#endif
}

double serial_inner_prod (const boost::numeric::ublas::vector <double> &lhs,
                          const boost::numeric::ublas::vector <double> &rhs,
                          int length) {
#ifdef HAVE_BLAS
    double ret = from_blas::ddot (length, &lhs(0), 1, &rhs(0), 1);
#else
    double ret = inner_prod (subrange (lhs, 0, length),
                             subrange (rhs, 0, length));
#endif
    return ret;
}

double serial_norm_inf (const boost::numeric::ublas::vector <double> &lhs,
                        int length) {
#ifdef HAVE_BLAS
    int i = from_blas::idamax (length, &lhs(0), 1);
    double ret = lhs(i);
#else
    double ret = norm_inf (subrange (lhs, 0, length));
#endif
    return ret;
}

#endif // PERF_HPP_INC
