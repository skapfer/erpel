#ifndef UBLAS_HPP_INCLUDED
#define UBLAS_HPP_INCLUDED

// some sensible uBLAS environment

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
namespace ublas = boost::numeric::ublas;
typedef ublas::matrix <double> matrix;
typedef ublas::vector <double> vector;
typedef ublas::zero_vector <double> zero_vector;
typedef ublas::zero_matrix <double> zero_matrix;

void load_matrix (matrix *x, const char *filename);

void set_up_zero (vector *, size_t);
void set_up_zero (std::vector <double> *, size_t);

#endif // UBLAS_HPP_INCLUDED
