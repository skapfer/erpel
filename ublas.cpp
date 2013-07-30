#include "ublas.hpp"
#include <fstream>

void load_matrix (matrix *x, const char *filename) {
    std::ifstream in (filename);
    double tmp;
    for (int i = 0; i != (int)x->size1 (); ++i)
    for (int j = 0; j != (int)x->size2 (); ++j) {
        in >> tmp;
        (*x)(i,j) = tmp;
        if (!in) abort ();
    }
    in >> tmp;
    if (in) {
        fprintf (stderr, "excess data in %s\n", filename);
        abort ();
    }
}

void set_up_zero (vector *dsy, size_t l) {
    *dsy = zero_vector (l);
}

void set_up_zero (std::vector <double> *dst, size_t l) {
    dst->resize (0u);
    dst->resize (l, 0.);
}
