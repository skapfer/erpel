
// dump a multi_array <double> into a VTK file

#include "config.hpp"
#include "writevtk.hpp"
#include <fstream>
#include <string>
#include <stdlib.h>
#include <math.h>

static
void write_vtk_header (std::ostream &os,
                       int nx, int ny, int nz,
                       const std::string &comment) {
    double dx = 1.;
    double dy = 1.;
    double dz = 1.;
    const std::string c_datatype = "double";
    os << "# vtk DataFile Version 2.0\n"
       << comment << "\nASCII\n"
       << "DATASET RECTILINEAR_GRID\nDIMENSIONS "
       << nx << " " << ny << " " << nz;
    os << "\nX_COORDINATES " << nx << " " << c_datatype << "\n";
    for (int i = 0; i != nx; ++i)
        os << (i+.5) * dx << " ";
    os << "\nY_COORDINATES " << ny << " " << c_datatype << "\n";
    for (int i = 0; i != ny; ++i)
        os << (i+.5) * dy << " ";
    os << "\nZ_COORDINATES " << nz << " " << c_datatype << "\n";
    for (int i = 0; i != nz; ++i)
        os << (i+.5) * dz << " ";
    os << "\n";
}

void write_to_vtk (const char *name, const char *filedesc,
                   const boost::multi_array <double, 3> &array) {

    std::ofstream ofs (name);
    int nx = array.shape ()[0];
    int ny = array.shape ()[1];
    int nz = array.shape ()[2];

    write_vtk_header (ofs, nx, ny, nz, filedesc);
    const std::string dataname = "default";
    const std::string datatype = "double";
    ofs << "POINT_DATA " << (nx*ny*nz)
        << "\nSCALARS " << dataname << " " << datatype << "\n"
        << "LOOKUP_TABLE default\n";
         
    for (int i = 0; i != nx; ++i)
    for (int j = 0; j != ny; ++j)
    for (int k = 0; k != nz; ++k) {
        if (isnan (array[i][j][k]))
            ofs << "0.0 ";
        else
            ofs << array[i][j][k] << " ";
    }
}

void write_to_vtk (std::string name, std::string filedesc,
                   const boost::multi_array <double, 3> &array) {
    write_to_vtk (name.c_str (), filedesc.c_str (), array);
}
