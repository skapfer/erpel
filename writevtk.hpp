#ifndef WRITE_VTK_ARRAY
#define WRITE_VTK_ARRAY

// dump a multi_array <double> into a VTK file

#include <boost/multi_array.hpp>
#include <string>

void write_to_vtk (std::string name, std::string filedesc,
                   const boost::multi_array <double, 3> &array);
void write_to_vtk (const char *name, const char *filedesc,
                   const boost::multi_array <double, 3> &array);

#endif // WRITE_VTK_ARRAY
