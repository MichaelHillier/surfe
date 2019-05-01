#ifndef EXPORT_TO_VTK_H
#define EXPORT_TO_VTK_H
#include <modelling_input.h>
#include <string>
#include <vector>

namespace Surfe {
void SURFE_LIB_EXPORT write_to_vtk(Basic_input input, std::string file_name);
}
#endif /* EXPORT_TO_VTK_H */
