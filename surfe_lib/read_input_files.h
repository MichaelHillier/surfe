#ifndef READ_INPUT_FILES_H
#define READ_INPUT_FILES_H

#include <csv.h>
#include <vector>
#include <modelling_input.h>
#include <grbf_exceptions.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataReader.h>
#include <vtkNew.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

std::string get_file_extension(const char *filename);

std::vector<Interface> build_interface_constraints(const char *interface_file);
std::vector<Planar> build_planar_constraints(const char *planar_file);
std::vector<Tangent> build_tangent_constraints(const char *tangent_file);
std::vector<Inequality> build_inequality_constraints(const char *inequality_file);

#endif
