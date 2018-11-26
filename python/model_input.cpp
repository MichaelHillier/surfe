#include <pybind11/pybind11.h>
#include <pybind11/stl.h> //auto conversion b/w stl and python types

#include <Eigen/LU>
#include <modelling_input.h>
#include <continuous_property.h>
#include <modeling_methods.h>

// ----------------
// Python interface
// ----------------

namespace py = pybind11;
using namespace Surfe;
PYBIND11_MODULE(surfe,m){
    //setup bindings for the input objects
    py::class_<Point>(m,"Point")
    .def(py::init<const double&,const double&,const double&>())
    .def("x",&Point::x)
    .def("y",&Point::y)
    .def("z",&Point::z)
    .def("c",&Point::c)
    .def("set_x",&Point::set_x)
    .def("set_y",&Point::set_y)
    .def("set_z",&Point::set_z)
    .def("set_c",&Point::set_c)
    //.def("scalar_field",&Point::scalar_field)
    .def("set_scalar_field",&Point::set_scalar_field)
    //.def("scalar_field",&Point::scalar_field)
    .def("get_field_list_size",&Point::get_field_list_size)
    .def("set_vector_field",&Point::set_vector_field)
    .def("nx_interp",&Point::nx_interp)
    .def("ny_interp",&Point::ny_interp)
    .def("nz_interp",&Point::nz_interp);

    py::class_<Interface, Point>(m, "Interface")
    .def(py::init<const double&,const double&,const double&,const double&>())
    .def("level",&Interface::level)
    .def("residual",&Interface::residual)
    .def("level_lower_bound",&Interface::level_lower_bound)
    .def("level_upper_bound",&Interface::level_upper_bound)
    .def("setResidual",&Interface::setResidual)
    .def("setLevel",&Interface::setLevel)
    .def("setLevelBounds",&Interface::setLevelBounds);
   py::class_<Inequality, Point>(m, "Inequality")
    .def(py::init<const double&,const double&,const double&,const double&>())
    .def("level",&Inequality::level)
    .def("residual",&Inequality::residual)
    .def("setResidual",&Inequality::setResidual);
   py::class_<Evaluation_Point, Point>(m,"Evaluation_Point")
   .def(py::init<const double&,const double&,const double&>());
   //.def("scalar_field",&Evaluation_Point::scalar_field);
   py::class_<Basic_input>(m,"Basic_input")
   .def(py::init<>())
   .def_readwrite("interfaces",&Basic_input::itrface)
   .def_readwrite("evaluation_pts",&Basic_input::evaluation_pts)
   .def_readwrite("inequality",&Basic_input::inequality)
   .def_readwrite("planar",&Basic_input::planar)
   .def_readwrite("tangent",&Basic_input::tangent);
  py::class_<Continuous_Property>(m,"Continuous_Property")
  .def(py::init<const model_parameters&, const Basic_input&>())
  .def("run",&Continuous_Property::run_algorithm);
  py::class_<model_parameters>(m,"model_parameters")
  .def(py::init<>())
  .def_readwrite("model_type",&model_parameters::model_type)
  .def_readwrite("min_stratigraphic_thickness",&model_parameters::min_stratigraphic_thickness)
  .def_readwrite("use_interface_data",&model_parameters::use_interface_data)
  .def_readwrite("use_planar_data",&model_parameters::use_planar_data)
  .def_readwrite("use_tangent",&model_parameters::use_tangent)
  .def_readwrite("use_inequality",&model_parameters::use_inequality)
  .def_readwrite("basis_type",&model_parameters::basis_type)
  .def_readwrite("shape_parameter",&model_parameters::shape_parameter)
  .def_readwrite("polynomial_order",&model_parameters::polynomial_order)
  .def_readwrite("shape_parameter",&model_parameters::shape_parameter)
  .def_readwrite("advanced_parameters",&model_parameters::advanced_parameters)
  .def_readwrite("model_global_anisotropy",&model_parameters::model_global_anisotropy)
  .def_readwrite("use_greedy",&model_parameters::use_greedy)
  .def_readwrite("use_restricted_range",&model_parameters::use_restricted_range)
  .def_readwrite("use_greedy",&model_parameters::use_greedy)
  .def_readwrite("interface_uncertainty",&model_parameters::interface_uncertainty)
  .def_readwrite("angular_uncertainty",&model_parameters::angular_uncertainty);
}
