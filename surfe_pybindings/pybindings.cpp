#include <pybind11/pybind11.h>
#include <pybind11/stl.h>       //auto conversion b/w stl and python types
#include <pybind11/iostream.h>  //auto conversion b/w stl and python types

#include <continuous_property.h>
#include <modeling_methods.h>
#include <modelling_input.h>
#include <modelling_methods_builder.h>
#include <surfe_api.h>

#include <Eigen/LU>

// ----------------
// Python interface
// ----------------

namespace py = pybind11;
PYBIND11_MODULE(surfe, m) {
	// setup bindings for Surfe_API
	py::class_<Surfe_API>(m, "Surfe_API")
		.def(py::init<const model_parameters &>())
		.def("ComputeInterpolant", &Surfe_API::ComputeInterpolant)
		.def("EvaluateInterpolantAtPoint", &Surfe_API::EvaluateInterpolantAtPoint)
		.def("EvaluateVectorInterpolantAtPoint", &Surfe_API::EvaluateVectorInterpolantAtPoint);

    // setup bindings for the input objects
    py::class_<Point>(m, "Point")
        .def(py::init<const double &, const double &, const double &>())
        .def("x", &Point::x)
        .def("y", &Point::y)
        .def("z", &Point::z)
        .def("c", &Point::c)
        .def("set_x", &Point::set_x)
        .def("set_y", &Point::set_y)
        .def("set_z", &Point::set_z)
        .def("set_c", &Point::set_c)
        .def("scalar_field", (double (Point::*)() const) & Point::scalar_field,
             "get scalar field value of point")
        .def("set_scalar_field", &Point::set_scalar_field)
        .def("scalar_field",
             (double (Point::*)(const int &)const) & Point::scalar_field,
             "get scalar field using index")
        .def("set_vector_field", &Point::set_vector_field)
        .def("nx_interp", &Point::nx_interp)
        .def("ny_interp", &Point::ny_interp)
        .def("nz_interp", &Point::nz_interp);

    py::class_<Interface, Point>(m, "Interface")
        .def(py::init<const double &, const double &, const double &,
                      const double &>())
        .def("level", &Interface::level)
        .def("residual", &Interface::residual)
        .def("level_lower_bound", &Interface::level_lower_bound)
        .def("level_upper_bound", &Interface::level_upper_bound)
        .def("setResidual", &Interface::setResidual)
        .def("setLevel", &Interface::setLevel)
        .def("setLevelBounds", &Interface::setLevelBounds);
    py::class_<Inequality, Point>(m, "Inequality")
        .def(py::init<const double &, const double &, const double &,
                      const double &>())
        .def("level", &Inequality::level)
        .def("residual", &Inequality::residual)
        .def("setResidual", &Inequality::setResidual);

    py::class_<Planar, Point>(m, "Planar")
        .def(py::init<const double&, const double&, const double&, const double&,
                      const double&, const double&>())
		.def(py::init<const double&, const double&, const double&, const double&,
			const double&, const int&>())
        //.def("getDipVector",&Planar::getDipVector)
        //.def("getStrikeVector",&Planar::getStrikeVector)
        .def("dip", &Planar::dip)
        .def("strike", &Planar::strike)
        .def("polarity", &Planar::polarity)
        .def("nx", &Planar::nx)
        .def("ny", &Planar::ny)
        .def("nz", &Planar::nz)
        .def("nx_lower_bound", &Planar::nx_lower_bound)
        .def("nx_upper_bound", &Planar::nx_upper_bound)
        .def("ny_lower_bound", &Planar::ny_lower_bound)
        .def("ny_upper_bound", &Planar::ny_upper_bound)
        .def("nz_lower_bound", &Planar::nz_lower_bound)
        .def("nz_upper_bound", &Planar::nz_upper_bound)
        .def("setNormalBounds", &Planar::setNormalBounds)
        .def("residual", &Planar::residual)
        .def("setResidual", &Planar::setResidual)
        .def("setNormal", &Planar::setNormal);
    py::class_<Tangent, Point>(m, "Tangent")
        .def(py::init<const double&, const double&, const double&, const double&,
                      const double&, const double&>())
        .def("tx", &Tangent::tx)
        .def("ty", &Tangent::ty)
        .def("tz", &Tangent::tz)
        .def("residual", &Tangent::residual)
        .def("angle_lower_bound", &Tangent::angle_lower_bound)
        .def("angle_upper_bound", &Tangent::angle_upper_bound)
        .def("inner_product_constraint", &Tangent::inner_product_constraint)
        .def("setResidual", &Tangent::setResidual)
        .def("setAngleBounds", &Tangent::setAngleBounds)
        .def("setInnerProductConstraint", &Tangent::setInnerProductConstraint);

    py::class_<model_parameters>(m, "model_parameters")
        .def(py::init<>())
        .def_readwrite("model_type", &model_parameters::model_type)
        .def_readwrite("min_stratigraphic_thickness",
                       &model_parameters::min_stratigraphic_thickness)
        .def_readwrite("use_interface_data",
                       &model_parameters::use_interface_data)
        .def_readwrite("use_planar_data", &model_parameters::use_planar_data)
        .def_readwrite("use_tangent", &model_parameters::use_tangent)
        .def_readwrite("use_inequality", &model_parameters::use_inequality)
        .def_readwrite("basis_type", &model_parameters::basis_type)
        .def_readwrite("shape_parameter", &model_parameters::shape_parameter)
        .def_readwrite("polynomial_order", &model_parameters::polynomial_order)
        .def_readwrite("shape_parameter", &model_parameters::shape_parameter)
        .def_readwrite("advanced_parameters",
                       &model_parameters::advanced_parameters)
        .def_readwrite("model_global_anisotropy",
                       &model_parameters::model_global_anisotropy)
        .def_readwrite("use_greedy", &model_parameters::use_greedy)
        .def_readwrite("use_restricted_range",
                       &model_parameters::use_restricted_range)
        .def_readwrite("use_greedy", &model_parameters::use_greedy)
        .def_readwrite("interface_uncertainty",
                       &model_parameters::interface_uncertainty)
        .def_readwrite("angular_uncertainty",
                       &model_parameters::angular_uncertainty);

    // bind parameter type enums
    py::class_<Parameter_Types> Parameter_Types(m, "Parameter_Types");
    py::enum_<Parameter_Types::ModelType>(Parameter_Types, "ModelType")
        .value("Single_surface", Parameter_Types::ModelType::Single_surface)
        .value("Lajaunie_approach",
               Parameter_Types::ModelType::Lajaunie_approach)
        .value("Stratigraphic_horizons",
               Parameter_Types::ModelType::Stratigraphic_horizons)
        .value("Continuous_property",
               Parameter_Types::ModelType::Continuous_property)
        .value("Vector_field", Parameter_Types::ModelType::Vector_field)
        .export_values();
    py::enum_<Parameter_Types::DWRT>(Parameter_Types, "DWRT")
        .value("PT1", Parameter_Types::DWRT::PT1)
        .value("PT2", Parameter_Types::DWRT::PT2)
        .export_values();
    py::enum_<Parameter_Types::SecondDerivatives>(Parameter_Types,
                                                  "SecondDerivatives")
        .value("DXDX", Parameter_Types::SecondDerivatives::DXDX)
        .value("DXDY", Parameter_Types::SecondDerivatives::DXDY)
        .value("DXDZ", Parameter_Types::SecondDerivatives::DXDZ)
        .value("DYDX", Parameter_Types::SecondDerivatives::DYDX)
        .value("DYDY", Parameter_Types::SecondDerivatives::DYDY)
        .value("DYDZ", Parameter_Types::SecondDerivatives::DYDZ)
        .value("DZDX", Parameter_Types::SecondDerivatives::DZDX)
        .value("DZDY", Parameter_Types::SecondDerivatives::DZDY)
        .value("DZDZ", Parameter_Types::SecondDerivatives::DZDZ)
        .export_values();
    py::enum_<Parameter_Types::FirstDerivatives>(Parameter_Types,
                                                 "FirstDerivatives")
        .value("DX", Parameter_Types::FirstDerivatives::DX)
        .value("DY", Parameter_Types::FirstDerivatives::DY)
        .value("DZ", Parameter_Types::FirstDerivatives::DZ)
        .export_values();
    py::enum_<Parameter_Types::RBF>(Parameter_Types, "RBF")
        .value("Cubic", Parameter_Types::RBF::Cubic)
        .value("Gaussian", Parameter_Types::RBF::Gaussian)
        .value("MQ", Parameter_Types::RBF::MQ)
        .value("IMQ", Parameter_Types::RBF::IMQ)
        .value("TPS", Parameter_Types::RBF::TPS)
        .value("R", Parameter_Types::RBF::R)
        .export_values();
}
