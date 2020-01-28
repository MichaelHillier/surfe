#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>       //auto conversion b/w stl and python types
#include <pybind11/iostream.h>  //auto conversion b/w stl and python types

#include <continuous_property.h>
#include <modeling_methods.h>
#include <modelling_input.h>
#include <surfe_api.h>

#include <Eigen/LU>

// ----------------
// Python interface
// ----------------

namespace py = pybind11;
PYBIND11_MODULE(surfepy, m) {
	// setup bindings for Surfe_API
	py::class_<Surfe_API>(m, "Surfe_API")
		.def(py::init<const int>())
		.def(py::init<const Parameters &>())
		.def("AddInterfaceConstraint",
		(void (Surfe_API::*)(const double&, const double&, const double&, const double&))
			&Surfe_API::AddInterfaceConstraint, "Add an interface constraint")
		.def("AddPlanarConstraintwNormal", &Surfe_API::AddPlanarConstraintwNormal)
		.def("AddPlanarConstraintwStrikeDipPolarity", &Surfe_API::AddPlanarConstraintwStrikeDipPolarity)
		.def("AddPlanarConstraintwAzimuthDipPolarity", &Surfe_API::AddPlanarConstraintwAzimuthDipPolarity)
		.def("AddTangentConstraint",
		(void (Surfe_API::*)(const double&, const double&, const double&, const double&, const double&, const double&))
			&Surfe_API::AddTangentConstraint, "Add a tangent constraint")
		.def("AddInequalityConstraint",
		(void (Surfe_API::*)(const double&, const double&, const double&, const double&))
			&Surfe_API::AddInequalityConstraint, "Add an inequality constraint")
		.def("ComputeInterpolant", &Surfe_API::ComputeInterpolant)
		.def("SetRegressionSmoothing", &Surfe_API::SetRegressionSmoothing)
		.def("SetGreedyAlgorithm", &Surfe_API::SetGreedyAlgorithm)
		.def("SetRestrictedRange", &Surfe_API::SetRestrictedRange)
		.def("SetRBFKernel",
		(void (Surfe_API::*)(const char *))
			&Surfe_API::SetRBFKernel, "Set the RBF kernel using a name")
		.def("SetRBFShapeParameter", &Surfe_API::SetRBFShapeParameter)
		.def("SetPolynomialOrder", &Surfe_API::SetPolynomialOrder)
		.def("SetGlobalAnisotropy", &Surfe_API::SetGlobalAnisotropy)
		.def("EvaluateInterpolantAtPoint", &Surfe_API::EvaluateInterpolantAtPoint)
		.def("EvaluateVectorInterpolantAtPoint", &Surfe_API::EvaluateVectorInterpolantAtPoint)
		.def("GetDataBoundsAndResolution", &Surfe_API::GetDataBoundsAndResolution)
		.def("GetInterfaceReferencePoints", &Surfe_API::GetInterfaceReferencePoints)
		.def("GetInterfaceConstraints", &Surfe_API::GetInterfaceConstraints)
		.def("SetInterfaceConstraints", &Surfe_API::GetInterfaceConstraints)
		.def("GetPlanarConstraints", &Surfe_API::GetPlanarConstraints)
		.def("SetPlanarConstraints", &Surfe_API::SetPlanarConstraints)
		.def("GetTangentConstraints", &Surfe_API::GetTangentConstraints)
		.def("SetTangentConstraints", &Surfe_API::SetTangentConstraints)
		.def("GetInequalityConstraints", &Surfe_API::GetInequalityConstraints)
		.def("SetInequalityConstraints", &Surfe_API::SetInequalityConstraints);
		
}