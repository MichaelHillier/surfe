#include <pybind11/pybind11.h>
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
		.def(py::init<>())
		.def(py::init<const Parameters &>())
		.def("GetParametersAndConstraints", &Surfe_API::GetParametersAndConstraints)
		.def("LoadConstraintsFromFiles", &Surfe_API::LoadConstraintsFromFiles)
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
		.def("SetModellingMode", &Surfe_API::SetModellingMode)
		.def("SetRegressionSmoothing", &Surfe_API::SetRegressionSmoothing)
		.def("SetGreedyAlgorithm", &Surfe_API::SetGreedyAlgorithm)
		.def("SetRestrictedRange", &Surfe_API::SetRestrictedRange)
		.def("SetRBFKernel",
		(void (Surfe_API::*)(const char *))
			&Surfe_API::SetRBFKernel, "Set the RBF kernel using a name")
		.def("SetRBFShapeParameter", &Surfe_API::SetRBFShapeParameter)
		.def("SetPolynomialOrder", &Surfe_API::SetPolynomialOrder)
		.def("SetGlobalAnisotropy", &Surfe_API::SetGlobalAnisotropy)
		.def("SetInterfaceUncertainty", &Surfe_API::SetInterfaceUncertainty)
		.def("SetAngularUncertainty", &Surfe_API::SetAngularUncertainty)
		.def("SetInterfaceDataFile", &Surfe_API::SetInterfaceDataFile)
		.def("SetPlanarDataFile", &Surfe_API::SetPlanarDataFile)
		.def("SetTangentDataFile", &Surfe_API::SetTangentDataFile)
		.def("SetInequalityDataFile", &Surfe_API::SetInequalityDataFile)
		.def("EvaluateInterpolantAtPoint", &Surfe_API::EvaluateInterpolantAtPoint)
		.def("EvaluateVectorInterpolantAtPoint",
			&Surfe_API::EvaluateVectorInterpolantAtPoint,
			py::return_value_policy::copy)
		.def("BuildRegularGrid",
		(void (Surfe_API::*)(const double&, const double &, const double &, const double &))
			&Surfe_API::BuildRegularGrid, "Build vtkImageData grid from zmin/zmax and resolution")
		.def("BuildRegularGrid",
		(void (Surfe_API::*)(const double&, const double &, const double &, const double &, const double &, const double&, const double &))
			&Surfe_API::BuildRegularGrid, "Build tkImageData grid from bounds and resolution")
		.def("GetEvaluatedVTKGridString", &Surfe_API::GetEvaluatedVTKGridString)
		.def("GetVTKIsosurfacesString", &Surfe_API::GetVTKIsosurfacesString)
		.def("GetVTKInterfaceConstraintsString", &Surfe_API::GetVTKInterfaceConstraintsString)
		.def("GetVTKTangentConstraintsString", &Surfe_API::GetVTKTangentConstraintsString)
		.def("GetVTKPlanarConstraintsString", &Surfe_API::GetVTKPlanarConstraintsString)
		.def("GetVTKInequalityConstraintString", &Surfe_API::GetVTKInequalityConstraintString)
		.def("WriteVTKInterfaceConstraints", &Surfe_API::WriteVTKInterfaceConstraints)
		.def("WriteVTKPlanarConstraints", &Surfe_API::WriteVTKPlanarConstraints)
		.def("WriteVTKTangentConstraints", &Surfe_API::WriteVTKTangentConstraints)
		.def("WriteVTKInequalityConstraints", &Surfe_API::WriteVTKInequalityConstraints)
		.def("WriteVTKEvaluationGrid", &Surfe_API::WriteVTKEvaluationGrid)
		.def("WriteVTKIsoSurfaces", &Surfe_API::WriteVTKIsoSurfaces)
		.def("VisualizeVTKData", &Surfe_API::VisualizeVTKData);
}