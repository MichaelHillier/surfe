#ifndef SURFE_API
#define SURFE_API

#include <surfe_lib_module.h>  // macro for importing / exporting dll

#include <modeling_methods.h>
#include <continuous_property.h>
#include <lajaunie.h>
#include <single_surface.h>
#include <stratigraphic_surfaces.h>
#include <vector_field.h>

#include <inputImpl.h>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkStructuredGrid.h>
#include <vtkImageData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkVertex.h>
#include <vtkMarchingCubes.h>
#include <vtkNew.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetMapper.h>
#include <vtkPointGaussianMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkImageSliceMapper.h>
#include <vtkImageSlice.h>
#include <vtkImageMapToColors.h>
#include <vtkImageProperty.h>
#include <vtkImageActor.h>
#include <vtkLookupTable.h>
#include <vtkArrowSource.h>
#include <vtkAssignAttribute.h>
#include <vtkGlyph3D.h>
#include <vtkActor.h>
#include <vtkProperty.h>

class SURFE_LIB_EXPORT Surfe_API {
private:
	// members
	GRBF_Modelling_Methods *method_;
	UI_Parameters params_;
	vtkSmartPointer<vtkImageData> grid_;

	bool have_interpolant_;
	bool evaluation_completed_;
	bool have_method_;
	bool parameters_changed_;
	bool constraint_files_changed_;
	bool constraints_changed_;

	std::string vtk_grid_string_;
	std::string vtk_interface_string_;
	std::string vtk_planar_string_;
	std::string vtk_tangent_string_;
	std::string vtk_inequality_string_;
	std::string vtk_isosurfaces_string_;

	// methods
	GRBF_Modelling_Methods* get_method(const UI_Parameters& params);
	void build_constraints_from_input_files();
public:
	Surfe_API();
	Surfe_API(const UI_Parameters& params);
	void GetUIParametersAndConstraints();
	void LoadConstraintsFromFiles();
	void AddInterfaceConstraint(
		const double &x, const double &y, const double &z,
		const double &level
	);
	void AddPlanarConstraintwNormal(
		const double &x, const double &y, const double &z,
		const double &nx, const double &ny, const double &nz
	);
	void AddPlanarConstraintwStrikeDipPolarity(
		const double &x, const double &y, const double &z,
		const double &strike, const double &dip, const double &polarity
	);
	void AddPlanarConstraintwAzimuthDipPolarity(
		const double &x, const double &y, const double &z,
		const double &azimuth, const double &dip, const double &polarity
	);
	void AddTangentConstraint(
		const double &x, const double &y, const double &z,
		const double &tx, const double &ty, const double &tz
	);
	void AddInequalityConstraint(
		const double &x, const double &y, const double &z,
		const double &level
	);
	void ComputeInterpolant();
	void SetRBFKernel(const Parameter_Types::RBF &rbf);
	void SetRBFKernel(const char *rbf_name);
	void SetRBFShapeParameter(const double &shape_param);
	void SetPolynomialOrder(const int &poly_order);
	void SetGlobalAnisotropy(const bool &g_anisotropy);
	void SetGreedy(const bool &greedy);
	void SetRestrictedRange(const bool &rr);
	void SetRegressionSmoothing(const bool &rs);
	void SetInterfaceUncertainty(const double &interface_uncertainty);
	void SetAngularUncertainty(const double &angular_uncertainty);
	void SetInterfaceDataFile(const char *interface_file);
	void SetPlanarDataFile(const char *planar_file);
	void SetTangentDataFile(const char *tangent_file);
	void SetInequalityDataFile(const char *inequality_file);
	double EvaluateInterpolantAtPoint(
		const double &x, const double &y, const double &z
	);
	double *EvaluateVectorInterpolantAtPoint(
		const double &x, const double &y, const double &z
	); // client responsible for deleting dynamically allocated array vector[3]
	void BuildRegularGrid(
		const double &zmin, const double &zmax,
		const double &resolution, const double &xy_percent_padding = 0
	);
	void BuildRegularGrid(
		const double &xmin, const double &xmax,
		const double &ymin, const double &ymax,
		const double &zmin, const double &zmax,
		const double &resolution);
	vtkSmartPointer<vtkImageData> GetEvaluatedGrid();
	vtkSmartPointer<vtkPolyData> GetIsoSurfaces();
	vtkSmartPointer<vtkPolyData> GetInterfaceConstraints();
	vtkSmartPointer<vtkPolyData> GetPlanarConstraints();
	vtkSmartPointer<vtkPolyData> GetTangentConstraints();
	vtkSmartPointer<vtkPolyData> GetInequalityConstraints();
	const char *GetEvaluatedVTKGridString();
	const char *GetVTKIsosurfacesString();
	const char *GetVTKInterfaceConstraintsString();
	const char *GetVTKPlanarConstraintsString();
	const char *GetVTKTangentConstraintsString();
	const char *GetVTKInequalityConstraintString();
	void WriteVTKInterfaceConstraints(const char *filename);
	void WriteVTKPlanarConstraints(const char *filename);
	void WriteVTKTangentConstraints(const char *filename);
	void WriteVTKInequalityConstraints(const char *filename);
	void WriteVTKEvaluationGrid(const char *filename);
	void WriteVTKIsoSurfaces(const char *filename);
	void VisualizeVTKData();
};

#endif // SURFE_API
