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
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkVertex.h>
#include <vtkDataObjectCollection.h>
#include <vtkMarchingCubes.h>
#include <vtkNew.h>

class SURFE_LIB_EXPORT Surfe_API {
private:
	// members
	GRBF_Modelling_Methods *method_;
	UI_Parameters params_;
	vtkSmartPointer<vtkStructuredGrid> sgrid;
	vtkSmartPointer<vtkDataObjectCollection> model_collection_;
	vtkSmartPointer<vtkPolyData> iso_surfaces_;

	bool have_interpolant_;
	bool evaluation_completed_;
	bool have_method_;
	bool parameters_changed_;
	bool constraint_files_changed_;
	bool constraints_changed_;

	// methods
	GRBF_Modelling_Methods* get_method(const UI_Parameters& params);
	vtkDataObjectCollection *convert_constraints_to_vtk();
	void build_constraints_from_csv_files();
public:
	Surfe_API();
	Surfe_API(const UI_Parameters& params);
	void GetUIParameters();
	void AddInterfaceConstraint(const Interface& pt);
	void AddInterfaceConstraint(
		const double &x,const double &y,const double &z, 
		const double &level
	);
	void AddPlanarConstraint(const Planar& planar_pt);
	void AddPlanarConstraintwNormal(
		const double &x,const double &y,const double &z,
		const double &nx,const double &ny,const double &nz
	);
	void AddPlanarConstraintwStrikeDipPolarity(
		const double &x,const double &y,const double &z,
		const double &strike, const double &dip,const double &polarity
	);
	void AddPlanarConstraintwAzimuthDipPolarity(
		const double &x, const double &y, const double &z,
		const double &azimuth, const double &dip, const double &polarity
	);
	void AddTangentConstraint(const Tangent& tangent_pt);
	void AddTangentConstraint(
		const double &x, const double &y, const double &z,
		const double &tx, const double &ty,	const double &tz
	);
	void AddInequalityConstraint(const Inequality& inequality_pt);
	void AddInequalityConstraint(
		const double &x, const double &y,const double &z,
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
	void SetInterfaceDataCSVFile(const char *interface_file);
	void SetPlanarDataCSVFile(const char *planar_file);
	void SetTangentDataCSVFile(const char *tangent_file);
	void SetInequalityDataCSVFile(const char *inequality_file);
	double EvaluateInterpolantAtPoint(
		const double &x, const double &y, const double &z
	);
	double *EvaluateVectorInterpolantAtPoint(
		const double &x,const double &y, const double &z
	);
	void ConstructRegularGridOutput(
		const double &zmin,	const double &zmax, 
		const double &resolution, const double &xy_percent_padding = 0
	);
	void ConstructRegularGridOutput(
		const double &xmin, const double &xmax,
		const double &ymin, const double &ymax,
		const double &zmin, const double &zmax,
		const double &resolution);
	vtkStructuredGrid *GetEvaluatedvtkStructuredGrid();
	vtkDataObjectCollection *GetConstraintsAndOutputAsVTKObjects();
	vtkPolyData *GetIsoSurfacesAsvtkPolyData();
	void WriteVTKConstraints(const char *filename);
	void WriteVTKEvaluationGrid(const char *filename);
	void WriteVTKIsoSurfaces(const char *filename);
};

#endif // SURFE_API
