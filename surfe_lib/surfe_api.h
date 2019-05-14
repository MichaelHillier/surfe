#ifndef SURFE_API
#define SURFE_API

#include <surfe_lib_module.h>  // macro for importing / exporting dll

#include <modeling_methods.h>
#include <continuous_property.h>
#include <lajaunie.h>
#include <single_surface.h>
#include <stratigraphic_surfaces.h>
#include <vector_field.h>

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
	GRBF_Modelling_Methods *method_;
	UI_Parameters params_;
	GRBF_Modelling_Methods* get_method(const UI_Parameters& params);
	vtkStructuredGrid *sgrid;

	bool have_interpolant_;
	bool evaluation_completed_;
	vtkDataObjectCollection *convert_constraints_to_vtk();
public:
	Surfe_API(const UI_Parameters& params) : params_(params)
	{
		sgrid = nullptr;
		have_interpolant_ = false;
		evaluation_completed_ = false;
		method_ = get_method(params_);
	}
	void AddInterfaceConstraint(const Interface& pt);
	void AddPlanarConstraint(const Planar& planar_pt);
	void AddTangentConstraint(const Tangent& tangent_pt);
	void AddInequalityConstraint(const Inequality& inequality_pt);
	void ComputeInterpolant();
	double EvaluateInterpolantAtPoint(
		const double &x, 
		const double &y, 
		const double &z);
	double *EvaluateVectorInterpolantAtPoint(
		const double &x,
		const double &y, 
		const double &z);
	void ConstructRegularGridOutput(
		const double &zmin,
		const double &zmax, 
		const double &resolution,
		const double &xy_percent_padding = 0);
	void ConstructRegularGridOutput(
		const double &xmin, const double &xmax,
		const double &ymin, const double &ymax,
		const double &zmin, const double &zmax,
		const double &resolution);
	vtkStructuredGrid *GetEvaluatedvtkStructuredGrid();
	vtkDataObjectCollection *GetConstraintsAndOutputAsVTKObjects();
	vtkPolyData *GetIsoSurfacesAsvtkPolyData();
};

#endif // SURFE_API
