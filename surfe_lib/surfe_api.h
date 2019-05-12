#ifndef SURFE_API
#define SURFE_API

#include <surfe_lib_module.h>  // macro for importing / exporting dll

#include <modeling_methods.h>
#include <continuous_property.h>
#include <lajaunie.h>
#include <single_surface.h>
#include <stratigraphic_surfaces.h>
#include <vector_field.h>

class SURFE_LIB_EXPORT Surfe_API {
private:
	GRBF_Modelling_Methods *method_;
	model_parameters params_;
	GRBF_Modelling_Methods* get_method(const model_parameters& params);
public:
	Surfe_API(const model_parameters& params) : params_(params)
	{
		method_ = get_method(params_);
	}
	void AddInterfaceConstraint(const Interface& pt);
	void AddPlanarConstraint(const Planar& planar_pt);
	void AddTangentConstraint(const Tangent& tangent_pt);
	void AddInequalityConstraint(const Inequality& inequality_pt);
	void ComputeInterpolant();
	double EvaluateInterpolantAtPoint(const double &x, const double &y, const double &z);
	std::vector<double> EvaluateVectorInterpolantAtPoint(const double &x, const double &y, const double &z);
	void ConstructRegularGridOutput(const double &zmin, const double &zmax, const double &resolution);
	void ConstructRegularGridOutput(
		const double &xmin, const double &xmax,
		const double &ymin, const double &ymax,
		const double &zmin, const double &zmax);
};

#endif // SURFE_API
