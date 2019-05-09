#ifndef SURFE_API
#define SURFE_API

#include <modelling_input.h>
#include <modelling_parameters.h>
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
	Constraints constraints_;
	GRBF_Modelling_Methods* get_method(const model_parameters& params, const Constraints& constraints);
public:
	Surfe_API(const model_parameters& params, const Constraints& constraints) : params_(params), constraints_(constraints)
	{
		method_ = get_method(params_, constraints_);
	}
	bool ComputeInterpolant();
	double EvaluateInterpolantAtPoint(const double &x, const double &y, const double &z);
	std::vector<double> EvaluateVectorInterpolantAtPoint(const double &x, const double &y, const double &z);
};

#endif // SURFE_API
