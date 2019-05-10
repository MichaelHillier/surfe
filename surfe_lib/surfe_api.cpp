#include <surfe_api.h>

GRBF_Modelling_Methods* Surfe_API::get_method(const model_parameters& params)
{
	// TODO
// 	if (params.model_type == Parameter_Types::Single_surface)
// 		return new Single_Surface(params, constraints);
// 	else if (params.model_type == Parameter_Types::Lajaunie_approach)
// 		return new Lajaunie_Approach(params, constraints);
// 	else if (params.model_type == Parameter_Types::Stratigraphic_horizons)
// 		return new Stratigraphic_Surfaces(params, constraints);
// 	else
// 		return new Continuous_Property(params, constraints);

	return nullptr;
}

bool Surfe_API::ComputeInterpolant()
{
	std::cout << " Starting SURFE algorithm " << std::endl;
	std::cout << " Processing input data...";
	if (!method_->process_input_data())
		return false;
	std::cout << "done!" << std::endl;
	std::cout << " Get method parameters...";
	if (!method_->get_method_parameters())
		return false;
	std::cout << "done!" << std::endl;
	std::cout << " Setup basis functions...";
	if (!method_->setup_basis_functions()) 
		return false;
	std::cout << "done!" << std::endl;
	std::cout << " Solve mathematical problem...";
	if (!method_->setup_system_solver()) {
		std::cout << "failed" << std::endl;
		return false;
	}
	std::cout << "done!" << std::endl;
}

double Surfe_API::EvaluateInterpolantAtPoint(const double &x, const double &y, const double &z)
{
	// convert x,y,z to Point
	Point pt(x, y, z);
	// evaluate scalar field at point
	method_->eval_scalar_interpolant_at_point(pt);
	return pt.scalar_field();
}

std::vector<double> Surfe_API::EvaluateVectorInterpolantAtPoint(const double &x, const double &y, const double &z)
{
	// convert x,y,z to Point
	Point pt(x, y, z);
	// evaluate vector field at point
	method_->eval_vector_interpolant_at_point(pt);
	std::vector<double> gradient;
	gradient.push_back(pt.nx_interp());
	gradient.push_back(pt.ny_interp());
	gradient.push_back(pt.nz_interp());
	return gradient;
}
