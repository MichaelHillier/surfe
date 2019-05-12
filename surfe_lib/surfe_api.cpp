#include <surfe_api.h>

GRBF_Modelling_Methods* Surfe_API::get_method(const model_parameters& params)
{
	if (params.model_type == Parameter_Types::Single_surface)
		return new Single_Surface(params);
	else if (params.model_type == Parameter_Types::Lajaunie_approach)
		return new Lajaunie_Approach(params);
	else if (params.model_type == Parameter_Types::Stratigraphic_horizons)
		return new Stratigraphic_Surfaces(params);
	else
		return new Continuous_Property(params);
}

void Surfe_API::AddInterfaceConstraint(const Interface& pt)
{
	method_->constraints.itrface.push_back(pt);
}

void Surfe_API::AddPlanarConstraint(const Planar& planar_pt)
{
	method_->constraints.planar.push_back(planar_pt);
}

void Surfe_API::AddTangentConstraint(const Tangent& tangent_pt)
{
	method_->constraints.tangent.push_back(tangent_pt);
}

void Surfe_API::AddInequalityConstraint(const Inequality& inequality_pt)
{
	method_->constraints.inequality.push_back(inequality_pt);
}

void Surfe_API::ComputeInterpolant()
{
	method_->remove_collocated_constraints();

	try
	{
		method_->process_input_data();
	}
	catch (std::exception& e)
	{
		SurfeExceptions exceptions(e);
		throw exceptions.what();
	}

	method_->get_method_parameters();

	try
	{
		method_->setup_basis_functions();
	}
	catch (std::exception& e)
	{
		SurfeExceptions exceptions(e);
		throw exceptions.what();
	}

	try
	{
		method_->setup_system_solver();
	}
	catch (std::exception& e)
	{
		SurfeExceptions exceptions(e);
		throw exceptions.what();
	}
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
