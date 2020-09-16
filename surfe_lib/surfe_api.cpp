#include <surfe_api.h>

#include <algorithm>
#include <time.h>
#include <vector>

void Surfe_API::progress(const float& precent_value)
{
	std::cout << int(precent_value * 100.0) << " %\r";
	std::cout.flush();
}

GRBF_Modelling_Methods* Surfe_API::get_method_from_parameters(const Parameters& params)
{
	if (params.model_type == Parameter_Types::Single_surface)
		return new Single_Surface(params);
	else if (params.model_type == Parameter_Types::Lajaunie_approach)
		return new Lajaunie_Approach(params);
	else if (params.model_type == Parameter_Types::Stratigraphic_horizons)
		return new Stratigraphic_Surfaces(params);
	else if (params.model_type == Parameter_Types::Vector_field)
		return new Vector_Field(params);
	else if (params.model_type == Parameter_Types::Continuous_property)
		return new Continuous_Property(params);
	else
		std::throw_with_nested(GRBF_Exceptions::unknown_modelling_mode);
}

SpatialParameters Surfe_API::GetDataBoundsAndResolution()
{
	// collect all constraints
	std::vector<Point> points;
	for (const auto &constraint : method_->constraints.inequality)
		points.emplace_back(Point(constraint.x(),constraint.y(),constraint.z()));
	for (const auto &constraint : method_->constraints.itrface)
		points.emplace_back(Point(constraint.x(), constraint.y(), constraint.z()));
	for (const auto &constraint : method_->constraints.planar)
		points.emplace_back(Point(constraint.x(), constraint.y(), constraint.z()));
	for (const auto &constraint : method_->constraints.tangent)
		points.emplace_back(Point(constraint.x(), constraint.y(), constraint.z()));

	// compute bounds and resolutions
	double xmin, xmax,ymin, ymax, zmin, zmax, resolution;
	if (!spatial_metrics(points, resolution, xmin, xmax, ymin, ymax, zmin, zmax))
		throw GRBF_Exceptions::problem_computing_spatial_parameters;
	SpatialParameters spatial;
	spatial.resolution = resolution;
	spatial.xmin = xmin;
	spatial.xmax = xmax;
	spatial.ymin = ymin;
	spatial.ymax = ymax;
	spatial.zmin = zmin;
	spatial.zmax = zmax;
	return spatial;
}

MatrixXd Surfe_API::GetInterfaceReferencePoints()
{
	std::vector<Interface> ref_pts = method_->interface_test_points;
	int n = ref_pts.size();
	MatrixXd ptarray(n, 3);
	for (int j = 0; j < n; j++) {
		ptarray(j, 0) = ref_pts[j].x();
		ptarray(j, 1) = ref_pts[j].y();
		ptarray(j, 2) = ref_pts[j].z();
	}
	return ptarray;
}

MatrixXd Surfe_API::GetInterfaceConstraints()
{
	std::vector<Interface> interface = method_->constraints.itrface;
	int n = interface.size();
	MatrixXd interface_constraints(n, 4);
	for (int j = 0; j < n; j++) {
		interface_constraints(j, 0) = interface[j].x();
		interface_constraints(j, 1) = interface[j].y();
		interface_constraints(j, 2) = interface[j].z();
		interface_constraints(j, 3) = interface[j].level();
	}
	//std::cout << " Interface Constraints Array:\n" << interface_constraints << std::endl;
	return interface_constraints;
}


void Surfe_API::SetInterfaceConstraints(const MatrixXd &interface_constraints)
{
	// does the interpolant already have interface constraints ?
	// if so, erase
	if (!method_->constraints.itrface.empty())
		method_->constraints.itrface.clear();

	int n = interface_constraints.rows();
	if (n != 0 && interface_constraints.cols() == 4)
	{
		for (int j = 0; j < n; j++)
			AddInterfaceConstraint(
				interface_constraints(j, 0), // x
				interface_constraints(j, 1), // y
				interface_constraints(j, 2), // z
				interface_constraints(j, 3)); // level
	}
	else
		throw GRBF_Exceptions::array_has_incorrect_dimensions;
}

MatrixXd Surfe_API::GetPlanarConstraints()
{
	std::vector<Planar> planar = method_->constraints.planar;
	int n = planar.size();
	MatrixXd planar_constraints(n, 6);
	for (int j = 0; j < n; j++) {
		planar_constraints(j, 0) = planar[j].x();
		planar_constraints(j, 1) = planar[j].y();
		planar_constraints(j, 2) = planar[j].z();
		planar_constraints(j, 3) = planar[j].nx();
		planar_constraints(j, 4) = planar[j].ny();
		planar_constraints(j, 5) = planar[j].nz();
	}
	return planar_constraints;
}

void Surfe_API::SetPlanarConstraints(const MatrixXd &planar_constraints)
{
	// does the interpolant already have planar constraints ?
	// if so, erase
	if (!method_->constraints.planar.empty())
		method_->constraints.planar.clear();

	int n = planar_constraints.rows();
	if (n != 0 && planar_constraints.cols() == 6)
	{
		for (int j = 0; j < n; j++)
			AddPlanarConstraintwNormal(
				planar_constraints(j, 0), // x
				planar_constraints(j, 1), // y
				planar_constraints(j, 2), // z
				planar_constraints(j, 3), // nx
				planar_constraints(j, 4), // ny
				planar_constraints(j, 5)); // nz
	}
	else
		throw GRBF_Exceptions::array_has_incorrect_dimensions;
}

MatrixXd Surfe_API::GetTangentConstraints()
{
	std::vector<Tangent> tangent = method_->constraints.tangent;
	int n = tangent.size();
	MatrixXd tangent_constraints(n, 6);
	for (int j = 0; j < n; j++) {
		tangent_constraints(j, 0) = tangent[j].x();
		tangent_constraints(j, 1) = tangent[j].y();
		tangent_constraints(j, 2) = tangent[j].z();
		tangent_constraints(j, 3) = tangent[j].tx();
		tangent_constraints(j, 4) = tangent[j].ty();
		tangent_constraints(j, 5) = tangent[j].tz();
	}
	return tangent_constraints;
}

void Surfe_API::SetTangentConstraints(const MatrixXd &tangent_constraints)
{
	// does the interpolant already have tangent constraints ?
	// if so, erase
	if (!method_->constraints.tangent.empty())
		method_->constraints.tangent.clear();

	int n = tangent_constraints.rows();
	if (n != 0 && tangent_constraints.cols() == 6)
	{
		for (int j = 0; j < n; j++)
			AddPlanarConstraintwNormal(
				tangent_constraints(j, 0), // x
				tangent_constraints(j, 1), // y
				tangent_constraints(j, 2), // z
				tangent_constraints(j, 3), // vx
				tangent_constraints(j, 4), // vy
				tangent_constraints(j, 5)); // vz
	}
	else
		throw GRBF_Exceptions::array_has_incorrect_dimensions;
}

MatrixXd Surfe_API::GetInequalityConstraints()
{
	std::vector<Inequality> ie = method_->constraints.inequality;
	int n = ie.size();
	MatrixXd inequality_constraints(n, 4);
	for (int j = 0; j < n; j++) {
		inequality_constraints(j, 0) = ie[j].x();
		inequality_constraints(j, 1) = ie[j].y();
		inequality_constraints(j, 2) = ie[j].z();
		inequality_constraints(j, 3) = ie[j].level();
	}
	return inequality_constraints;
}

void Surfe_API::SetInequalityConstraints(const MatrixXd &inequality_constraints)
{
	// does the interpolant already have inequality constraints ?
	// if so, erase
	if (!method_->constraints.inequality.empty())
		method_->constraints.inequality.clear();

	int n = inequality_constraints.rows();
	if (n != 0 && inequality_constraints.cols() == 4)
	{
		for (int j = 0; j < n; j++)
			AddInequalityConstraint(
				inequality_constraints(j, 0), // x
				inequality_constraints(j, 1), // y
				inequality_constraints(j, 2), // z
				inequality_constraints(j, 3)); // level
	}
	else
		throw GRBF_Exceptions::array_has_incorrect_dimensions;
}

int Surfe_API::GetNumberOfInterfaces()
{
	return (int)method_->interface_test_points.size();
}

Surfe_API::Surfe_API(const Parameters& params)
{
	method_ = nullptr;
	have_interpolant_ = false;
	parameters_changed_ = true;
	constraints_changed_ = false;

	method_ = get_method_from_parameters(params);
}

Surfe_API::Surfe_API(const int &modelling_method)
{
	if (modelling_method == 1) // Parameter_Types::Single_surface;
		method_ = new Single_Surface();
	else if (modelling_method == 2) // Parameter_Types::Lajaunie_approach
		method_ = new Lajaunie_Approach();
	else if (modelling_method == 3) // Parameter_Types::Vector_field
		method_ = new Vector_Field();
	else if (modelling_method == 4) // Parameter_Types::Stratigraphic_horizons
		method_ = new Stratigraphic_Surfaces();
	else if (modelling_method == 5) // Parameter_Types::Continuous_property
		method_ = new Continuous_Property();
	else
		throw GRBF_Exceptions::unknown_modelling_mode;

	parameters_changed_ = true;
}

void Surfe_API::AddInterfaceConstraint(const double &x, const double &y, const double &z, const double &level)
{
	Interface interface_constraint(x, y, z, level);
	method_->constraints.itrface.push_back(interface_constraint);

	method_->parameters.use_interface = true;
	constraints_changed_ = true;
}

void Surfe_API::AddPlanarConstraintwNormal(const double &x, const double &y, const double &z, const double &nx, const double &ny, const double &nz)
{
	Planar planar_constraint(x, y, z, nx, ny, nz);
	method_->constraints.planar.push_back(planar_constraint);

	method_->parameters.use_planar = true;
	constraints_changed_ = true;
}

void Surfe_API::AddPlanarConstraintwStrikeDipPolarity(const double &x, const double &y, const double &z, const double &strike, const double &dip, const int &polarity)
{
	Planar planar_constraint(x, y, z, dip, strike, polarity);
	method_->constraints.planar.push_back(planar_constraint);

	method_->parameters.use_planar = true;
	constraints_changed_ = true;
}

void Surfe_API::AddPlanarConstraintwAzimuthDipPolarity(const double &x, const double &y, const double &z, const double &azimuth, const double &dip, const int &polarity)
{
	double strike = 0.0;
	// convert azimuth to strike
	if (azimuth >= 90.0)
		strike = azimuth - 90.0;
	else
		strike = azimuth + 270.0;
	Planar planar_constraint(x, y, z, dip, strike, polarity);
	method_->constraints.planar.push_back(planar_constraint);

	method_->parameters.use_planar = true;
	constraints_changed_ = true;
}

void Surfe_API::AddTangentConstraint(const double &x, const double &y, const double &z, const double &tx, const double &ty, const double &tz)
{
	Tangent tangent_constraint(x, y, z, tx, ty, tz);
	method_->constraints.tangent.push_back(tangent_constraint);

	method_->parameters.use_tangent = true;
	constraints_changed_ = true;
}

void Surfe_API::AddInequalityConstraint(const double &x, const double &y, const double &z, const double &level)
{
	Inequality inequality_constraint(x, y, z, level);
	method_->constraints.inequality.push_back(inequality_constraint);

	method_->parameters.use_inequality = true;
	constraints_changed_ = true;
}

void Surfe_API::ComputeInterpolant()
{
	method_->remove_collocated_constraints();

	try
	{
		method_->process_input_data();
	}
	catch (const std::exception& e)
	{
		SurfeExceptions exceptions(e);
		throw exceptions;
	}

	method_->get_method_parameters();

	try
	{
		method_->setup_basis_functions();
	}
	catch (const std::exception& e)
	{
		SurfeExceptions exceptions(e);
		throw exceptions;
	}

	try
	{
		method_->setup_system_solver();
	}
	catch (const std::exception& e)
	{
		SurfeExceptions exceptions(e);
		throw exceptions;
	}

	std::cout << "Interpolant has been computed" << std::endl;

	have_interpolant_ = true;
	constraints_changed_ = false;
	parameters_changed_ = false;
}

void Surfe_API::SetRegressionSmoothing(const bool &use_regression_smoothing, const double &amount /*= 0*/)
{
	method_->parameters.use_regression_smoothing = true;
	method_->parameters.smoothing_amount = amount;

	parameters_changed_ = true;
}

void Surfe_API::SetGreedyAlgorithm(const bool &use_greedy, const double &interface_uncertainty /*= 0*/, const double &angular_uncertainty /*= 0*/)
{
	method_->parameters.use_greedy = true;
	method_->parameters.interface_uncertainty = interface_uncertainty;
	method_->parameters.angular_uncertainty = angular_uncertainty;

	parameters_changed_ = true;
}

void Surfe_API::SetRBFKernel(const Parameter_Types::RBF &rbf)
{
	method_->parameters.basis_type = rbf;

	parameters_changed_ = true;
}

void Surfe_API::SetRBFKernel(const char *rbf_name)
{
	if (strcmp(rbf_name, "r3") == 0)
		method_->parameters.basis_type = Parameter_Types::RBF::Cubic;
	else if (strcmp(rbf_name, "WendlandC2") == 0)
		method_->parameters.basis_type = Parameter_Types::RBF::WendlandC2;
	else if (strcmp(rbf_name, "r") == 0)
		method_->parameters.basis_type = Parameter_Types::RBF::R;
	else if (strcmp(rbf_name, "Gaussian") == 0)
		method_->parameters.basis_type = Parameter_Types::RBF::Gaussian;
	else if (strcmp(rbf_name, "Multiquadratics") == 0)
		method_->parameters.basis_type = Parameter_Types::RBF::MQ;
	else if (strcmp(rbf_name, "Thin Plate Spline") == 0)
		method_->parameters.basis_type = Parameter_Types::RBF::TPS;
	else if (strcmp(rbf_name, "Inverse Multiquadratics") == 0)
		method_->parameters.basis_type = Parameter_Types::RBF::IMQ;
	else if (strcmp(rbf_name, "MaternC4") == 0)
		method_->parameters.basis_type = Parameter_Types::MaternC4;
	else
		throw GRBF_Exceptions::unknown_rbf;

	parameters_changed_ = true;
}

void Surfe_API::SetRBFShapeParameter(const double &shape_param)
{
	method_->parameters.shape_parameter = shape_param;

	parameters_changed_ = true;
}

void Surfe_API::SetPolynomialOrder(const int &poly_order)
{
	method_->parameters.polynomial_order = poly_order;

	parameters_changed_ = true;
}

void Surfe_API::SetGlobalAnisotropy(const bool &g_anisotropy)
{
	method_->parameters.model_global_anisotropy = g_anisotropy;

	parameters_changed_ = true;
}


void Surfe_API::SetRestrictedRange(const bool &use_restricted_range, const double &interface_uncertainty /*= 0*/, const double &angular_uncertainty /*= 0*/)
{
	method_->parameters.use_restricted_range = use_restricted_range;
	method_->parameters.interface_uncertainty = interface_uncertainty;
	method_->parameters.angular_uncertainty = angular_uncertainty;

	parameters_changed_ = true;
}

double Surfe_API::EvaluateInterpolantAtPoint(const double &x, const double &y, const double &z)
{
	if (have_interpolant_)
	{
		if (constraints_changed_ || parameters_changed_)
			throw GRBF_Exceptions::interpolant_needs_update;

		// convert x,y,z to Point
		Point pt(x, y, z);
		// evaluate scalar field at point
		method_->eval_scalar_interpolant_at_point(pt);
		return pt.scalar_field();
	}
	else
		throw GRBF_Exceptions::missing_interpolant;
}
VectorXd Surfe_API::EvaluateInterpolantAtPoints(const MatrixXd &locations)
{
	std::cout << " Evaluating interpolant at list of points..." << std::endl;
	if (have_interpolant_)
	{
		int n = locations.rows();
		int evaluations_completed = 0;
		int percent_increment = 0;
		VectorXd interpolant(n);
		if (n != 0 && locations.cols() == 3)
		{
#pragma omp parallel for
			for (int j = 0; j < n; j++)
			{
				// convert x,y,z to Point
				Point pt(locations(j,0),locations(j,1),locations(j,2));

				// evaluate scalar field at this point
				method_->eval_scalar_interpolant_at_point(pt);

				// set scalar field value for this point in vector
				interpolant(j) = pt.scalar_field();

				evaluations_completed++;
				float percent_completed = ((float)evaluations_completed / (float)n);
				int percent_integer = std::round(percent_completed * 100.0);

				if (percent_integer > percent_increment)
				{
					percent_increment = percent_integer;
					progress(percent_completed);
				}

			}

			progress(1.0);
			std::cout << "" << std::endl;
		
			return interpolant;
		}
		else
			throw GRBF_Exceptions::array_has_incorrect_dimensions;
		
		
	}
	else
		throw GRBF_Exceptions::missing_interpolant;

	
}
Vector3d Surfe_API::EvaluateVectorInterpolantAtPoint(const double &x, const double &y, const double &z)
{
	if (have_interpolant_)
	{
		if (constraints_changed_ || parameters_changed_)
			throw GRBF_Exceptions::interpolant_needs_update;

		// convert x,y,z to Point
		Point pt(x, y, z);
		// evaluate vector field at point
		double gradient[3];
		method_->eval_vector_interpolant_at_point(pt);
		gradient[0] = pt.nx_interp();
		gradient[1] = pt.ny_interp();
		gradient[2] = pt.nz_interp();
		Vector3d grad(gradient[0], gradient[1], gradient[2]);
		return grad;
	}
	else
		throw GRBF_Exceptions::missing_interpolant;
}

MatrixXd Surfe_API::EvaluateVectorInterpolantAtPoints(const MatrixXd &locations)
{
	std::cout << " Evaluating vector interpolant at list of points..." << std::endl;
	if (have_interpolant_)
	{
		int n = locations.rows();
		int evaluations_completed = 0;
		int percent_increment = 0;

		MatrixXd interpolant(n,3);
		if (n != 0 && locations.cols() == 3)
		{
			for (int j = 0; j < n; j++)
			{
				// convert x,y,z to Point
				Point pt(locations(j,0),locations(j,1),locations(j,2));

				// evaluate scalar field at this point
				method_->eval_vector_interpolant_at_point(pt);

				// set vector components field value for this point in vector

				interpolant(j,0) = pt.nx_interp();
				interpolant(j,1) = pt.ny_interp();
				interpolant(j,2) = pt.nz_interp();

				evaluations_completed++;
				float percent_completed = ((float)evaluations_completed / (float)n);
				int percent_integer = std::round(percent_completed * 100.0);

				if (percent_integer > percent_increment)
				{
					percent_increment = percent_integer;
					progress(percent_completed);
				}
			}

			progress(1.0);
			std::cout << "" << std::endl;
		
			return interpolant;
		}
		else
			throw GRBF_Exceptions::array_has_incorrect_dimensions;
		
		
	}
	else
		throw GRBF_Exceptions::missing_interpolant;

	
}