#include <surfe_api.h>

#include <read_input_files.h>

#include <vtkNew.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImagePlaneWidget.h>
#include <vtkCellPicker.h>

#include <algorithm>
#include <time.h>

GRBF_Modelling_Methods* Surfe_API::get_method(const Parameters& params)
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

void Surfe_API::build_constraints_from_input_files()
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;
	try
	{
		if (!input_.interface_file.empty()) {
			std::vector<Interface> interface_constraints;
			std::string extension = get_file_extension(input_.interface_file.c_str());
			if (extension == "csv")
			{
				CSVInterfaceConstraintFileReader reader = 
					CSVInterfaceConstraintFileReader::CreateUsingDefaultPropertyNames(input_.interface_file.c_str());
				interface_constraints = reader.GetConstraints();
			}
			else if (extension == "vtp" || extension == "vtk")
			{
				VTKInterfaceConstraintFileReader reader = 
					VTKInterfaceConstraintFileReader::CreateUsingDefaultPropertyNames(input_.interface_file.c_str());
				interface_constraints = reader.GetConstraints();
			}
			method_->ui_parameters.use_interface = true;
			method_->constraints.itrface = interface_constraints;
		}
		if (!input_.inequality_file.empty()) {
			std::vector<Inequality> inequality_constraints;
			std::string extension = get_file_extension(input_.inequality_file.c_str());
			if (extension == "csv")
			{
				CSVInequalityConstraintFileReader reader =
					CSVInequalityConstraintFileReader::CreateUsingDefaultPropertyNames(input_.inequality_file.c_str());
				inequality_constraints = reader.GetConstraints();
			}
			else if (extension == "vtp" || extension == "vtk")
			{
				VTKInequalityConstraintFileReader reader =
					VTKInequalityConstraintFileReader::CreateUsingDefaultPropertyNames(input_.inequality_file.c_str());
				inequality_constraints = reader.GetConstraints();
			}
			method_->ui_parameters.use_inequality = true;
			method_->constraints.inequality = inequality_constraints;
		}

		if (!input_.planar_file.empty()) {
			std::vector<Planar> planar_constraints;
			std::string extension = get_file_extension(input_.planar_file.c_str());
			if (extension == "csv")
			{
				CSVPlanarConstraintFileReader reader =
					CSVPlanarConstraintFileReader::CreateUsingDefaultPropertyNames(input_.planar_file.c_str());
				planar_constraints = reader.GetConstraints();
			}
			else if (extension == "vtp" || extension == "vtk")
			{
				VTKPlanarConstraintFileReader reader =
					VTKPlanarConstraintFileReader::CreateUsingDefaultPropertyNames(input_.planar_file.c_str());
				planar_constraints = reader.GetConstraints();
			}
			method_->ui_parameters.use_planar = true;
			method_->constraints.planar = planar_constraints;
		}
		if (!input_.tangent_file.empty()) {
			std::vector<Tangent> tangent_constraints;
			std::string extension = get_file_extension(input_.planar_file.c_str());
			if (extension == "csv")
			{
				CSVTangentConstraintFileReader reader =
					CSVTangentConstraintFileReader::CreateUsingDefaultPropertyNames(input_.tangent_file.c_str());
				tangent_constraints = reader.GetConstraints();
			}
			else if (extension == "vtp" || extension == "vtk")
			{
				VTKTangentConstraintFileReader reader =
					VTKTangentConstraintFileReader::CreateUsingDefaultPropertyNames(input_.tangent_file.c_str());
				tangent_constraints = reader.GetConstraints();
			}
			method_->ui_parameters.use_tangent = true;
			method_->constraints.tangent = tangent_constraints;
		}
	}
	catch (const std::exception&e)
	{
		std::rethrow_if_nested(e);
	}

	constraint_files_changed_ = false; // since they have been loaded
	constraints_changed_ = true;
}

SpatialParameters Surfe_API::GetDataBoundsAndResolution()
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

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

int Surfe_API::GetNumberOfInterfaces()
{
	return (int)method_->interface_test_points.size();
}

double Surfe_API::GetScalarFieldValueAtInterfaceX(const int &i)
{
	if (i <= method_->interface_test_points.size())
	{
		Interface interface_pt = method_->interface_test_points[i];
		// evaluate interpolant at this point
		Point point(interface_pt.x(), interface_pt.y(), interface_pt.z());
		method_->eval_scalar_interpolant_at_point(point);
		return point.scalar_field();
	}
	else
		return -99999;
}

Surfe_API::Surfe_API()
{
	method_ = nullptr;
	have_interpolant_ = false;
	have_method_ = false;
	parameters_changed_ = false;
	constraint_files_changed_ = true;
	constraints_changed_ = false;
}

Surfe_API::Surfe_API(const Parameters& params)
{
	input_.parameters = params;

	method_ = nullptr;
	have_interpolant_ = false;
	have_method_ = true;
	parameters_changed_ = true;
	constraint_files_changed_ = true;
	constraints_changed_ = false;

	method_ = get_method(input_.parameters);

	try
	{
		build_constraints_from_input_files();
	}
	catch (const std::exception&e)
	{
		std::cout << "Exception: " << e.what() << " occurred. " << std::endl;
		throw;
	}
	evaluation_completed_ = false;
}

void Surfe_API::GetParametersAndConstraints()
{
	input_ = InputImpl::GetDialogParameters();

	parameters_changed_ = true;
	constraint_files_changed_ = true;

	// if method exists erase it, user could have changed parameters. 
	// leaving existing method not valid. delete it to not have 
	// memory leak
	if (method_) {
		delete method_;
		method_ = nullptr;
	}

	method_ = get_method(input_.parameters);
	try
	{
		build_constraints_from_input_files();
	}
	catch (const std::exception&e)
	{
		SurfeExceptions exceptions(e);
		throw exceptions;
	}
}

void Surfe_API::LoadConstraintsFromFiles()
{
	// if method exists erase it, user could have changed parameters. 
	// leaving existing method not valid. delete it to not have 
	// memory leak
	if (method_) {
		delete method_;
		method_ = nullptr;
	}

	method_ = get_method(input_.parameters);

	try
	{
		build_constraints_from_input_files();
	}
	catch (const std::exception&e)
	{
 		SurfeExceptions exceptions(e);
		throw exceptions;
	}
}

void Surfe_API::AddInterfaceConstraint(const double &x, const double &y, const double &z, const double &level)
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	Interface interface_constraint(x, y, z, level);
	method_->constraints.itrface.push_back(interface_constraint);

	method_->ui_parameters.use_interface = true;
	constraints_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::AddPlanarConstraintwNormal(const double &x, const double &y, const double &z, const double &nx, const double &ny, const double &nz)
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	Planar planar_constraint(x, y, z, nx, ny, nz);
	method_->constraints.planar.push_back(planar_constraint);

	method_->ui_parameters.use_planar = true;
	constraints_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::AddPlanarConstraintwStrikeDipPolarity(const double &x, const double &y, const double &z, const double &strike, const double &dip, const double &polarity)
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	Planar planar_constraint(x, y, z, dip, strike, polarity);
	method_->constraints.planar.push_back(planar_constraint);

	method_->ui_parameters.use_planar = true;
	constraints_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::AddPlanarConstraintwAzimuthDipPolarity(const double &x, const double &y, const double &z, const double &azimuth, const double &dip, const double &polarity)
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	double strike = 0.0;
	// convert azimuth to strike
	if (azimuth >= 90.0)
		strike = azimuth - 90.0;
	else
		strike = azimuth + 270.0;
	Planar planar_constraint(x, y, z, dip, strike, polarity);
	method_->constraints.planar.push_back(planar_constraint);

	method_->ui_parameters.use_planar = true;
	constraints_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::AddTangentConstraint(const double &x, const double &y, const double &z, const double &tx, const double &ty, const double &tz)
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	Tangent tangent_constraint(x, y, z, tx, ty, tz);
	method_->constraints.tangent.push_back(tangent_constraint);

	method_->ui_parameters.use_tangent = true;
	constraints_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::AddInequalityConstraint(const double &x, const double &y, const double &z, const double &level)
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	Inequality inequality_constraint(x, y, z, level);
	method_->constraints.inequality.push_back(inequality_constraint);

	method_->ui_parameters.use_inequality = true;
	constraints_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::ComputeInterpolant()
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	if (constraint_files_changed_)
	{
		try
		{
			build_constraints_from_input_files();
		}
		catch (const std::exception&e)
		{
			SurfeExceptions exceptions(e);
			throw exceptions;
		}
	}

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
	evaluation_completed_ = false;

}

void Surfe_API::SetModellingMode(const int &mode)
{

	if (mode == 1)
		input_.parameters.model_type = Parameter_Types::Single_surface;
	else if (mode == 2) 
		input_.parameters.model_type = Parameter_Types::Lajaunie_approach;
	else if (mode == 3)
		input_.parameters.model_type = Parameter_Types::Vector_field;
	else if (mode == 4)
		input_.parameters.model_type = Parameter_Types::Stratigraphic_horizons;
	else if (mode == 5)
		input_.parameters.model_type = Parameter_Types::Continuous_property;
	else
		throw GRBF_Exceptions::unknown_modelling_mode;

	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetRegressionSmoothing(const bool &use_regression_smoothing, const double &amount /*= 0*/)
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	method_->ui_parameters.use_regression_smoothing = true;
	method_->ui_parameters.smoothing_amount = amount;

	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetGreedyAlgorithm(const bool &use_greedy, const double &interface_uncertainty /*= 0*/, const double &angular_uncertainty /*= 0*/)
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	method_->ui_parameters.use_greedy = true;
	method_->ui_parameters.interface_uncertainty = interface_uncertainty;
	method_->ui_parameters.angular_uncertainty = angular_uncertainty;

	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetRBFKernel(const Parameter_Types::RBF &rbf)
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	method_->ui_parameters.basis_type = rbf;

	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetRBFKernel(const char *rbf_name)
{
	if (strcmp(rbf_name, "r3") == 0)
		input_.parameters.basis_type = Parameter_Types::RBF::Cubic;
	else if (strcmp(rbf_name, "WendlandC2") == 0)
		input_.parameters.basis_type = Parameter_Types::RBF::WendlandC2;
	else if (strcmp(rbf_name, "r") == 0)
		input_.parameters.basis_type = Parameter_Types::RBF::R;
	else if (strcmp(rbf_name, "Gaussian") == 0)
		input_.parameters.basis_type = Parameter_Types::RBF::Gaussian;
	else if (strcmp(rbf_name, "Multiquadratics") == 0)
		input_.parameters.basis_type = Parameter_Types::RBF::MQ;
	else if (strcmp(rbf_name, "Thin Plate Spline") == 0)
		input_.parameters.basis_type = Parameter_Types::RBF::TPS;
	else if (strcmp(rbf_name, "Inverse Multiquadratics") == 0)
		input_.parameters.basis_type = Parameter_Types::RBF::IMQ;
	else if (strcmp(rbf_name, "MaternC4") == 0)
		input_.parameters.basis_type = Parameter_Types::MaternC4;
	else
		throw GRBF_Exceptions::unknown_rbf;

	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetRBFShapeParameter(const double &shape_param)
{
	input_.parameters.shape_parameter = shape_param;

	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetPolynomialOrder(const int &poly_order)
{
	input_.parameters.polynomial_order = poly_order;

	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetGlobalAnisotropy(const bool &g_anisotropy)
{
	input_.parameters.model_global_anisotropy = g_anisotropy;

	parameters_changed_ = true;
	evaluation_completed_ = false;
}


void Surfe_API::SetRestrictedRange(const bool &use_restricted_range, const double &interface_uncertainty /*= 0*/, const double &angular_uncertainty /*= 0*/)
{
	input_.parameters.use_regression_smoothing = use_restricted_range;
	input_.parameters.interface_uncertainty = interface_uncertainty;
	input_.parameters.angular_uncertainty = angular_uncertainty;

	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetInterfaceUncertainty(const double &interface_uncertainty)
{
	input_.parameters.interface_uncertainty = interface_uncertainty;

	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetAngularUncertainty(const double &angular_uncertainty)
{
	input_.parameters.angular_uncertainty = angular_uncertainty;

	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetInterfaceDataFile(const char *interfaceFile)
{
	input_.interface_file = interfaceFile;

	constraint_files_changed_ = true;
	constraints_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetPlanarDataFile(const char *planarFile)
{
	input_.planar_file = planarFile;

	constraint_files_changed_ = true;
	constraints_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetTangentDataFile(const char *tangentFile)
{
	input_.tangent_file = tangentFile;

	constraint_files_changed_ = true;
	constraints_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetInequalityDataFile(const char *inequalityFile)
{
	input_.inequality_file = inequalityFile;

	constraint_files_changed_ = true;
	constraints_changed_ = true;
	evaluation_completed_ = false;
}

double Surfe_API::EvaluateInterpolantAtPoint(const double &x, const double &y, const double &z)
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

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

double *Surfe_API::EvaluateVectorInterpolantAtPoint(const double &x, const double &y, const double &z)
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	if (have_interpolant_)
	{
		if (constraints_changed_ || parameters_changed_)
			throw GRBF_Exceptions::interpolant_needs_update;

		double *gradient = new double[3];
		// convert x,y,z to Point
		Point pt(x, y, z);
		// evaluate vector field at point
		method_->eval_vector_interpolant_at_point(pt);
		gradient[0] = pt.nx_interp();
		gradient[1] = pt.ny_interp();
		gradient[2] = pt.nz_interp();
		return gradient;
	}
	else
		throw GRBF_Exceptions::missing_interpolant;
}