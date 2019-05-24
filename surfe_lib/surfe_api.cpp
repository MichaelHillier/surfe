#include <surfe_api.h>

#include <read_input_files.h>

#include <vtkNew.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLImageDataWriter.h>

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
	else
		return new Continuous_Property(params);
}

void Surfe_API::build_constraints_from_input_files()
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;
	try
	{
		if (!input_.interface_file.empty()) {
			std::vector<Interface> interface_constraints;
			interface_constraints = build_interface_constraints(input_.interface_file.c_str());
			method_->ui_parameters.use_interface = true;
			method_->constraints.itrface = interface_constraints;
		}
		if (!input_.inequality_file.empty()) {
			std::vector<Inequality> inequality_constraints;
			inequality_constraints = build_inequality_constraints(input_.inequality_file.c_str());
			method_->ui_parameters.use_inequality = true;
			method_->constraints.inequality = inequality_constraints;
		}

		if (!input_.planar_file.empty()) {
			std::vector<Planar> planar_constraints;
			planar_constraints = build_planar_constraints(input_.planar_file.c_str());
			method_->ui_parameters.use_planar = true;
			method_->constraints.planar = planar_constraints;
		}
		if (!input_.tangent_file.empty()) {
			std::vector<Tangent> tangent_constraints;
			tangent_constraints = build_tangent_constraints(input_.tangent_file.c_str());
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

void Surfe_API::progress(const float &progress_value)
{
	int barWidth = 70;

	std::cout << "[";
	int pos = barWidth * progress_value;
	for (int i = 0; i < barWidth; ++i) {
		if (i < pos) std::cout << "=";
		else if (i == pos) std::cout << ">";
		else std::cout << " ";
	}
	std::cout << "] " << int(progress_value * 100.0) << " %\r";
	std::cout.flush();
}

Surfe_API::Surfe_API()
{
	method_ = nullptr;
	grid_ = nullptr;
	have_interpolant_ = false;
	have_method_ = false;
	evaluation_completed_ = false;
	parameters_changed_ = false;
	constraint_files_changed_ = true;
	constraints_changed_ = false;
	evaluation_completed_ = false;
}

Surfe_API::Surfe_API(const Parameters& params)
{
	input_.parameters = params;

	method_ = nullptr;
	grid_ = nullptr;
	have_interpolant_ = false;
	have_method_ = true;
	evaluation_completed_ = false;
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
	catch (std::exception& e)
	{
		SurfeExceptions exceptions(e);
		throw exceptions;
	}

	method_->get_method_parameters();

	try
	{
		method_->setup_basis_functions();
	}
	catch (std::exception& e)
	{
		SurfeExceptions exceptions(e);
		throw exceptions;
	}

	try
	{
		method_->setup_system_solver();
	}
	catch (std::exception& e)
	{
		SurfeExceptions exceptions(e);
		throw exceptions;
	}

	std::cout << "Interpolant has been computed" << std::endl;

	have_interpolant_ = true;
	constraints_changed_ = false;
	parameters_changed_ = false;

}

void Surfe_API::SetModellingMode(const int &mode)
{
	// if method exists erase it, user could have changed parameters. 
	// leaving existing method not valid. delete it to not have 
	// memory leak

	if (mode == 1) {
		if (method_) {
			delete method_;
			method_ = nullptr;
		}
		method_ = new Single_Surface;
		method_->ui_parameters.model_type = Parameter_Types::Single_surface;
	}
	else if (mode == 2) {
		if (method_) {
			delete method_;
			method_ = nullptr;
		}
		method_ = new Lajaunie_Approach;
		method_->ui_parameters.model_type = Parameter_Types::Lajaunie_approach;
	}
	else if (mode == 3) {
		if (method_) {
			delete method_;
			method_ = nullptr;
		}
		method_ = new Vector_Field;
		method_->ui_parameters.model_type = Parameter_Types::Vector_field;
	}
	else if (mode == 4) {
		if (method_) {
			delete method_;
			method_ = nullptr;
		}
		method_ = new Stratigraphic_Surfaces;
		method_->ui_parameters.model_type = Parameter_Types::Stratigraphic_horizons;
	}
	else if (mode == 5) {
		if (method_) {
			delete method_;
			method_ = nullptr;
		}
		method_ = new Continuous_Property;
		method_->ui_parameters.model_type = Parameter_Types::Continuous_property;
	}
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
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	if (rbf_name == "r3")
		method_->ui_parameters.basis_type = Parameter_Types::RBF::Cubic;
	else if (rbf_name == "WendlandC2")
		method_->ui_parameters.basis_type = Parameter_Types::RBF::WendlandC2;
	else if (rbf_name == "r")
		method_->ui_parameters.basis_type = Parameter_Types::RBF::R;
	else if (rbf_name == "Gaussian")
		method_->ui_parameters.basis_type = Parameter_Types::RBF::Gaussian;
	else if (rbf_name == "Multiquadratics")
		method_->ui_parameters.basis_type = Parameter_Types::RBF::MQ;
	else if (rbf_name == "Thin Plate Spline")
		method_->ui_parameters.basis_type = Parameter_Types::RBF::TPS;
	else if (rbf_name = "Inverse Multiquadratics")
		method_->ui_parameters.basis_type = Parameter_Types::RBF::IMQ; // Inverse Multiquadratics
	else
		throw GRBF_Exceptions::unknown_rbf;

	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetRBFShapeParameter(const double &shape_param)
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	method_->ui_parameters.shape_parameter = shape_param;

	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetPolynomialOrder(const int &poly_order)
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	method_->ui_parameters.polynomial_order = poly_order;

	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetGlobalAnisotropy(const bool &g_anisotropy)
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	method_->ui_parameters.model_global_anisotropy = g_anisotropy;

	parameters_changed_ = true;
	evaluation_completed_ = false;
}


void Surfe_API::SetRestrictedRange(const bool &use_restricted_range, const double &interface_uncertainty /*= 0*/, const double &angular_uncertainty /*= 0*/)
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	method_->ui_parameters.use_regression_smoothing = use_restricted_range;
	method_->ui_parameters.interface_uncertainty = interface_uncertainty;
	method_->ui_parameters.angular_uncertainty = angular_uncertainty;

	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetInterfaceUncertainty(const double &interface_uncertainty)
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	method_->ui_parameters.interface_uncertainty = interface_uncertainty;

	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetAngularUncertainty(const double &angular_uncertainty)
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	method_->ui_parameters.angular_uncertainty = angular_uncertainty;

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

void Surfe_API::BuildRegularGrid(const double &zmin, const double &zmax, const double &resolution, const double &xy_percent_padding /*= 0*/)
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	std::vector<Point> agg_pts = convert_constraints_to_points(method_->constraints);

	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	for (const auto &point : agg_pts)
		pts->InsertNextPoint(point.x(), point.y(), point.z());

	double bounds[6];
	pts->GetBounds(bounds);

	int nx = (bounds[1] - bounds[0]) / resolution;
	int ny = (bounds[3] - bounds[2]) / resolution;
	int nz = (zmax - zmin) / resolution;

	double origin[3] = { bounds[0],bounds[2], zmin };

	vtkSmartPointer<vtkImageData> constructed_grid = vtkSmartPointer<vtkImageData>::New();
	constructed_grid->SetDimensions(nx + 1, ny + 1, nz + 1);
	constructed_grid->SetOrigin(origin);
	constructed_grid->SetSpacing(resolution, resolution, resolution);

	grid_ = constructed_grid;
}

void Surfe_API::BuildRegularGrid(
	const double &xmin, const double &xmax,
	const double &ymin, const double &ymax,
	const double &zmin, const double &zmax,
	const double &resolution
)
{
	int nx = (xmax - xmin) / resolution;
	int ny = (ymax - ymin) / resolution;
	int nz = (zmax - zmin) / resolution;

	double origin[3] = { xmin, ymin, zmin };

	vtkSmartPointer<vtkImageData> constructed_grid = vtkSmartPointer<vtkImageData>::New();
	constructed_grid->SetDimensions(nx + 1, ny + 1, nz + 1);
	constructed_grid->SetOrigin(origin);
	constructed_grid->SetSpacing(resolution, resolution, resolution);

	grid_ = constructed_grid;
}

vtkSmartPointer<vtkImageData> Surfe_API::GetEvaluatedGrid()
{
	if (!grid_)
		throw GRBF_Exceptions::no_sgrid_exists;

	if (!have_interpolant_ || parameters_changed_ || constraints_changed_)
	{
		try
		{
			ComputeInterpolant();
		}
		catch (std::exception& e)
		{
			std::cout << "Exception: " << e.what() << " occurred. " << std::endl;
			throw;
		}
	}

	vtkSmartPointer<vtkDoubleArray> sfield = vtkSmartPointer<vtkDoubleArray>::New();
	sfield->SetName("Scalar Field");
	sfield->SetNumberOfComponents(1);
	sfield->SetNumberOfTuples(grid_->GetNumberOfPoints());

	clock_t this_time = clock();
	clock_t last_time = this_time;
	double time_counter = 0.0;

	int N = grid_->GetNumberOfPoints();
	int evaluations_completed = 0;
	std::cout << "Evaluating interpolant in grid: " << std::endl;
	#pragma omp parallel for
	for (int j = 0; j < N; j++) {
		double point[3];
		grid_->GetPoint(j, point);
		Point pt(point[0], point[1], point[2]);
		// evaluate scalar field at point
		method_->eval_scalar_interpolant_at_point(pt);
		double scalar_field = pt.scalar_field();
		sfield->SetTuple1(j, scalar_field);
		// Print progress every 1s
		evaluations_completed++;
		this_time = clock();
		time_counter += (double)(this_time - last_time);
		last_time = this_time;
		if (time_counter > (double)(1 * CLOCKS_PER_SEC))
		{
			time_counter -= (double)(1 * CLOCKS_PER_SEC);
			float percent_completed = ((float)evaluations_completed / (float)N);
			progress(percent_completed);
		}		
	}
	progress(1);
	std::cout<<std::endl;
	grid_->GetPointData()->SetScalars(sfield);

	evaluation_completed_ = true;

	cout << "Finished evaluating interpolant in grid" << endl;

	return grid_;
}

vtkSmartPointer<vtkPolyData> Surfe_API::GetIsoSurfaces()
{
	if (!grid_)
		throw GRBF_Exceptions::no_sgrid_exists;

	if (!evaluation_completed_ || parameters_changed_ || constraints_changed_) {
		try
		{
			GetEvaluatedGrid();
		}
		catch (const std::exception&)
		{
			throw;
		}
	}

	vtkNew<vtkMarchingCubes> mcube;
	mcube->SetInputData(grid_);
	mcube->ComputeScalarsOn();
	for (int j = 0; j < method_->interface_test_points.size(); j++) {
		Interface interface_pt = method_->interface_test_points[j];
		// evaluate interpolant at this point
		Point point(interface_pt.x(), interface_pt.y(), interface_pt.z());
		method_->eval_scalar_interpolant_at_point(point);
		double iso_value = point.scalar_field();
		mcube->SetValue(j, iso_value);
	}
	mcube->Update();

	vtkSmartPointer<vtkPolyData> iso_surfaces = vtkSmartPointer<vtkPolyData>::New();
	iso_surfaces = mcube->GetOutput();

	cout << "Finished marching cubes" << endl;

	return iso_surfaces;
}

vtkSmartPointer<vtkPolyData> Surfe_API::GetInterfaceConstraints()
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	if (!method_->constraints.itrface.empty())
	{
		// interface
		vtkSmartPointer<vtkPolyData> interface_constraints = vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPoints> interface_points = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkDoubleArray> interface_scalar_field = vtkSmartPointer<vtkDoubleArray>::New();
		interface_scalar_field->SetName("level");
		interface_scalar_field->SetNumberOfComponents(1);
		interface_scalar_field->SetNumberOfTuples(grid_->GetNumberOfPoints());
		for (int j = 0; j < (int)method_->constraints.itrface.size(); j++) {
			Interface *interface_pt = &method_->constraints.itrface[j];
			interface_points->InsertNextPoint(interface_pt->x(), interface_pt->y(), interface_pt->z());
			interface_scalar_field->SetTuple1(j, interface_pt->scalar_field());
			// 			vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
			// 			vertex->GetPointIds()->SetId(0, j);
			// 			interface_constraints->Allocate(1, 1);
			// 			interface_constraints->InsertNextCell(vertex->GetCellType(), vertex->GetPointIds());
		}
		interface_constraints->SetPoints(interface_points);
		interface_constraints->GetPointData()->AddArray(interface_scalar_field);
		return interface_constraints;
	}
	else
		return nullptr;
}

vtkSmartPointer<vtkPolyData> Surfe_API::GetPlanarConstraints()
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	if (!method_->constraints.planar.empty())
	{
		// planar
		vtkSmartPointer<vtkPolyData> planar_constraints = vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPoints> planar_points = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkDoubleArray> planar_gradient_field = vtkSmartPointer<vtkDoubleArray>::New();
		int n_tuples = (int)method_->constraints.planar.size();
		planar_gradient_field->SetName("normal");
		planar_gradient_field->SetNumberOfTuples(n_tuples);
		planar_gradient_field->SetNumberOfComponents(3);
		planar_gradient_field->SetComponentName(0, "nx");
		planar_gradient_field->SetComponentName(1, "ny");
		planar_gradient_field->SetComponentName(2, "nz");
		// initalization for vector data
		for (int k = 0; k < n_tuples; k++) {
			for (int l = 0; l < 3; l++)
				planar_gradient_field->InsertComponent(k, l, 0.0);
		}
		for (int j = 0; j < (int)method_->constraints.planar.size(); j++) {
			Planar *planar_pt = &method_->constraints.planar[j];
			planar_points->InsertNextPoint(planar_pt->x(), planar_pt->y(), planar_pt->z());
			double gradient[3] = { planar_pt->nx(),planar_pt->ny(),planar_pt->nz() };
			planar_gradient_field->SetTuple(j, gradient);
			// 			vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
			// 			vertex->GetPointIds()->SetId(0, j);
			// 			planar_constraints->Allocate(1, 1);
			// 			planar_constraints->InsertNextCell(vertex->GetCellType(), vertex->GetPointIds());
		}
		planar_constraints->SetPoints(planar_points);
		planar_constraints->GetPointData()->AddArray(planar_gradient_field);
		return planar_constraints;
	}
	else
		return nullptr;
}

vtkSmartPointer<vtkPolyData> Surfe_API::GetTangentConstraints()
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	if (!method_->constraints.tangent.empty())
	{
		// tangent
		vtkSmartPointer<vtkPolyData> tangent_constraints = vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPoints> tangent_points = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkDoubleArray> tangent_vector = vtkSmartPointer<vtkDoubleArray>::New();
		tangent_vector->SetName("tangent");
		tangent_vector->SetNumberOfComponents(3);
		tangent_vector->SetComponentName(0, "tx");
		tangent_vector->SetComponentName(1, "ty");
		tangent_vector->SetComponentName(2, "tz");
		for (int j = 0; j < (int)method_->constraints.tangent.size(); j++) {
			Tangent*tangent_pt = &method_->constraints.tangent[j];
			tangent_points->InsertNextPoint(tangent_pt->x(), tangent_pt->y(), tangent_pt->z());
			double tangent[3] = { tangent_pt->tx(),tangent_pt->ty(),tangent_pt->tz() };
			tangent_vector->SetTuple(j, tangent);
			// 			vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
			// 			vertex->GetPointIds()->SetId(0, j);
			// 			tangent_constraints->Allocate(1, 1);
			// 			tangent_constraints->InsertNextCell(vertex->GetCellType(), vertex->GetPointIds());
		}
		tangent_constraints->SetPoints(tangent_points);
		tangent_constraints->GetPointData()->AddArray(tangent_vector);
		return tangent_constraints;
	}
	else
		return nullptr;
}

vtkSmartPointer<vtkPolyData> Surfe_API::GetInequalityConstraints()
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;

	if (!method_->constraints.inequality.empty())
	{
		// inequality
		vtkSmartPointer<vtkPolyData> inequality_constraints = vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPoints> inequality_points = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkDoubleArray> inequality_scalar_field = vtkSmartPointer<vtkDoubleArray>::New();
		inequality_scalar_field->SetName("level");
		inequality_scalar_field->SetNumberOfComponents(1);
		inequality_scalar_field->SetNumberOfTuples(grid_->GetNumberOfPoints());
		for (int j = 0; j < (int)method_->constraints.inequality.size(); j++) {
			Inequality *inequality_pt = &method_->constraints.inequality[j];
			inequality_points->InsertNextPoint(inequality_pt->x(), inequality_pt->y(), inequality_pt->z());
			inequality_scalar_field->SetTuple1(j, inequality_pt->scalar_field());
			// 			vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
			// 			vertex->GetPointIds()->SetId(0, j);
			// 			inequality_constraints->Allocate(1, 1);
			// 			inequality_constraints->InsertNextCell(vertex->GetCellType(), vertex->GetPointIds());
		}
		inequality_constraints->SetPoints(inequality_points);
		inequality_constraints->GetPointData()->AddArray(inequality_scalar_field);
		return inequality_constraints;
	}
	else
		return nullptr;
}

const char *Surfe_API::GetEvaluatedVTKGridString()
{
	if (!grid_)
		throw GRBF_Exceptions::no_sgrid_exists;

	if (!evaluation_completed_ || parameters_changed_ || constraints_changed_) {
		try
		{
			GetEvaluatedGrid();
		}
		catch (const std::exception&)
		{
			throw;
		}
	}

	vtkNew<vtkXMLImageDataWriter> writer;
	writer->SetInputData(grid_);
	writer->WriteToOutputStringOn();
	writer->Write();

	vtk_grid_string_ = writer->GetOutputString();

	return vtk_grid_string_.c_str();
}

const char * Surfe_API::GetVTKIsosurfacesString()
{
	try
	{
		vtkSmartPointer<vtkPolyData> isosurfaces = GetIsoSurfaces();

		if (isosurfaces)
		{
			vtkNew<vtkXMLPolyDataWriter> writer;
			writer->SetInputData(isosurfaces);
			writer->WriteToOutputStringOn();
			writer->Write();

			vtk_isosurfaces_string_ = writer->GetOutputString();
			return vtk_isosurfaces_string_.c_str();
		}
	}
	catch (std::exception& e)
	{
		std::cout << "Exception: " << e.what() << " occurred. " << std::endl;
		throw;
	}
}

const char * Surfe_API::GetVTKInterfaceConstraintsString()
{
	vtkSmartPointer<vtkPolyData> poly = GetInterfaceConstraints();

	if (poly) {
		vtkNew<vtkXMLPolyDataWriter> writer;
		writer->SetInputData(poly);
		writer->WriteToOutputStringOn();
		writer->Write();

		vtk_interface_string_ = writer->GetOutputString();
		return vtk_interface_string_.c_str();
	}
	else
		return nullptr;
}

const char * Surfe_API::GetVTKPlanarConstraintsString()
{
	vtkSmartPointer<vtkPolyData> poly = GetPlanarConstraints();

	if (poly) {
		vtkNew<vtkXMLPolyDataWriter> writer;
		writer->SetInputData(poly);
		writer->WriteToOutputStringOn();
		writer->Write();

		vtk_planar_string_ = writer->GetOutputString();
		return vtk_planar_string_.c_str();
	}
	else
		return nullptr;
}

const char * Surfe_API::GetVTKTangentConstraintsString()
{
	vtkSmartPointer<vtkPolyData> poly = GetTangentConstraints();

	if (poly) {
		vtkNew<vtkXMLPolyDataWriter> writer;
		writer->SetInputData(poly);
		writer->WriteToOutputStringOn();
		writer->Write();

		vtk_tangent_string_ = writer->GetOutputString();
		return vtk_tangent_string_.c_str();
	}
	else
		return nullptr;
}

const char * Surfe_API::GetVTKInequalityConstraintString()
{
	vtkSmartPointer<vtkPolyData> poly = GetInequalityConstraints();

	if (poly) {
		vtkNew<vtkXMLPolyDataWriter> writer;
		writer->SetInputData(poly);
		writer->WriteToOutputStringOn();
		writer->Write();

		vtk_inequality_string_ = writer->GetOutputString();
		return vtk_inequality_string_.c_str();
	}
	else
		return nullptr;
}

void Surfe_API::WriteVTKInterfaceConstraints(const char *filename)
{
	vtkSmartPointer<vtkPolyData> poly = GetInterfaceConstraints();

	if (poly) {
		vtkNew<vtkXMLPolyDataWriter> writer;
		writer->SetInputData(poly);
		writer->SetFileName(filename);
		writer->Write();
	}
}

void Surfe_API::WriteVTKPlanarConstraints(const char *filename)
{
	vtkSmartPointer<vtkPolyData> poly = GetPlanarConstraints();

	if (poly) {
		vtkNew<vtkXMLPolyDataWriter> writer;
		writer->SetInputData(poly);
		writer->SetFileName(filename);
		writer->Write();
	}
}

void Surfe_API::WriteVTKTangentConstraints(const char *filename)
{
	vtkSmartPointer<vtkPolyData> poly = GetTangentConstraints();

	if (poly) {
		vtkNew<vtkXMLPolyDataWriter> writer;
		writer->SetInputData(poly);
		writer->SetFileName(filename);
		writer->Write();
	}
}

void Surfe_API::WriteVTKInequalityConstraints(const char *filename)
{
	vtkSmartPointer<vtkPolyData> poly = GetInequalityConstraints();

	if (poly) {
		vtkNew<vtkXMLPolyDataWriter> writer;
		writer->SetInputData(poly);
		writer->SetFileName(filename);
		writer->Write();
	}
}

void Surfe_API::WriteVTKEvaluationGrid(const char *filename)
{
	if (evaluation_completed_ && grid_) {
		vtkNew<vtkXMLImageDataWriter> writer;
		writer->SetInputData(grid_);
		writer->SetFileName(filename);
		writer->SetDataModeToBinary();
		writer->Write();
	}
}

void Surfe_API::WriteVTKIsoSurfaces(const char *filename)
{
	vtkPolyData *isosurfaces = GetIsoSurfaces();

	if (isosurfaces) {
		vtkNew<vtkXMLPolyDataWriter> writer;
		writer->SetInputData(isosurfaces);
		writer->SetFileName(filename);
		writer->SetDataModeToBinary();
		writer->Write();
	}
}

void Surfe_API::VisualizeVTKData()
{
	// get grid
	vtkImageData *grid = GetEvaluatedGrid();
	double spacing[3];
	grid->GetSpacing(spacing);
	double min_scale = *std::max_element(spacing, spacing + 3);

	// get constraints
	vtkSmartPointer<vtkPolyData> interface = GetInterfaceConstraints();
	vtkSmartPointer<vtkPolyData> planar = GetPlanarConstraints();
	vtkSmartPointer<vtkPolyData> tangent = GetTangentConstraints();
	vtkSmartPointer<vtkPolyData> inequality = GetInequalityConstraints();

	// get iso surfaces
	vtkSmartPointer<vtkPolyData> isosurfaces = GetIsoSurfaces();

	vtkNew<vtkRenderer> ren;
	ren->SetBackground(0.1, 0.15, 0.3);
	vtkNew<vtkRenderWindow> renWin;
	renWin->SetSize(1000, 1000);
	renWin->AddRenderer(ren);

	vtkNew<vtkRenderWindowInteractor> iren;
	iren->SetRenderWindow(renWin);

	vtkNew<vtkLookupTable> lut;
	lut->SetNumberOfColors(256);
	lut->SetHueRange(0.0, 1.0);
	lut->SetRange(grid->GetScalarRange());
	lut->Build();

	// Grid Objects
	int dimensions[3];
	grid->GetDimensions(dimensions);
	// x slice
	vtkNew<vtkImageSliceMapper> grid_xslice_mapper;
	grid_xslice_mapper->SetInputData(grid);
	grid_xslice_mapper->SetOrientationToX();
	grid_xslice_mapper->SetSliceNumber(dimensions[0] - 1);

	vtkNew<vtkImageSlice> grid_xslice;
	grid_xslice->SetMapper(grid_xslice_mapper);
	grid_xslice->GetProperty()->SetLookupTable(lut);
	grid_xslice->GetProperty()->UseLookupTableScalarRangeOn();
	grid_xslice->GetProperty()->SetInterpolationTypeToLinear();
	grid_xslice->GetProperty()->SetOpacity(0.35);

	ren->AddViewProp(grid_xslice);

	// y slice
	vtkNew<vtkImageSliceMapper> grid_yslice_mapper;
	grid_yslice_mapper->SetInputData(grid);
	grid_yslice_mapper->SetOrientationToY();
	grid_yslice_mapper->SetSliceNumber(dimensions[1] - 1);

	vtkNew<vtkImageSlice> grid_yslice;
	grid_yslice->SetMapper(grid_yslice_mapper);
	grid_yslice->GetProperty()->SetLookupTable(lut);
	grid_yslice->GetProperty()->UseLookupTableScalarRangeOn();
	grid_yslice->GetProperty()->SetInterpolationTypeToLinear();
	grid_yslice->GetProperty()->SetOpacity(0.35);

	ren->AddViewProp(grid_yslice);

	// z slice
	vtkNew<vtkImageSliceMapper> grid_zslice_mapper;
	grid_zslice_mapper->SetInputData(grid);
	grid_zslice_mapper->SetOrientationToZ();
	grid_zslice_mapper->SetSliceNumber(0);

	vtkNew<vtkImageSlice> grid_zslice;
	grid_zslice->SetMapper(grid_zslice_mapper);
	grid_zslice->GetProperty()->SetLookupTable(lut);
	grid_zslice->GetProperty()->UseLookupTableScalarRangeOn();
	grid_zslice->GetProperty()->SetInterpolationTypeToLinear();
	grid_zslice->GetProperty()->SetOpacity(0.35);

	ren->AddViewProp(grid_zslice);

	// isosurface mapper
	vtkNew<vtkPolyDataMapper> isosurface_mapper;
	isosurface_mapper->SetInputData(isosurfaces);
	// isosurface actor
	vtkNew<vtkActor> isosurface_actor;
	isosurface_actor->SetMapper(isosurface_mapper);
	ren->AddActor(isosurface_actor);

	if (interface)
	{
		// interface mapper
		vtkNew<vtkPointGaussianMapper> interface_mapper;
		interface_mapper->SetInputData(interface);
		interface_mapper->SetScaleFactor(0.0);
		// interface actor
		vtkNew<vtkActor> interface_actor;
		interface_actor->SetMapper(interface_mapper);
		interface_actor->GetProperty()->SetColor(57, 152, 0);
		interface_actor->GetProperty()->SetPointSize(5);
		ren->AddActor(interface_actor);
	}
	if (planar)
	{
		vtkNew<vtkArrowSource> arrow;

		vtkNew<vtkAssignAttribute> vector;
		vector->SetInputData(planar);
		vector->Assign("normal", vtkDataSetAttributes::VECTORS, vtkAssignAttribute::POINT_DATA);
		vector->Update();

		vtkNew<vtkGlyph3D> glyph;
		glyph->SetInputConnection(0, vector->GetOutputPort());
		glyph->SetInputConnection(1, arrow->GetOutputPort());
		glyph->SetVectorModeToUseVector();
		glyph->SetScaleFactor(7448);
		glyph->OrientOn();
		glyph->Update();
	

		vtkNew<vtkPolyDataMapper> planar_mapper;
		planar_mapper->SetInputConnection(glyph->GetOutputPort());
		planar_mapper->ScalarVisibilityOff();

		vtkNew<vtkActor> planar_actor;
		planar_actor->SetMapper(planar_mapper);
		planar_actor->GetProperty()->SetColor(0.6902, 0.7686, 0.8706);
		ren->AddActor(planar_actor);
	}
	if (tangent)
	{
		vtkNew<vtkArrowSource> arrow;

		vtkNew<vtkAssignAttribute> vector;
		vector->SetInputData(tangent);
		vector->Assign("tangent", vtkDataSetAttributes::VECTORS, vtkAssignAttribute::POINT_DATA);
		vector->Update();

		vtkNew<vtkGlyph3D> glyph;
		glyph->SetInputConnection(0, vector->GetOutputPort());
		glyph->SetInputConnection(1, arrow->GetOutputPort());
		glyph->SetVectorModeToUseVector();
		glyph->SetScaleFactor(7448);
		glyph->OrientOn();
		glyph->Update();


		vtkNew<vtkPolyDataMapper> tangent_mapper;
		tangent_mapper->SetInputConnection(glyph->GetOutputPort());
		tangent_mapper->ScalarVisibilityOff();

		vtkNew<vtkActor> tangent_actor;
		tangent_actor->SetMapper(tangent_mapper);
		tangent_actor->GetProperty()->SetColor(0.6902, 0.7686, 0.8706);
		ren->AddActor(tangent_actor);
	}
	if (inequality)
	{
		// interface mapper
		vtkNew<vtkPointGaussianMapper> inequality_mapper;
		inequality_mapper->SetInputData(inequality);
		inequality_mapper->SetScaleFactor(0.0);
		// interface actor
		vtkNew<vtkActor> inequality_actor;
		inequality_actor->SetMapper(inequality_mapper);
		inequality_actor->GetProperty()->SetColor(57, 152, 0);
		inequality_actor->GetProperty()->SetPointSize(5);
		ren->AddActor(inequality_actor);
	}

	iren->Initialize();
	renWin->Render();
	iren->Start();
}
