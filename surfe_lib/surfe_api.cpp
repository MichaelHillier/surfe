#include <surfe_api.h>

#include <read_csv_files.h>

#include <vtkNew.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLImageDataWriter.h>

GRBF_Modelling_Methods* Surfe_API::get_method(const UI_Parameters& params)
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

void Surfe_API::build_constraints_from_csv_files()
{
	if (!method_)
		throw GRBF_Exceptions::grbf_method_is_null;
	try
	{
		if (strlen(params_.interface_file) != 0) {
			std::vector<Interface> interface_constraints;
			interface_constraints = build_interface_constraints(params_.interface_file);
			method_->constraints.itrface = interface_constraints;
		}
		if (strlen(params_.inequality_file) != 0) {
			std::vector<Inequality> inequality_constraints;
			inequality_constraints = build_inequality_constraints(params_.inequality_file);
			method_->constraints.inequality = inequality_constraints;
		}

		if (strlen(params_.planar_file) != 0) {
			std::vector<Planar> planar_constraints;
			planar_constraints = build_planar_constraints(params_.planar_file);
			method_->constraints.planar = planar_constraints;
		}
		if (strlen(params_.tangent_file) != 0) {
			std::vector<Tangent> tangent_constraints;
			tangent_constraints = build_tangent_constraints(params_.tangent_file);
			method_->constraints.tangent = tangent_constraints;
		}
	}
	catch (const std::exception&e)
	{
		SurfeExceptions exceptions(e);
		throw exceptions.what();
	}

	constraint_files_changed_ = false; // since they have been loaded
	constraints_changed_ = true;
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

Surfe_API::Surfe_API(const UI_Parameters& params)
{
	params_ = params;

	method_ = nullptr;
	grid_ = nullptr;
	have_interpolant_ = false;
	have_method_ = true;
	evaluation_completed_ = false;
	parameters_changed_ = true;
	constraint_files_changed_ = true;
	constraints_changed_ = false;

	method_ = get_method(params_);

	try
	{
		build_constraints_from_csv_files();
	}
	catch (const std::exception&e)
	{
		throw;
	}
	evaluation_completed_ = false;
}

void Surfe_API::GetUIParametersAndConstraints()
{
	params_ = InputImpl::GetDialogParameters();
	parameters_changed_ = true;
	constraint_files_changed_ = true;
	method_ = get_method(params_);
	try
	{
		build_constraints_from_csv_files();
	}
	catch (const std::exception&e)
	{
		throw;
	}
}

void Surfe_API::LoadConstraintsFromFiles()
{
	method_ = get_method(params_);
	try
	{
		build_constraints_from_csv_files();
	}
	catch (const std::exception&e)
	{
		throw;
	}
}

void Surfe_API::AddInterfaceConstraint(const double &x, const double &y, const double &z, const double &level)
{
	Interface interface_constraint(x, y, z, level);
	method_->constraints.itrface.push_back(interface_constraint);
	constraints_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::AddPlanarConstraintwNormal(const double &x, const double &y, const double &z, const double &nx, const double &ny, const double &nz)
{
	Planar planar_constraint(x, y, z, nx, ny, nz);
	method_->constraints.planar.push_back(planar_constraint);
	constraints_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::AddPlanarConstraintwStrikeDipPolarity(const double &x, const double &y, const double &z, const double &strike, const double &dip, const double &polarity)
{
	Planar planar_constraint(x, y, z, dip, strike, polarity);
	method_->constraints.planar.push_back(planar_constraint);
	constraints_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::AddPlanarConstraintwAzimuthDipPolarity(const double &x, const double &y, const double &z, const double &azimuth, const double &dip, const double &polarity)
{
	double strike = 0.0;
	// convert azimuth to strike
	if (azimuth >= 90.0)
		strike = azimuth - 90.0;
	else
		strike = azimuth + 270.0;
	Planar planar_constraint(x, y, z, dip, strike, polarity);
	method_->constraints.planar.push_back(planar_constraint);
	constraints_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::AddTangentConstraint(const double &x, const double &y, const double &z, const double &tx, const double &ty, const double &tz)
{
	Tangent tangent_constraint(x, y, z, tx, ty, tz);
	method_->constraints.tangent.push_back(tangent_constraint);
	constraints_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::AddInequalityConstraint(const double &x, const double &y, const double &z, const double &level)
{
	Inequality inequality_constraint(x, y, z, level);
	method_->constraints.inequality.push_back(inequality_constraint);
	constraints_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::ComputeInterpolant()
{
	if (constraint_files_changed_)
	{
		try
		{
			build_constraints_from_csv_files();
		}
		catch (const std::exception&e)
		{
			throw;
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

	have_interpolant_ = true;
	constraints_changed_ = false;
	parameters_changed_ = false;
}

void Surfe_API::SetRBFKernel(const Parameter_Types::RBF &rbf)
{
	params_.basis_type = rbf;
	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetRBFKernel(const char *rbf_name)
{
	if (rbf_name == "r3")
		params_.basis_type = Parameter_Types::RBF::Cubic;
	else if (rbf_name == "WendlandC2")
		params_.basis_type = Parameter_Types::RBF::WendlandC2;
	else if (rbf_name == "r")
		params_.basis_type = Parameter_Types::RBF::R;
	else if (rbf_name == "Gaussian")
		params_.basis_type = Parameter_Types::RBF::Gaussian;
	else if (rbf_name == "Multiquadratics")
		params_.basis_type = Parameter_Types::RBF::MQ;
	else if (rbf_name == "Thin Plate Spline")
		params_.basis_type = Parameter_Types::RBF::TPS;
	else if (rbf_name = "Inverse Multiquadratics")
		params_.basis_type = Parameter_Types::RBF::IMQ; // Inverse Multiquadratics
	else
		throw GRBF_Exceptions::unknown_rbf;

	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetRBFShapeParameter(const double &shape_param)
{
	params_.shape_parameter = shape_param;
	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetPolynomialOrder(const int &poly_order)
{
	params_.polynomial_order = poly_order;
	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetGlobalAnisotropy(const bool &g_anisotropy)
{
	params_.model_global_anisotropy = g_anisotropy;
	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetGreedy(const bool &greedy)
{
	params_.use_greedy = greedy;
	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetRestrictedRange(const bool &rr)
{
	params_.use_restricted_range = rr;
	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetRegressionSmoothing(const bool &rs)
{
	params_.use_regression_smoothing = rs;
	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetInterfaceUncertainty(const double &interface_uncertainty)
{
	params_.interface_uncertainty = interface_uncertainty;
	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetAngularUncertainty(const double &angular_uncertainty)
{
	params_.angular_uncertainty = angular_uncertainty;
	parameters_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetInterfaceDataCSVFile(const char *interface_file)
{
	auto params_interface_file_length = strlen(params_.interface_file);
	if (params_interface_file_length > 0) {
		delete[] params_.interface_file;
		params_.interface_file = nullptr;
	}

	auto interface_file_length = strlen(interface_file);
	char *temp_interface_filename = new char[interface_file_length + 1];
	strncpy_s(
		temp_interface_filename,
		interface_file_length + 1,
		interface_file,
		interface_file_length
	);
	temp_interface_filename[interface_file_length] = '\0';

	params_.interface_file = temp_interface_filename;
	constraint_files_changed_ = true;
	constraints_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetPlanarDataCSVFile(const char *planar_file)
{
	auto params_planar_file_length = strlen(params_.planar_file);
	if (params_planar_file_length > 0) {
		delete[] params_.planar_file;
		params_.planar_file = nullptr;
	}

	auto planar_file_length = strlen(planar_file);
	char *temp_planar_filename = new char[planar_file_length + 1];
	strncpy_s(
		temp_planar_filename,
		planar_file_length + 1,
		planar_file,
		planar_file_length
	);
	temp_planar_filename[planar_file_length] = '\0';

	params_.planar_file = temp_planar_filename;
	constraint_files_changed_ = true;
	constraints_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetTangentDataCSVFile(const char *tangent_file)
{
	auto params_tangent_file_length = strlen(params_.tangent_file);
	if (params_tangent_file_length > 0) {
		delete[] params_.tangent_file;
		params_.tangent_file = nullptr;
	}

	auto tangent_file_length = strlen(tangent_file);
	char *temp_tangent_filename = new char[tangent_file_length + 1];
	strncpy_s(
		temp_tangent_filename,
		tangent_file_length + 1,
		tangent_file,
		tangent_file_length
	);
	temp_tangent_filename[tangent_file_length] = '\0';

	params_.tangent_file = temp_tangent_filename;
	constraint_files_changed_ = true;
	constraints_changed_ = true;
	evaluation_completed_ = false;
}

void Surfe_API::SetInequalityDataCSVFile(const char *inequality_file)
{
	auto params_inequality_file_length = strlen(params_.inequality_file);
	if (params_inequality_file_length > 0) {
		delete[] params_.inequality_file;
		params_.inequality_file = nullptr;
	}

	auto inequality_file_length = strlen(inequality_file);
	char *temp_inequality_filename = new char[inequality_file_length + 1];
	strncpy_s(
		temp_inequality_filename,
		inequality_file_length + 1,
		inequality_file,
		inequality_file_length
	);
	temp_inequality_filename[inequality_file_length] = '\0';

	params_.inequality_file = temp_inequality_filename;
	constraint_files_changed_ = true;
	constraints_changed_ = true;
	evaluation_completed_ = false;
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

double *Surfe_API::EvaluateVectorInterpolantAtPoint(const double &x, const double &y, const double &z)
{
	if (have_interpolant_)
	{
		if (constraints_changed_ || parameters_changed_)
			throw GRBF_Exceptions::interpolant_needs_update;

		double gradient[3] = { 0,0,0 };
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
	constructed_grid->SetSpacing(resolution,resolution,resolution);

	grid_ = constructed_grid;
}

vtkImageData * Surfe_API::GetEvaluatedGrid()
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
			throw;
		}
	}

	vtkSmartPointer<vtkDoubleArray> sfield = vtkSmartPointer<vtkDoubleArray>::New();
	sfield->SetName("Scalar Field");
	sfield->SetNumberOfComponents(1);
	sfield->SetNumberOfTuples(grid_->GetNumberOfPoints());
	
	#pragma omp parallel for schedule(dynamic)
	for (int j = 0; j < grid_->GetNumberOfPoints(); j++) {
		double point[3];
		grid_->GetPoint(j, point);
		Point pt(point[0], point[1], point[2]);
		// evaluate scalar field at point
		method_->eval_scalar_interpolant_at_point(pt);
		double scalar_field = pt.scalar_field();
		sfield->SetTuple1(j, scalar_field);
	}
	grid_->GetPointData()->SetScalars(sfield);

	evaluation_completed_ = true;

	cout << " Finished evaluating interpolant in grid" << endl;

	return grid_;
}

vtkPolyData * Surfe_API::GetIsoSurfaces()
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

	return iso_surfaces;
}

vtkPolyData * Surfe_API::GetInterfaceConstraints()
{
	if (!method_->constraints.itrface.empty())
	{
		// interface
		vtkSmartPointer<vtkPolyData> interface_constraints = vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPoints> interface_points = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkDoubleArray> interface_scalar_field = vtkSmartPointer<vtkDoubleArray>::New();
		interface_scalar_field->SetName("Scalar Field");
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

vtkPolyData	* Surfe_API::GetPlanarConstraints()
{
	if (!method_->constraints.planar.empty())
	{
		// planar
		vtkSmartPointer<vtkPolyData> planar_constraints = vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPoints> planar_points = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkDoubleArray> planar_gradient_field = vtkSmartPointer<vtkDoubleArray>::New();
		int n_tuples = (int)method_->constraints.planar.size();
		planar_gradient_field->SetName("Gradient Field");
		planar_gradient_field->SetNumberOfTuples(n_tuples);
		planar_gradient_field->SetNumberOfComponents(3);
		planar_gradient_field->SetComponentName(0, "Gx");
		planar_gradient_field->SetComponentName(1, "Gy");
		planar_gradient_field->SetComponentName(2, "Gz");
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

vtkPolyData * Surfe_API::GetTangentConstraints()
{
	if (!method_->constraints.tangent.empty())
	{
		// tangent
		vtkSmartPointer<vtkPolyData> tangent_constraints = vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPoints> tangent_points = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkDoubleArray> tangent_vector = vtkSmartPointer<vtkDoubleArray>::New();
		tangent_vector->SetName("Tangent Vector");
		tangent_vector->SetNumberOfComponents(3);
		tangent_vector->SetComponentName(0, "Tx");
		tangent_vector->SetComponentName(1, "Ty");
		tangent_vector->SetComponentName(2, "Tz");
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

vtkPolyData * Surfe_API::GetInequalityConstraints()
{
	if (!method_->constraints.inequality.empty())
	{
		// inequality
		vtkSmartPointer<vtkPolyData> inequality_constraints = vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPoints> inequality_points = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkDoubleArray> inequality_scalar_field = vtkSmartPointer<vtkDoubleArray>::New();
		inequality_scalar_field->SetName("Scalar Field");
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
		vtkPolyData *isosurfaces = GetIsoSurfaces();
		vtkNew<vtkXMLPolyDataWriter> writer;
		writer->SetInputData(isosurfaces);
		writer->WriteToOutputStringOn();
		writer->Write();

		vtk_isosurfaces_string_ = writer->GetOutputString();
		return vtk_isosurfaces_string_.c_str();
	}
	catch (std::exception& e)
	{
		throw;
	}
}

const char * Surfe_API::GetVTKInterfaceConstraintsString()
{
	vtkPolyData *poly = GetInterfaceConstraints();
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
	vtkPolyData *poly = GetPlanarConstraints();
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
	vtkPolyData *poly = GetTangentConstraints();
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
	vtkPolyData *poly = GetInequalityConstraints();
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
	vtkPolyData *poly = GetInterfaceConstraints();
	if (poly) {
		vtkNew<vtkXMLPolyDataWriter> writer;
		writer->SetInputData(poly);
		writer->SetFileName(filename);
		writer->Write();
	}
}

void Surfe_API::WriteVTKPlanarConstraints(const char *filename)
{
	vtkPolyData *poly = GetPlanarConstraints();
	if (poly) {
		vtkNew<vtkXMLPolyDataWriter> writer;
		writer->SetInputData(poly);
		writer->SetFileName(filename);
		writer->Write();
	}
}

void Surfe_API::WriteVTKTangentConstraints(const char *filename)
{
	vtkPolyData *poly = GetTangentConstraints();
	if (poly) {
		vtkNew<vtkXMLPolyDataWriter> writer;
		writer->SetInputData(poly);
		writer->SetFileName(filename);
		writer->Write();
	}
}

void Surfe_API::WriteVTKInequalityConstraints(const char *filename)
{
	vtkPolyData *poly = GetInequalityConstraints();
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
