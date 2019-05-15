#include <surfe_api.h>

#include <read_csv_files.h>

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

vtkDataObjectCollection * Surfe_API::convert_constraints_to_vtk()
{
	vtkSmartPointer<vtkDataObjectCollection> collection = vtkSmartPointer<vtkDataObjectCollection>::New();

	// inequality
	vtkSmartPointer<vtkPolyData> inequality_constraints = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> inequality_points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkDoubleArray> inequality_scalar_field = vtkSmartPointer<vtkDoubleArray>::New();
	inequality_scalar_field->SetName("Scalar Field");
	inequality_scalar_field->SetNumberOfComponents(1);
	inequality_scalar_field->SetNumberOfTuples(sgrid->GetNumberOfPoints());
	for (int j = 0; j < (int)method_->constraints.inequality.size(); j++ ){
		Inequality *inequality_pt = &method_->constraints.inequality[j];
		inequality_points->InsertNextPoint(inequality_pt->x(), inequality_pt->y(), inequality_pt->z());
		inequality_scalar_field->SetTuple1(j, inequality_pt->scalar_field());
		vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
		vertex->GetPointIds()->SetId(0, j);
		inequality_constraints->InsertNextCell(vertex->GetCellType(), vertex->GetPointIds());
	}
	inequality_constraints->SetPoints(inequality_points);
	inequality_constraints->GetPointData()->AddArray(inequality_scalar_field);
	collection->AddItem(inequality_constraints);
	
	// interface
	vtkSmartPointer<vtkPolyData> interface_constraints = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> interface_points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkDoubleArray> interface_scalar_field = vtkSmartPointer<vtkDoubleArray>::New();
	interface_scalar_field->SetName("Scalar Field");
	interface_scalar_field->SetNumberOfComponents(1);
	interface_scalar_field->SetNumberOfTuples(sgrid->GetNumberOfPoints());
	for (int j = 0; j < (int)method_->constraints.itrface.size(); j++) {
		Interface *interface_pt = &method_->constraints.itrface[j];
		interface_points->InsertNextPoint(interface_pt->x(), interface_pt->y(), interface_pt->z());
		interface_scalar_field->SetTuple1(j, interface_pt->scalar_field());
		vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
		vertex->GetPointIds()->SetId(0, j);
		interface_constraints->InsertNextCell(vertex->GetCellType(), vertex->GetPointIds());
	}
	interface_constraints->SetPoints(inequality_points);
	interface_constraints->GetPointData()->AddArray(interface_scalar_field);
	collection->AddItem(interface_constraints);

	// planar
	vtkSmartPointer<vtkPolyData> planar_constraints = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> planar_points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkDoubleArray> planar_gradient_field = vtkSmartPointer<vtkDoubleArray>::New();
	planar_gradient_field->SetName("Gradient Field");
	planar_gradient_field->SetNumberOfComponents(3);
	planar_gradient_field->SetComponentName(0, "Gx");
	planar_gradient_field->SetComponentName(1, "Gy");
	planar_gradient_field->SetComponentName(2, "Gz");
	for (int j = 0; j < (int)method_->constraints.planar.size(); j++) {
		Planar *planar_pt = &method_->constraints.planar[j];
		planar_points->InsertNextPoint(planar_pt->x(), planar_pt->y(), planar_pt->z());
		double gradient[3] = { planar_pt->nx(),planar_pt->ny(),planar_pt->nz() };
		planar_gradient_field->SetTuple(j, gradient);
		vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
		vertex->GetPointIds()->SetId(0, j);
		planar_constraints->InsertNextCell(vertex->GetCellType(), vertex->GetPointIds());
	}
	planar_constraints->SetPoints(planar_points);
	planar_constraints->GetPointData()->AddArray(planar_gradient_field);
	collection->AddItem(planar_constraints);

	// planar
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
		vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
		vertex->GetPointIds()->SetId(0, j);
		tangent_constraints->InsertNextCell(vertex->GetCellType(), vertex->GetPointIds());
	}
	tangent_constraints->SetPoints(tangent_points);
	tangent_constraints->GetPointData()->AddArray(tangent_vector);
	collection->AddItem(tangent_constraints);

	return collection;
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

void Surfe_API::GetUIParameters()
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

void Surfe_API::AddInterfaceConstraint(const Interface& pt)
{
	method_->constraints.itrface.push_back(pt);
	constraints_changed_ = true;
}

void Surfe_API::AddPlanarConstraint(const Planar& planar_pt)
{
	method_->constraints.planar.push_back(planar_pt);
	constraints_changed_ = true;
}

void Surfe_API::AddTangentConstraint(const Tangent& tangent_pt)
{
	method_->constraints.tangent.push_back(tangent_pt);
	constraints_changed_ = true;
}

void Surfe_API::AddInequalityConstraint(const Inequality& inequality_pt)
{
	method_->constraints.inequality.push_back(inequality_pt);
	constraints_changed_ = true;
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
}

void Surfe_API::SetRBFShapeParameter(const double &shape_param)
{
	params_.shape_parameter = shape_param;
	parameters_changed_ = true;
}

void Surfe_API::SetPolynomialOrder(const int &poly_order)
{
	params_.polynomial_order = poly_order;
	parameters_changed_ = true;
}

void Surfe_API::SetGlobalAnisotropy(const bool &g_anisotropy)
{
	params_.model_global_anisotropy = g_anisotropy;
	parameters_changed_ = true;
}

void Surfe_API::SetGreedy(const bool &greedy)
{
	params_.use_greedy = greedy;
	parameters_changed_ = true;
}

void Surfe_API::SetRestrictedRange(const bool &rr)
{
	params_.use_restricted_range = rr;
	parameters_changed_ = true;
}

void Surfe_API::SetRegressionSmoothing(const bool &rs)
{
	params_.use_regression_smoothing = rs;
	parameters_changed_ = true;
}

void Surfe_API::SetInterfaceUncertainty(const double &interface_uncertainty)
{
	params_.interface_uncertainty = interface_uncertainty;
	parameters_changed_ = true;
}

void Surfe_API::SetAngularUncertainty(const double &angular_uncertainty)
{
	params_.angular_uncertainty = angular_uncertainty;
	parameters_changed_ = true;
}

void Surfe_API::SetInterfaceDataCSVFile(const char *interface_file)
{
	params_.interface_file = interface_file;
	constraint_files_changed_ = true;
}

void Surfe_API::SetPlanarDataCSVFile(const char *planar_file)
{
	params_.planar_file = planar_file;
	constraint_files_changed_ = true;
}

void Surfe_API::SetTangentDataCSVFile(const char *tangent_file)
{
	params_.tangent_file = tangent_file;
	constraint_files_changed_ = true;
}

void Surfe_API::SetInequalityDataCSVFile(const char *inequality_file)
{
	params_.inequality_file = inequality_file;
	constraint_files_changed_ = true;
}

double Surfe_API::EvaluateInterpolantAtPoint(const double &x, const double &y, const double &z)
{
	if (have_interpolant_)
	{
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

void Surfe_API::ConstructRegularGridOutput(const double &zmin, const double &zmax, const double &resolution, const double &xy_percent_padding /*= 0*/)
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

	vtkSmartPointer<vtkStructuredGrid> s_grid = vtkSmartPointer<vtkStructuredGrid>::New();
	vtkSmartPointer<vtkPoints> grid_points = vtkSmartPointer<vtkPoints>::New();
	s_grid->SetDimensions(nx + 1, ny + 1, nz + 1);

	// create grid points
	for (int j = 0; j < (nz + 1); j++) {
		for (int k = 0; k < (ny + 1); k++) {
			for (int l = 0; l < (nx + 1); l++) {
				double grid_pt[3];
				grid_pt[0] = origin[0] + l * resolution;
				grid_pt[1] = origin[1] + k * resolution;
				grid_pt[2] = origin[2] + j * resolution;
				grid_points->InsertNextPoint(grid_pt);
			}
		}
	}
	s_grid->SetPoints(grid_points);

	sgrid = s_grid;
}

void Surfe_API::ConstructRegularGridOutput(
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

	vtkSmartPointer<vtkStructuredGrid> s_grid = vtkSmartPointer<vtkStructuredGrid>::New();
	vtkSmartPointer<vtkPoints> grid_points = vtkSmartPointer<vtkPoints>::New();
	s_grid->SetDimensions(nx + 1, ny + 1, nz + 1);

	// create grid points
	for (int j = 0; j < (nz + 1); j++) {
		for (int k = 0; k < (ny + 1); k++) {
			for (int l = 0; l < (nx + 1); l++) {
				double grid_pt[3];
				grid_pt[0] = origin[0] + l * resolution;
				grid_pt[1] = origin[1] + k * resolution;
				grid_pt[2] = origin[2] + j * resolution;
				grid_points->InsertNextPoint(grid_pt);
			}
		}
	}
	s_grid->SetPoints(grid_points);

	sgrid = s_grid;
}

vtkStructuredGrid * Surfe_API::GetEvaluatedvtkStructuredGrid()
{
	if (!sgrid)
		throw GRBF_Exceptions::no_sgrid_exists;

	if (!have_interpolant_)
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
	sfield->SetNumberOfTuples(sgrid->GetNumberOfPoints());
	
	#pragma omp parallel for schedule(dynamic)
	for (int j = 0; j < sgrid->GetNumberOfPoints(); j++) {
		double point[3];
		sgrid->GetPoint(j, point);
		Point pt(point[0], point[1], point[2]);
		// evaluate scalar field at point
		method_->eval_scalar_interpolant_at_point(pt);
		double scalar_field = pt.scalar_field();
		sfield->SetTuple1(j, scalar_field);
	}

	sgrid->GetPointData()->AddArray(sfield);

	evaluation_completed_ = true;

	return sgrid;
}

vtkDataObjectCollection * Surfe_API::GetConstraintsAndOutputAsVTKObjects()
{
	if (!sgrid)
		throw GRBF_Exceptions::no_sgrid_exists;

	// Get Evaluated structured grid using the interpolant
	try
	{
		sgrid = GetEvaluatedvtkStructuredGrid();
	}
	catch (std::exception& e) 
	{
		throw;
	}

	// Get Constraints as vtk objects
	vtkDataObjectCollection *constraints = convert_constraints_to_vtk();

	// Get Iso surfaces
	vtkPolyData *iso_surfaces = GetIsoSurfacesAsvtkPolyData();

	constraints->AddItem(iso_surfaces);
	constraints->AddItem(sgrid);

	return constraints;
}

vtkPolyData * Surfe_API::GetIsoSurfacesAsvtkPolyData()
{
	if (!sgrid)
		throw GRBF_Exceptions::no_sgrid_exists;

	if (!evaluation_completed_) {
		try
		{
			sgrid = GetEvaluatedvtkStructuredGrid();
		}
		catch (const std::exception&)
		{
			throw;
		}
	}

	vtkNew<vtkMarchingCubes> mcube;
	mcube->SetInputData(sgrid);
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
