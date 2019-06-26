#include <geo_builder.h>

InputParameters Geo_Builder::getGUIParameters()
{
	return InputImpl::GetDialogParameters();
}

void Geo_Builder::progress(const float &progress_value)
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

void Geo_Builder::build_constraints_from_input_files()
{
	try
	{
		if (!input_.interface_file.empty()) {
			std::string extension = get_file_extension(input_.interface_file.c_str());
			if (extension == "csv")
			{
				CSVInterfaceConstraintFileReader reader =
					CSVInterfaceConstraintFileReader::CreateUsingDefaultPropertyNames(surfe, input_.interface_file.c_str());
				reader.GetConstraints();
			}
			else if (extension == "vtp" || extension == "vtk")
			{
				VTKInterfaceConstraintFileReader reader =
					VTKInterfaceConstraintFileReader::CreateUsingDefaultPropertyNames(surfe, input_.interface_file.c_str());
				reader.GetConstraints();
			}
		}
		if (!input_.inequality_file.empty()) {
			std::vector<Inequality> inequality_constraints;
			std::string extension = get_file_extension(input_.inequality_file.c_str());
			if (extension == "csv")
			{
				CSVInequalityConstraintFileReader reader =
					CSVInequalityConstraintFileReader::CreateUsingDefaultPropertyNames(surfe, input_.inequality_file.c_str());
				reader.GetConstraints();
			}
			else if (extension == "vtp" || extension == "vtk")
			{
				VTKInequalityConstraintFileReader reader =
					VTKInequalityConstraintFileReader::CreateUsingDefaultPropertyNames(surfe, input_.inequality_file.c_str());
				reader.GetConstraints();
			}
		}

		if (!input_.planar_file.empty()) {
			std::vector<Planar> planar_constraints;
			std::string extension = get_file_extension(input_.planar_file.c_str());
			if (extension == "csv")
			{
				CSVPlanarConstraintFileReader reader =
					CSVPlanarConstraintFileReader::CreateUsingDefaultPropertyNames(surfe, input_.planar_file.c_str());
				reader.GetConstraints();
			}
			else if (extension == "vtp" || extension == "vtk")
			{
				VTKPlanarConstraintFileReader reader =
					VTKPlanarConstraintFileReader::CreateUsingDefaultPropertyNames(surfe, input_.planar_file.c_str());
				reader.GetConstraints();
			}
		}
		if (!input_.tangent_file.empty()) {
			std::vector<Tangent> tangent_constraints;
			std::string extension = get_file_extension(input_.planar_file.c_str());
			if (extension == "csv")
			{
				CSVTangentConstraintFileReader reader =
					CSVTangentConstraintFileReader::CreateUsingDefaultPropertyNames(surfe, input_.tangent_file.c_str());
				reader.GetConstraints();
			}
			else if (extension == "vtp" || extension == "vtk")
			{
				VTKTangentConstraintFileReader reader =
					VTKTangentConstraintFileReader::CreateUsingDefaultPropertyNames(surfe, input_.tangent_file.c_str());
				reader.GetConstraints();
			}
		}
	}
	catch (const std::exception&e)
	{
		std::rethrow_if_nested(e);
	}
}

void Geo_Builder::CreateGRBFInterpolantFromGUIParameters()
{
	input_ = getGUIParameters();

	surfe = new Surfe_API(input_.parameters);

	LoadConstraintsFromFiles();
}

void Geo_Builder::LoadConstraintsFromFiles()
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

void Geo_Builder::SetInterfaceDataFile(const char *interfaceFile)
{
	input_.interface_file = interfaceFile;

	evaluation_completed_ = false;
}

void Geo_Builder::SetPlanarDataFile(const char *planarFile)
{
	input_.planar_file = planarFile;

	evaluation_completed_ = false;

}

void Geo_Builder::SetTangentDataFile(const char *tangentFile)
{
	input_.tangent_file = tangentFile;

	evaluation_completed_ = false;
}

void Geo_Builder::SetInequalityDataFile(const char *inequalityFile)
{
	input_.inequality_file = inequalityFile;

	evaluation_completed_ = false;
}

void Geo_Builder::BuildRegularGrid(const double &zmin, const double &zmax, const double &resolution, const double &xy_percent_padding /*= 0*/)
{
	SpatialParameters spatial = surfe->GetDataBoundsAndResolution();
	double xmin = spatial.xmin;
	double xmax = spatial.xmax;
	double ymin = spatial.ymin;
	double ymax = spatial.ymax;

	if (xy_percent_padding != 0 && xy_percent_padding < 100 && xy_percent_padding > 0)
	{
		double dx = xmax - xmin;
		double dy = ymax - ymin;
		xmin -= dx * (xy_percent_padding / 100.0);
		xmax += dx * (xy_percent_padding / 100.0);
		ymin -= dy * (xy_percent_padding / 100.0);
		ymax += dy * (xy_percent_padding / 100.0);
	}

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

void Geo_Builder::BuildRegularGrid(
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

void Geo_Builder::BuildRegularGrid(const double &xy_percent_padding /*= 0*/)
{
	try
	{
		// get spatial parameters
		SpatialParameters spatial = surfe->GetDataBoundsAndResolution();

		double bounds[6];
		bounds[0] = spatial.xmin;
		bounds[1] = spatial.xmax;
		bounds[2] = spatial.ymin;
		bounds[3] = spatial.ymax;
		bounds[4] = spatial.zmin;
		bounds[5] = spatial.zmax;
		if (xy_percent_padding != 0 && xy_percent_padding < 100 && xy_percent_padding > 0)
		{
			double dx = spatial.xmax - spatial.xmin;
			double dy = spatial.ymax - spatial.ymin;
			double dz = spatial.zmax - spatial.zmin;
			bounds[0] -= dx * (xy_percent_padding / 100.0);
			bounds[1] += dx * (xy_percent_padding / 100.0);
			bounds[2] -= dy * (xy_percent_padding / 100.0);
			bounds[3] += dy * (xy_percent_padding / 100.0);
			bounds[4] -= dy * (xy_percent_padding / 100.0);
			bounds[5] += dy * (xy_percent_padding / 100.0);
		}
		if (spatial.resolution == 0)
			throw "problem_computing_grid";

		int nx = (bounds[1] - bounds[0]) / spatial.resolution;
		int ny = (bounds[3] - bounds[2]) / spatial.resolution;
		int nz = (bounds[5] - bounds[4]) / spatial.resolution;

		if (nx == 0 || ny == 0 || nz == 0)
			throw "problem_computing_grid";

		double origin[3] = { bounds[0],bounds[2], bounds[4] };

		vtkSmartPointer<vtkImageData> constructed_grid = vtkSmartPointer<vtkImageData>::New();
		constructed_grid->SetDimensions(nx + 1, ny + 1, nz + 1);
		constructed_grid->SetOrigin(origin);
		constructed_grid->SetSpacing(spatial.resolution, spatial.resolution, spatial.resolution);

		grid_ = constructed_grid;

	}
	catch (const std::exception&)
	{
		throw;
	}
}

void Geo_Builder::BuildRegularGrid(const double &resolution, const double &xy_percent_padding /*= 0*/)
{
	try
	{
		// get spatial parameters
		SpatialParameters spatial = surfe->GetDataBoundsAndResolution();

		double bounds[6];
		bounds[0] = spatial.xmin;
		bounds[1] = spatial.xmax;
		bounds[2] = spatial.ymin;
		bounds[3] = spatial.ymax;
		bounds[4] = spatial.zmin;
		bounds[5] = spatial.zmax;
		if (xy_percent_padding != 0 && xy_percent_padding < 100 && xy_percent_padding > 0)
		{
			double dx = spatial.xmax - spatial.xmin;
			double dy = spatial.ymax - spatial.ymin;
			double dz = spatial.zmax - spatial.zmin;
			bounds[0] -= dx * (xy_percent_padding / 100.0);
			bounds[1] += dx * (xy_percent_padding / 100.0);
			bounds[2] -= dy * (xy_percent_padding / 100.0);
			bounds[3] += dy * (xy_percent_padding / 100.0);
			bounds[4] -= dy * (xy_percent_padding / 100.0);
			bounds[5] += dy * (xy_percent_padding / 100.0);
		}
		if (spatial.resolution == 0)
			throw "problem_computing_grid";

		int nx = (bounds[1] - bounds[0]) / resolution;
		int ny = (bounds[3] - bounds[2]) / resolution;
		int nz = (bounds[5] - bounds[4]) / resolution;

		if (nx == 0 || ny == 0 || nz == 0)
			throw "problem_computing_grid";

		double origin[3] = { bounds[0],bounds[2], bounds[4] };

		vtkSmartPointer<vtkImageData> constructed_grid = vtkSmartPointer<vtkImageData>::New();
		constructed_grid->SetDimensions(nx + 1, ny + 1, nz + 1);
		constructed_grid->SetOrigin(origin);
		constructed_grid->SetSpacing(resolution, resolution, resolution);

		grid_ = constructed_grid;

	}
	catch (const std::exception&)
	{
		throw;
	}
}

vtkSmartPointer<vtkImageData> Geo_Builder::GetEvaluatedGrid()
{
	if (!grid_)
		throw "No vtkImageData grid exists";

	if (evaluation_completed_)
		return grid_;

	if (!surfe->InterpolantComputed())
	{
		try
		{
			surfe->ComputeInterpolant();
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

	// 	vtkSmartPointer<vtkDoubleArray> normfield = vtkSmartPointer<vtkDoubleArray>::New();
	// 	normfield->SetName("Gradient Norm");
	// 	normfield->SetNumberOfComponents(1);
	// 	normfield->SetNumberOfTuples(grid_->GetNumberOfPoints());

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
		// evaluate scalar field at point
		double scalar_field = surfe->EvaluateInterpolantAtPoint(point[0], point[1], point[2]);

		sfield->SetTuple1(j, scalar_field);
		// 		method_->eval_vector_interpolant_at_point(pt);
		// 		double nx = pt.nx_interp();
		// 		double ny = pt.ny_interp();
		// 		double nz = pt.nz_interp();
		// 		double norm = sqrt(nx*nx + ny * ny + nz * nz);
		// 		normfield->SetTuple1(j, norm);

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
	std::cout << std::endl;
	grid_->GetPointData()->SetScalars(sfield);
	//grid_->GetPointData()->AddArray(normfield);

	evaluation_completed_ = true;

	cout << "Finished evaluating interpolant in grid" << endl;

	return grid_;
}

vtkSmartPointer<vtkPolyData> Geo_Builder::GetIsoSurfaces()
{
	if (!grid_)
		throw "No vtkImageData grid exists";

	try
	{
		GetEvaluatedGrid();
	}
	catch (const std::exception&)
	{
		throw;
	}

	vtkNew<vtkMarchingCubes> mcube;
	mcube->SetInputData(grid_);
	mcube->ComputeScalarsOn();

	MatrixXd interface_ref_points = surfe->GetInterfaceReferencePoints();
	for (int j = 0; j < interface_ref_points.rows(); j++) {
		// evaluate interpolant at point found in row j
		double iso_value = surfe->EvaluateInterpolantAtPoint(interface_ref_points(j,0), interface_ref_points(j, 1), interface_ref_points(j, 2));
		mcube->SetValue(j, iso_value);
	}
	mcube->Update();

	vtkSmartPointer<vtkPolyData> iso_surfaces = vtkSmartPointer<vtkPolyData>::New();
	iso_surfaces = mcube->GetOutput();

	cout << "Finished marching cubes" << endl;

	return iso_surfaces;
}

vtkSmartPointer<vtkPolyData> Geo_Builder::GetInterfaceConstraints()
{
	MatrixXd interface = surfe->GetInterfaceConstraints();

	if (interface.rows() != 0)
	{
		// interface
		vtkSmartPointer<vtkPolyData> interface_constraints = vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPoints> interface_points = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkDoubleArray> interface_levels_array = vtkSmartPointer<vtkDoubleArray>::New();
		interface_levels_array->SetName("level");
		interface_levels_array->SetNumberOfComponents(1);
		interface_levels_array->SetNumberOfTuples(interface.rows());
		for (int j = 0; j < interface.rows(); j++) {
			interface_points->InsertNextPoint(interface(j, 0), interface(j, 1), interface(j, 2));
			interface_levels_array->SetTuple1(j, interface(j, 3));
			// 			vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
			// 			vertex->GetPointIds()->SetId(0, j);
			// 			interface_constraints->Allocate(1, 1);
			// 			interface_constraints->InsertNextCell(vertex->GetCellType(), vertex->GetPointIds());
		}
		interface_constraints->SetPoints(interface_points);
		interface_constraints->GetPointData()->AddArray(interface_levels_array);
		return interface_constraints;
	}
	else
		return nullptr;
}

vtkSmartPointer<vtkPolyData> Geo_Builder::GetPlanarConstraints()
{
	MatrixXd planar = surfe->GetPlanarConstraints();

	if (planar.rows() != 0)
	{
		// planar
		vtkSmartPointer<vtkPolyData> planar_constraints = vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPoints> planar_points = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkDoubleArray> planar_normal_array = vtkSmartPointer<vtkDoubleArray>::New();
		int n_tuples = planar.rows();
		planar_normal_array->SetName("normal");
		planar_normal_array->SetNumberOfTuples(n_tuples);
		planar_normal_array->SetNumberOfComponents(3);
		planar_normal_array->SetComponentName(0, "nx");
		planar_normal_array->SetComponentName(1, "ny");
		planar_normal_array->SetComponentName(2, "nz");
		// initalization for vector data
		for (int k = 0; k < n_tuples; k++) {
			for (int l = 0; l < 3; l++)
				planar_normal_array->InsertComponent(k, l, 0.0);
		}
		for (int j = 0; j < n_tuples; j++) {
			planar_points->InsertNextPoint(planar(j, 0), planar(j, 1), planar(j, 2));
			double normal[3] = { planar(j, 3), planar(j, 4), planar(j, 5) };
			planar_normal_array->SetTuple(j, normal);
			// 			vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
			// 			vertex->GetPointIds()->SetId(0, j);
			// 			planar_constraints->Allocate(1, 1);
			// 			planar_constraints->InsertNextCell(vertex->GetCellType(), vertex->GetPointIds());
		}
		planar_constraints->SetPoints(planar_points);
		planar_constraints->GetPointData()->AddArray(planar_normal_array);
		return planar_constraints;
	}
	else
		return nullptr;
}

vtkSmartPointer<vtkPolyData> Geo_Builder::GetTangentConstraints()
{
	MatrixXd tangent = surfe->GetTangentConstraints();

	if (tangent.rows() != 0)
	{
		// tangent
		vtkSmartPointer<vtkPolyData> tangent_constraints = vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPoints> tangent_points = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkDoubleArray> vector_array = vtkSmartPointer<vtkDoubleArray>::New();
		int n_tuples = tangent.rows();
		vector_array->SetName("tangent");
		vector_array->SetNumberOfComponents(3);
		vector_array->SetComponentName(0, "tx");
		vector_array->SetComponentName(1, "ty");
		vector_array->SetComponentName(2, "tz");
		// initalization for vector data
		for (int k = 0; k < n_tuples; k++) {
			for (int l = 0; l < 3; l++)
				vector_array->InsertComponent(k, l, 0.0);
		}
		for (int j = 0; j < n_tuples; j++) {
			tangent_points->InsertNextPoint(tangent(j, 0), tangent(j, 1), tangent(j, 2));
			double vector[3] = { tangent(j, 3), tangent(j, 4), tangent(j, 5) };
			vector_array->SetTuple(j, vector);
			// 			vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
			// 			vertex->GetPointIds()->SetId(0, j);
			// 			tangent_constraints->Allocate(1, 1);
			// 			tangent_constraints->InsertNextCell(vertex->GetCellType(), vertex->GetPointIds());
		}
		tangent_constraints->SetPoints(tangent_points);
		tangent_constraints->GetPointData()->AddArray(vector_array);
		return tangent_constraints;
	}
	else
		return nullptr;
}

vtkSmartPointer<vtkPolyData> Geo_Builder::GetInequalityConstraints()
{
	MatrixXd inequality = surfe->GetInequalityConstraints();

	if (inequality.rows() != 0)
	{
		// inequality
		vtkSmartPointer<vtkPolyData> inequality_constraints = vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPoints> inequality_points = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkDoubleArray> inequality_level_array = vtkSmartPointer<vtkDoubleArray>::New();
		inequality_level_array->SetName("level");
		inequality_level_array->SetNumberOfComponents(1);
		inequality_level_array->SetNumberOfTuples(inequality.rows());
		for (int j = 0; j < inequality.rows(); j++) {
			inequality_points->InsertNextPoint(inequality(j, 0), inequality(j, 1), inequality(j, 2));
			inequality_level_array->SetTuple1(j, inequality(j, 3));
			// 			vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
			// 			vertex->GetPointIds()->SetId(0, j);
			// 			inequality_constraints->Allocate(1, 1);
			// 			inequality_constraints->InsertNextCell(vertex->GetCellType(), vertex->GetPointIds());
		}
		inequality_constraints->SetPoints(inequality_points);
		inequality_constraints->GetPointData()->AddArray(inequality_level_array);
		return inequality_constraints;
	}
	else
		return nullptr;
}

void Geo_Builder::WriteVTKInterfaceConstraints(const char *filename)
{
	vtkSmartPointer<vtkPolyData> poly = GetInterfaceConstraints();

	if (poly) {
		vtkNew<vtkXMLPolyDataWriter> writer;
		writer->SetInputData(poly);
		writer->SetFileName(filename);
		writer->Write();
	}
}

void Geo_Builder::WriteVTKPlanarConstraints(const char *filename)
{
	vtkSmartPointer<vtkPolyData> poly = GetPlanarConstraints();

	if (poly) {
		vtkNew<vtkXMLPolyDataWriter> writer;
		writer->SetInputData(poly);
		writer->SetFileName(filename);
		writer->Write();
	}
}

void Geo_Builder::WriteVTKTangentConstraints(const char *filename)
{
	vtkSmartPointer<vtkPolyData> poly = GetTangentConstraints();

	if (poly) {
		vtkNew<vtkXMLPolyDataWriter> writer;
		writer->SetInputData(poly);
		writer->SetFileName(filename);
		writer->Write();
	}
}

void Geo_Builder::WriteVTKInequalityConstraints(const char *filename)
{
	vtkSmartPointer<vtkPolyData> poly = GetInequalityConstraints();

	if (poly) {
		vtkNew<vtkXMLPolyDataWriter> writer;
		writer->SetInputData(poly);
		writer->SetFileName(filename);
		writer->Write();
	}
}

void Geo_Builder::WriteVTKEvaluationGrid(const char *filename)
{

	if (!grid_)
		throw "No vtkImageData grid exists";

	if (!evaluation_completed_) {
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
	writer->SetFileName(filename);
	writer->SetDataModeToBinary();
	writer->Write();
}

void Geo_Builder::WriteVTKIsoSurfaces(const char *filename)
{
	vtkSmartPointer<vtkPolyData> isosurfaces = GetIsoSurfaces();

	if (isosurfaces) {
		vtkNew<vtkXMLPolyDataWriter> writer;
		writer->SetInputData(isosurfaces);
		writer->SetFileName(filename);
		writer->SetDataModeToBinary();
		writer->Write();
	}
}

void Geo_Builder::VisualizeVTKData()
{
	// get grid
	vtkImageData *grid = GetEvaluatedGrid();
	double spacing[3];
	grid->GetSpacing(spacing);
	double min_scale = *std::max_element(spacing, spacing + 3);
	// 	SpatialParameters spatial = compute_constraint_bounds_and_resolution();
	// 	min_scale = spatial.resolution;

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

	vtkNew<vtkCellPicker> picker;
	picker->SetTolerance(0.005);

	vtkNew<vtkImagePlaneWidget> x_plane_widget;
	x_plane_widget->SetInputData(grid);
	x_plane_widget->SetInteractor(iren);
	x_plane_widget->SetPicker(picker);
	x_plane_widget->RestrictPlaneToVolumeOn();
	x_plane_widget->SetLookupTable(lut);
	x_plane_widget->SetResliceInterpolateToLinear();
	x_plane_widget->SetPlaneOrientation(0);
	x_plane_widget->SetSliceIndex(0);
	x_plane_widget->SetDefaultRenderer(ren);
	x_plane_widget->On();
	x_plane_widget->InteractionOn();

	vtkNew<vtkImagePlaneWidget> y_plane_widget;
	y_plane_widget->SetInputData(grid);
	y_plane_widget->SetInteractor(iren);
	y_plane_widget->SetPicker(picker);
	y_plane_widget->RestrictPlaneToVolumeOn();
	y_plane_widget->SetLookupTable(lut);
	y_plane_widget->SetResliceInterpolateToLinear();
	y_plane_widget->SetPlaneOrientation(1);
	y_plane_widget->SetSliceIndex(dimensions[1] / 2);
	y_plane_widget->SetDefaultRenderer(ren);
	y_plane_widget->On();
	y_plane_widget->InteractionOn();

	vtkNew<vtkImagePlaneWidget> z_plane_widget;
	z_plane_widget->SetInputData(grid);
	z_plane_widget->SetInteractor(iren);
	z_plane_widget->SetPicker(picker);
	z_plane_widget->RestrictPlaneToVolumeOn();
	z_plane_widget->SetLookupTable(lut);
	z_plane_widget->SetResliceInterpolateToLinear();
	z_plane_widget->SetPlaneOrientation(2);
	z_plane_widget->SetSliceIndex(0);
	z_plane_widget->SetDefaultRenderer(ren);
	z_plane_widget->On();
	z_plane_widget->InteractionOn();

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
		glyph->SetScaleFactor(min_scale * 3);
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
		glyph->SetScaleFactor(min_scale * 3);
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
