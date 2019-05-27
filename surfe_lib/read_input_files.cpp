#include <read_input_files.h>

std::string get_file_extension(const char *filename)
{
	std::string str(filename);
	return str.substr(str.find_last_of(".") + 1);
}


std::vector<Interface> build_interface_constraints(const char *interface_file)
{
	std::vector<Interface> interface_constraints;
	std::string filename_extension = get_file_extension(interface_file);

	if (filename_extension == "csv")
	{
		try
		{
			io::CSVReader<4> in(interface_file);
			in.read_header(io::ignore_missing_column, "x", "y", "z", "level");
			if (!in.has_column("x") || !in.has_column("y") || !in.has_column("z"))
				std::throw_with_nested(GRBF_Exceptions::missing_coords_in_file);
			if (!in.has_column("level"))
				std::throw_with_nested(GRBF_Exceptions::missing_property_info_in_file);
			double x;
			double y;
			double z;
			double level;
			while (in.read_row(x, y, z, level)) {
				interface_constraints.emplace_back(Interface(x, y, z, level));
			}
			return interface_constraints;
		}
		catch (std::exception &e)
		{
			std::throw_with_nested(e);
		}
	}
	else if (filename_extension == "vtk" || filename_extension == "vtp")
	{
		vtkSmartPointer<vtkPolyData> input = vtkSmartPointer<vtkPolyData>::New();
		if (filename_extension == "vtk")
		{
			vtkNew<vtkPolyDataReader> reader;
			reader->SetFileName(interface_file);
			reader->Update();
			input = reader->GetOutput();
		}
		if (filename_extension == "vtp")
		{
			vtkNew<vtkXMLPolyDataReader> reader;
			reader->SetFileName(interface_file);
			reader->Update();
			input = reader->GetOutput();
		}
		if (input->GetNumberOfPoints() != 0)
		{
			// search the geometry for level data
			vtkPointData *point_data = input->GetPointData();
			vtkDoubleArray *level_data = nullptr;
			for (int j = 0; j < point_data->GetNumberOfArrays(); j++) {
				std::string property_name = point_data->GetArrayName(j);
				if (property_name == "level" || property_name == "scalar field")
				{
					level_data = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(j));
					break;
				}
			}
			if (level_data)
			{
				for (int j = 0; j < input->GetNumberOfPoints(); j++) {
					double location[3];
					input->GetPoint(j, location);
					interface_constraints.emplace_back(Interface(location[0], location[1], location[2], level_data->GetTuple1(j)));
				}
				return interface_constraints;
			}
			else // no level property. if mode is single surface we set level = 0
			{
				for (int j = 0; j < input->GetNumberOfPoints(); j++) {
					double location[3];
					input->GetPoint(j, location);
					interface_constraints.emplace_back(Interface(location[0], location[1], location[2], 0));
				}
				return interface_constraints;
			}
		}
		else
			std::throw_with_nested(GRBF_Exceptions::input_file_problem);
	}
	else
		std::throw_with_nested(GRBF_Exceptions::unsupported_file);
}

std::vector<Planar> build_planar_constraints(const char *planar_file)
{
	std::vector<Planar> planar_constraints;
	std::string filename_extension = get_file_extension(planar_file);


	bool have_normal = false;
	bool have_strike_dip_polarity = false;
	bool have_azimuth_dip_polarity = false;

	if (filename_extension == "csv")
	{
		try
		{
			io::CSVReader<10> in(planar_file);
			in.read_header(io::ignore_missing_column, "x", "y", "z", "nx", "ny", "nz", "dip", "strike", "azimuth", "polarity");
			if (!in.has_column("x") || !in.has_column("y") || !in.has_column("z"))
				std::throw_with_nested(GRBF_Exceptions::missing_coords_in_file);
			if (in.has_column("nx") && in.has_column("ny") && in.has_column("nz"))
				have_normal = true;
			else if (in.has_column("strike") && in.has_column("dip") && in.has_column("polarity"))
				have_strike_dip_polarity = true;
			else if (in.has_column("azimuth") && in.has_column("dip") && in.has_column("polarity"))
				have_azimuth_dip_polarity = true;
			else
			{
				std::throw_with_nested(GRBF_Exceptions::missing_property_info_in_file);
			}
			double x;
			double y;
			double z;
			double nx;
			double ny;
			double nz;
			double dip;
			double azimuth;
			double strike;
			int polarity;
			while (in.read_row(x, y, z, nx, ny, nz, dip, strike, azimuth, polarity)) {
				if (have_normal)
					planar_constraints.emplace_back(Planar(x, y, z, nx, ny, nz));
				else if (have_strike_dip_polarity)
					planar_constraints.emplace_back(Planar(x, y, z, dip, strike, polarity));
				else if (have_azimuth_dip_polarity) {
					if (azimuth >= 90.0)
						strike = azimuth - 90.0;
					else
						strike = azimuth + 270.0;
					planar_constraints.emplace_back(Planar(x, y, z, dip, strike, polarity));
				}
			}
			return planar_constraints;
		}
		catch (std::exception &e)
		{
			std::throw_with_nested(e);
		}
	}
	else if (filename_extension == "vtk" || filename_extension == "vtp")
	{
		vtkSmartPointer<vtkPolyData> input = vtkSmartPointer<vtkPolyData>::New();
		if (filename_extension == "vtk")
		{
			vtkNew<vtkPolyDataReader> reader;
			reader->SetFileName(planar_file);
			reader->Update();
			if (reader->GetOutput()->GetNumberOfPoints() != 0)
				input = reader->GetOutput();
		}
		if (filename_extension == "vtp")
		{
			vtkNew<vtkXMLPolyDataReader> reader;
			reader->SetFileName(planar_file);
			reader->Update();
			if (reader->GetOutput()->GetNumberOfPoints() != 0)
				input = reader->GetOutput();
		}
		if (input->GetNumberOfPoints() != 0)
		{
			// search the geometry for level data
			vtkPointData *point_data = input->GetPointData();
			vtkDoubleArray *normal = nullptr;
			vtkDoubleArray *nx = nullptr;
			vtkDoubleArray *ny = nullptr;
			vtkDoubleArray *nz = nullptr;
			vtkDoubleArray *dip = nullptr;
			vtkDoubleArray *azimuth = nullptr;
			vtkDoubleArray *strike = nullptr;
			vtkDoubleArray *polarity = nullptr;
			for (int j = 0; j < point_data->GetNumberOfArrays(); j++) {
				std::string property_name = point_data->GetArrayName(j);
				if (property_name == "nx")
					nx = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(j));
				if (property_name == "ny")
					ny = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(j));
				if (property_name == "nz")
					nz = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(j));
				if (property_name == "dip")
					dip = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(j));
				if (property_name == "strike")
					strike = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(j));
				if (property_name == "azimuth")
					azimuth = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(j));
				if (property_name == "polarity")
					polarity = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(j));
				if (property_name == "normal")
					normal = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(j));
			}

			if (dip && strike && polarity) {
				for (int j = 0; j < input->GetNumberOfPoints(); j++) {
					double location[3];
					input->GetPoint(j, location);
					double d = dip->GetTuple1(j);
					double s = strike->GetTuple1(j);
					int p = (int)polarity->GetTuple1(j);
					planar_constraints.emplace_back(Planar(location[0], location[1], location[2], d, s, p));
				}
				return planar_constraints;
			}
			if (dip && azimuth && polarity) {
				for (int j = 0; j < input->GetNumberOfPoints(); j++) {
					double location[3];
					input->GetPoint(j, location);
					double d = dip->GetTuple1(j);
					double a = azimuth->GetTuple1(j);
					double s;
					if (a >= 90.0)
						s = a - 90.0;
					else
						s = a + 270.0;
					int p = (int)polarity->GetTuple1(j);
					planar_constraints.emplace_back(Planar(location[0], location[1], location[2], d, s, p));
				}
				return planar_constraints;
			}
			if (nx && ny && nz) {
				for (int j = 0; j < input->GetNumberOfPoints(); j++) {
					double location[3];
					input->GetPoint(j, location);
					double n_x = nx->GetTuple1(j);
					double n_y = ny->GetTuple1(j);
					double n_z = nz->GetTuple1(j);
					planar_constraints.emplace_back(Planar(location[0], location[1], location[2], n_x, n_y, n_z));
				}
				return planar_constraints;
			}
			if (normal)
			{
				if (normal->GetNumberOfComponents() != 3)
					std::throw_with_nested(GRBF_Exceptions::input_file_problem);
				else {
					for (int j = 0; j < input->GetNumberOfPoints(); j++) {
						double location[3];
						input->GetPoint(j, location);
						double n_x = normal->GetComponent(j, 0);
						double n_y = normal->GetComponent(j, 1);
						double n_z = normal->GetComponent(j, 2);
						planar_constraints.emplace_back(Planar(location[0], location[1], location[2], n_x, n_y, n_z));
					}
				}
				return planar_constraints;
			}
			std::throw_with_nested(GRBF_Exceptions::missing_property_info_in_file);
		}
		else
			std::throw_with_nested(GRBF_Exceptions::input_file_problem);
	}
	else
		std::throw_with_nested(GRBF_Exceptions::unsupported_file);
}

std::vector<Tangent> build_tangent_constraints(const char *tangent_file)
{
	std::vector<Tangent> tangent_constraints;
	std::string filename_extension = get_file_extension(tangent_file);

	if (filename_extension == "csv")
	{
		try
		{
			bool have_normal = false;
			bool have_strike_dip_polarity = false;
			bool have_azimuth_dip_polarity = false;

			io::CSVReader<6> in(tangent_file);
			in.read_header(io::ignore_missing_column, "x", "y", "z", "tx", "ty", "tz");
			if (!in.has_column("x") || !in.has_column("y") || !in.has_column("z"))
				std::throw_with_nested(GRBF_Exceptions::missing_coords_in_file);
			if (!in.has_column("tx") || !in.has_column("ty") || !in.has_column("tz"))
				std::throw_with_nested(GRBF_Exceptions::missing_property_info_in_file);
			double x;
			double y;
			double z;
			double tx;
			double ty;
			double tz;
			while (in.read_row(x, y, z, tx, ty, tz)) {
				tangent_constraints.emplace_back(Tangent(x, y, z, tx, ty, tz));
			}
			return tangent_constraints;
		}
		catch (std::exception &e)
		{
			std::throw_with_nested(e);
		}
	}
	else if (filename_extension == "vtk" || filename_extension == "vtp")
	{
		vtkSmartPointer<vtkPolyData> input = vtkSmartPointer<vtkPolyData>::New();
		if (filename_extension == "vtk")
		{
			vtkNew<vtkPolyDataReader> reader;
			reader->SetFileName(tangent_file);
			reader->Update();
			if (reader->GetOutput()->GetNumberOfPoints() != 0)
				input = reader->GetOutput();
		}
		if (filename_extension == "vtp")
		{
			vtkNew<vtkXMLPolyDataReader> reader;
			reader->SetFileName(tangent_file);
			reader->Update();
			if (reader->GetOutput()->GetNumberOfPoints() != 0)
				input = reader->GetOutput();
		}
		if (input->GetNumberOfPoints() != 0)
		{
			// search the geometry for level data
			vtkPointData *point_data = input->GetPointData();
			vtkDoubleArray *tangent = nullptr;
			vtkDoubleArray *tx = nullptr;
			vtkDoubleArray *ty = nullptr;
			vtkDoubleArray *tz = nullptr;
			for (int j = 0; j < point_data->GetNumberOfArrays(); j++) {
				std::string property_name = point_data->GetArrayName(j);
				if (property_name == "tx" || property_name == "vx")
					tx = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(j));
				if (property_name == "ty" || property_name == "vy")
					ty = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(j));
				if (property_name == "tz" || property_name == "vz")
					tz = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(j));
				if (property_name == "tangent" || property_name == "vector" || property_name == "linear")
					tangent = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(j));
			}
			if (tx && ty && tz) {
				for (int j = 0; j < input->GetNumberOfPoints(); j++) {
					double location[3];
					input->GetPoint(j, location);
					double t_x = tx->GetTuple1(j);
					double t_y = ty->GetTuple1(j);
					double t_z = tz->GetTuple1(j);
					tangent_constraints.emplace_back(Tangent(location[0], location[1], location[2], t_x, t_y, t_z));
				}
				return tangent_constraints;
			}
			if (tangent)
			{
				if (tangent->GetNumberOfComponents() != 3)
					std::throw_with_nested(GRBF_Exceptions::input_file_problem);
				else {
					for (int j = 0; j < input->GetNumberOfPoints(); j++) {
						double location[3];
						input->GetPoint(j, location);
						double t_x = tangent->GetComponent(j, 0);
						double t_y = tangent->GetComponent(j, 1);
						double t_z = tangent->GetComponent(j, 2);
						tangent_constraints.emplace_back(Tangent(location[0], location[1], location[2], t_x, t_y, t_z));
					}
				}
				return tangent_constraints;
			}
			std::throw_with_nested(GRBF_Exceptions::missing_property_info_in_file);
		}
		else
			std::throw_with_nested(GRBF_Exceptions::input_file_problem);
	}
	else
		std::throw_with_nested(GRBF_Exceptions::unsupported_file);
}

std::vector<Inequality> build_inequality_constraints(const char *inequality_file)
{
	std::vector<Inequality> inequality_constraints;
	std::string filename_extension = get_file_extension(inequality_file);

	if (filename_extension == "csv")
	{
		try
		{
			io::CSVReader<4> in(inequality_file);
			in.read_header(io::ignore_missing_column, "x", "y", "z", "level");
			if (!in.has_column("x") || !in.has_column("y") || !in.has_column("z"))
				std::throw_with_nested(GRBF_Exceptions::missing_coords_in_file);
			if (!in.has_column("level"))
				std::throw_with_nested(GRBF_Exceptions::missing_property_info_in_file);
			double x;
			double y;
			double z;
			double level;
			while (in.read_row(x, y, z, level)) {
				inequality_constraints.emplace_back(Inequality(x, y, z, level));
			}
			return inequality_constraints;
		}
		catch (std::exception &e)
		{
			std::throw_with_nested(e);
		}
	}
	else if (filename_extension == "vtk" || filename_extension == "vtp")
	{
		vtkSmartPointer<vtkPolyData> input = vtkSmartPointer<vtkPolyData>::New();
		if (filename_extension == "vtk")
		{
			vtkNew<vtkPolyDataReader> reader;
			reader->SetFileName(inequality_file);
			reader->Update();
			if (reader->GetOutput()->GetNumberOfPoints() != 0)
				input = reader->GetOutput();
		}
		if (filename_extension == "vtp")
		{
			vtkNew<vtkXMLPolyDataReader> reader;
			reader->SetFileName(inequality_file);
			reader->Update();
			if (reader->GetOutput()->GetNumberOfPoints() != 0)
				input = reader->GetOutput();
		}
		if (input->GetNumberOfPoints() != 0)
		{
			// search the geometry for level data
			vtkPointData *point_data = input->GetPointData();
			vtkDoubleArray *level_data = nullptr;
			for (int j = 0; j < point_data->GetNumberOfArrays(); j++) {
				std::string property_name = point_data->GetArrayName(j);
				if (property_name == "level")
				{
					level_data = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(j));
					break;
				}
			}
			if (level_data)
			{
				for (int j = 0; j < input->GetNumberOfPoints(); j++) {
					double location[3];
					input->GetPoint(j, location);
					inequality_constraints.emplace_back(Inequality(location[0], location[1], location[2], level_data->GetTuple1(j)));
				}
				return inequality_constraints;
			}
			else
				std::throw_with_nested(GRBF_Exceptions::missing_property_info_in_file);
		}
		else
			std::throw_with_nested(GRBF_Exceptions::input_file_problem);
	}
	else
		std::throw_with_nested(GRBF_Exceptions::unsupported_file);
}

void ConstraintFileReader::SetFilenameAndExtension(const char *filename)
{
	filename_ = filename;
	extension_ = get_file_extension(filename);
}

std::vector<Interface> CSVInterfaceConstraintFileReader::GetConstraints()
{
	std::vector<Interface> interface_constraints;

	if (GetXName().empty() && GetYName().empty() && GetZName().empty())
		std::throw_with_nested(GRBF_Exceptions::missing_coords_in_file);

	rapidcsv::Document file(filename_);
	std::vector<double> x = file.GetColumn<double>(GetXName());
	std::vector<double> y = file.GetColumn<double>(GetYName());
	std::vector<double> z = file.GetColumn<double>(GetZName());
	std::vector<double> level;
	if (!GetLevelPropertyName().empty())
		level = file.GetColumn<double>(GetLevelPropertyName());

	if (x.size() == y.size() && x.size() == z.size())
	{
		for (int j = 0; j < (int)x.size(); j++) {
			if (!GetLevelPropertyName().empty())
				if (level.size() == x.size())
					interface_constraints.emplace_back(Interface(x[j], y[j], z[j], level[j]));
				else
					std::throw_with_nested(GRBF_Exceptions::input_file_problem);
			else
				interface_constraints.emplace_back(Interface(x[j], y[j], z[j], 0));
		}
	}
	else
		std::throw_with_nested(GRBF_Exceptions::input_file_problem);

	return interface_constraints;
}

void CSVInterfaceConstraintFileReader::SearchForDefaultPropertyNames()
{
	// find default property names
	for (const auto &prop_name : property_names_) {
		if (prop_name == "x" || prop_name == "X")
			SetXName(prop_name.c_str());
		if (prop_name == "y" || prop_name == "Y")
			SetYName(prop_name.c_str());
		if (prop_name == "z" || prop_name == "Z")
			SetZName(prop_name.c_str());
		if (prop_name == "level" || prop_name == "Level")
			SetLevelPropertyName(prop_name.c_str());
	}
}

std::vector<Interface> VTKInterfaceConstraintFileReader::GetConstraints()
{
	std::vector<Interface> interface_constraints;

	vtkPointData *point_data = constraints_vtk->GetPointData();
	vtkDoubleArray *level = nullptr;
	if (!GetLevelPropertyName().empty()) {
		std::string level_name = GetLevelPropertyName();
		level = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(level_name.c_str()));
	}

	for (int j = 0; j < constraints_vtk->GetNumberOfPoints(); j++) {
		double location[3];
		constraints_vtk->GetPoint(j, location);
		if (level)
			interface_constraints.emplace_back(Interface(location[0], location[1], location[2], level->GetTuple1(j)));
		else
			interface_constraints.emplace_back(Interface(location[0], location[1], location[2], 0));
	}

	return interface_constraints;
}

void VTKInterfaceConstraintFileReader::SearchForDefaultPropertyNames()
{
	// find default property names
	for (const auto &prop_name : property_names_) {
		if (prop_name == "level" || prop_name == "Level")
			SetLevelPropertyName(prop_name.c_str());
	}
}

void CSVConstraintFileReader::ReadPropertyNames()
{
	rapidcsv::Document file(filename_);
	property_names_ = file.GetColumnNames();
}

void CSVConstraintFileReader::SearchForDefaultCoordinatePropertyNames()
{
	for (const auto &prop_name : property_names_) {
		if (prop_name == "x" || prop_name == "X")
			SetXName(prop_name.c_str());
		if (prop_name == "y" || prop_name == "Y")
			SetYName(prop_name.c_str());
		if (prop_name == "z" || prop_name == "Z")
			SetZName(prop_name.c_str());
	}
}

void VTKConstraintFileReader::ReadPropertyNames()
{
	constraints_vtk = vtkSmartPointer<vtkPolyData>::New();

	if (extension_ == "vtk")
	{
		vtkNew<vtkPolyDataReader> reader;
		reader->SetFileName(filename_.c_str());
		reader->Update();
		if (reader->GetOutput()->GetNumberOfPoints() != 0)
			constraints_vtk = reader->GetOutput();
	}
	if (extension_ == "vtp")
	{
		vtkNew<vtkXMLPolyDataReader> reader;
		reader->SetFileName(filename_.c_str());
		reader->Update();
		if (reader->GetOutput()->GetNumberOfPoints() != 0)
			constraints_vtk = reader->GetOutput();
	}
	if (constraints_vtk->GetNumberOfPoints() > 0)
	{
		vtkPointData *pd = constraints_vtk->GetPointData();
		for (int j = 0; j < pd->GetNumberOfArrays(); j++)
			property_names_.emplace_back(pd->GetAbstractArray(j)->GetName());
	}
}

void PlanarConstraintFileReader::SearchForDefaultPropertyNames()
{
// 	for (const auto &prop_name : property_names_) {
// 		if (prop_name == "x")
// 			SetXName(prop_name.c_str());
// 		if (prop_name == "X")
// 			SetXName(prop_name.c_str());
// 		if (prop_name == "y")
// 			SetYName(prop_name.c_str());
// 		if (prop_name == "Y")
// 			SetYName(prop_name.c_str());
// 		if (prop_name == "z")
// 			SetZName(prop_name.c_str());
// 		if (prop_name == "Z")
// 			SetZName(prop_name.c_str());
// 	}
}

void CSVPlanarConstraintFileReader::SearchForDefaultPropertyNames()
{
// 	for (const auto &prop_name : property_names_) {
//  		if (prop_name == "x")
// 			SetXName(prop_name.c_str());
// 		if (prop_name == "X")
//  			SetXName(prop_name.c_str());
//  		if (prop_name == "y")
//  			SetYName(prop_name.c_str());
//  		if (prop_name == "Y")
//  			SetYName(prop_name.c_str());
//  		if (prop_name == "z")
//  			SetZName(prop_name.c_str());
//  		if (prop_name == "Z")
//  			SetZName(prop_name.c_str());
// 		if ()
//  	}
}

void InterfaceConstraintFileReader::SearchForDefaultPropertyNames()
{
	// find default property names
	for (const auto &prop_name : property_names_) {
		if (prop_name == "level" || prop_name == "Level")
			SetLevelPropertyName(prop_name.c_str());
	}
}
