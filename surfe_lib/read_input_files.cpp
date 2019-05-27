#include <read_input_files.h>

std::string get_file_extension(const char *filename)
{
	std::string str(filename);
	return str.substr(str.find_last_of(".") + 1);
}

vtkSmartPointer<vtkPolyData> GetvtkPolyDataGeometryFromFile(const char *filename)
{
	vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
	std::string extension = get_file_extension(filename);

	if (extension == "vtk")
	{
		vtkNew<vtkPolyDataReader> reader;
		reader->SetFileName(filename);
		reader->Update();
		if (reader->GetOutput()->GetNumberOfPoints() != 0)
			poly = reader->GetOutput();
	}
	if (extension == "vtp")
	{
		vtkNew<vtkXMLPolyDataReader> reader;
		reader->SetFileName(filename);
		reader->Update();
		if (reader->GetOutput()->GetNumberOfPoints() != 0)
			poly = reader->GetOutput();
	}
	if (poly->GetNumberOfPoints() > 0)
		return poly;
	else
		return nullptr;
}

std::vector<std::string> ReadPropertyNamesFromCSV(const char *filename)
{
	std::vector<std::string> property_names;

	ifstream file(filename);

	if (file)
	{
		std::string header;
		std::getline(file, header);

		std::string token;
		std::istringstream tokenStream(header);
		while (std::getline(tokenStream, token, ','))
		{
			property_names.push_back(token);
		}
	}
	else
	{
		std::throw_with_nested(GRBF_Exceptions::error_parsing_csv_header);
	}

	return property_names;
}

std::vector<std::string> ReadPropertyNamesFromVTK(const char *filename)
{
	std::vector<std::string> property_names;

	vtkSmartPointer<vtkPolyData> poly = GetvtkPolyDataGeometryFromFile(filename);

	if (poly)
	{
		if (poly->GetNumberOfPoints() > 0)
		{
			vtkPointData *pd = poly->GetPointData();
			for (int j = 0; j < pd->GetNumberOfArrays(); j++)
				property_names.emplace_back(pd->GetAbstractArray(j)->GetName());
			return property_names;
		}
	}
	
	return property_names;
}

void ConstraintFileReader::SetFilenameAndExtension(const char *filename)
{
	filename_ = filename;
	extension_ = get_file_extension(filename);
}

std::vector<Interface> CSVInterfaceConstraintFileReader::GetConstraints()
{
	std::vector<Interface> interface_constraints;

	std::string xname = GetXName();
	std::string yname = GetYName();
	std::string zname = GetZName();

	if (xname.empty() && yname.empty() && zname.empty())
		std::throw_with_nested(GRBF_Exceptions::missing_coords_in_file);

	if (GetLevelPropertyName().empty())
	{
		io::CSVReader<3> in(filename_);
		in.read_header(io::ignore_missing_column, xname, yname, zname);
		if (!in.has_column(xname) || !in.has_column(yname) || !in.has_column(zname))
			std::throw_with_nested(GRBF_Exceptions::missing_coords_in_file);
		double x;
		double y;
		double z;
		while (in.read_row(x, y, z)) {
			interface_constraints.emplace_back(Interface(x, y, z, 0));
		}
		return interface_constraints;
	}
	else
	{
		std::string levelname = GetLevelPropertyName();

		io::CSVReader<4> in(filename_);
		in.read_header(io::ignore_missing_column, xname, yname, zname,levelname);
		if (!in.has_column(xname) || !in.has_column(yname) || !in.has_column(zname))
			std::throw_with_nested(GRBF_Exceptions::missing_coords_in_file);
		if (!in.has_column(levelname))
			std::throw_with_nested(GRBF_Exceptions::missing_property_info_in_file);
		double x;
		double y;
		double z;
		double level;
		while (in.read_row(x, y, z,level)) {
			interface_constraints.emplace_back(Interface(x, y, z, level));
		}
		return interface_constraints;
	}
}

std::vector<Interface> VTKInterfaceConstraintFileReader::GetConstraints()
{
	std::vector<Interface> interface_constraints;

	vtkSmartPointer<vtkPolyData> poly = GetvtkPolyDataGeometryFromFile(filename_.c_str());
	vtkPointData *point_data = poly->GetPointData();
	vtkDoubleArray *level = nullptr;
	std::string level_name = GetLevelPropertyName();
	if (level_name.empty())
		level = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(level_name.c_str()));

	for (int j = 0; j < poly->GetNumberOfPoints(); j++) {
		double location[3];
		poly->GetPoint(j, location);
		if (level)
			interface_constraints.emplace_back(Interface(location[0], location[1], location[2], level->GetTuple1(j)));
		else
			interface_constraints.emplace_back(Interface(location[0], location[1], location[2], 0));
	}

	return interface_constraints;
}

void InterfaceConstraintFileReader::SearchForDefaultPropertyNames()
{
	// Search for x,y,z, level
	for (const auto &prop_name : property_names_) {
		std::string lowercase_string = prop_name;
		std::transform(lowercase_string.begin(), lowercase_string.end(), lowercase_string.begin(), ::tolower);
		if (lowercase_string == "x")
			SetXName(prop_name.c_str());
		if (lowercase_string == "y")
			SetYName(prop_name.c_str());
		if (lowercase_string == "z")
			SetZName(prop_name.c_str());
		if (lowercase_string == "level")
			SetLevelPropertyName(prop_name);
	}
}

void PlanarConstraintFileReader::SearchForDefaultPropertyNames()
{
	// Search for x,y,z,dip,strike,azimuth,polarity,normal,nx,ny,nz
	for (const auto &prop_name : property_names_) {
		std::string lowercase_string = prop_name;
		std::transform(lowercase_string.begin(), lowercase_string.end(), lowercase_string.begin(), ::tolower);
		if (lowercase_string == "x")
			SetXName(prop_name.c_str());
		if (lowercase_string == "y")
			SetYName(prop_name.c_str());
		if (lowercase_string == "z")
			SetZName(prop_name.c_str());
		if (lowercase_string == "dip")
			SetDipPropertyName(prop_name);
		if (lowercase_string == "strike")
			SetStrikePropertyName(prop_name);
		if (lowercase_string == "azimuth" || lowercase_string == "dip_direction")
			SetAzimuthPropertyName(prop_name);
		if (lowercase_string == "polarity")
			SetPolarityPropertyName(prop_name);
		if (lowercase_string == "normal")
			SetNormalVecPropertyName(prop_name);
		if (lowercase_string == "nx")
			SetNxPropertyName(prop_name);
		if (lowercase_string == "ny")
			SetNyPropertyName(prop_name);
		if (lowercase_string == "nz")
			SetNzPropertyName(prop_name);
	}
}

void PlanarConstraintFileReader::setCaseFromDefaultSearch()
{
	if (!GetNormalVecPropertyName().empty())
	{
		have_dip_strike_polarity = false;
		have_dip_azimuth_polarity = false;
		have_normal = true;
		have_normal_components = false;
		return;
	}
	else if (!GetNxPropertyName().empty() && !GetNyPropertyName().empty() && !GetNzPropertyName().empty())
	{
		have_dip_strike_polarity = false;
		have_dip_azimuth_polarity = false;
		have_normal = false;
		have_normal_components = true;
		return;
	}
	else if (!GetDipPropertyName().empty() && !GetStrikePropertyName().empty() && !GetPolarityPropertyName().empty())
	{
		have_dip_strike_polarity = true;
		have_dip_azimuth_polarity = false;
		have_normal = false;
		have_normal_components = false;
		return;
	}
	else if (!GetDipPropertyName().empty() && !GetAzimuthPropertyName().empty() && !GetPolarityPropertyName().empty())
	{
		have_dip_strike_polarity = false;
		have_dip_azimuth_polarity = true;
		have_normal = false;
		have_normal_components = false;
		return;
	}
}

void TangentConstraintFileReader::setCaseFromDefaultSearch()
{
	if (!GetVectorPropertyName().empty())
	{
		have_vector = true;
		have_vector_components = false;
	}
	else if (!GetVxPropertyName().empty() && !GetVyPropertyName().empty() && !GetVzPropertyName().empty())
	{
		have_vector = false;
		have_vector_components = true;
	}

}

void TangentConstraintFileReader::SearchForDefaultPropertyNames()
{
	// Search for x,y,z,vector,vx,vy,vz
	for (const auto &prop_name : property_names_) {
		std::string lowercase_string = prop_name;
		std::transform(lowercase_string.begin(), lowercase_string.end(),lowercase_string.begin(), ::tolower);
		if (lowercase_string == "x")
			SetXName(prop_name.c_str());
		if (lowercase_string == "y")
			SetYName(prop_name.c_str());
		if (lowercase_string == "z")
			SetZName(prop_name.c_str());
		if (lowercase_string == "vector" || lowercase_string == "tangent")
			SetVectorPropertyName(prop_name);
		if (lowercase_string == "vx" || lowercase_string == "tx")
			SetVxPropertyName(prop_name);
		if (lowercase_string == "vy" || lowercase_string == "ty")
			SetVyPropertyName(prop_name);
		if (lowercase_string == "vz" || lowercase_string == "tz")
			SetVzPropertyName(prop_name);
	}
}

void InequalityConstraintFileReader::SearchForDefaultPropertyNames()
{
	// Search for x,y,z, level
	for (const auto &prop_name : property_names_) {
		std::string lowercase_string = prop_name;
		std::transform(lowercase_string.begin(), lowercase_string.end(), lowercase_string.begin(), ::tolower);
		if (lowercase_string == "x")
			SetXName(prop_name.c_str());
		if (lowercase_string == "y")
			SetYName(prop_name.c_str());
		if (lowercase_string == "z")
			SetZName(prop_name.c_str());
		if (lowercase_string == "level")
			SetLevelPropertyName(prop_name);
	}
}

std::vector<Planar> CSVPlanarConstraintFileReader::GetConstraints()
{
	std::vector<Planar> planar_constraints;

	std::string xname = GetXName();
	std::string yname = GetYName();
	std::string zname = GetZName();

	if (xname.empty() && yname.empty() && zname.empty())
		std::throw_with_nested(GRBF_Exceptions::missing_coords_in_file);

	if (have_normal_components)
	{
		std::string nx = GetNxPropertyName();
		std::string ny = GetNyPropertyName();
		std::string nz = GetNzPropertyName();
		try
		{
			io::CSVReader<6> in(filename_);
			in.read_header(io::ignore_missing_column, xname, yname, zname, nx, ny, nz);
			if (!in.has_column(xname) || !in.has_column(yname) || !in.has_column(zname))
				std::throw_with_nested(GRBF_Exceptions::missing_coords_in_file);
			if (!in.has_column(nx) || !in.has_column(ny) || !in.has_column(nz))
				std::throw_with_nested(GRBF_Exceptions::missing_property_info_in_file);

			double x;
			double y;
			double z;
			double n_x;
			double n_y;
			double n_z;
			while (in.read_row(x, y, z, n_x, n_y, n_z)) {
				planar_constraints.emplace_back(Planar(x, y, z, n_x, n_y, n_z));
			}
		}
		catch (std::exception &e)
		{
			std::throw_with_nested(e);
		}
	}
	else if (have_dip_strike_polarity)
	{
		std::string dipname = GetDipPropertyName();
		std::string strikename = GetStrikePropertyName();
		std::string polarityname = GetPolarityPropertyName();
		try
		{
			io::CSVReader<6> in(filename_);
			in.read_header(io::ignore_missing_column, xname, yname, zname, dipname, strikename, polarityname);
			if (!in.has_column(xname) || !in.has_column(yname) || !in.has_column(zname))
				std::throw_with_nested(GRBF_Exceptions::missing_coords_in_file);
			if (!in.has_column(dipname) || !in.has_column(strikename) || !in.has_column(polarityname))
				std::throw_with_nested(GRBF_Exceptions::missing_property_info_in_file);

			double x;
			double y;
			double z;
			double dip;
			double strike;
			int polariy;
			while (in.read_row(x, y, z, dip, strike, polariy)) {
				planar_constraints.emplace_back(Planar(x, y, z, dip, strike, polariy));
			}
		}
		catch (std::exception &e)
		{
			std::throw_with_nested(e);
		}
	}
	else if (have_dip_azimuth_polarity)
	{
		std::string dipname = GetDipPropertyName();
		std::string azimuthname = GetAzimuthPropertyName();
		std::string polarityname = GetPolarityPropertyName();
		try
		{
			io::CSVReader<6> in(filename_);
			in.read_header(io::ignore_missing_column, xname, yname, zname, dipname, azimuthname, polarityname);
			if (!in.has_column(xname) || !in.has_column(yname) || !in.has_column(zname))
				std::throw_with_nested(GRBF_Exceptions::missing_coords_in_file);
			if (!in.has_column(dipname) || !in.has_column(azimuthname) || !in.has_column(polarityname))
				std::throw_with_nested(GRBF_Exceptions::missing_property_info_in_file);

			double x;
			double y;
			double z;
			double dip;
			double azimuth;
			int polarity;
			while (in.read_row(x, y, z, dip, azimuth, polarity)) {
				double strike;
				if (azimuth >= 90.0)
					strike = azimuth - 90.0;
				else
					strike = azimuth + 270.0;
				planar_constraints.emplace_back(Planar(x, y, z, dip, strike, polarity));
			}
		}
		catch (std::exception &e)
		{
			std::throw_with_nested(e);
		}
	}

	return planar_constraints;
}

std::vector<Planar> VTKPlanarConstraintFileReader::GetConstraints()
{
	std::vector<Planar> planar_constraints;

	// search the geometry for level data
	vtkSmartPointer<vtkPolyData> poly = GetvtkPolyDataGeometryFromFile(filename_.c_str());

	vtkPointData *point_data = poly->GetPointData();
	vtkDoubleArray *normal = nullptr;
	vtkDoubleArray *nx = nullptr;
	vtkDoubleArray *ny = nullptr;
	vtkDoubleArray *nz = nullptr;
	vtkDoubleArray *dip = nullptr;
	vtkDoubleArray *azimuth = nullptr;
	vtkDoubleArray *strike = nullptr;
	vtkDoubleArray *polarity = nullptr;
	if (have_normal)
		normal = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(GetNormalVecPropertyName().c_str()));
	else if (have_normal_components)
	{
		nx = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(GetNxPropertyName().c_str()));
		ny = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(GetNyPropertyName().c_str()));
		nz = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(GetNzPropertyName().c_str()));
	}
	else if (have_dip_strike_polarity)
	{
		dip = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(GetDipPropertyName().c_str()));
		strike = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(GetStrikePropertyName().c_str()));
		polarity = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(GetPolarityPropertyName().c_str()));
	}
	else if (have_dip_azimuth_polarity)
	{
		dip = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(GetDipPropertyName().c_str()));
		azimuth = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(GetAzimuthPropertyName().c_str()));
		polarity = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(GetPolarityPropertyName().c_str()));
	}

	if (normal)
	{
		if (normal->GetNumberOfComponents() != 3)
			std::throw_with_nested(GRBF_Exceptions::input_file_problem);
		else {
			for (int j = 0; j < poly->GetNumberOfPoints(); j++) {
				double location[3];
				poly->GetPoint(j, location);
				double n_x = normal->GetComponent(j, 0);
				double n_y = normal->GetComponent(j, 1);
				double n_z = normal->GetComponent(j, 2);
				planar_constraints.emplace_back(Planar(location[0], location[1], location[2], n_x, n_y, n_z));
			}
		}
		return planar_constraints;
	}
	else if (nx && ny && nz) {
		for (int j = 0; j < poly->GetNumberOfPoints(); j++) {
			double location[3];
			poly->GetPoint(j, location);
			double n_x = nx->GetTuple1(j);
			double n_y = ny->GetTuple1(j);
			double n_z = nz->GetTuple1(j);
			planar_constraints.emplace_back(Planar(location[0], location[1], location[2], n_x, n_y, n_z));
		}
		return planar_constraints;
	}
	else if (dip && strike && polarity) {
		for (int j = 0; j < poly->GetNumberOfPoints(); j++) {
			double location[3];
			poly->GetPoint(j, location);
			double d = dip->GetTuple1(j);
			double s = strike->GetTuple1(j);
			int p = (int)polarity->GetTuple1(j);
			planar_constraints.emplace_back(Planar(location[0], location[1], location[2], d, s, p));
		}
		return planar_constraints;
	}
	else if (dip && azimuth && polarity) {
		for (int j = 0; j < poly->GetNumberOfPoints(); j++) {
			double location[3];
			poly->GetPoint(j, location);
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
	std::throw_with_nested(GRBF_Exceptions::missing_property_info_in_file);
}

std::vector<Tangent> CSVTangentConstraintFileReader::GetConstraints()
{
	std::vector<Tangent> tangent_constraints;

	std::string xname = GetXName();
	std::string yname = GetYName();
	std::string zname = GetZName();

	if (xname.empty() && yname.empty() && zname.empty())
		std::throw_with_nested(GRBF_Exceptions::missing_coords_in_file);

	if (have_vector_components)
	{
		std::string vx_name = GetVxPropertyName();
		std::string vy_name = GetVyPropertyName();
		std::string vz_name = GetVzPropertyName();
		try
		{
			io::CSVReader<6> in(filename_);
			in.read_header(io::ignore_missing_column, xname, yname, zname, vx_name, vy_name, vz_name);
			if (!in.has_column(xname) || !in.has_column(yname) || !in.has_column(zname))
				std::throw_with_nested(GRBF_Exceptions::missing_coords_in_file);
			if (!in.has_column(vx_name) || !in.has_column(vy_name) || !in.has_column(vz_name))
				std::throw_with_nested(GRBF_Exceptions::missing_property_info_in_file);

			double x;
			double y;
			double z;
			double vx;
			double vy;
			double vz;
			while (in.read_row(x, y, z, vx, vy, vz)) {
				tangent_constraints.emplace_back(Tangent(x, y, z, vx, vy, vz));
			}
		}
		catch (std::exception &e)
		{
			std::throw_with_nested(e);
		}
	}

	return tangent_constraints;
}

std::vector<Tangent> VTKTangentConstraintFileReader::GetConstraints()
{
	std::vector<Tangent> tangent_constraints;

	// search the geometry for level data
	vtkSmartPointer<vtkPolyData> poly = GetvtkPolyDataGeometryFromFile(filename_.c_str());

	vtkPointData *point_data = poly->GetPointData();
	vtkDoubleArray *tangent = nullptr;
	vtkDoubleArray *tx = nullptr;
	vtkDoubleArray *ty = nullptr;
	vtkDoubleArray *tz = nullptr;
	if (have_vector)
	{
		tangent = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(GetVectorPropertyName().c_str()));
		if (tangent->GetNumberOfComponents() != 3)
			std::throw_with_nested(GRBF_Exceptions::input_file_problem);
		else {
			for (int j = 0; j < poly->GetNumberOfPoints(); j++) {
				double location[3];
				poly->GetPoint(j, location);
				double t_x = tangent->GetComponent(j, 0);
				double t_y = tangent->GetComponent(j, 1);
				double t_z = tangent->GetComponent(j, 2);
				tangent_constraints.emplace_back(Tangent(location[0], location[1], location[2], t_x, t_y, t_z));
			}
		}
		return tangent_constraints;
	}
	else if (have_vector_components)
	{
		tx = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(GetVxPropertyName().c_str()));
		ty = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(GetVyPropertyName().c_str()));
		tz = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(GetVzPropertyName().c_str()));
		for (int j = 0; j < poly->GetNumberOfPoints(); j++) {
			double location[3];
			poly->GetPoint(j, location);
			tangent_constraints.emplace_back(Tangent(location[0], location[1], location[2], tx->GetTuple1(j), ty->GetTuple1(j), tz->GetTuple1(j)));
		}
		return tangent_constraints;
	}	

	return tangent_constraints;
}

std::vector<Inequality> CSVInequalityConstraintFileReader::GetConstraints()
{
	std::vector<Inequality> inequality_constraints;

	std::string xname = GetXName();
	std::string yname = GetYName();
	std::string zname = GetZName();
	std::string levelname = GetLevelPropertyName();

	if (xname.empty() && yname.empty() && zname.empty())
		std::throw_with_nested(GRBF_Exceptions::missing_coords_in_file);

	try
	{

		io::CSVReader<4> in(filename_);
		in.read_header(io::ignore_missing_column, xname, yname, zname, levelname);
		if (!in.has_column(xname) || !in.has_column(yname) || !in.has_column(zname))
			std::throw_with_nested(GRBF_Exceptions::missing_coords_in_file);
		if (!in.has_column(levelname))
			std::throw_with_nested(GRBF_Exceptions::missing_property_info_in_file);
		double x;
		double y;
		double z;
		double level;
		while (in.read_row(x, y, z, level)) {
			inequality_constraints.emplace_back(Inequality(x, y, z, level));
		}
	}
	catch (std::exception &e)
	{
		std::throw_with_nested(e);
	}

	return  inequality_constraints;
}

std::vector<Inequality> VTKInequalityConstraintFileReader::GetConstraints()
{
	std::vector<Inequality> inequality_constraints;

	// search the geometry for level data
	vtkSmartPointer<vtkPolyData> poly = GetvtkPolyDataGeometryFromFile(filename_.c_str());

	if (!poly)
		std::throw_with_nested(GRBF_Exceptions::input_file_problem);

	vtkPointData *point_data = poly->GetPointData();
	vtkDoubleArray *level = vtkDoubleArray::SafeDownCast(point_data->GetAbstractArray(GetLevelPropertyName().c_str()));
	if (!level)
		std::throw_with_nested(GRBF_Exceptions::missing_property_info_in_file);

	for (int j = 0; j < poly->GetNumberOfPoints(); j++) {
		double location[3];
		poly->GetPoint(j, location);
		inequality_constraints.emplace_back(Inequality(location[0], location[1], location[2], level->GetTuple1(j)));
	}
	return inequality_constraints;
}
