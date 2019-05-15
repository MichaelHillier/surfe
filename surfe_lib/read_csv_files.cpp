#include <read_csv_files.h>

std::vector<Interface> build_interface_constraints(const char *interface_file)
{
	std::vector<Interface> interface_constraints;

	try
	{
		io::CSVReader<4> in(interface_file);
		in.read_header(io::ignore_extra_column, "x", "y", "z", "level");
		if (!in.has_column("x") || !in.has_column("y") || !in.has_column("z"))
			std::throw_with_nested(GRBF_Exceptions::missing_coords_in_file);
		if (!in.has_column("level"))
			std::throw_with_nested(GRBF_Exceptions::missing_level_info_in_file);
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

std::vector<Planar> build_planar_constraints(const char *planar_file)
{
	std::vector<Planar> planar_constraints;

	try
	{
		bool have_normal = false;
		bool have_strike_dip_polarity = false;
		bool have_azimuth_dip_polarity = false;

		io::CSVReader<10> in(planar_file);
		in.read_header(io::ignore_extra_column, "x", "y", "z", "nx", "ny", "nz", "dip", "strike", "azimuth", "polarity");
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
			std::throw_with_nested(GRBF_Exceptions::missing_orientation_in_file);
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
		while (in.read_row(x, y, z, nx, ny, nz, dip, azimuth, strike, polarity)) {
			if (have_normal)
				planar_constraints.emplace_back(Planar(x, y, z, nx, ny, nz));
			else if (have_strike_dip_polarity)
				planar_constraints.emplace_back(Planar(x, y, z, dip, strike, polarity));
			else {
				// convert azimuth to strike
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

std::vector<Tangent> build_tangent_constraints(const char *tangent_file)
{
	std::vector<Tangent> tangent_constraints;

	try
	{
		bool have_normal = false;
		bool have_strike_dip_polarity = false;
		bool have_azimuth_dip_polarity = false;

		io::CSVReader<6> in(tangent_file);
		in.read_header(io::ignore_extra_column, "x", "y", "z", "tx", "ty", "tz");
		if (!in.has_column("x") || !in.has_column("y") || !in.has_column("z"))
			std::throw_with_nested(GRBF_Exceptions::missing_coords_in_file);
		if (!in.has_column("tx") || !in.has_column("ty") || !in.has_column("tz"))
			std::throw_with_nested(GRBF_Exceptions::missing_orientation_in_file);
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

std::vector<Inequality> build_inequality_constraints(const char *inequality_file)
{
	std::vector<Inequality> inequality_constraints;

	try
	{
		io::CSVReader<4> in(inequality_file);
		in.read_header(io::ignore_extra_column, "x", "y", "z", "level");
		if (!in.has_column("x") || !in.has_column("y") || !in.has_column("z"))
			std::throw_with_nested(GRBF_Exceptions::missing_coords_in_file);
		if (!in.has_column("level"))
			std::throw_with_nested(GRBF_Exceptions::missing_level_info_in_file);
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
