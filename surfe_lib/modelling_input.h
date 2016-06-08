#ifndef _modelling_input_h
#define _modelling_input_h

#include <surfe_lib_module.h>

#include <vector>

class SURFE_LIB_EXPORT Point{
private:
	double _x;
	double _y;
	double _z;
	double _c;
	double _scalar_field;
	double _field_normal[3];
public:
	Point(const double& x_coord,const double& y_coord, const double &z_coord, const double &c_coord = 0)
		: _x(x_coord), _y(y_coord), _z(z_coord), _c(c_coord)
	{
		_scalar_field = NULL;
		for (int j = 0; j < 3; j++ ) _field_normal[j] = NULL; 
	}
	double x() const { return _x; }
	double y() const { return _y; }
	double z() const { return _z; }
	double c() const { return _c; }
	void set_x(const double &x_coord) { _x = x_coord; }
	void set_y(const double &y_coord) { _y = y_coord; }
	void set_z(const double &z_coord) { _z = z_coord; }
	void set_c(const double &c_coord) { _c = c_coord; }
	double scalar_field() const { return _scalar_field; }
	void set_scalar_field(const double &scalar_field_value) { _scalar_field = scalar_field_value; }
	void set_vector_field(const double &nx, const double &ny, const double &nz) { _field_normal[0] = nx; _field_normal[1] = ny; _field_normal[2] = nz; }
	double nx_interp() const { return _field_normal[0]; }
	double ny_interp() const { return _field_normal[1]; }
	double nz_interp() const { return _field_normal[2]; }
};

class SURFE_LIB_EXPORT Evaluation_Point : public Point{
public:
	Evaluation_Point(const double &x_coord,
		const double &y_coord,
		const double &z_coord,
		const double &c_coord = NULL)
		: Point(x_coord,y_coord,z_coord,c_coord)
	{ 
	} 
};

class SURFE_LIB_EXPORT Interface : public Point {
private:
	double _level;
	double _residual;
public:
	Interface(const double &x_coord,
		const double &y_coord,
		const double &z_coord,
		const double &lvl,
		const double &c_coord = NULL)
		: Point(x_coord,y_coord,z_coord,c_coord), _level(lvl) { _residual = 0.0; } 
	double level() const { return _level; }
	double residual() const { return _residual; }
	void setResidual(const double &res) { _residual = res; }
	void setLevel(const double &v) { _level = v; }
};

class SURFE_LIB_EXPORT Inequality : public Point {
private:
	double _inequality_level;
	bool _residual;
public:
	Inequality(const double &x_coord,
		const double &y_coord,
		const double &z_coord,
		const double &lvl,
		const double &c_coord = NULL)
		: Point(x_coord,y_coord,z_coord,c_coord), _inequality_level(lvl) { _residual = true; } 
	double level() const { return _inequality_level; }
	bool residual() const { return _residual; }
	void setResidual(const bool &res) { _residual = res; }
};

class SURFE_LIB_EXPORT Planar : public Point {
private:
	double _dip;
	double _strike;
	int _polarity;
	double _normal[3];
	double _residual;
	double _normal_error[3][2];
	bool _compute_strike_dip_polarity_from_normal();
	bool _compute_normal_from_strike_dip_polarity();
public:
	Planar (const double &x_coord,
		const double &y_coord,
		const double &z_coord,
		const double &nx,
		const double &ny,
		const double &nz,
		const double &c_coord = NULL)
		: Point(x_coord,y_coord,z_coord,c_coord)
	{
		_normal[0] = nx;
		_normal[1] = ny;
		_normal[2] = nz;
		_residual = 0.0;
		// compute strike, dip, and polarity
		_compute_strike_dip_polarity_from_normal();
	}
	Planar (const double &x_coord,
		const double &y_coord,
		const double &z_coord,
		const double &dip,
		const double &strike,
		const int &polarity,
		const double &c_coord = NULL)
		: Point(x_coord,y_coord,z_coord,c_coord), _dip(dip), _strike(strike), _polarity(polarity)
	{
		// compute normal
		_compute_normal_from_strike_dip_polarity();
		// residual default
		_residual = 0.0;
	}
	bool getDipVector(double (&vector)[3]);
	bool getStrikeVector(double (&vector)[3]);
	double dip() const { return _dip; }
	double strike() const { return _strike; }
	int polarity() const { return _polarity; }
	double nx() const { return _normal[0]; }
	double ny() const { return _normal[1]; }
	double nz() const { return _normal[2]; }
	bool getNormalError(double (&matrix)[3][2]);
	double residual() const { return _residual; }
	void setResidual(const double &res) { _residual = res; }
}; 

class SURFE_LIB_EXPORT Tangent : public Point {
private:
	double _tangent[3];
	double _residual;
public:
	Tangent (const double &x_coord,
		const double &y_coord,
		const double &z_coord,
		const double &tx,
		const double &ty,
		const double &tz,
		const double &c_coord = NULL)
		: Point(x_coord,y_coord,z_coord,c_coord)
	{
		_tangent[0] = tx;
		_tangent[1] = ty;
		_tangent[2] = tz;
		_residual = 0.0;
	}
	double tx() const { return _tangent[0]; }
	double ty() const { return _tangent[1]; }
	double tz() const { return _tangent[2]; }
	double residual() const { return _residual; }
	void setResidual(const double &res) { _residual = res; }
}; 


class SURFE_LIB_EXPORT Basic_input{
private:
	// interface
	void _get_distinct_interface_iso_values();
	void _get_interface_points();
	bool _interface_points_are_coplanar(){ return true; } // Not implemented yet. should be tested when 2nd order polynomials are used. Also when unisolvent points are used this should be called.

public:
	Basic_input()
	{
		inequality = new std::vector<Inequality>();
		interface = new std::vector<Interface>();
		planar = new std::vector<Planar>();
		tangent = new std::vector<Tangent>();

		evaluation_pts = new std::vector<Evaluation_Point>();

		interface_iso_values = new std::vector<double>();
		interface_point_lists = new std::vector< std::vector < Interface > >();
		interface_test_points = new std::vector< Interface >();
	}
	~Basic_input()
	{
		delete inequality;
		delete interface;
		delete planar;
		delete tangent;
		delete evaluation_pts;
		delete interface_iso_values;
		delete interface_point_lists;
		delete interface_test_points;
	}
	// input data 
	std::vector< Inequality > *inequality;
	std::vector< Interface > *interface;
	std::vector< Planar > *planar;
	std::vector< Tangent > *tangent;

	// evaluation sites in grid
	std::vector< Evaluation_Point > *evaluation_pts;

	// for interface data
	std::vector < double > *interface_iso_values;
	std::vector < std::vector < Interface > > *interface_point_lists;
	std::vector < Interface > *interface_test_points;
	bool get_interface_data();

};

inline double distance_btw_pts(const Point &p1, const Point &p2);
int nearest_neighbour_index(const Point &p, const std::vector < Point > &pts);
int furtherest_neighbour_index(const Point &p, const std::vector < Point > &pts);
int furtherest_neighbour_index(const std::vector < Point > &pts1, const std::vector < Point > &pts2);
void calculate_bounds(const std::vector< Point > &pts, double (&bounds)[6]);
std::vector<int> get_extremal_point_data_indices_from_points(const std::vector< Point > &pts);
bool is_index_in_list(const int &index, const std::vector < int > &list);
bool find_fill_distance(const Basic_input &input, double &fill_distance);

#endif