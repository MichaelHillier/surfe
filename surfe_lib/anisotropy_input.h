#ifndef _anisotropy_input_h
#define _anisotropy_input_h

#include <surfe_lib_module.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include <vector>

#define D2R 0.01745329251994329576923690768489 // degrees to radians conversion factor
#define R2D 57.295779513082320876798154814105  // radians to degrees conversion factor

using namespace Eigen;

class SURFE_LIB_EXPORT Orientation{
private:
	double _x;
	double _y;
	double _z;
	double _dip;
	double _strike;
	double _normal[3];
	Vector3d _normal;
	Vector3d _dipvector;
	Vector3d _strikevector;
	bool _compute_strike_dip_polarity_from_normal();
	bool _compute_normal_from_strike_dip_polarity();
	void _set_eigensystem();
public:
	Orientation(const double &x_coord,
		const double &y_coord,
		const double &z_coord,
		const double &nx,
		const double &ny,
		const double &nz)
	{
		_x = x_coord;
		_y = y_coord;
		_z = z_coord;
		_normal[0] = nx;
		_normal[1] = ny;
		_normal[2] = nz;
		// compute strike, dip, and polarity
		_compute_strike_dip_polarity_from_normal();
	}
	Orientation(const double &x_coord,
		const double &y_coord,
		const double &z_coord,
		const double &dip,
		const double &strike)
	{
		_x = x_coord;
		_y = y_coord;
		_z = z_coord;
		_dip = dip;
		_strike = strike;
		// compute normal
		_compute_normal_from_strike_dip_polarity();
		_set_eigensystem();
	}
	bool getDipVector(double(&vector)[3]);
	bool getStrikeVector(double(&vector)[3]);
	double dip() const { return _dip; }
	double strike() const { return _strike; }
	double nx() const { return _normal[0]; }
	double ny() const { return _normal[1]; }
	double nz() const { return _normal[2]; }
	void setNormal(const double &nx, const double &ny, const double &nz) { _normal[0] = nx; _normal[1] = ny; _normal[2] = nz; }
	Matrix3d eigenvectors; // principal directions of anisotropy
	Vector3d eigenvalues;
};

class SURFE_LIB_EXPORT local_anisotropy_input{
private:
	// Attributes
	Matrix3d _eigenvectors; // principal directions of anisotropy
	Vector3d _eigenvalues;
public:
	local_anisotropy_input()
	{
		orientation = new std::vector<Orientation>;
	}
	~local_anisotropy_input()
	{
		delete orientation;
	}
	// input data 
	std::vector< Orientation > *orientation;
	void getMeanTensor(Matrix3d &y,Matrix3d &d);

#endif