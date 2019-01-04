#ifndef _anisotropy_input_h
#define _anisotropy_input_h

#include <surfe_lib_module.h>
#include <modelling_parameters.h>
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
	double _dip_direction; // for euro folk
	Vector3d _normal;
	Vector3d _dipvector;
	Vector3d _strikevector;
	bool _compute_strike_dip_polarity_from_normal();
	bool _compute_a_normal_from_strike_dip();
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
		const double &strike,
		bool euro_convention = false)
	{
		_x = x_coord;
		_y = y_coord;
		_z = z_coord;
		_dip = dip;
		if (euro_convention) // in this case the 'strike' is not the strike but the dip direction
		{
			if (strike >= 90) _strike = strike - 90.0;
			else _strike = strike + 270.0;
		}
		else _strike = strike;
		// compute normal
		_compute_a_normal_from_strike_dip();
		_set_eigensystem();
	}
	bool getDipVector(double(&vector)[3]);
	bool getStrikeVector(double(&vector)[3]);
	double dip() const { return _dip; }
	double strike() const { return _strike; }
	double nx() const { return _normal[0]; }
	double ny() const { return _normal[1]; }
	double nz() const { return _normal[2]; }
	double x() const { return _x; }
	double y() const { return _y; }
	double z() const { return _z; }
	void setNormal(const double &nx, const double &ny, const double &nz) { _normal[0] = nx; _normal[1] = ny; _normal[2] = nz; }
	Matrix3d eigenvectors; // principal directions of anisotropy
	Vector3d eigenvalues;
	Matrix3d U;
	Matrix3d D;
	Matrix3d S;
	Matrix3d Tensor;
	Matrix3d Transform;
};

class SURFE_LIB_EXPORT TensorEvaluationPoints{
private:
	double _x;
	double _y;
	double _z;
public:
	TensorEvaluationPoints(const double &x_coord,
		const double &y_coord,
		const double &z_coord)
	{
		_x = x_coord;
		_y = y_coord;
		_z = z_coord;
	}
 	double x() const { return _x; }
 	double y() const { return _y; }
 	double z() const { return _z; }
	Matrix3d eigenvectors; // principal directions of anisotropy
	Vector3d eigenvalues;
	Vector3d GetDipVector();
	Vector3d GetStrikeVector();
	Vector3d GetNormalVector();
};

class SURFE_LIB_EXPORT TensorInput{
private:
	// Attributes
	std::vector < std::vector < Orientation > > neighbourhoods;
	void _sort_eigensystem(Matrix3d &evectors, Vector3d &evalues);
	bool _get_neighbourhoods(const int &n_neighbors);
	void _get_local_anisotropy_from_neighbourhoods(std::vector < std::vector < Orientation > > &nh);
public:
	TensorInput()
	{
		orientation = new std::vector<Orientation>;
		evalpts = new std::vector<TensorEvaluationPoints>;
	}
	~TensorInput(){}
	bool GetLocalAnisotropy(const model_parameters &parameters);
	// input data 
	std::vector< Orientation > *orientation;
	std::vector< TensorEvaluationPoints > *evalpts;
};
#endif