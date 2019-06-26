#ifndef SURFE_API
#define SURFE_API

#include <surfe_lib_module.h>  // macro for importing / exporting dll

#include <modeling_methods.h>
#include <continuous_property.h>
#include <lajaunie.h>
#include <single_surface.h>
#include <stratigraphic_surfaces.h>
#include <vector_field.h>

class SURFE_LIB_EXPORT Surfe_API {
private:
	// members
	GRBF_Modelling_Methods *method_;

	bool have_interpolant_;
	bool parameters_changed_;
	bool constraints_changed_;

	// methods
	GRBF_Modelling_Methods* get_method_from_parameters(const Parameters& params);

public:
	Surfe_API(const int &modelling_method);
	Surfe_API(const Parameters& params);
	
	// Manually (no input files) adding constraints to interpolant methods
	void AddInterfaceConstraint(
		const double &x, const double &y, const double &z,
		const double &level
	);
	void AddPlanarConstraintwNormal(
		const double &x, const double &y, const double &z,
		const double &nx, const double &ny, const double &nz
	);
	void AddPlanarConstraintwStrikeDipPolarity(
		const double &x, const double &y, const double &z,
		const double &strike, const double &dip, const int &polarity
	);
	void AddPlanarConstraintwAzimuthDipPolarity(
		const double &x, const double &y, const double &z,
		const double &azimuth, const double &dip, const int &polarity
	);
	void AddTangentConstraint(
		const double &x, const double &y, const double &z,
		const double &tx, const double &ty, const double &tz
	);
	void AddInequalityConstraint(
		const double &x, const double &y, const double &z,
		const double &level
	);
	void ComputeInterpolant();
	void SetRegressionSmoothing(const bool &use_regression_smoothing, const double &amount);
	void SetGreedyAlgorithm(const bool &use_greedy, const double &interface_uncertainty, const double &angular_uncertainty);
	void SetRestrictedRange(const bool &use_restricted_range, const double &interface_uncertainty = 0, const double &angular_uncertainty = 0);
	void SetRBFKernel(const Parameter_Types::RBF &rbf);
	void SetRBFKernel(const char *rbf_name);
	void SetRBFShapeParameter(const double &shape_param);
	void SetPolynomialOrder(const int &poly_order);
	void SetGlobalAnisotropy(const bool &g_anisotropy);
	double EvaluateInterpolantAtPoint(
		const double &x, const double &y, const double &z
	);
	Vector3d EvaluateVectorInterpolantAtPoint(
		const double &x, const double &y, const double &z
	);

	SpatialParameters GetDataBoundsAndResolution();

	// Array of interface reference points: 1 per interface. 
	// n x 3 matrix, n = number of interfaces 
	MatrixXd GetInterfaceReferencePoints();

	// Array of interface constraints
	// n x 4 matrix, n = number of interface points
	// columns: x, y, z, level
	// level: structural levels (interface codes)
	MatrixXd GetInterfaceConstraints();
	void SetInterfaceConstraints(const MatrixXd &interface_constraints);

	// Array of planar constraints
	// n x 6 matrix, n = number of planar points
	// columns: x, y, z, nx, ny, nz
	// nx, ny, nz : components of the normal vector
	MatrixXd GetPlanarConstraints();
	void SetPlanarConstraints(const MatrixXd &planar_constraints);

	// Array of tangent constraints
	// n x 6 matrix, n = number of tangent points
	// columns: x, y, z, vx, vy, vz
	// vx, vy, vz : components of the tangent vector
	MatrixXd GetTangentConstraints();
	void SetTangentConstraints(const MatrixXd &tangent_constraints);

	// Array of inequality constraints
	// n x 4 matrix, n = number of inequality points
	// columns: x, y, z, level
	// level: structural level (compatible with interface level codes)
	MatrixXd GetInequalityConstraints();
	void SetInequalityConstraints(const MatrixXd &inequality_constraints);

	int GetNumberOfInterfaces();

	bool InterpolantComputed() const { return have_interpolant_; }

};

#endif // SURFE_API
