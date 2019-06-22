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
	InputParameters input_;

	bool have_interpolant_;
	bool evaluation_completed_;
	bool have_method_;
	bool parameters_changed_;
	bool constraint_files_changed_;
	bool constraints_changed_;

	// methods
	GRBF_Modelling_Methods* get_method(const Parameters& params);
	void build_constraints_from_input_files();
public:
	Surfe_API();
	Surfe_API(const Parameters& params);
	void LoadConstraintsFromFiles();
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
		const double &strike, const double &dip, const double &polarity
	);
	void AddPlanarConstraintwAzimuthDipPolarity(
		const double &x, const double &y, const double &z,
		const double &azimuth, const double &dip, const double &polarity
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
	void SetModellingMode(const int &mode);
	void SetRegressionSmoothing(const bool &use_regression_smoothing, const double &amount = 0);
	void SetGreedyAlgorithm(const bool &use_greedy, const double &interface_uncertainty = 0, const double &angular_uncertainty = 0);
	void SetRestrictedRange(const bool &use_restricted_range, const double &interface_uncertainty = 0, const double &angular_uncertainty = 0);
	void SetRBFKernel(const Parameter_Types::RBF &rbf);
	void SetRBFKernel(const char *rbf_name);
	void SetRBFShapeParameter(const double &shape_param);
	void SetPolynomialOrder(const int &poly_order);
	void SetGlobalAnisotropy(const bool &g_anisotropy);
	void SetInterfaceUncertainty(const double &interface_uncertainty);
	void SetAngularUncertainty(const double &angular_uncertainty);
	void SetInterfaceDataFile(const char *interfaceFile);
	void SetPlanarDataFile(const char *planarFile);
	void SetTangentDataFile(const char *tangentFile);
	void SetInequalityDataFile(const char *inequalityFile);
	double EvaluateInterpolantAtPoint(
		const double &x, const double &y, const double &z
	);
	double *EvaluateVectorInterpolantAtPoint(
		const double &x, const double &y, const double &z
	); // client responsible for deleting dynamically allocated array vector[3]
	SpatialParameters GetDataBoundsAndResolution();
	int GetNumberOfInterfaces();
	double GetScalarFieldValueAtInterfaceX(const int &i);
	GRBF_Modelling_Methods *GetMethod() { return method_; }
};

#endif // SURFE_API
