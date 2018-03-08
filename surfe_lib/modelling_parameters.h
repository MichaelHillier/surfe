#ifndef modelling_parameters_h
#define modelling_parameters_h

#include <surfe_lib_module.h>

#define D2R 0.01745329251994329576923690768489 // degrees to radians conversion factor
#define R2D 57.295779513082320876798154814105  // radians to degrees conversion factor

struct SURFE_LIB_EXPORT Parameter_Types{
	enum DWRT {PT1,PT2};
	enum SecondDerivatives {DXDX,DXDY,DXDZ,DYDX,DYDY,DYDZ,DZDX,DZDY,DZDZ};
	enum FirstDerivatives {DX,DY,DZ};
	enum RBF {Cubic,Gaussian,MQ,IMQ,TPS,R};
	enum SolverType {Linear,Quadratic};
	enum ModelType {Single_surface,Lajaunie_approach,Stratigraphic_horizons,Continuous_property,Vector_field};
	enum AXIS {Xaxis,Yaxis,Zaxis};
};

struct SURFE_LIB_EXPORT model_parameters{
	////////////////////////////////
	//        UI parameters       //
	////////////////////////////////

	// model type
	Parameter_Types::ModelType model_type;
	double min_stratigraphic_thickness;
	// interface input
	bool use_interface_data;
	bool use_planar_data;
	bool use_tangent;
	bool use_inequality;
	// basis parameters
	Parameter_Types::RBF basis_type;
	double shape_parameter;
	int polynomial_order;

	bool advanced_parameters;
	bool model_global_anisotropy;
	bool use_greedy;
	bool use_restricted_range;
	double interface_uncertainty;
	double angular_uncertainty;
	bool use_regression_smoothing;
	double smoothing_amount;
	int nearest_neighbours;
	int idw_power;

	// initialization ...
	model_parameters() : model_type(Parameter_Types::Single_surface), min_stratigraphic_thickness(0),
		use_interface_data(true), use_planar_data(true), use_tangent(false), use_inequality(false),
		basis_type(Parameter_Types::Cubic), shape_parameter(100), polynomial_order(1),
		advanced_parameters(false), model_global_anisotropy(false), use_greedy(false), use_restricted_range(false), use_regression_smoothing(false), interface_uncertainty(0), angular_uncertainty(0), nearest_neighbours(2),idw_power(4) {}
};

struct SURFE_LIB_EXPORT basic_parameters{
	// number of constraints, for each basic constraint type
	unsigned int n_interface;
	unsigned int n_planar;
	unsigned int n_inequality;
	unsigned int n_tangent;
	unsigned int n_constraints;
	unsigned int n_equality;
	// basis function parameters
	bool modified_basis;
	// polynomial parameters
	bool poly_term;
	unsigned int n_poly_terms;
	// type of problem
	Parameter_Types::SolverType problem_type;
	//bool greedy; // greedy algorithm
	bool restricted_range; // restricted range constraints - bounded inequalities
	// initialization 
	basic_parameters() : n_interface(0), n_planar(0), n_inequality(0), n_tangent(0), n_constraints(0), n_equality(0),
		modified_basis(false), poly_term(true), n_poly_terms(4), problem_type(Parameter_Types::Linear),restricted_range(false){}
};

#endif