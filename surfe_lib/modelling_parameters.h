#ifndef base_definitions_h
#define base_definitions_h

#include <vector>

struct Parameter_Types{
	enum DWRT {PT1,PT2};
	enum SecondDerivatives {DXDX,DXDY,DXDZ,DYDX,DYDY,DYDZ,DZDX,DZDY,DZDZ};
	enum FirstDerivatives {DX,DY,DZ};
	enum RBF {Cubic,Gaussian,MQ,IMQ,TPS,R};
	enum SolverType {Linear,Quadratic};
	enum ModelType {Single_surface,Lajaunie_approach,Stratigraphic_horizons,Continuous_property};
	enum AXIS {Xaxis,Yaxis,Zaxis};
};

struct model_parameters{
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
	bool use_smoothing;
	double interface_slack; // interface slack
	double gradient_slack; // gradient slack

	// initialization ...
	model_parameters() : model_type(Parameter_Types::Single_surface), min_stratigraphic_thickness(0),
		use_interface_data(true), use_planar_data(true), use_tangent(false), use_inequality(false),
		basis_type(Parameter_Types::Cubic), shape_parameter(100), polynomial_order(1),
		advanced_parameters(false), model_global_anisotropy(false), use_smoothing(false), interface_slack(0), gradient_slack(0) {}
};

struct basic_parameters{
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
	// initialization 
	basic_parameters() : n_interface(0), n_planar(0), n_inequality(0), n_tangent(0), n_constraints(0), n_equality(0),
		modified_basis(false), poly_term(true), n_poly_terms(4), problem_type(Parameter_Types::Linear){}
};

#endif