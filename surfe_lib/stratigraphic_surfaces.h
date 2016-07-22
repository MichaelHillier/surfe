#ifndef stratigraphic_surfaces_h
#define stratigraphic_surfaces_h

#include <surfe_lib_module.h> // macro for importing / exporting dll

#include <modeling_methods.h>

class SURFE_LIB_EXPORT Stratigraphic_Surfaces : public GRBF_Modelling_Methods {
private:
	// Methods
	bool _get_increment_pairs();
	std::vector < std::vector < Point > > _get_lithostratigraphic_increment_pairs_for_inequality_point(const Inequality &ie_pt);
	double _get_closest_horizon_level_above_given_level(const double &given_level,const std::vector<double> &horizon_levels);
	double _get_closest_horizon_level_below_given_level(const double &given_level,const std::vector<double> &horizon_levels);
	// Attributes
	int _n_increment_pairs;
	int _n_sequenced_interface_pairs;
	int _n_sequenced_inequality_pairs;
	int _n_interface_pairs;
	std::vector < std::vector < Point > > *_increment_pairs;
public:
	// Constructor/Destructor
	Stratigraphic_Surfaces(const model_parameters& m_p, const Basic_input& basic_i);
	~Stratigraphic_Surfaces(){};
	// Methods
	bool get_interpolation_matrix(MatrixXd &interpolation_matrix);
	bool get_equality_values(VectorXd &equality_values);
	bool get_inequality_matrix(const MatrixXd &interpolation_matrix, MatrixXd &inequality_matrix);
	bool get_inequality_values(VectorXd &inequality_values);
	void eval_scalar_interpolant_at_point(Point &p);
	void eval_vector_interpolant_at_point(Point &p);
	bool get_method_parameters();
	bool process_input_data();
	bool setup_system_solver();
	bool get_minimial_and_excluded_input(Basic_input &greedy_input, Basic_input &excluded_input) { return true; } // TO implement
	bool measure_residuals(Basic_input &input) { return true; } // TO implement
	bool append_greedy_input(Basic_input &input) { return true; } // TO implement
	GRBF_Modelling_Methods *clone() { return new Stratigraphic_Surfaces(*this); }
	// Attributes
	Polynomial_Basis *p_basis;
};

#endif