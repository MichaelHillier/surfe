#ifndef lajaunie_h
#define lajaunie_h

#include <surfe_lib_module.h> // macro for importing / exporting dll

#include "modeling_methods.h"

class SURFE_LIB_EXPORT Lajaunie_Approach : public GRBF_Modelling_Methods {
private:
	int _n_increment_pair;
	bool _get_polynomial_matrix_block(MatrixXd &poly_matrix);
	bool _insert_polynomial_matrix_blocks_in_interpolation_matrix(const MatrixXd &poly_matrix, MatrixXd &interpolation_matrix);
	std::vector < std::vector < Interface > > *_increment_pairs;
	bool _get_increment_pairs();
public:
	// Constructor/Destructor
	Lajaunie_Approach(const model_parameters& m_p, const Basic_input& basic_i);
	~Lajaunie_Approach(){};
	// Methods
	Polynomial_Basis *create_polynomial_basis(const int &poly_order);
	bool get_interpolation_matrix(MatrixXd &interpolation_matrix);
	bool get_equality_values(VectorXd &equality_values);
	bool get_inequality_values(VectorXd &b, VectorXd &r);
	void eval_scalar_interpolant_at_point(Point &p);
	void eval_vector_interpolant_at_point(Point &p);
	bool get_method_parameters();
	bool process_input_data();
	bool setup_system_solver();
	bool get_minimial_and_excluded_input(Basic_input &greedy_input, Basic_input &excluded_input);
	bool measure_residuals(Basic_input &input);
	bool append_greedy_input(Basic_input &input);
	bool convert_modified_kernel_to_rbf_kernel();
	GRBF_Modelling_Methods *clone() { return new Lajaunie_Approach(*this); }
	// Attributes
	Polynomial_Basis *p_basis;
};

#endif