#ifndef continuous_property_h
#define continuous_property_h

#include <surfe_lib_module.h> // macro for importing / exporting dll

#include <modeling_methods.h>

class SURFE_LIB_EXPORT Continuous_Property : public GRBF_Modelling_Methods {
private:
	bool _get_polynomial_matrix_block(MatrixXd &poly_matrix);
	bool _insert_polynomial_matrix_blocks_in_interpolation_matrix(const MatrixXd &poly_matrix, MatrixXd &interpolation_matrix);
public:
	// Constructor/Destructor
	Continuous_Property(const model_parameters& m_p, const Basic_input& basic_i);
	Continuous_Property(){};
	// Methods
	Polynomial_Basis *create_polynomial_basis(const int &poly_order);
	bool get_interpolation_matrix(MatrixXd &interpolation_matrix);
	bool get_equality_values(VectorXd &equality_values);
	void eval_scalar_interpolant_at_point(Point &p);
	void eval_vector_interpolant_at_point(Point &p);
	bool get_method_parameters();
	bool process_input_data();
	bool setup_system_solver();
	bool get_minimial_and_excluded_input(Basic_input &greedy_input, Basic_input &excluded_input) {return true;}
	bool measure_residuals(Basic_input &input);
	bool append_greedy_input(Basic_input &input);
	bool convert_modified_kernel_to_rbf_kernel() { return true; } // To IMPLEMENT
	GRBF_Modelling_Methods *clone() { return new Continuous_Property(*this); }
	// Attributes
	Polynomial_Basis *p_basis;
};

#endif
