#ifndef vector_field_h
#define vector_field_h

#include <surfe_lib_module.h> // macro for importing / exporting dll

#include <modeling_methods.h>


class SURFE_LIB_EXPORT Vector_Field : public Greedy_Method {
public:
	// Constructor/Destructor
	Vector_Field(const model_parameters& m_p, const Basic_input& basic_i) 
	{
		m_parameters = m_p; 
		b_input = basic_i; 

		_avg_nn_dist_ie = 0;
		_avg_nn_dist_p = 0;
		_avg_nn_dist_itr = 0;
		_avg_nn_dist_t = 0;
		_iteration = 0;
	}
	~Vector_Field(){};
	// Methods
	bool get_interpolation_matrix(std::vector< std::vector <double> > &interpolation_matrix);
	bool get_equality_values(std::vector<double> &equality_values);
	void eval_scalar_interpolant_at_point(Point &p);
	void eval_vector_interpolant_at_point(Point &p);
	bool get_method_parameters();
	bool process_input_data();
	bool setup_system_solver();
	bool get_minimial_and_excluded_input(Basic_input &greedy_input, Basic_input &excluded_input) { return true; } // TO implement
	bool measure_residuals(Basic_input &input) { return true; } // TO implement
	bool append_greedy_input(Basic_input &input) { return true; } // TO implement
	Greedy_Method *clone() { return new Vector_Field(*this); }
};

#endif