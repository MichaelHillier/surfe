#ifndef matrix_solver_h
#define matrix_solver_h

#include <gmpxx.h>

#include <vector>

class System_Solver {
public:
	System_Solver(){}
	virtual ~System_Solver(){}
	std::vector<double> weights;
	virtual bool solve() = 0;
 	virtual bool validate_matrix_systems() = 0;
 	//virtual bool validate_constraint_vectors() = 0;
 	//virtual bool check_solution() = 0;
};

class Linear_LU_decomposition : public System_Solver {
private:
	std::vector < std::vector < double > > _interpolation_matrix;
	std::vector<double> _constraint_values;
public:
	Linear_LU_decomposition(const std::vector < std::vector < double > > &matrix, const std::vector<double> &vector)
	{ _interpolation_matrix = matrix; _constraint_values = vector; }
	Linear_LU_decomposition(){}
	virtual ~Linear_LU_decomposition() {}
 	bool solve();
 	bool validate_matrix_systems();
 	//bool validate_constraint_vectors();
 	bool check_solution();
};

class Quadratic_Predictor_Corrector : public System_Solver {
private:
	mpf_class _largest_element;
	std::vector < std::vector < double > > _interpolation_matrixD;
	std::vector < std::vector < mpf_class > > _hessian_matrix;
	std::vector < std::vector < mpf_class > > _interpolation_matrix;
	std::vector < std::vector < mpf_class > > _equality_matrix;
	std::vector < std::vector < mpf_class > > _inequality_matrix;
	std::vector < mpf_class > _equality_vector;
	std::vector < mpf_class > _inequality_vector;
	std::vector < std::vector < mpf_class > > _convert_double_matrix_2_mpf(const std::vector < std::vector < double > > &matrix);
	std::vector < mpf_class > _convert_double_vector_2_mpf(const std::vector <double > &vector);
	std::vector < double> _convert_mpf_vector_2_double(const std::vector < mpf_class > &vector);
	std::vector < std::vector <mpf_class > > _get_hessian_matrix(const std::vector < std::vector < mpf_class > > &matrix);
public:
	Quadratic_Predictor_Corrector(const std::vector < std::vector < double > > &interpolation_matrix,
			                        const std::vector < std::vector < double > > &equality_matrix,
									const std::vector < std::vector < double > > &inequality_matrix,
									const std::vector < double > equality_vector,
									const std::vector < double > inequality_vector)
	{
		_interpolation_matrixD = interpolation_matrix;
		_interpolation_matrix = _convert_double_matrix_2_mpf(interpolation_matrix);
		_equality_matrix = _convert_double_matrix_2_mpf(equality_matrix);
		_inequality_matrix = _convert_double_matrix_2_mpf(inequality_matrix);
		_equality_vector = _convert_double_vector_2_mpf(equality_vector);
		_inequality_vector = _convert_double_vector_2_mpf(inequality_vector);

		_hessian_matrix = _get_hessian_matrix(_interpolation_matrix);
	}
	bool solve();
	bool validate_matrix_systems();

};

#endif

