#ifndef matrix_solver_h
#define matrix_solver_h

#include <gmpxx.h>

#include <vector>
#include <Eigen/Core>

using namespace Eigen;

class System_Solver {
public:
	System_Solver(){}
	virtual ~System_Solver(){}
	VectorXd weights;
	virtual bool solve() = 0;
 	virtual bool validate_matrix_systems() = 0;
 	//virtual bool validate_constraint_vectors() = 0;
 	//virtual bool check_solution() = 0;
};

class Linear_LU_decomposition : public System_Solver {
private:
	MatrixXd _interpolation_matrix;
	VectorXd _constraint_values;
public:
	Linear_LU_decomposition(const MatrixXd &matrix, const VectorXd &vector)
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
	MatrixXd _interpolation_matrixD;
	Matrix <mpf_class, Dynamic, Dynamic> _hessian_matrix;
	Matrix <mpf_class, Dynamic, Dynamic> _interpolation_matrix;
	Matrix <mpf_class, Dynamic, Dynamic> _equality_matrix;
	Matrix <mpf_class, Dynamic, Dynamic> _inequality_matrix;
	Matrix <mpf_class, Dynamic, 1> _equality_vector;
	Matrix <mpf_class, Dynamic, 1> _inequality_vector;
	Matrix <mpf_class, Dynamic, Dynamic> _convert_double_matrix_2_mpf(const MatrixXd &matrix);
	Matrix <mpf_class, Dynamic, 1> _convert_double_vector_2_mpf(const VectorXd &vector);
	VectorXd _convert_mpf_vector_2_double(const Matrix <mpf_class, Dynamic, 1> &vector);
	Matrix <mpf_class, Dynamic, Dynamic> _get_hessian_matrix(const Matrix <mpf_class, Dynamic, Dynamic> &matrix);
public:
	Quadratic_Predictor_Corrector(const MatrixXd &interpolation_matrix,
		                          const MatrixXd &equality_matrix,
								  const MatrixXd &inequality_matrix,
								  const VectorXd equality_vector,
								  const VectorXd inequality_vector)
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

