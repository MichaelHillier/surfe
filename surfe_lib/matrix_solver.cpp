#include <matrix_solver.h>
#include <math_methods.h>

#include <iostream>
#include <iomanip>
#include <fstream>

bool Linear_LU_decomposition::solve()
{

	if (_constraint_values.rows() != _interpolation_matrix.rows()) return false;

	weights = _interpolation_matrix.partialPivLu().solve(_constraint_values);

	if ( !weights.allFinite()) return false; 

	return true;
}

bool Linear_LU_decomposition::validate_matrix_systems()
{
	// check if there are any NaN or INF values
	if (_interpolation_matrix.allFinite()) return true;
	else return false;

	// the other potential problem could be the input has un-initialized values in it.
}

bool Linear_LU_decomposition::check_solution()
{
	double relative_error = (_interpolation_matrix*weights - _constraint_values).norm() / _constraint_values.norm(); // norm() is L2 norm
	std::cout << "The relative error is:\n" << relative_error << std::endl;
	return true;
}

Matrix <mpf_class, Dynamic, Dynamic> Quadratic_Predictor_Corrector::_convert_double_matrix_2_mpf(const MatrixXd &matrix)
{
	int nrows = (int)matrix.rows();
	int ncols = (int)matrix.cols();

	Matrix <mpf_class, Dynamic, Dynamic> mpf_matrix(nrows, ncols);

	for (int j = 0; j < nrows; j++ ){
		for (int k = 0; k < ncols; k++ ){
			mpf_matrix(j,k) = matrix(j,k);
		}
	}

	return mpf_matrix;
}

Matrix <mpf_class, Dynamic, 1> Quadratic_Predictor_Corrector::_convert_double_vector_2_mpf(const VectorXd &vector)
{
	int nrows = (int)vector.rows();
	Matrix <mpf_class, Dynamic, 1>  mpf_vector(nrows);

	for (int j = 0; j < nrows ; j++ ) mpf_vector(j) = vector(j);
	return mpf_vector;
}

VectorXd Quadratic_Predictor_Corrector::_convert_mpf_vector_2_double(const Matrix <mpf_class, Dynamic, 1> &vector)
{
	int nrows = (int)vector.rows();
	VectorXd v(nrows);

	for (int j = 0; j < nrows; j++ ) v(j) = vector(j).get_d();
	return v;
}

Matrix <mpf_class, Dynamic, Dynamic> Quadratic_Predictor_Corrector::_get_hessian_matrix(const Matrix <mpf_class, Dynamic, Dynamic> &matrix)
{
	int nrows = (int)matrix.rows();
	int ncols = (int)matrix.cols();
	Matrix <mpf_class, Dynamic, Dynamic> hessian(nrows, ncols);

	hessian = 2.0*matrix;

	return hessian;
}

bool Quadratic_Predictor_Corrector::solve()
{
	int n = (int)_hessian_matrix.rows();

	Matrix <mpf_class, Dynamic, 1> fvalues(n);

	//if (!validate_matrix_systems()) return false;

	if (!Math_methods::quadratic_solver(_hessian_matrix,_equality_matrix,_inequality_matrix,_equality_vector,_inequality_vector,fvalues)) return false;

	weights = _convert_mpf_vector_2_double(fvalues);

	return true;
}

bool Quadratic_Predictor_Corrector::validate_matrix_systems()
{
	if (!_interpolation_matrixD.allFinite()) return false;

	LLT<MatrixXd> lltofMatrix(_interpolation_matrixD);
	if (lltofMatrix.info() == NumericalIssue) return false;

	return true;
}
