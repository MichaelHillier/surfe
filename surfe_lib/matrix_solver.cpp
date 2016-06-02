#include <matrix_solver.h>
#include <math_methods.h>

#include <iostream>
#include <iomanip>
#include <fstream>

bool Linear_LU_decomposition::solve()
{
	// create placeholder matrix/vector for LU decomposition algorithm
	std::vector< std::vector < double > > matrix = _interpolation_matrix;

	if (_constraint_values.size() != _interpolation_matrix.size()) return false;
	std::vector<double> constraint_values = _constraint_values;

	//if (!validate_matrix_systems()) return false;

	if (!Math_methods::solve_sym_linear_system_via_LU_dcmp<double>(matrix,constraint_values)) return false;
	// assign weight vector to solution vector
	weights = constraint_values;
	return true;
}

bool Linear_LU_decomposition::validate_matrix_systems()
{
	// check if matrix is square
	int n = (int)_interpolation_matrix.size();
	for (int j = 0; j < n; j++  ){
		if ( (int)_interpolation_matrix[j].size() != n ) return false;
	}

	// check is matrix is symmetric 
	std::vector< std::vector < double > > transpose_matrix = Math_methods::matrix_transpose(_interpolation_matrix);
	for (int j = 0; j < n; j++ ){
		for (int k = 0; k < n; k++ ){
			if ( abs(_interpolation_matrix[j][k] - transpose_matrix[j][k]) > 1e-6 ) return false; // found e.g. in lajaunie where diff was 1e-8. SO changed to 1e-6 from 1e-9. Should re-look at this situation
		}
	}

	return true;
}

bool Linear_LU_decomposition::check_solution()
{
	std::vector<double> result;
	result.resize((int)weights.size());

	Math_methods::matrix_vector_multiply(_interpolation_matrix,weights,result);
	
	std::vector<double> residual;
	residual.resize((int)weights.size());
	for (int j = 0; j < (int)result.size();j++ ){
		residual[j] = _constraint_values[j] - result[j];
	}

	return true;
}

std::vector < std::vector < mpf_class > > Quadratic_Predictor_Corrector::_convert_double_matrix_2_mpf( const std::vector < std::vector < double > > &matrix )
{
	int nrows = (int)matrix.size();
	int ncols = (int)matrix[0].size();

	std::vector < std::vector < mpf_class > > mpf_matrix = Math_methods::make_std_matrix<mpf_class>(nrows,ncols);

	for (int j = 0; j < nrows; j++ ){
		for (int k = 0; k < ncols; k++ ){
			mpf_matrix[j][k] = matrix[j][k];
		}
	}

	return mpf_matrix;
}

std::vector < mpf_class > Quadratic_Predictor_Corrector::_convert_double_vector_2_mpf( const std::vector <double > &vector )
{
	std::vector <mpf_class>  mpf_vector;

	for (int j = 0; j < (int)vector.size(); j++ ) mpf_vector.push_back(vector[j]);
	return mpf_vector;
}

std::vector < double> Quadratic_Predictor_Corrector::_convert_mpf_vector_2_double( const std::vector < mpf_class > &vector )
{
	std::vector <double> v;
	v.resize((int)vector.size());

	for (int j = 0; j < (int)v.size(); j++ ) v[j] = vector[j].get_d();
	return v;
}

std::vector < std::vector <mpf_class > > Quadratic_Predictor_Corrector::_get_hessian_matrix( const std::vector < std::vector < mpf_class > > &matrix )
{
	int nrows = (int)matrix.size();
	int ncols = (int)matrix[0].size();
	std::vector < std::vector <mpf_class> > hessian = Math_methods::make_std_matrix<mpf_class>(nrows,ncols);

	for (int j = 0; j < nrows; j++ ){
		for (int k = 0; k < ncols; k++ ){
			hessian[j][k] = 2.0*matrix[j][k];
		}
	}

	return hessian;
}

bool Quadratic_Predictor_Corrector::solve()
{
	std::vector< mpf_class > fvalues;
	fvalues.resize((int)_hessian_matrix.size());

	//if (!validate_matrix_systems()) return false;

	if (!Math_methods::quadratic_solver(_hessian_matrix,_equality_matrix,_inequality_matrix,_equality_vector,_inequality_vector,fvalues)) return false;

	weights = _convert_mpf_vector_2_double(fvalues);

	return true;
}

bool Quadratic_Predictor_Corrector::validate_matrix_systems()
{
	// determine if interpolation matrix is positive definite
	std::vector< std::vector < double > > Eigenvector_matrix = Math_methods::make_std_matrix<double>((int)_interpolation_matrixD.size(),(int)_interpolation_matrixD.size());
	std::vector< double > Eigenvalues;
	Eigenvalues.resize((int)_interpolation_matrixD.size());
	Math_methods::eigenanalysis(_interpolation_matrixD,Eigenvector_matrix,Eigenvalues);
	for (int j = 0; j < (int)Eigenvalues.size(); j++ ){
		if ( Eigenvalues[j] < 0 ) return false; // matrix is not positive definite
	}
	return true;
}
