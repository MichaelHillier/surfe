// SURFace Estimator(SURFE) - Terms and Conditions of Use

// Unless otherwise noted, computer program source code of the SURFace Estimator(SURFE)
// is covered under Crown Copyright, Government of Canada, and is distributed under the MIT License.

// The Canada wordmark and related graphics associated with this distribution are protected under 
// trademark law and copyright law.No permission is granted to use them outside the parameters of 
// the Government of Canada's corporate identity program. For more information, see
// http://www.tbs-sct.gc.ca/fip-pcim/index-eng.asp

// Copyright title to all 3rd party software distributed with the SURFace Estimator(SURFE) is held 
// by the respective copyright holders as noted in those files.Users are asked to read the 3rd Party
// Licenses referenced with those assets.

// MIT License

// Copyright(c) 2017 Government of Canada

// Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
// and associated documentation files(the "Software"), to deal in the Software without restriction,
// including without limitation the rights to use, copy, modify, merge, publish, distribute, 
// sublicense, and / or sell copies of the Software, and to permit persons to whom the Software is 
// furnished to do so, subject to the following conditions :

// The above copyright notice and this permission notice shall be included in all copies or
// substantial portions of the Software.

//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
// NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef matrix_solver_h
#define matrix_solver_h

#include <gmpxx.h>

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <Eigen/Core>

using namespace Eigen;
using namespace std;

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
	{ 
		_interpolation_matrix = matrix; 
		_constraint_values = vector; 
// 		std::ofstream file1("interpM.txt");
// 		std::ofstream file2("consV.txt");
// 		if (file1)
// 		{
// 			file1 << matrix <<"\n";
// 			file1.close();
// 		}
// 		if (file2)
// 		{
// 			file2 << vector <<"\n";
// 			file2.close();
// 		}
	}
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
	MatrixXd _hessian_matrixD;
	MatrixXd _equality_matrixD;
	MatrixXd _inequality_matrixD;
	VectorXd _equality_vectorD;
	VectorXd _inequality_vectorD;
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
		_hessian_matrixD = 2.0*interpolation_matrix;
		_equality_matrixD = equality_matrix;
		_inequality_matrixD = inequality_matrix;
		_equality_vectorD = equality_vector;
		_inequality_vectorD = inequality_vector;
// 		_interpolation_matrix = _convert_double_matrix_2_mpf(interpolation_matrix);
// 		_equality_matrix = _convert_double_matrix_2_mpf(equality_matrix);
// 		_inequality_matrix = _convert_double_matrix_2_mpf(inequality_matrix);
// 		_equality_vector = _convert_double_vector_2_mpf(equality_vector);
// 		_inequality_vector = _convert_double_vector_2_mpf(inequality_vector);

// 		std::ofstream file1("interpM.txt");
// 		std::ofstream file2("ieM.txt");
// 		std::ofstream file3("ieV.txt");
// 		std::ofstream file4("eM.txt");
// 		std::ofstream file5("eV.txt");
// 		if (file1)
// 		{
// 			file1 << interpolation_matrix <<"\n";
// 			file1.close();
// 		}
// 		if (file2)
// 		{
// 			file2 << inequality_matrix <<"\n";
// 			file2.close();
// 		}
// 		if (file3)
// 		{
// 			file3 << inequality_vector <<"\n";
// 			file3.close();
// 		}
// 		if (file4)
// 		{
// 			file4 << equality_matrix <<"\n";
// 			file4.close();
// 		}
// 		if (file5)
// 		{
// 			file5 << equality_vector <<"\n";
// 			file5.close();
// 		}

		//_hessian_matrix = _get_hessian_matrix(_interpolation_matrix);
	}
	bool solve();
	bool validate_matrix_systems();

};

class Quadratic_Predictor_Corrector_LOQO : public System_Solver {
private:
	MatrixXd _H;
	MatrixXd _A;
	VectorXd _b;
	VectorXd _r;
public:
	Quadratic_Predictor_Corrector_LOQO(
		const MatrixXd &interpolation_matrix,
		const MatrixXd &inequality_matrix,
		const VectorXd &constraints,
		const VectorXd &constraints_ranges)
	{
		_H = 2.0*interpolation_matrix;
		_A = inequality_matrix;
		_b = constraints;
		_r = constraints_ranges;

		// 		// Debug
		//  		cout<<" Hessian matrix:\n"<< _H << endl;
		//  		cout<<" Inequality matrix:\n"<< _A << endl;
		// 		cout<<" Constraints :\n"<< _b << endl;
		//  		cout<<" Constraint ranges :\n"<< _r <<endl;
		// end debug
	}
	bool solve();
	bool validate_matrix_systems();
};

#endif

