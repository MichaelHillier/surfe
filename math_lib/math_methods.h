// SURFace Estimator(SURFE) - Terms and Conditions of Use

// Unless otherwise noted, computer program source code of the SURFace
// Estimator(SURFE) is covered under Crown Copyright, Government of Canada, and
// is distributed under the MIT License.

// The Canada wordmark and related graphics associated with this distribution
// are protected under trademark law and copyright law.No permission is granted
// to use them outside the parameters of the Government of Canada's corporate
// identity program. For more information, see
// http://www.tbs-sct.gc.ca/fip-pcim/index-eng.asp

// Copyright title to all 3rd party software distributed with the SURFace
// Estimator(SURFE) is held by the respective copyright holders as noted in
// those files.Users are asked to read the 3rd Party Licenses referenced with
// those assets.

// MIT License

// Copyright(c) 2017 Government of Canada

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files(the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and / or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions :

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT
// NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
// PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef math_methods_h
#define math_methods_h

#include <math_lib_module.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
using namespace Eigen;

class MATH_LIB_EXPORT Math_methods {
private:
	static double _find_step_length(
		const VectorXd &a,
		const VectorXd &da,
		const VectorXd &b,
		const VectorXd &db);
	static void _rot(std::vector<std::vector<double> > &a, const double &s, const double &tau,
		const int &i, const int &j, const int &k, const int &l);
	static double _get_double(const double &d) { return d; }
	static double _find_step(const VectorXd &da, const VectorXd &a);
	static double _find_positivity_step(
		const VectorXd &da, const VectorXd &a,
		const VectorXd &db, const VectorXd &b,
		const VectorXd &dc, const VectorXd &c,
		const VectorXd &dd, const VectorXd &d);

public:
	static bool sort_vector_w_index(std::vector<double> &arr, std::vector<int> &brr);
	static double max_element_wrt_zero(const double &a, const double &b);
	template <class T>
	static void SWAP(T &a, T &b) {
		T x = a;
		a = b;
		b = x;
	}
	static bool angle_btw_2_vectors(
		const std::vector<double> &v1, const std::vector<double> &v2, double &angle);
	static double RandomDouble(const double &min, const double &max);
	static bool quadratic_solver(
		const MatrixXd &H,
		const MatrixXd &A,
		const MatrixXd &C,
		const VectorXd &b,
		const VectorXd &d,
		VectorXd &fvalues);
	static bool quadratic_solver_loqo(
		const MatrixXd &H,
		const MatrixXd &A,
		const VectorXd &b,
		const VectorXd &r,
		VectorXd &fvalues);
};


#endif
