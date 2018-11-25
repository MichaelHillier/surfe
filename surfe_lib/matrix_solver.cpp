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

#include <math_methods.h>
#include <matrix_solver.h>

#include <fstream>
#include <iomanip>
#include <iostream>

using namespace Surfe;
bool Linear_LU_decomposition::solve() {

  if (_constraint_values.rows() != _interpolation_matrix.rows())
    return false;

  weights = _interpolation_matrix.partialPivLu().solve(_constraint_values);

  if (!weights.allFinite())
    return false;

  return true;
}

bool Linear_LU_decomposition::validate_matrix_systems() {
  // check if there are any NaN or INF values
  if (_interpolation_matrix.allFinite())
    return true;
  else
    return false;

  // the other potential problem could be the input has un-initialized values in
  // it.
}

bool Linear_LU_decomposition::check_solution() {
  double relative_error =
      (_interpolation_matrix * weights - _constraint_values).norm() /
      _constraint_values.norm(); // norm() is L2 norm
  std::cout << "The relative error is:\n" << relative_error << std::endl;
  return true;
}

Matrix<mpf_class, Dynamic, Dynamic>
Quadratic_Predictor_Corrector::_convert_double_matrix_2_mpf(
    const MatrixXd &matrix) {
  int nrows = (int)matrix.rows();
  int ncols = (int)matrix.cols();

  Matrix<mpf_class, Dynamic, Dynamic> mpf_matrix(nrows, ncols);

  for (int j = 0; j < nrows; j++) {
    for (int k = 0; k < ncols; k++) {
      mpf_matrix(j, k) = matrix(j, k);
    }
  }

  return mpf_matrix;
}

Matrix<mpf_class, Dynamic, 1>
Quadratic_Predictor_Corrector::_convert_double_vector_2_mpf(
    const VectorXd &vector) {
  int nrows = (int)vector.rows();
  Matrix<mpf_class, Dynamic, 1> mpf_vector(nrows);

  for (int j = 0; j < nrows; j++)
    mpf_vector(j) = vector(j);
  return mpf_vector;
}

VectorXd Quadratic_Predictor_Corrector::_convert_mpf_vector_2_double(
    const Matrix<mpf_class, Dynamic, 1> &vector) {
  int nrows = (int)vector.rows();
  VectorXd v(nrows);

  for (int j = 0; j < nrows; j++)
    v(j) = vector(j).get_d();
  return v;
}

Matrix<mpf_class, Dynamic, Dynamic>
Quadratic_Predictor_Corrector::_get_hessian_matrix(
    const Matrix<mpf_class, Dynamic, Dynamic> &matrix) {
  int nrows = (int)matrix.rows();
  int ncols = (int)matrix.cols();
  Matrix<mpf_class, Dynamic, Dynamic> hessian(nrows, ncols);

  hessian = 2.0 * matrix;

  return hessian;
}

bool Quadratic_Predictor_Corrector::solve() {
  int n = (int)_hessian_matrixD.rows();

  // Matrix <mpf_class, Dynamic, 1> fvalues(n);
  Matrix<double, Dynamic, 1> fvalues(n);

  // if (!validate_matrix_systems()) return false;

  // if
  // (!Math_methods::quadratic_solver(_hessian_matrix,_equality_matrix,_inequality_matrix,_equality_vector,_inequality_vector,fvalues))
  // return false;

  if (!Math_methods::quadratic_solver(_hessian_matrixD, _equality_matrixD,
                                      _inequality_matrixD, _equality_vectorD,
                                      _inequality_vectorD, fvalues))
    return false;

  // weights = _convert_mpf_vector_2_double(fvalues);

  weights = fvalues;

  return true;
}

bool Quadratic_Predictor_Corrector::validate_matrix_systems() {
  if (!_interpolation_matrixD.allFinite())
    return false;

  LLT<MatrixXd> lltofMatrix(_interpolation_matrixD);
  if (lltofMatrix.info() == NumericalIssue)
    return false;

  return true;
}

bool Quadratic_Predictor_Corrector_LOQO::solve() {
  int n = (int)_H.rows();

  VectorXd w(n);
  if (!Math_methods::quadratic_solver_loqo(_H, _A, _b, _r, w))
    return false;

  weights = w;

  return true;
}

bool Quadratic_Predictor_Corrector_LOQO::validate_matrix_systems() {
  if (!_H.allFinite())
    return false;

  LLT<MatrixXd> lltofMatrix(_H);
  if (lltofMatrix.info() == NumericalIssue)
    return false;

  return true;
}
