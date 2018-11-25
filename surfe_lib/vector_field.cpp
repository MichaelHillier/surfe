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

#include <algorithm>
#include <basis.h>
#include <math_methods.h>
#include <matrix_solver.h>
#include <vector>
#include <vector_field.h>

#include <fstream>
#include <iomanip>
#include <iostream>

using namespace Surfe;
bool Vector_Field::get_method_parameters() {
  // # of constraints for each constraint type ...
  b_parameters.n_interface = 0;
  b_parameters.n_inequality = 0;
  b_parameters.n_planar = (int)b_input.planar->size();
  b_parameters.n_tangent = 0;
  // Total number of constraints ...
  b_parameters.n_constraints =
      b_parameters.n_interface + b_parameters.n_inequality +
      3 * b_parameters.n_planar + b_parameters.n_tangent;
  // Total number of equality constraints
  b_parameters.n_equality = b_parameters.n_interface +
                            3 * b_parameters.n_planar + b_parameters.n_tangent;

  // polynomial parameters ...

  b_parameters.poly_term = false;
  b_parameters.modified_basis = false;
  b_parameters.problem_type = Parameter_Types::Linear;

  b_parameters.n_poly_terms = 0;

  return true;
}

bool Vector_Field::process_input_data() { return true; }

bool Vector_Field::get_equality_values(VectorXd &equality_values) {
  for (int j = 0; j < (int)b_input.planar->size(); j++) {
    equality_values(3 * j) = b_input.planar->at(j).nx();
    equality_values(3 * j + 1) = b_input.planar->at(j).ny();
    equality_values(3 * j + 2) = b_input.planar->at(j).nz();
  }

  return true;
}

bool Vector_Field::get_interpolation_matrix(MatrixXd &interpolation_matrix) {
  int n_p = b_parameters.n_planar;

  // Base Matrix Structure
  // | p_x/p_x p_x/p_y p_x/p_z |
  // | p_y/p_x p_y/p_y p_y/p_z |
  // | p_z/p_x p_z/p_y p_z/p_z |

  // Planar Constraints
  for (int j = 0; j < n_p; j++) {
    // Row:planar/Column:planar block
    for (int k = 0; k < n_p; k++) {
      kernel->set_points(b_input.planar->at(j), b_input.planar->at(k));
      interpolation_matrix(3 * j, 3 * k) =
          kernel->basis_planar_planar(Parameter_Types::DXDX);
      interpolation_matrix(3 * j, 3 * k + 1) =
          kernel->basis_planar_planar(Parameter_Types::DXDY);
      interpolation_matrix(3 * j, 3 * k + 2) =
          kernel->basis_planar_planar(Parameter_Types::DXDZ);
      interpolation_matrix(3 * j + 1, 3 * k) =
          kernel->basis_planar_planar(Parameter_Types::DYDX);
      interpolation_matrix(3 * j + 1, 3 * k + 1) =
          kernel->basis_planar_planar(Parameter_Types::DYDY);
      interpolation_matrix(3 * j + 1, 3 * k + 2) =
          kernel->basis_planar_planar(Parameter_Types::DYDZ);
      interpolation_matrix(3 * j + 2, 3 * k) =
          kernel->basis_planar_planar(Parameter_Types::DZDX);
      interpolation_matrix(3 * j + 2, 3 * k + 1) =
          kernel->basis_planar_planar(Parameter_Types::DZDY);
      interpolation_matrix(3 * j + 2, 3 * k + 2) =
          kernel->basis_planar_planar(Parameter_Types::DZDZ);
    }
  }

  return true;
}

bool Vector_Field::setup_system_solver() {

  int n = b_parameters.n_equality + b_parameters.n_poly_terms;

  VectorXd equality_values;
  get_equality_values(equality_values);

  MatrixXd interpolation_matrix(n, n);
  if (!get_interpolation_matrix(interpolation_matrix))
    return false;

  Linear_LU_decomposition llu(interpolation_matrix, equality_values);
  if (!llu.solve())
    return false;
  solver = &llu;

  return true;
}

void Vector_Field::eval_scalar_interpolant_at_point(Point &p) {
  int n_p = b_parameters.n_planar;

  Kernel *kernel_j = kernel->clone();
  double elemsum_x = 0.0;
  double elemsum_y = 0.0;
  double elemsum_z = 0.0;
  for (int k = 0; k < n_p; k++) {
    kernel_j->set_points(p, b_input.planar->at(k));
    elemsum_x += solver->weights[3 * k] *
                 kernel_j->basis_planar_planar(Parameter_Types::DXDX);
    elemsum_x += solver->weights[3 * k + 1] *
                 kernel_j->basis_planar_planar(Parameter_Types::DXDY);
    elemsum_x += solver->weights[3 * k + 2] *
                 kernel_j->basis_planar_planar(Parameter_Types::DXDZ);

    elemsum_y += solver->weights[3 * k] *
                 kernel_j->basis_planar_planar(Parameter_Types::DYDX);
    elemsum_y += solver->weights[3 * k + 1] *
                 kernel_j->basis_planar_planar(Parameter_Types::DYDY);
    elemsum_y += solver->weights[3 * k + 2] *
                 kernel_j->basis_planar_planar(Parameter_Types::DYDZ);

    elemsum_z += solver->weights[3 * k] *
                 kernel_j->basis_planar_planar(Parameter_Types::DZDX);
    elemsum_z += solver->weights[3 * k + 1] *
                 kernel_j->basis_planar_planar(Parameter_Types::DZDY);
    elemsum_z += solver->weights[3 * k + 2] *
                 kernel_j->basis_planar_planar(Parameter_Types::DZDZ);
  }
  p.set_vector_field(elemsum_x, elemsum_y, elemsum_z);
  delete kernel_j;
}

void Vector_Field::eval_vector_interpolant_at_point(Point &p) {
  // this method is a copy of the eval_scalar_interpolant_at_points() method
  // TO DO: Fix eval_scalar_interpolant_at_points() method to actually computer
  // the scalar field and not the vector field.
  int n_p = b_parameters.n_planar;

  Kernel *kernel_j = kernel->clone();
  double elemsum_x = 0.0;
  double elemsum_y = 0.0;
  double elemsum_z = 0.0;
  for (int k = 0; k < n_p; k++) {
    kernel_j->set_points(p, b_input.planar->at(k));
    elemsum_x += solver->weights[3 * k] *
                 kernel_j->basis_planar_planar(Parameter_Types::DXDX);
    elemsum_x += solver->weights[3 * k + 1] *
                 kernel_j->basis_planar_planar(Parameter_Types::DXDY);
    elemsum_x += solver->weights[3 * k + 2] *
                 kernel_j->basis_planar_planar(Parameter_Types::DXDZ);

    elemsum_y += solver->weights[3 * k] *
                 kernel_j->basis_planar_planar(Parameter_Types::DYDX);
    elemsum_y += solver->weights[3 * k + 1] *
                 kernel_j->basis_planar_planar(Parameter_Types::DYDY);
    elemsum_y += solver->weights[3 * k + 2] *
                 kernel_j->basis_planar_planar(Parameter_Types::DYDZ);

    elemsum_z += solver->weights[3 * k] *
                 kernel_j->basis_planar_planar(Parameter_Types::DZDX);
    elemsum_z += solver->weights[3 * k + 1] *
                 kernel_j->basis_planar_planar(Parameter_Types::DZDY);
    elemsum_z += solver->weights[3 * k + 2] *
                 kernel_j->basis_planar_planar(Parameter_Types::DZDZ);
  }
  p.set_vector_field(elemsum_x, elemsum_y, elemsum_z);
  delete kernel_j;
}
