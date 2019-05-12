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
#include <continuous_property.h>
#include <math_methods.h>
#include <matrix_solver.h>
#include <vector>

#include <fstream>
#include <iomanip>
#include <iostream>

bool Continuous_Property::_get_polynomial_matrix_block(MatrixXd &poly_matrix) {
    int n_ie = b_parameters.n_inequality;
    int n_i = b_parameters.n_interface;
    int n_p = b_parameters.n_planar;
    int n_t = b_parameters.n_tangent;

    int nl = n_i + n_ie;  // n_ie should always be zero.

    p_basis = create_polynomial_basis(m_parameters.polynomial_order);

    if ((int)poly_matrix.size() != b_parameters.n_poly_terms) return false;
    // for interface points ...
    for (int j = 0; j < nl; j++) {
        p_basis->set_point(constraints.itrface[j]);
        VectorXd b = p_basis->basis();
        if ((int)b.rows() != b_parameters.n_poly_terms) return false;
        for (int k = 0; k < (int)b.rows(); k++) poly_matrix(k, j) = b(k);
    }
    // for planar points ...
    for (int j = 0; j < n_p; j++) {
        p_basis->set_point(constraints.planar[j]);
        VectorXd bx = p_basis->dx();
        VectorXd by = p_basis->dy();
        VectorXd bz = p_basis->dz();
        if ((int)bx.rows() != b_parameters.n_poly_terms) return false;
        for (int k = 0; k < (int)bx.rows(); k++) {
            poly_matrix(k, nl + 3 *j) = bx(k);
            poly_matrix(k, nl + 3 *j + 1) = by(k);
            poly_matrix(k, nl + 3 *j + 2) = bz(k);
        }
    }
    // for tangent points ...
    for (int j = 0; j < n_t; j++) {
        p_basis->set_point(constraints.tangent[j]);
        VectorXd bx = p_basis->dx();
        VectorXd by = p_basis->dy();
        VectorXd bz = p_basis->dz();
        if ((int)bx.rows() != b_parameters.n_poly_terms) return false;
        for (int k = 0; k < (int)bx.rows(); k++) {
            poly_matrix(k, nl + 3 *n_p + j) =
                constraints.tangent[j].tx() * bx(k) +
                constraints.tangent[j].ty() * by(k) +
                constraints.tangent[j].tz() * bz(k);
        }
    }

    return true;
}

bool Continuous_Property::_insert_polynomial_matrix_blocks_in_interpolation_matrix(
    const MatrixXd &poly_matrix, MatrixXd &interpolation_matrix) {
    int n_ie = b_parameters.n_inequality;
    int n_i = b_parameters.n_interface;
    int n_p = b_parameters.n_planar;
    int n_t = b_parameters.n_tangent;

    // build polynomial blocks
    // | A PT |
    // | P 0  |
    // start with P
    for (int j = 0; j < (int)poly_matrix.rows(); j++) {
        for (int k = 0; k < (int)poly_matrix.cols(); k++) {
            interpolation_matrix(n_ie + n_i + 3 * n_p + n_t + j, k) = poly_matrix(j, k);
            interpolation_matrix(k, n_ie + n_i + 3 *n_p + n_t + j) = interpolation_matrix(n_ie + n_i + 3 * n_p + n_t + j, k);
        }
    }

    for (int j = 0; j < (int)poly_matrix.size(); j++) {
        for (int k = 0; k < (int)poly_matrix.size(); k++) {
            interpolation_matrix(n_ie + n_i + 3 * n_p + n_t + j, n_ie + n_i + 3 *n_p + n_t + k) = 0;
        }
    }

    return true;
}

Continuous_Property::Continuous_Property(const model_parameters& mparams)
{
	// set GUI parameters
	m_parameters = mparams;

	_iteration = 0;
}

Continuous_Property::~Continuous_Property() {
    std::cout << "dest" << std::endl;
}
Polynomial_Basis *Continuous_Property::create_polynomial_basis(
    const int &poly_order) {
    if (poly_order == 0)
        return new Poly_Zero;
    else if (poly_order == 1)
        return new Poly_First;
    else
        return new Poly_Second;
}

void Continuous_Property::get_method_parameters() {
    // # of constraints for each constraint type ...
    b_parameters.n_interface = (int)constraints.itrface.size();
    b_parameters.n_inequality = 0;
    b_parameters.n_planar = 0;
    b_parameters.n_tangent = 0;
    // Total number of constraints ...
    b_parameters.n_constraints = b_parameters.n_interface + b_parameters.n_inequality +
        3 * b_parameters.n_planar + b_parameters.n_tangent;
    // Total number of equality constraints
    b_parameters.n_equality = b_parameters.n_interface + 3 * b_parameters.n_planar + b_parameters.n_tangent;

    // polynomial parameters ...
    if (b_parameters.n_inequality == 0) 
	{
        b_parameters.poly_term = false;  // NOTE: May want to have this as an option when using SPD functions
        b_parameters.modified_basis = false;
        b_parameters.problem_type = Parameter_Types::Linear;
    } 
	else
	{
        b_parameters.poly_term = false;
        b_parameters.modified_basis = true;
        b_parameters.problem_type = Parameter_Types::Quadratic;
    }

    int m = m_parameters.polynomial_order + 1;
    b_parameters.n_poly_terms = 0;  // for 3D only...
}

void Continuous_Property::setup_system_solver() {
    int n = b_parameters.n_equality + b_parameters.n_poly_terms;

    VectorXd equality_values(n);
    get_equality_values(equality_values);
    MatrixXd interpolation_matrix(n, n);

	if (!get_interpolation_matrix(interpolation_matrix))
		std::throw_with_nested(GRBF_Exceptions::error_computing_interpolation_matrix);

    Linear_LU_decomposition *llu = new Linear_LU_decomposition(interpolation_matrix, equality_values);
	if (!llu->solve())
		std::throw_with_nested(GRBF_Exceptions::linear_solver_failure);
    solver = llu;
}

// bool Continuous_Property::get_minimial_greedy_input( Basic_input
// &greedy_input )
// {
// 	// get extremal interface points
// 	// Note:
// 	// if poly order = 1 find 4 interface points nicely sampling the volume
// 	// if poly order = 2 find 6 interface points nicely sampling the volume
// 	std::vector < Interface > extremal_interace_pts;
// 	std::vector < Point >
// pts(b_input.interface.begin(),b_input.interface.end()); 	std::vector <
// int >
// interface_indices = get_extremal_point_data_indices_from_points(pts); 	if
// (
// (int)interface_indices.size() < b_parameters.n_poly_terms) return false;
// for
// (int j = 0; j < b_parameters.n_poly_terms; j++ )
// extremal_interace_pts.push_back(b_input.interface[interface_indices[j]]);
//
// 	greedy_input.interface = extremal_interace_pts;
//
// 	return true;
// }

bool Continuous_Property::measure_residuals(Constraints &input)
{
    if (solver == nullptr) return false;

    // inequalities points
	for (auto &inequality_pt: input.inequality){
        eval_scalar_interpolant_at_point(inequality_pt);
        if (inequality_pt.level() >= 0) {
            if (inequality_pt.scalar_field() >= 0)
				inequality_pt.setResidual(true);
            else
				inequality_pt.setResidual(false);
        } else {
            if (inequality_pt.scalar_field() < 0)
				inequality_pt.setResidual(true);
            else
				inequality_pt.setResidual(false);
        }
    }
    // interface points
	for (auto &interface_pt : input.itrface){
        eval_scalar_interpolant_at_point(interface_pt);
		interface_pt.setResidual(abs(interface_pt.scalar_field() - interface_pt.level())); // THIS DOES NOT WORK FOR INCREMENT POINTS
    }
    // planar points
	for (auto &planar_pt : input.planar){
        eval_vector_interpolant_at_point(planar_pt);
        double angle = 0.0;
        std::vector<double> v1;
        std::vector<double> v2;
        v1.push_back(planar_pt.nx());
        v1.push_back(planar_pt.ny());
        v1.push_back(planar_pt.nz());
        v2.push_back(planar_pt.nx_interp());
        v2.push_back(planar_pt.ny_interp());
        v2.push_back(planar_pt.nz_interp());
        Math_methods::angle_btw_2_vectors(v1, v2, angle);
		planar_pt.setResidual(angle);
    }
    // tangent points
	for (auto &tangent_pt: input.tangent){
        eval_vector_interpolant_at_point(tangent_pt);
        double angle = 0.0;
        std::vector<double> v1;
        std::vector<double> v2;
        v1.push_back(tangent_pt.tx());
        v1.push_back(tangent_pt.ty());
        v1.push_back(tangent_pt.tz());
        v2.push_back(tangent_pt.nx_interp());
        v2.push_back(tangent_pt.ny_interp());
        v2.push_back(tangent_pt.nz_interp());
        Math_methods::angle_btw_2_vectors<double>(v1, v2, angle);
		tangent_pt.setResidual(angle);
    }

    return true;
}

bool Continuous_Property::append_greedy_input(Constraints &input) {
    // planar > tangent > interface > inequalities

    double r2d = 57.29577951308232;

    // PLANAR Observations
    std::vector<double> large_planar_residuals;
    std::vector<int> large_planar_residuals_indices;
    for (int j = 0; j < (int)input.planar.size(); j++) {
        double grad_err = input.planar[j].residual() * r2d;
        if (grad_err > m_parameters.angular_uncertainty) {
            large_planar_residuals.push_back(grad_err);
            large_planar_residuals_indices.push_back(j);
        }
    }
    if (!large_planar_residuals.empty()) {
        Math_methods::sort_vector_w_index(large_planar_residuals, large_planar_residuals_indices);
        constraints.planar.emplace_back(input.planar[large_planar_residuals_indices[large_planar_residuals.size() - 1]]);
        return true;
    }
    // TANGENT Observations
	for (auto &tangent_pt : input.tangent){
        if (tangent_pt.residual() * r2d > m_parameters.angular_uncertainty) {
            this->constraints.tangent.emplace_back(tangent_pt);
            return true;
        }
    }
    // INTERFACE Observations
    std::vector<double> large_interface_residuals;
    std::vector<int> large_interface_residuals_indices;
    for (int j = 0; j < (int)input.itrface.size(); j++) {
        double interface_err = input.itrface.at(j).residual();
        if (interface_err > m_parameters.interface_uncertainty) {
            large_interface_residuals.push_back(interface_err);
            large_interface_residuals_indices.push_back(j);
        }
    }
    if (!large_interface_residuals.empty()) {
        Math_methods::sort_vector_w_index(large_interface_residuals, large_interface_residuals_indices);
		constraints.itrface.emplace_back(input.itrface[large_interface_residuals_indices[large_interface_residuals.size() - 1]]);
        return true;
    }
    // INEQUALITIES Observations
	for (auto &inequality_pt : input.inequality) {
        if (!inequality_pt.residual()) {
            constraints.inequality.emplace_back(inequality_pt);
            return true;
        }
    }

    return false;
}

void Continuous_Property::eval_scalar_interpolant_at_point(Point &p) {
    int n_i = b_parameters.n_interface;
    int n_p = b_parameters.n_planar;
    int n_t = b_parameters.n_tangent;

    Kernel *kernel_j = kernel->clone();
    double elemsum_1 = 0.0;
    double elemsum_2 = 0.0;
    double elemsum_3 = 0.0;
    double elemsum_4 = 0.0;
    double poly = 0.0;
    for (int k = 0; k < n_i; k++) {
        kernel_j->set_points(p, constraints.itrface[k]);
        elemsum_2 += solver->weights[k] * kernel_j->basis_pt_pt();
    }

    for (int k = 0; k < n_p; k++) {
        kernel_j->set_points(p, constraints.planar[k]);
        elemsum_3 += solver->weights[n_i + 3 * k] * kernel_j->basis_pt_planar_x();
        elemsum_3 += solver->weights[n_i + 3 * k + 1] * kernel_j->basis_pt_planar_y();
        elemsum_3 += solver->weights[n_i + 3 * k + 2] * kernel_j->basis_pt_planar_z();
    }

    for (int k = 0; k < n_t; k++) {
        kernel_j->set_points(p, constraints.tangent[k]);
        elemsum_4 += solver->weights[n_i + 3 * n_p + k] * kernel_j->basis_pt_tangent();
    }

    if (b_parameters.poly_term) {
        Polynomial_Basis *p_basis_j = p_basis->clone();
        p_basis_j->set_point(p);
        VectorXd b = p_basis_j->basis();
        for (int k = 0; k < (int)b.size(); k++)
            poly += b(k) * solver->weights[n_i + 3 * n_p + n_t + k];
        delete p_basis_j;
    }
    p.set_scalar_field(elemsum_1 + elemsum_2 + elemsum_3 + elemsum_4 + poly);
    delete kernel_j;
}

void Continuous_Property::eval_vector_interpolant_at_point(Point &p) {
    int n_i = b_parameters.n_interface;
    int n_p = b_parameters.n_planar;
    int n_t = b_parameters.n_tangent;

    Kernel *kernel_j = kernel->clone();
    double elemsum_1_x = 0.0;
    double elemsum_1_y = 0.0;
    double elemsum_1_z = 0.0;
    double elemsum_2_x = 0.0;
    double elemsum_2_y = 0.0;
    double elemsum_2_z = 0.0;
    double elemsum_3_x = 0.0;
    double elemsum_3_y = 0.0;
    double elemsum_3_z = 0.0;
    double poly_x = 0.0;
    double poly_y = 0.0;
    double poly_z = 0.0;
    // interface constraints
    for (int k = 0; k < n_i; k++) {
        kernel->set_points(p, constraints.itrface[k]);
        elemsum_1_x += solver->weights[k] * kernel->basis_planar_x_pt();
        elemsum_1_y += solver->weights[k] * kernel->basis_planar_y_pt();
        elemsum_1_z += solver->weights[k] * kernel->basis_planar_z_pt();
    }
    // normal constraints
    for (int k = 0; k < n_p; k++) {
        kernel->set_points(p, constraints.planar[k]);
        elemsum_2_x += solver->weights[n_i + 0 + 3 * k] * kernel->basis_planar_planar(Parameter_Types::DXDX);
        elemsum_2_x += solver->weights[n_i + 1 + 3 * k] * kernel->basis_planar_planar(Parameter_Types::DXDY);
        elemsum_2_x += solver->weights[n_i + 2 + 3 * k] * kernel->basis_planar_planar(Parameter_Types::DXDZ);
        elemsum_2_y += solver->weights[n_i + 0 + 3 * k] * kernel->basis_planar_planar(Parameter_Types::DYDX);
        elemsum_2_y += solver->weights[n_i + 1 + 3 * k] * kernel->basis_planar_planar(Parameter_Types::DYDY);
        elemsum_2_y += solver->weights[n_i + 2 + 3 * k] * kernel->basis_planar_planar(Parameter_Types::DYDZ);
        elemsum_2_z += solver->weights[n_i + 0 + 3 * k] * kernel->basis_planar_planar(Parameter_Types::DZDX);
        elemsum_2_z += solver->weights[n_i + 1 + 3 * k] * kernel->basis_planar_planar(Parameter_Types::DZDY);
        elemsum_2_z += solver->weights[n_i + 2 + 3 * k] * kernel->basis_planar_planar(Parameter_Types::DZDZ);
    }
    // tangent constraints
    for (int k = 0; k < n_t; k++) {
        kernel->set_points(p, constraints.tangent[k]);
        elemsum_3_x += solver->weights[n_i + 3 * n_p + k] * kernel->basis_planar_tangent(Parameter_Types::DX);
        elemsum_3_y += solver->weights[n_i + 3 * n_p + k] * kernel->basis_planar_tangent(Parameter_Types::DY);
        elemsum_3_z += solver->weights[n_i + 3 * n_p + k] * kernel->basis_planar_tangent(Parameter_Types::DZ);
    }
    if (b_parameters.poly_term) {
        Polynomial_Basis *p_basis_j = p_basis->clone();
        p_basis_j->set_point(p);
        VectorXd bx = p_basis_j->dx();
        VectorXd by = p_basis_j->dy();
        VectorXd bz = p_basis_j->dz();
        for (int k = 0; k < (int)bx.size(); k++) {
            poly_x += bx(k) * solver->weights[n_i + 3 * n_p + n_t + k];
            poly_y += by(k) * solver->weights[n_i + 3 * n_p + n_t + k];
            poly_z += bz(k) * solver->weights[n_i + 3 * n_p + n_t + k];
        }
        delete p_basis_j;
    }
    double nx = elemsum_1_x + elemsum_2_x + elemsum_3_x + poly_x;
    double ny = elemsum_1_y + elemsum_2_y + elemsum_3_y + poly_y;
    double nz = elemsum_1_z + elemsum_2_z + elemsum_3_z + poly_z;
    p.set_vector_field(nx, ny, nz);
    delete kernel_j;
}

bool Continuous_Property::get_equality_values(VectorXd &equality_values) {
    int j = 0;
    int k = 0;
    int l = 0;
    int m = 0;
    for (j = 0; j < (int)constraints.itrface.size(); j++)
        equality_values(j) = constraints.itrface[j].level();
    for (k = 0; k < (int)constraints.planar.size(); k++) {
        equality_values(3 *k + j) = constraints.planar[k].nx();
        equality_values(3 *k + j + 1) = constraints.planar[k].ny();
        equality_values(3 *k + j + 2) = constraints.planar[k].nz();
    }
    for (l = 0; l < (int)constraints.tangent.size(); l++)
        equality_values(l + 3 *k + j) = 0.0;
    if (b_parameters.poly_term)
        for (m = 0; m < (int)b_parameters.n_poly_terms; m++)
            equality_values(m + l + 3 *k + j) = 0.0;

    return true;
}

void Continuous_Property::process_input_data() {
	if (constraints.itrface.empty())
		std::throw_with_nested(GRBF_Exceptions::no_iterface_data);
}

bool Continuous_Property::get_interpolation_matrix(
    MatrixXd &interpolation_matrix) {
    int n_i = b_parameters.n_interface;
    int n_p = b_parameters.n_planar;
    int n_t = b_parameters.n_tangent;

    // Row and Column constraint order : interface (itr) -> planar (p_x,p_y,p_z)
    // -> tangent (t)

    // Base Matrix Structure
    // | itr/itr itr/p_y itr/p_y itr/p_z itr/t |
    // | p_x/itr p_x/p_x p_x/p_y p_x/p_z p_x/t |
    // | p_y/itr p_y/p_x p_y/p_y p_y/p_z p_y/t |
    // | p_z/itr p_z/p_x p_z/p_y p_z/p_z p_z/t |
    // |   t/itr   t/p_x   t/p_y   t/p_z   t/t |

    // Interface Constraints:
    for (int j = 0; j < n_i; j++) {
        // Row:interface/Column:interface block
        for (int k = 0; k < n_i; k++) {
            kernel->set_points(constraints.itrface[j], constraints.itrface[k]);
            interpolation_matrix(j, k) = kernel->basis_pt_pt();
        }
        // Row:interface/Column:planar block
        for (int k = 0; k < n_p; k++) {
            kernel->set_points(constraints.itrface[j], constraints.planar[k]);
            interpolation_matrix(j, 3 *k + n_i) = kernel->basis_pt_planar_x();
            interpolation_matrix(j, 3 *k + n_i + 1) = kernel->basis_pt_planar_y();
            interpolation_matrix(j, 3 *k + n_i + 2) = kernel->basis_pt_planar_z();
        }
        // Row:interface/Column:tangent block
        for (int k = 0; k < n_t; k++) {
            kernel->set_points(constraints.itrface[j], constraints.tangent[k]);
            interpolation_matrix(j, n_i + 3 *n_p + k) = kernel->basis_pt_tangent();
        }
    }
    // Planar Constraints
    for (int j = 0; j < n_p; j++) {
        // Row:planar/Column:interface block
        for (int k = 0; k < n_i; k++) {
            kernel->set_points(constraints.planar[j], constraints.itrface[k]);
            interpolation_matrix(3 * j + n_i, k) = kernel->basis_planar_x_pt();
            interpolation_matrix(3 * j + n_i + 1, k) = kernel->basis_planar_y_pt();
            interpolation_matrix(3 * j + n_i + 2, k) = kernel->basis_planar_z_pt();
        }
        // Row:planar/Column:planar block
        for (int k = 0; k < n_p; k++) {
            kernel->set_points(constraints.planar[j], constraints.planar[k]);
            interpolation_matrix(3 * j + n_i, 3 *k + n_i) = kernel->basis_planar_planar(Parameter_Types::DXDX);
            interpolation_matrix(3 * j + n_i, 3 *k + n_i + 1) = kernel->basis_planar_planar(Parameter_Types::DXDY);
            interpolation_matrix(3 * j + n_i, 3 *k + n_i + 2) = kernel->basis_planar_planar(Parameter_Types::DXDZ);
            interpolation_matrix(3 * j + n_i + 1, 3 *k + n_i) = kernel->basis_planar_planar(Parameter_Types::DYDX);
            interpolation_matrix(3 * j + n_i + 1, 3 *k + n_i + 1) = kernel->basis_planar_planar(Parameter_Types::DYDY);
            interpolation_matrix(3 * j + n_i + 1, 3 *k + n_i + 2) = kernel->basis_planar_planar(Parameter_Types::DYDZ);
            interpolation_matrix(3 * j + n_i + 2, 3 *k + n_i) = kernel->basis_planar_planar(Parameter_Types::DZDX);
            interpolation_matrix(3 * j + n_i + 2, 3 *k + n_i + 1) = kernel->basis_planar_planar(Parameter_Types::DZDY);
            interpolation_matrix(3 * j + n_i + 2, 3 *k + n_i + 2) = kernel->basis_planar_planar(Parameter_Types::DZDZ);
        }
        // Row:planar/Column:tangent block
        for (int k = 0; k < n_t; k++) {
            kernel->set_points(constraints.planar[j], constraints.tangent[k]);
            interpolation_matrix(3 * j + n_i, n_i + 3 *n_p + k) = kernel->basis_planar_tangent(Parameter_Types::DX);
            interpolation_matrix(3 * j + n_i + 1, n_i + 3 *n_p + k) = kernel->basis_planar_tangent(Parameter_Types::DY);
            interpolation_matrix(3 * j + n_i + 2, n_i + 3 *n_p + k) = kernel->basis_planar_tangent(Parameter_Types::DZ);
        }
    }
    // Tangent Constraints
    for (int j = 0; j < n_t; j++) {
        // Row:tangent/Column:interface block
        for (int k = 0; k < n_i; k++) {
            kernel->set_points(constraints.tangent[j], constraints.itrface[k]);
            interpolation_matrix(j + n_i + 3 * n_p, k) = kernel->basis_tangent_pt();
        }
        // Row:tangent/Column:planar block
        for (int k = 0; k < n_p; k++) {
            kernel->set_points(constraints.tangent[j], constraints.planar[k]);
            interpolation_matrix(j + n_i + 3 * n_p, 3 *k + n_i) = kernel->basis_tangent_planar(Parameter_Types::DX);
            interpolation_matrix(j + n_i + 3 * n_p, 3 *k + n_i + 1) = kernel->basis_tangent_planar(Parameter_Types::DY);
            interpolation_matrix(j + n_i + 3 * n_p, 3 *k + n_i + 2) = kernel->basis_tangent_planar(Parameter_Types::DZ);
        }
        // Row:tangent/Column:tangent block
        for (int k = 0; k < n_t; k++) {
            kernel->set_points(constraints.tangent[j], constraints.tangent[k]);
            interpolation_matrix(j + n_i + 3 * n_p, n_i + 3 *n_p + k) = kernel->basis_tangent_tangent();
        }
    }

    // build polynomial blocks if required
    // | A PT |
    // | P 0  |
    if (b_parameters.poly_term) {
        MatrixXd poly_matrix(b_parameters.n_poly_terms, b_parameters.n_constraints);
        if (!_get_polynomial_matrix_block(poly_matrix)) 
			return false;
        // fill remaining matrix blocks (P, PT, 0)
        if (!_insert_polynomial_matrix_blocks_in_interpolation_matrix(poly_matrix, interpolation_matrix))
            return false;
    }

    return true;
}
