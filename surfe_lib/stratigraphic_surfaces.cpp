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
#include <stratigraphic_surfaces.h>
#include <vector>

#include <fstream>
#include <iomanip>
#include <iostream>

bool Stratigraphic_Surfaces::_get_increment_pairs() {
	// 1) sequenced contacts
	_n_sequenced_interface_pairs = (int)interface_test_points.size() - 1;
	_increment_pairs.resize(_n_sequenced_interface_pairs);
	for (int j = 0; j < _n_sequenced_interface_pairs; j++) {
		_increment_pairs[j].push_back(interface_test_points[j]);
		_increment_pairs[j].push_back(interface_test_points[j + 1]);
	}
	// 2) lithostratigraphic inequalities
	_n_sequenced_inequality_pairs = 0;
	int idx = 0;
	for (const auto &inequality : constraints.inequality) {
		std::vector<std::vector<Point> > strati_seq_ine_pairs =
			_get_lithostratigraphic_increment_pairs_for_inequality_point(inequality);
		for (const auto &seq_ineq_pair : strati_seq_ine_pairs) {
			Interface a_itr_pt1(seq_ineq_pair[0].x(), seq_ineq_pair[0].y(), seq_ineq_pair[0].z(), 0.0);
			Interface a_itr_pt2(seq_ineq_pair[1].x(), seq_ineq_pair[1].y(), seq_ineq_pair[1].z(), 0.0);
			std::vector<Interface> pair;
			pair.push_back(a_itr_pt1);
			pair.push_back(a_itr_pt2);
			_increment_pairs.push_back(pair);
		}
		int j_size = (int)strati_seq_ine_pairs.size();
		_n_sequenced_inequality_pairs += j_size;
		idx++;
	}
	// 3) the interface increment pairs
	_n_interface_pairs = 0;
	for (const auto &point_list : interface_point_lists) {
		_n_interface_pairs += ((int)point_list.size() - 1);
		for (int k = 0; k < ((int)point_list.size() - 1); k++) {
			std::vector<Interface> interface_incr_p;
			interface_incr_p.push_back(point_list[0]);
			interface_incr_p.push_back(point_list[k + 1]);
			_increment_pairs.push_back(interface_incr_p);
		}
	}
	_n_increment_pairs = _n_sequenced_interface_pairs + _n_sequenced_inequality_pairs + _n_interface_pairs;

	return true;
}

std::vector<std::vector<Point> > Stratigraphic_Surfaces::
_get_lithostratigraphic_increment_pairs_for_inequality_point(const Inequality &ie_pt)
{
	std::vector<std::vector<Point> > strati_incr_p;
	// if horizon is above current pt (ie_pt)
	// strati_incr_p[0][0] = above_horizon_pt (interface)
	// strati_incr_p[0][1] = ie_pt
	// if horizon is below current pt
	// strati_incr_p[1][0] = ie_pt
	// strati_incr_p[1][1] = below_horizon_pt (interface)
	// the above structure is needed to ensure the inequality constraint is
	// s(p1) > s (p2) e.g. where p1 = strati_incr_p[1][0] p2 = strati_incr_p[1][1]
	int n_pair = 0;

	std::vector<double> horizon_levels;
	for (const auto &test_point : interface_test_points)
		horizon_levels.push_back(test_point.level());

	// are there horizons above ie_pt ?
	double level_above = 0;
	if (_get_closest_horizon_level_above_given_level(ie_pt.level(), horizon_levels, level_above)) {
		strati_incr_p.resize(n_pair + 1);
		for (const auto &test_point : interface_test_points) {
			if (test_point.level() == level_above) {
				strati_incr_p[n_pair].push_back(test_point);
				strati_incr_p[n_pair].push_back(ie_pt);
				break;
			}
		}
		n_pair++;
	}
	// are there horizons below ie_pt ?
	double level_below = 0;
	if (_get_closest_horizon_level_below_given_level(ie_pt.level(), horizon_levels, level_below)) {
		strati_incr_p.resize(n_pair + 1);
		for (const auto &test_point : interface_test_points) {
			if (test_point.level() == level_below) {
				strati_incr_p[n_pair].push_back(ie_pt);
				strati_incr_p[n_pair].push_back(test_point);
				break;
			}
		}
		n_pair++;
	}

	return strati_incr_p;
}

bool Stratigraphic_Surfaces::_get_closest_horizon_level_above_given_level(
	const double &given_level, const std::vector<double> &horizon_levels, double &above_level)
{
	if (horizon_levels.empty()) return NULL;
	std::vector<double> diff;
	std::vector<int> indx;

	for (int j = 0; j < (int)horizon_levels.size(); j++) {
		diff.push_back(horizon_levels[j] - given_level);
		indx.push_back(j);
	}
	// sort these level differences so that we can find/determine which horizon/interface 
	// is ABOVE the given inequality level
	Math_methods::sort_vector_w_index(diff, indx);
	for (int j = 0; j < diff.size(); j++) {
		if (diff[j] > 0) {
			above_level = horizon_levels[indx[j]];
			return true;
		}
	}
	// if we got here than there are no horizons/interfaces ABOVE the current
	// inequality level (given_level)
	return false;
}

bool Stratigraphic_Surfaces::_get_closest_horizon_level_below_given_level(
	const double &given_level, const std::vector<double> &horizon_levels, double &below_level)
{
	if (horizon_levels.empty()) return NULL;
	std::vector<double> diff;
	std::vector<int> indx;

	for (int j = 0; j < (int)horizon_levels.size(); j++) {
		diff.push_back(given_level - horizon_levels[j]);
		indx.push_back(j);
	}
	// sort these level differences so that we can find/determine which horizon/interface 
	// is BELOW the given inequality level
	Math_methods::sort_vector_w_index(diff, indx);
	for (int j = 0; j < diff.size(); j++) {
		if (diff[j] > 0) {
			below_level = horizon_levels[indx[j]];
			return true;
		}
	}
	// if we got here than there are no horizons/interfaces BELOW the current
	// inequality level (given_level)
	return false;
}

bool Stratigraphic_Surfaces::_get_polynomial_matrix_block(
	MatrixXd &poly_matrix) {
	int n_ip = _n_increment_pairs;
	int n_p = intern_params.n_planar;
	int n_t = intern_params.n_tangent;

	p_basis = create_polynomial_basis(1);

	if ((int)poly_matrix.rows() != intern_params.n_poly_terms) return false;

	// for Interface Increment Pair Constraints:
	for (int j = 0; j < (int)_increment_pairs.size(); j++) {
		p_basis->set_point(_increment_pairs[j][0]);
		VectorXd b1 = p_basis->basis();
		p_basis->set_point(_increment_pairs[j][1]);
		VectorXd b2 = p_basis->basis();
		if ((int)b1.rows() != intern_params.n_poly_terms ||
			(int)b2.rows() != intern_params.n_poly_terms)
			return false;
		for (int k = 0; k < (int)b1.rows(); k++)
			poly_matrix(k, j) = b1(k) - b2(k);
	}
	// for planar points ...
	for (int j = 0; j < n_p; j++) {
		p_basis->set_point(constraints.planar[j]);
		VectorXd bx = p_basis->dx();
		VectorXd by = p_basis->dy();
		VectorXd bz = p_basis->dz();
		if ((int)bx.rows() != intern_params.n_poly_terms ||
			(int)by.rows() != intern_params.n_poly_terms ||
			(int)bz.rows() != intern_params.n_poly_terms)
			return false;
		for (int k = 0; k < (int)bx.rows(); k++) {
			poly_matrix(k, n_ip + 3 * j) = bx(k);
			poly_matrix(k, n_ip + 3 * j + 1) = by(k);
			poly_matrix(k, n_ip + 3 * j + 2) = bz(k);
		}
	}
	// for tangent points ...
	for (int j = 0; j < n_t; j++) {
		p_basis->set_point(constraints.tangent[j]);
		VectorXd bx = p_basis->dx();
		VectorXd by = p_basis->dy();
		VectorXd bz = p_basis->dz();
		if ((int)bx.rows() != intern_params.n_poly_terms ||
			(int)by.rows() != intern_params.n_poly_terms ||
			(int)bz.rows() != intern_params.n_poly_terms)
			return false;
		for (int k = 0; k < (int)bx.rows(); k++) {
			poly_matrix(k, n_ip + 3 * n_p + j) =
				constraints.tangent[j].tx() * bx(k) +
				constraints.tangent[j].ty() * by(k) +
				constraints.tangent[j].tz() * bz(k);
		}
	}

	return true;
}

bool Stratigraphic_Surfaces::_insert_polynomial_matrix_blocks_in_interpolation_matrix(
	const MatrixXd &poly_matrix, MatrixXd &interpolation_matrix)
{
	int n_ip = _n_increment_pairs;
	int n_p = intern_params.n_planar;
	int n_t = intern_params.n_tangent;

	if ((n_ip + 3 * n_p + n_t + poly_matrix.rows()) > interpolation_matrix.rows() ||
		(n_ip + 3 * n_p + n_t + poly_matrix.rows()) > interpolation_matrix.cols())
		return false;
	// build polynomial blocks
	// | A PT |
	// | P 0  |
	// start with P
	for (int j = 0; j < (int)poly_matrix.rows(); j++) {
		for (int k = 0; k < (int)poly_matrix.cols(); k++) {
			interpolation_matrix(n_ip + 3 * n_p + n_t + j, k) = poly_matrix(j, k);
			interpolation_matrix(k, n_ip + 3 * n_p + n_t + j) =
				interpolation_matrix(n_ip + 3 * n_p + n_t + j, k);
		}
	}

	for (int j = 0; j < (int)poly_matrix.rows(); j++) {
		for (int k = 0; k < (int)poly_matrix.rows(); k++) {
			interpolation_matrix(n_ip + 3 * n_p + n_t + j, n_ip + 3 * n_p + n_t + k) = 0;
		}
	}

	return true;
}

Stratigraphic_Surfaces::Stratigraphic_Surfaces(const Parameters& m_params)
{
	solver = nullptr;
	kernel = nullptr;
	rbf_kernel = nullptr;
	p_basis = nullptr;

	// set GUI parameters
	parameters = m_params;

	_n_increment_pairs = 0;
	_n_sequenced_interface_pairs = 0;
	_n_sequenced_inequality_pairs = 0;
	_n_interface_pairs = 0;

	_iteration = 0;
}

Stratigraphic_Surfaces::Stratigraphic_Surfaces()
{
	solver = nullptr;
	kernel = nullptr;
	rbf_kernel = nullptr;
	p_basis = nullptr;

	_n_increment_pairs = 0;
	_n_sequenced_interface_pairs = 0;
	_n_sequenced_inequality_pairs = 0;
	_n_interface_pairs = 0;

	_iteration = 0;
}

void Stratigraphic_Surfaces::get_method_parameters() {
	// # of constraints for each constraint type ...
	intern_params.n_interface = (int)constraints.itrface.size();
	intern_params.n_inequality = (int)constraints.inequality.size();
	intern_params.n_planar = (int)constraints.planar.size();
	intern_params.n_tangent = (int)constraints.tangent.size();
	// Total number of constraints ...
	intern_params.n_constraints =
		_n_increment_pairs + 3 * intern_params.n_planar + intern_params.n_tangent;
	// Total number of equality constraints
	if (parameters.use_restricted_range)
		intern_params.restricted_range = true;
	else {
		intern_params.n_equality = _n_interface_pairs +
			3 * intern_params.n_planar +
			intern_params.n_tangent;
		intern_params.n_inequality =
			_n_sequenced_interface_pairs + _n_sequenced_inequality_pairs;
	}

	// polynomial parameters ...
	intern_params.poly_term = false;
	intern_params.modified_basis = true;
	intern_params.problem_type = Parameter_Types::Quadratic;

	int m = parameters.polynomial_order + 1;
	intern_params.n_poly_terms = (int)(m * (m + 1) * (m + 2) / 6);
}

void Stratigraphic_Surfaces::process_input_data() {
	if (!get_interface_data())
		throw GRBF_Exceptions::no_iterface_data;
	if (!_get_increment_pairs())
		throw GRBF_Exceptions::no_interface_increment_pairs;
	if (!check_input_data())
		throw GRBF_Exceptions::invalid_input_data;

	if (parameters.use_restricted_range) {
		for (int j = 0; j < (int)constraints.planar.size(); j++) {
			constraints.planar[j].setNormalBounds(
				parameters.angular_uncertainty,
				parameters.angular_uncertainty / 2);
			// Need more ROBUST METHOD. Try large statistical sampling
			// from von Mises spherical distribution
			std::cout << " Planar[" << j << "] Bounds: " << std::endl;
			std::cout << "	nx: " << constraints.planar[j].nx_lower_bound()
				<< " <= " << constraints.planar[j].nx()
				<< " <= " << constraints.planar[j].nx_upper_bound() << std::endl;
			std::cout << "	ny: " << constraints.planar[j].ny_lower_bound()
				<< " <= " << constraints.planar[j].ny()
				<< " <= " << constraints.planar[j].ny_upper_bound() << std::endl;
			std::cout << "	nz: " << constraints.planar[j].nz_lower_bound()
				<< " <= " << constraints.planar[j].nz()
				<< " <= " << constraints.planar[j].nz_upper_bound() << std::endl;
		}
		for (int j = 0; j < (int)constraints.tangent.size(); j++) {
			constraints.tangent[j].setAngleBounds(parameters.angular_uncertainty);
			std::cout << " Tangent[" << j << "] Bounds: " << std::endl;
			std::cout << "	" << constraints.tangent[j].angle_lower_bound()
				<< " <= " << constraints.tangent[j].inner_product_constraint()
				<< " <= " << constraints.tangent[j].angle_upper_bound()
				<< std::endl;
		}
	}
}

bool Stratigraphic_Surfaces::get_equality_values(VectorXd &equality_values) {
	int j = 0;
	int k = 0;
	int l = 0;
	int m = 0;

	// B/c the _increment_pairs list contains pairs for the sequence of the
	// stratigraphy and the inequality points we need to know the offset where
	// is
	// oncontact increment pairs
	int offset = _n_sequenced_interface_pairs + _n_sequenced_inequality_pairs;

	for (j = 0; j < _n_interface_pairs; j++)
		equality_values(j) = _increment_pairs[j + offset][0].level() -
		_increment_pairs[j + offset][1].level();
	for (k = 0; k < (int)constraints.planar.size(); k++) {
		equality_values(3 * k + j) = constraints.planar[k].nx();
		equality_values(3 * k + j + 1) = constraints.planar[k].ny();
		equality_values(3 * k + j + 2) = constraints.planar[k].nz();
	}
	for (l = 0; l < (int)constraints.tangent.size(); l++)
		equality_values(l + 3 * k + j) = constraints.tangent[l].inner_product_constraint();
	if (intern_params.poly_term)
		for (m = 0; m < (int)intern_params.n_poly_terms; m++)
			equality_values(m + l + 3 * k + j) = 0.0;

	return true;
}

bool Stratigraphic_Surfaces::get_inequality_values(VectorXd &inequality_values) {
	int j = 0;
	int k = 0;
	for (j = 0; j < _n_sequenced_interface_pairs; j++)
		inequality_values(j) = parameters.min_stratigraphic_thickness;
	for (k = 0; k < _n_sequenced_inequality_pairs; k++)
		inequality_values(k + j) = 0.0;

	return true;
}

bool Stratigraphic_Surfaces::get_inequality_values(VectorXd &b, VectorXd &r) {
	// compute largest distance between points
	// this distance will represent the upper
	// bound for inequality constraints
	// NOTE: this notion depends on the norm of gradient of the scalar field
	// to be ~1. Which is not true in reality. Will affect results.
	// a possible direction for future work
	std::vector<Point> aggregated_pts = convert_constraints_to_points(constraints);
	double distance = get_largest_distance_between_points(aggregated_pts);
	// this could be dangerous when combined with greedy (expensive compute with dense grids)
	int j = 0;
	int k = 0;
	// Sequenced Interface Points 1st
	for (j = 0; j < _n_sequenced_interface_pairs; j++) {
		b(j) = parameters.min_stratigraphic_thickness;
		r(j) = distance - parameters.min_stratigraphic_thickness;
	}
	// Sequenced Inequality Points 2nd
	for (k = 0; k < _n_sequenced_inequality_pairs; k++) {
		b(k + j) = 0.0;
		r(k + j) = parameters.min_stratigraphic_thickness;
	}
	// Interface increment pairs
	for (int l = 0; l < _n_interface_pairs; l++) {
		b(k + j + l) = -parameters.interface_uncertainty;
		r(k + j + l) = 2 * parameters.interface_uncertainty;
	}

	int n_ip = _n_increment_pairs;
	int n_p = intern_params.n_planar;
	int n_t = intern_params.n_tangent;

	// planar data
	for (int j = 0; j < n_p; j++) {
		// x-component
		b(n_ip + 3 * j + 0) = constraints.planar[j].nx_lower_bound();
		r(n_ip + 3 * j + 0) = constraints.planar[j].nx_upper_bound() -
			constraints.planar[j].nx_lower_bound();
		// y-component
		b(n_ip + 3 * j + 1) = constraints.planar[j].ny_lower_bound();
		r(n_ip + 3 * j + 1) = constraints.planar[j].ny_upper_bound() -
			constraints.planar[j].ny_lower_bound();
		// z-component
		b(n_ip + 3 * j + 2) = constraints.planar[j].nz_lower_bound();
		r(n_ip + 3 * j + 2) = constraints.planar[j].nz_upper_bound() -
			constraints.planar[j].nz_lower_bound();
	}

	// tangent data
	for (int j = 0; j < n_t; j++) {
		b(n_ip + 3 * n_p + j) = constraints.tangent[j].angle_lower_bound();
		r(n_ip + 3 * n_p + j) = constraints.tangent[j].angle_upper_bound() -
			constraints.tangent[j].angle_lower_bound();
	}

	return true;
}

Polynomial_Basis *Stratigraphic_Surfaces::create_polynomial_basis(const int &poly_order) {
	return new Poly_First(true);
}

bool Stratigraphic_Surfaces::get_interpolation_matrix(MatrixXd &interpolation_matrix) {
	int n_ip = _n_increment_pairs;
	int n_p = intern_params.n_planar;
	int n_t = intern_params.n_tangent;

	// Row and Column constraint order : increment pair (ip) -> planar
	// (p_x,p_y,p_z) -> tangent (t) Note ip encapsulates 3 types of increment
	// constraints : 1) sequenced stratigraphic reference increments 2)
	// sequenced
	// stratigraphic inequality increments 3) interface increments
	// _increment_pair[] container must be built in the above order (e.g.
	// 1),2),3))

	// Base Matrix Structure
	// |  ip/ip  ip/p_y  ip/p_y  ip/p_z  ip/t |
	// | p_x/ip p_x/p_x p_x/p_y p_x/p_z p_x/t |
	// | p_y/ip p_y/p_x p_y/p_y p_y/p_z p_y/t |
	// | p_z/ip p_z/p_x p_z/p_y p_z/p_z p_z/t |
	// |   t/ip   t/p_x   t/p_y   t/p_z   t/t |

	// Interface Increment Pair Constraints:
	for (int j = 0; j < (int)_increment_pairs.size(); j++) {
		// Row:interface increment pair/Column:interface increment pair block
		for (int k = 0; k < (int)_increment_pairs.size(); k++) {
			kernel->set_points(_increment_pairs[j][0], _increment_pairs[k][0]);
			double v1 = kernel->basis_pt_pt();
			kernel->set_points(_increment_pairs[j][0], _increment_pairs[k][1]);
			double v2 = kernel->basis_pt_pt();
			kernel->set_points(_increment_pairs[j][1], _increment_pairs[k][0]);
			double v3 = kernel->basis_pt_pt();
			kernel->set_points(_increment_pairs[j][1], _increment_pairs[k][1]);
			double v4 = kernel->basis_pt_pt();
			interpolation_matrix(j, k) = (v1 - v2) - (v3 - v4);
		}
		// Row:interface increment pair/Column:planar block
		for (int k = 0; k < n_p; k++) {
			kernel->set_points(_increment_pairs[j][0], constraints.planar[k]);
			double v1x = kernel->basis_pt_planar_x();
			double v1y = kernel->basis_pt_planar_y();
			double v1z = kernel->basis_pt_planar_z();
			kernel->set_points(_increment_pairs[j][1], constraints.planar[k]);
			double v2x = kernel->basis_pt_planar_x();
			double v2y = kernel->basis_pt_planar_y();
			double v2z = kernel->basis_pt_planar_z();
			interpolation_matrix(j, 3 * k + n_ip) = v1x - v2x;
			interpolation_matrix(j, 3 * k + n_ip + 1) = v1y - v2y;
			interpolation_matrix(j, 3 * k + n_ip + 2) = v1z - v2z;
		}
		// Row:interface increment pair/Column:tangent block
		for (int k = 0; k < n_t; k++) {
			kernel->set_points(_increment_pairs[j][0], constraints.tangent[k]);
			double v1 = kernel->basis_pt_tangent();
			kernel->set_points(_increment_pairs[j][1], constraints.tangent[k]);
			double v2 = kernel->basis_pt_tangent();
			interpolation_matrix(j, n_ip + 3 * n_p + k) = v1 - v2;
		}
	}
	// Planar Constraints
	for (int j = 0; j < n_p; j++) {
		// Row:planar/Column:interface increment pair
		for (int k = 0; k < (int)_increment_pairs.size(); k++) {
			kernel->set_points(constraints.planar[j], _increment_pairs[k][0]);
			double v1x = kernel->basis_planar_x_pt();
			double v1y = kernel->basis_planar_y_pt();
			double v1z = kernel->basis_planar_z_pt();
			kernel->set_points(constraints.planar[j], _increment_pairs[k][1]);
			double v2x = kernel->basis_planar_x_pt();
			double v2y = kernel->basis_planar_y_pt();
			double v2z = kernel->basis_planar_z_pt();
			interpolation_matrix(3 * j + n_ip, k) = v1x - v2x;
			interpolation_matrix(3 * j + n_ip + 1, k) = v1y - v2y;
			interpolation_matrix(3 * j + n_ip + 2, k) = v1z - v2z;
		}
		// Row:planar/Column:planar block
		for (int k = 0; k < n_p; k++) {
			kernel->set_points(constraints.planar[j], constraints.planar[k]);
			interpolation_matrix(3 * j + n_ip, 3 * k + n_ip) =
				kernel->basis_planar_planar(Parameter_Types::DXDX);
			interpolation_matrix(3 * j + n_ip, 3 * k + n_ip + 1) =
				kernel->basis_planar_planar(Parameter_Types::DXDY);
			interpolation_matrix(3 * j + n_ip, 3 * k + n_ip + 2) =
				kernel->basis_planar_planar(Parameter_Types::DXDZ);
			interpolation_matrix(3 * j + n_ip + 1, 3 * k + n_ip) =
				kernel->basis_planar_planar(Parameter_Types::DYDX);
			interpolation_matrix(3 * j + n_ip + 1, 3 * k + n_ip + 1) =
				kernel->basis_planar_planar(Parameter_Types::DYDY);
			interpolation_matrix(3 * j + n_ip + 1, 3 * k + n_ip + 2) =
				kernel->basis_planar_planar(Parameter_Types::DYDZ);
			interpolation_matrix(3 * j + n_ip + 2, 3 * k + n_ip) =
				kernel->basis_planar_planar(Parameter_Types::DZDX);
			interpolation_matrix(3 * j + n_ip + 2, 3 * k + n_ip + 1) =
				kernel->basis_planar_planar(Parameter_Types::DZDY);
			interpolation_matrix(3 * j + n_ip + 2, 3 * k + n_ip + 2) =
				kernel->basis_planar_planar(Parameter_Types::DZDZ);
		}
		// Row:planar/Column:tangent block
		for (int k = 0; k < n_t; k++) {
			kernel->set_points(constraints.planar[j], constraints.tangent[k]);
			interpolation_matrix(3 * j + n_ip, n_ip + 3 * n_p + k) =
				kernel->basis_planar_tangent(Parameter_Types::DX);
			interpolation_matrix(3 * j + n_ip + 1, n_ip + 3 * n_p + k) =
				kernel->basis_planar_tangent(Parameter_Types::DY);
			interpolation_matrix(3 * j + n_ip + 2, n_ip + 3 * n_p + k) =
				kernel->basis_planar_tangent(Parameter_Types::DZ);
		}
	}
	// Tangent Constraints
	for (int j = 0; j < n_t; j++) {
		// Row:tangent/Column:interface increment pair block
		for (int k = 0; k < (int)_increment_pairs.size(); k++) {
			kernel->set_points(constraints.tangent[j], _increment_pairs[k][0]);
			double v1 = kernel->basis_tangent_pt();
			kernel->set_points(constraints.tangent[j], _increment_pairs[k][1]);
			double v2 = kernel->basis_tangent_pt();
			interpolation_matrix(j + n_ip + 3 * n_p, k) = v1 - v2;
		}
		// Row:tangent/Column:planar block
		for (int k = 0; k < n_p; k++) {
			kernel->set_points(constraints.tangent[j], constraints.planar[k]);
			interpolation_matrix(j + n_ip + 3 * n_p, 3 * k + n_ip) =
				kernel->basis_tangent_planar(Parameter_Types::DX);
			interpolation_matrix(j + n_ip + 3 * n_p, 3 * k + n_ip + 1) =
				kernel->basis_tangent_planar(Parameter_Types::DY);
			interpolation_matrix(j + n_ip + 3 * n_p, 3 * k + n_ip + 2) =
				kernel->basis_tangent_planar(Parameter_Types::DZ);
		}
		// Row:tangent/Column:tangent block
		for (int k = 0; k < n_t; k++) {
			kernel->set_points(constraints.tangent[j], constraints.tangent[k]);
			interpolation_matrix(j + n_ip + 3 * n_p, n_ip + 3 * n_p + k) =
				kernel->basis_tangent_tangent();
		}
	}

	// build polynomial blocks if required
	// | A PT |
	// | P 0  |
	if (intern_params.poly_term) {
		MatrixXd poly_matrix(intern_params.n_poly_terms, intern_params.n_constraints);
		if (!_get_polynomial_matrix_block(poly_matrix)) return false;
		// fill remaining matrix blocks (P, PT, 0)
		if (!_insert_polynomial_matrix_blocks_in_interpolation_matrix(
			poly_matrix, interpolation_matrix))
			return false;
	}

	return true;
}

bool Stratigraphic_Surfaces::get_inequality_matrix(
	const MatrixXd &interpolation_matrix, MatrixXd &inequality_matrix)
{
	if (inequality_matrix.rows() == 0 ||
		inequality_matrix.rows() > interpolation_matrix.rows() ||
		inequality_matrix.cols() != interpolation_matrix.cols())
		return false;

	for (int j = 0; j < (int)inequality_matrix.rows(); j++) {
		for (int k = 0; k < (int)inequality_matrix.cols(); k++) {
			inequality_matrix(j, k) = interpolation_matrix(j, k);
		}
	}

	return true;
}

void Stratigraphic_Surfaces::setup_system_solver() {
	// only way to solve this problem is via a quadratic optimization problem
	int n_ie = intern_params.n_inequality;
	int n_e = intern_params.n_equality;
	int n_c = intern_params.n_constraints;

	if (intern_params.restricted_range) {
		VectorXd b(n_c);
		VectorXd r(n_c);
		get_inequality_values(b, r);

		MatrixXd interpolation_matrix(n_c, n_c);
		if (!get_interpolation_matrix(interpolation_matrix))
			throw GRBF_Exceptions::error_computing_interpolation_matrix;

		MatrixXd inequality_matrix(n_c, n_c);
		inequality_matrix = interpolation_matrix;

		Quadratic_Predictor_Corrector_LOQO *qpc =
			new Quadratic_Predictor_Corrector_LOQO(interpolation_matrix, inequality_matrix, b, r);
		if (!qpc->solve())
			throw GRBF_Exceptions::loqo_quadratic_solver_failure;
		solver = qpc;
	}
	else {
		VectorXd inequality_values(n_ie);
		get_inequality_values(inequality_values);

		VectorXd equality_values(n_e);
		get_equality_values(equality_values);

		MatrixXd interpolation_matrix(n_c, n_c);
		if (!get_interpolation_matrix(interpolation_matrix))
			throw GRBF_Exceptions::error_computing_interpolation_matrix;

		MatrixXd inequality_matrix(n_ie, n_c);
		if (!get_inequality_matrix(interpolation_matrix, inequality_matrix))
			throw GRBF_Exceptions::error_computing_inequality_vector;

		MatrixXd equality_matrix(n_e, n_c);
		if (!get_equality_matrix(interpolation_matrix, equality_matrix))
			throw GRBF_Exceptions::error_computing_equality_vector;

		Quadratic_Predictor_Corrector *qpc = new Quadratic_Predictor_Corrector(
			interpolation_matrix, equality_matrix, inequality_matrix,
			equality_values, inequality_values);
		if (!qpc->solve())
			throw GRBF_Exceptions::pc_quadratic_solver_failure;
		solver = qpc;
	}

	if (!_update_interface_iso_values())
		throw GRBF_Exceptions::error_updating_interface_iso_values;
}

bool Stratigraphic_Surfaces::convert_modified_kernel_to_rbf_kernel() {
	if (!rbf_kernel || !kernel) return false;

	// prep for linear prob...
	// set the constraints...
	// for Lajaunie and Stratigraphic Surface methods we don't update the
	// itrface[].level()'s we update the _increment_pairs[][] level()'s
	for (auto &increment_pair : _increment_pairs) {
		eval_scalar_interpolant_at_point(increment_pair[0]);
		eval_scalar_interpolant_at_point(increment_pair[1]);
		increment_pair[0].setLevel(increment_pair[0].scalar_field());
		increment_pair[1].setLevel(increment_pair[1].scalar_field());
	}
	for (auto &planar_pt : constraints.planar) {
		eval_vector_interpolant_at_point(planar_pt);
		double normal[3] = { planar_pt.nx_interp(),
							planar_pt.ny_interp(),
							planar_pt.nz_interp() };
		planar_pt.setNormal(normal[0], normal[1], normal[2]);
	}
	for (auto &tangent_pt : constraints.tangent) {
		eval_vector_interpolant_at_point(tangent_pt);
		double vf[3] = { tangent_pt.nx_interp(),
						tangent_pt.ny_interp(),
						tangent_pt.nz_interp() };
		double inner_product = vf[0] * tangent_pt.tx() +
			vf[1] * tangent_pt.ty() +
			vf[2] * tangent_pt.tz();
		tangent_pt.setInnerProductConstraint(inner_product);
	}

	// switch from modified kernel to normal rbf kernel
	kernel = rbf_kernel;

	if (parameters.use_restricted_range)
		intern_params.restricted_range = false;
	intern_params.n_equality =
		_n_increment_pairs + 3 * intern_params.n_planar + intern_params.n_tangent;
	_n_interface_pairs = _n_increment_pairs;  // hack b/c defn of
	// Stratigraphic_Surfaces::get_equality_values()
	_n_sequenced_interface_pairs = 0;
	_n_sequenced_inequality_pairs = 0;
	intern_params.poly_term = true;
	intern_params.n_poly_terms = 3;
	intern_params.modified_basis = false;
	intern_params.problem_type = Parameter_Types::Linear;
	int n_e = intern_params.n_equality;
	int n_p = intern_params.n_poly_terms;
	VectorXd equality_values(n_e + n_p);
	get_equality_values(equality_values);

	MatrixXd interpolation_matrix(n_e + n_p, n_e + n_p);
	if (!get_interpolation_matrix(interpolation_matrix)) return false;

	Linear_LU_decomposition *llu =
		new Linear_LU_decomposition(interpolation_matrix, equality_values);
	if (!llu->solve()) return false;
	solver = llu;

	if (!_update_interface_iso_values()) return false;

	return true;
}

void Stratigraphic_Surfaces::eval_scalar_interpolant_at_point(Point &p) {
	int n_ip = _n_increment_pairs;
	int n_p = intern_params.n_planar;
	int n_t = intern_params.n_tangent;

	Kernel *kernel_j = kernel->clone();
	double elemsum_1 = 0.0;
	double elemsum_2 = 0.0;
	double elemsum_3 = 0.0;
	double poly = 0.0;
	for (int k = 0; k < (int)_increment_pairs.size(); k++) {
		kernel_j->set_points(p, _increment_pairs[k][0]);
		double v1 = kernel_j->basis_pt_pt();
		kernel_j->set_points(p, _increment_pairs[k][1]);
		double v2 = kernel_j->basis_pt_pt();
		elemsum_1 += solver->weights[k] * (v1 - v2);
	}
	for (int k = 0; k < n_p; k++) {
		kernel_j->set_points(p, constraints.planar[k]);
		elemsum_2 +=
			solver->weights[n_ip + 3 * k] * kernel_j->basis_pt_planar_x();
		elemsum_2 +=
			solver->weights[n_ip + 3 * k + 1] * kernel_j->basis_pt_planar_y();
		elemsum_2 +=
			solver->weights[n_ip + 3 * k + 2] * kernel_j->basis_pt_planar_z();
	}
	for (int k = 0; k < n_t; k++) {
		kernel_j->set_points(p, constraints.tangent[k]);
		elemsum_3 +=
			solver->weights[n_ip + 3 * n_p + k] * kernel_j->basis_pt_tangent();
	}
	if (intern_params.poly_term) {
		Polynomial_Basis *p_basis_j = p_basis->clone();
		p_basis_j->set_point(p);
		VectorXd b = p_basis_j->basis();
		for (int k = 0; k < (int)b.size(); k++)
			poly += b(k) * solver->weights[n_ip + 3 * n_p + n_t + k];
		delete p_basis_j;
	}
	p.set_scalar_field(elemsum_1 + elemsum_2 + elemsum_3 + poly);
	delete kernel_j;
}

void Stratigraphic_Surfaces::eval_vector_interpolant_at_point(Point &p) {
	int n_ip = _n_increment_pairs;
	int n_p = intern_params.n_planar;
	int n_t = intern_params.n_tangent;

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
	for (int k = 0; k < n_ip; k++) {
		kernel->set_points(p, _increment_pairs[k][0]);
		double v1x = kernel->basis_planar_x_pt();
		double v1y = kernel->basis_planar_y_pt();
		double v1z = kernel->basis_planar_z_pt();
		kernel->set_points(p, _increment_pairs[k][1]);
		double v2x = kernel->basis_planar_x_pt();
		double v2y = kernel->basis_planar_y_pt();
		double v2z = kernel->basis_planar_z_pt();
		elemsum_1_x += solver->weights[k] * (v1x - v2x);
		elemsum_1_y += solver->weights[k] * (v1y - v2y);
		elemsum_1_z += solver->weights[k] * (v1z - v2z);
	}
	// planar constraints
	for (int k = 0; k < n_p; k++) {
		kernel->set_points(p, constraints.planar[k]);
		elemsum_2_x += solver->weights[n_ip + 0 + 3 * k] *
			kernel->basis_planar_planar(Parameter_Types::DXDX);
		elemsum_2_x += solver->weights[n_ip + 1 + 3 * k] *
			kernel->basis_planar_planar(Parameter_Types::DXDY);
		elemsum_2_x += solver->weights[n_ip + 2 + 3 * k] *
			kernel->basis_planar_planar(Parameter_Types::DXDZ);
		elemsum_2_y += solver->weights[n_ip + 0 + 3 * k] *
			kernel->basis_planar_planar(Parameter_Types::DYDX);
		elemsum_2_y += solver->weights[n_ip + 1 + 3 * k] *
			kernel->basis_planar_planar(Parameter_Types::DYDY);
		elemsum_2_y += solver->weights[n_ip + 2 + 3 * k] *
			kernel->basis_planar_planar(Parameter_Types::DYDZ);
		elemsum_2_z += solver->weights[n_ip + 0 + 3 * k] *
			kernel->basis_planar_planar(Parameter_Types::DZDX);
		elemsum_2_z += solver->weights[n_ip + 1 + 3 * k] *
			kernel->basis_planar_planar(Parameter_Types::DZDY);
		elemsum_2_z += solver->weights[n_ip + 2 + 3 * k] *
			kernel->basis_planar_planar(Parameter_Types::DZDZ);
	}
	// Tangent constraints
	for (int k = 0; k < n_t; k++) {
		kernel->set_points(p, constraints.tangent[k]);
		elemsum_3_x += solver->weights[n_ip + 3 * n_p + k] *
			kernel->basis_planar_tangent(Parameter_Types::DX);
		elemsum_3_y += solver->weights[n_ip + 3 * n_p + k] *
			kernel->basis_planar_tangent(Parameter_Types::DY);
		elemsum_3_z += solver->weights[n_ip + 3 * n_p + k] *
			kernel->basis_planar_tangent(Parameter_Types::DZ);
	}
	if (intern_params.poly_term) {
		Polynomial_Basis *p_basis_j = p_basis->clone();
		p_basis_j->set_point(p);
		VectorXd bx = p_basis_j->dx();
		VectorXd by = p_basis_j->dy();
		VectorXd bz = p_basis_j->dz();
		for (int k = 0; k < (int)bx.size(); k++) {
			poly_x += bx(k) * solver->weights[n_ip + 3 * n_p + n_t + k];
			poly_y += by(k) * solver->weights[n_ip + 3 * n_p + n_t + k];
			poly_z += bz(k) * solver->weights[n_ip + 3 * n_p + n_t + k];
		}
		delete p_basis_j;
	}
	double nx = elemsum_1_x + elemsum_2_x + elemsum_3_x + poly_x;
	double ny = elemsum_1_y + elemsum_2_y + elemsum_3_y + poly_y;
	double nz = elemsum_1_z + elemsum_2_z + elemsum_3_z + poly_z;
	p.set_vector_field(nx, ny, nz);
	delete kernel_j;
}