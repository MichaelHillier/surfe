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
#include <single_surface.h>
#include <vector>

#include <fstream>
#include <iomanip>
#include <iostream>

bool Single_Surface::_get_polynomial_matrix_block(MatrixXd &poly_matrix) {
	int n_ie = intern_params.n_inequality;
	int n_i = intern_params.n_interface;
	int n_p = intern_params.n_planar;
	int n_t = intern_params.n_tangent;

	int nl = n_i + n_ie;  // n_ie should always be zero.

	p_basis = create_polynomial_basis(parameters.polynomial_order);

	if ((int)poly_matrix.rows() != intern_params.n_poly_terms) return false;
	// for interface points ...
	for (int j = 0; j < nl; j++) {
		p_basis->set_point(constraints.itrface[j]);
		VectorXd b = p_basis->basis();
		if ((int)b.rows() != intern_params.n_poly_terms) return false;
		for (int k = 0; k < (int)b.rows(); k++) poly_matrix(k, j) = b(k);
	}
	// for planar points ...
	for (int j = 0; j < n_p; j++) {
		p_basis->set_point(constraints.planar[j]);
		VectorXd bx = p_basis->dx();
		VectorXd by = p_basis->dy();
		VectorXd bz = p_basis->dz();
		if ((int)bx.rows() != intern_params.n_poly_terms) return false;
		for (int k = 0; k < (int)bx.rows(); k++) {
			poly_matrix(k, nl + 3 * j) = bx(k);
			poly_matrix(k, nl + 3 * j + 1) = by(k);
			poly_matrix(k, nl + 3 * j + 2) = bz(k);
		}
	}
	// for tangent points ...
	for (int j = 0; j < n_t; j++) {
		p_basis->set_point(constraints.tangent[j]);
		VectorXd bx = p_basis->dx();
		VectorXd by = p_basis->dy();
		VectorXd bz = p_basis->dz();
		if ((int)bx.rows() != intern_params.n_poly_terms) return false;
		for (int k = 0; k < (int)bx.rows(); k++) {
			poly_matrix(k, nl + 3 * n_p + j) =
				constraints.tangent[j].tx() * bx(k) +
				constraints.tangent[j].ty() * by(k) +
				constraints.tangent[j].tz() * bz(k);
		}
	}

	return true;
}

bool Single_Surface::_insert_polynomial_matrix_blocks_in_interpolation_matrix(
	const MatrixXd &poly_matrix, MatrixXd &interpolation_matrix) {
	int n_ie = intern_params.n_inequality;
	int n_i = intern_params.n_interface;
	int n_p = intern_params.n_planar;
	int n_t = intern_params.n_tangent;

	// build polynomial blocks
	// | A PT |
	// | P 0  |
	// start with P
	for (int j = 0; j < (int)poly_matrix.rows(); j++) {
		for (int k = 0; k < (int)poly_matrix.cols(); k++) {
			interpolation_matrix(n_ie + n_i + 3 * n_p + n_t + j, k) = poly_matrix(j, k);
			interpolation_matrix(k, n_ie + n_i + 3 * n_p + n_t + j) =
				interpolation_matrix(n_ie + n_i + 3 * n_p + n_t + j, k);
		}
	}

	for (int j = 0; j < (int)poly_matrix.rows(); j++) {
		for (int k = 0; k < (int)poly_matrix.rows(); k++) {
			interpolation_matrix(n_ie + n_i + 3 * n_p + n_t + j,
				n_ie + n_i + 3 * n_p + n_t + k) = 0;
		}
	}

	return true;
}

Single_Surface::Single_Surface(const Parameters& mparams)
{
	// set GUI parameters
	parameters = mparams;

	solver = nullptr;
	kernel = nullptr;
	rbf_kernel = nullptr;
	p_basis = nullptr;

	_iteration = 0;
}

Single_Surface::Single_Surface()
{
	solver = nullptr;
	kernel = nullptr;
	rbf_kernel = nullptr;
	p_basis = nullptr;

	_iteration = 0;
}

Polynomial_Basis *Single_Surface::create_polynomial_basis(
	const int &poly_order) {
	if (poly_order == 0)
		return new Poly_Zero;
	else if (poly_order == 1)
		return new Poly_First;
	else
		return new Poly_Second;
}

void Single_Surface::get_method_parameters() {
	// # of constraints for each constraint type ...
	intern_params.n_interface = (int)constraints.itrface.size();
	intern_params.n_inequality = (int)constraints.inequality.size();
	intern_params.n_planar = (int)constraints.planar.size();
	intern_params.n_tangent = (int)constraints.tangent.size();
	// Total number of constraints ...
	intern_params.n_constraints =
		intern_params.n_interface + intern_params.n_inequality +
		3 * intern_params.n_planar + intern_params.n_tangent;
	// Total number of equality constraints
	if (parameters.use_restricted_range)
		intern_params.restricted_range = true;
	else
		intern_params.n_equality = intern_params.n_interface +
		3 * intern_params.n_planar +
		intern_params.n_tangent;

	// polynomial parameters ...
	if (intern_params.n_inequality != 0 ||
		parameters.use_restricted_range != 0) {
		intern_params.poly_term = false;
		intern_params.modified_basis = true;
		intern_params.problem_type = Parameter_Types::Quadratic;
	}
	else {
		intern_params.poly_term = true;  // NOTE: May want to have this as an option when using SPD functions
		intern_params.modified_basis = false;
		intern_params.problem_type = Parameter_Types::Linear;
	}

	int m = parameters.polynomial_order + 1;
	intern_params.n_poly_terms = (int)(m * (m + 1) * (m + 2) / 6);  // for 3D only...
}

void Single_Surface::setup_system_solver() {
	// the type of mathematical solver to be used can depends on the following
	// 1) Linear or Quadratic system
	// 2) What type of RBF is used
	// 3) Smoothing -> maybe use least squares / right now we are doing matrix
	// smoothing

	int n_ie = intern_params.n_inequality;
	int n_e = intern_params.n_equality;
	int n_c = intern_params.n_constraints;
	if (intern_params.problem_type == Parameter_Types::Quadratic) {
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

			Quadratic_Predictor_Corrector *qpc =
				new Quadratic_Predictor_Corrector(
					interpolation_matrix, equality_matrix, inequality_matrix,
					equality_values, inequality_values);
			if (!qpc->solve())
				throw GRBF_Exceptions::pc_quadratic_solver_failure;
			solver = qpc;
		}
	}
	else  // Linear
	{
		int n_p = intern_params.n_poly_terms;
		VectorXd equality_values(n_e + n_p);
		get_equality_values(equality_values);

		MatrixXd interpolation_matrix(n_e + n_p, n_e + n_p);
		if (!get_interpolation_matrix(interpolation_matrix))
			throw GRBF_Exceptions::error_computing_interpolation_matrix;

		Linear_LU_decomposition *llu =
			new Linear_LU_decomposition(interpolation_matrix, equality_values);
		if (!llu->solve())
			throw GRBF_Exceptions::linear_solver_failure;
		solver = llu;
	}

	// check_interpolant();
}

bool Single_Surface::get_minimial_and_excluded_input(Constraints &greedy_input, Constraints &excluded_input)
{
	// get extremal interface points
	// Note:
	// if poly order = 1 find 4 interface points nicely sampling the volume
	// if poly order = 2 find 6 interface points nicely sampling the volume
	std::vector<Interface> extremal_interace_pts;
	std::vector<Point> pts(constraints.itrface.begin(), constraints.itrface.end());
	std::vector<int> interface_indices =
		get_extremal_point_data_indices_from_points(pts);
	int num_extremal_itr_pts = 4;
	if (parameters.polynomial_order == 2) num_extremal_itr_pts = 6;
	if ((int)interface_indices.size() < num_extremal_itr_pts) return false;
	for (const auto& inequality_pt : constraints.inequality)
		excluded_input.inequality.push_back(inequality_pt);
	for (int j = 0; j < (int)constraints.itrface.size(); j++) {
		bool exclude_index = true;
		for (int k = 0; k < num_extremal_itr_pts; k++) {
			if (interface_indices[k] == j) {
				exclude_index = false;
				greedy_input.itrface.push_back(constraints.itrface[interface_indices[k]]);
				break;
			}
		}
		if (exclude_index)
			excluded_input.itrface.push_back(constraints.itrface[j]);
	}
	greedy_input.planar.push_back(constraints.planar[0]);
	for (int j = 1; j < (int)constraints.planar.size(); j++)
		excluded_input.planar.push_back(constraints.planar[j]);
	for (const auto &tangent_pt : constraints.tangent)
		excluded_input.tangent.push_back(tangent_pt);

	return true;
}

bool Single_Surface::measure_residuals(Constraints &input) {
	if (solver == nullptr) return false;

#pragma omp parallel sections
	{
#pragma omp section
		{
			// inequalities points
			for (auto &inequality_pt : input.inequality) {
				eval_scalar_interpolant_at_point(inequality_pt);
				if (inequality_pt.level() >= 0) {
					if (inequality_pt.scalar_field() >= 0)
						inequality_pt.setResidual(true);
					else
						inequality_pt.setResidual(false);
				}
				else {
					if (inequality_pt.scalar_field() < 0)
						inequality_pt.setResidual(true);
					else
						inequality_pt.setResidual(false);
				}
			}
		}
#pragma omp section
		{
			// interface points
			for (auto &interface_pt : input.itrface) {
				eval_scalar_interpolant_at_point(interface_pt);
				interface_pt.setResidual(abs(interface_pt.scalar_field()));
			}
		}
#pragma omp section
		{
			// planar points
			for (auto &planar_pt : input.planar) {
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
		}
#pragma omp section
		{
			// tangent points
			for (auto &tangent_pt : input.tangent) {
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
				Math_methods::angle_btw_2_vectors(v1, v2, angle);
				tangent_pt.setResidual(angle);
			}
		}
	}

	return true;
}

bool Single_Surface::append_greedy_input(Constraints &input) {
	// This function can most likely be promoted the the Greedy parent class
	// Below section can be a lot of computations - leverage parallelism
	std::vector<int> planar_indices_to_include;      // PLANAR Observations
	std::vector<int> tangent_indices_to_include;     // TANGENT Observations
	std::vector<int> interface_indices_to_include;   // INTERFACE Observations
	std::vector<int> inequality_indices_to_include;  // INEQUALITIES
													 // Observations
	// For the first iteration only consider adding additional planar
	// constraints
	// These additional constraints could force large changes in the scalar
	// field
	// which could consequently pass closely to interface points not yet
	// included
	if (_iteration == 0)
		planar_indices_to_include =
		Get_Planar_STL_Vector_Indices_With_Large_Residuals(
			input.planar, parameters.angular_uncertainty,
			constraints.GetPlanarAvgNNDist());
	else {
#pragma omp parallel sections
		{
#pragma omp section
			{
				// PLANAR Observations
				planar_indices_to_include =
					Get_Planar_STL_Vector_Indices_With_Large_Residuals(
						input.planar, parameters.angular_uncertainty,
						constraints.GetPlanarAvgNNDist());
			}
#pragma omp section
			{
				// TANGENT Observations
				tangent_indices_to_include =
					Get_Tangent_STL_Vector_Indices_With_Large_Residuals(
						input.tangent, parameters.angular_uncertainty,
						constraints.GetPlanarAvgNNDist());
			}
#pragma omp section
			{
				// INTERFACE Observations
				interface_indices_to_include =
					Get_Interface_STL_Vector_Indices_With_Large_Residuals(
						input.itrface, parameters.interface_uncertainty,
						constraints.GetInterfaceAvgNNDist());
			}
#pragma omp section
			{
				// INEQUALITIES Observations
				inequality_indices_to_include =
					Get_Inequality_STL_Vector_Indices_With_Large_Residuals(
						input.inequality, constraints.GetInequalityAvgNNDist());
			}
		}
	}

	int pI2i = (int)planar_indices_to_include.size();
	int tI2i = (int)tangent_indices_to_include.size();
	int itrI2i = (int)interface_indices_to_include.size();
	int ieI2i = (int)inequality_indices_to_include.size();

	// Add to current greedy input
	for (int j = 0; j < pI2i; j++)
		this->constraints.planar.push_back(
			input.planar[planar_indices_to_include[j]]);
	for (int j = 0; j < tI2i; j++)
		this->constraints.tangent.push_back(
			input.tangent[tangent_indices_to_include[j]]);
	for (int j = 0; j < itrI2i; j++)
		this->constraints.itrface.push_back(
			input.itrface[interface_indices_to_include[j]]);
	for (int j = 0; j < ieI2i; j++)
		this->constraints.inequality.push_back(
			input.inequality[inequality_indices_to_include[j]]);

	// Remove data from input so that residuals do not have to be measured at
	// locations where the constraints are supplied
	for (int j = 0; j < pI2i; j++) {
		input.planar.erase(input.planar.begin() + planar_indices_to_include[j]);
		for (int k = j; k < pI2i; k++) planar_indices_to_include[k]--;
	}
	for (int j = 0; j < tI2i; j++) {
		input.tangent.erase(input.tangent.begin() + tangent_indices_to_include[j]);
		for (int k = j; k < tI2i; k++) tangent_indices_to_include[k]--;
	}
	for (int j = 0; j < itrI2i; j++) {
		input.itrface.erase(input.itrface.begin() + interface_indices_to_include[j]);
		for (int k = j; k < itrI2i; k++) interface_indices_to_include[k]--;
	}
	for (int j = 0; j < ieI2i; j++) {
		input.inequality.erase(input.inequality.begin() + inequality_indices_to_include[j]);
		for (int k = j; k < ieI2i; k++) inequality_indices_to_include[k]--;
	}

	if (pI2i != 0 || tI2i != 0 || itrI2i != 0 || ieI2i != 0)
		return true;
	else
		return false;
}

bool Single_Surface::convert_modified_kernel_to_rbf_kernel() {
	if (!rbf_kernel || !kernel) return false;

	// prep for linear prob...
	// set the constraints...
	// temp container for interface points needed for cases where inequalities
	// are
	// used b/c inequalities are turned into 'interface points' (they do not
	// actually represent horizon interfaces points though)
	std::vector<Interface> temp_itr;
	for (const auto& inequality_pt : constraints.inequality) {
		Point test_point(inequality_pt.x(), inequality_pt.y(), inequality_pt.z());
		// debug
		// cout<<" Inequality["<<j<<"]: "<<endl;
		eval_scalar_interpolant_at_point(test_point);
		// cout<<"	Scalar field = "<<test_point.scalar_field()<<endl;
		Interface itr_point(test_point.x(), test_point.y(), test_point.z(), test_point.scalar_field());
		temp_itr.push_back(itr_point);
	}
	for (auto &interface_pt : constraints.itrface) {
		Point test_point(interface_pt.x(), interface_pt.y(), interface_pt.z());
		// debug
		// cout<<" Interface["<<j<<"]: "<<endl;
		eval_scalar_interpolant_at_point(test_point);
		// cout<<"	Scalar field = "<<test_point.scalar_field()<<endl;
		interface_pt.setLevel(test_point.scalar_field());

		Interface itr_point(test_point.x(), test_point.y(), test_point.z(), test_point.scalar_field());
		temp_itr.push_back(itr_point);
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
		// 		cout<<" Tangent["<<j<<"]: "<<endl;
		// 		cout<<"	Tx = "<<b_input.tangent->at(j).tx()<<" Nx
		// interpolated =
		// "<<b_input.tangent->at(j).nx_interp()<<endl; 		cout<<"	Ty
		// =
		// "<<b_input.tangent->at(j).ty()<<" Ny interpolated =
		// "<<b_input.tangent->at(j).ny_interp()<<endl; 		cout<<"	Tz
		// =
		// "<<b_input.tangent->at(j).tz()<<" Nz interpolated =
		// "<<b_input.tangent->at(j).nz_interp()<<endl; 		cout<<" Tx*nx
		// + Ty*ny +
		// Tz*nz =
		// "<<b_input.tangent->at(j).tx()*b_input.tangent->at(j).nx_interp()
		// + b_input.tangent->at(j).ty()*b_input.tangent->at(j).ny_interp() +
		// 			b_input.tangent->at(j).tz()*b_input.tangent->at(j).nz_interp()<<endl;
		double vf[3] = { tangent_pt.nx_interp(),tangent_pt.ny_interp(),tangent_pt.nz_interp() };
		double inner_product = vf[0] * tangent_pt.tx() +
			vf[1] * tangent_pt.ty() +
			vf[2] * tangent_pt.tz();
		tangent_pt.setInnerProductConstraint(inner_product);
	}

	// switch from modified kernel to normal rbf kernel
	kernel = rbf_kernel;

	constraints.inequality.clear();
	constraints.itrface.clear();
	constraints.itrface = temp_itr;

	int n_p = intern_params.n_poly_terms;
	if (parameters.use_restricted_range)
		intern_params.restricted_range = false;
	intern_params.n_interface = (int)constraints.itrface.size();
	intern_params.n_inequality = (int)constraints.inequality.size();
	intern_params.n_equality = intern_params.n_interface +
		3 * intern_params.n_planar +
		intern_params.n_tangent;
	intern_params.poly_term = true;
	intern_params.modified_basis = false;
	intern_params.problem_type = Parameter_Types::Linear;
	int n_e = intern_params.n_equality;
	VectorXd equality_values(n_e + n_p);
	get_equality_values(equality_values);

	MatrixXd interpolation_matrix(n_e + n_p, n_e + n_p);
	if (!get_interpolation_matrix(interpolation_matrix)) return false;

	Linear_LU_decomposition *llu =
		new Linear_LU_decomposition(interpolation_matrix, equality_values);
	if (!llu->solve()) return false;
	solver = llu;

	return true;
}

void Single_Surface::eval_scalar_interpolant_at_point(Point &p) {
	int n_ie = intern_params.n_inequality;
	int n_i = intern_params.n_interface;
	int n_p = intern_params.n_planar;
	int n_t = intern_params.n_tangent;

	Kernel *kernel_j = kernel->clone();
	double elemsum_1 = 0.0;
	double elemsum_2 = 0.0;
	double elemsum_3 = 0.0;
	double elemsum_4 = 0.0;
	double poly = 0.0;
	for (int k = 0; k < n_ie; k++) {
		kernel_j->set_points(p, constraints.inequality[k]);
		elemsum_1 += solver->weights[k] * kernel_j->basis_pt_pt();
	}
	for (int k = 0; k < n_i; k++) {
		kernel_j->set_points(p, constraints.itrface[k]);
		elemsum_2 += solver->weights[n_ie + k] * kernel_j->basis_pt_pt();
	}
	for (int k = 0; k < n_p; k++) {
		kernel_j->set_points(p, constraints.planar[k]);
		elemsum_3 += solver->weights[n_ie + n_i + 3 * k] * kernel_j->basis_pt_planar_x();
		elemsum_3 += solver->weights[n_ie + n_i + 3 * k + 1] * kernel_j->basis_pt_planar_y();
		elemsum_3 += solver->weights[n_ie + n_i + 3 * k + 2] * kernel_j->basis_pt_planar_z();
	}
	for (int k = 0; k < n_t; k++) {
		kernel_j->set_points(p, constraints.tangent[k]);
		elemsum_4 += solver->weights[n_ie + n_i + 3 * n_p + k] * kernel_j->basis_pt_tangent();
	}
	if (intern_params.poly_term) {
		Polynomial_Basis *p_basis_j = p_basis->clone();
		p_basis_j->set_point(p);
		VectorXd b = p_basis_j->basis();
		for (int k = 0; k < (int)b.size(); k++)
			poly += b(k) * solver->weights[n_ie + n_i + 3 * n_p + n_t + k];
		delete p_basis_j;
	}
	p.set_scalar_field(elemsum_1 + elemsum_2 + elemsum_3 + elemsum_4 + poly);
	delete kernel_j;
}

void Single_Surface::eval_vector_interpolant_at_point(Point &p) {
	int n_ie = intern_params.n_inequality;
	int n_i = intern_params.n_interface;
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
	// inequality constraints
	for (int k = 0; k < n_ie; k++) {
		kernel->set_points(p, constraints.inequality[k]);
		elemsum_1_x += solver->weights[k] * kernel->basis_planar_x_pt();
		elemsum_1_y += solver->weights[k] * kernel->basis_planar_y_pt();
		elemsum_1_z += solver->weights[k] * kernel->basis_planar_z_pt();
	}
	// interface constraints
	for (int k = 0; k < n_i; k++) {
		kernel->set_points(p, constraints.itrface[k]);
		elemsum_1_x += solver->weights[n_ie + k] * kernel->basis_planar_x_pt();
		elemsum_1_y += solver->weights[n_ie + k] * kernel->basis_planar_y_pt();
		elemsum_1_z += solver->weights[n_ie + k] * kernel->basis_planar_z_pt();
	}
	// normal constraints
	for (int k = 0; k < n_p; k++) {
		kernel->set_points(p, constraints.planar[k]);
		elemsum_2_x += solver->weights[n_ie + n_i + 0 + 3 * k] *
			kernel->basis_planar_planar(Parameter_Types::DXDX);
		elemsum_2_x += solver->weights[n_ie + n_i + 1 + 3 * k] *
			kernel->basis_planar_planar(Parameter_Types::DXDY);
		elemsum_2_x += solver->weights[n_ie + n_i + 2 + 3 * k] *
			kernel->basis_planar_planar(Parameter_Types::DXDZ);
		elemsum_2_y += solver->weights[n_ie + n_i + 0 + 3 * k] *
			kernel->basis_planar_planar(Parameter_Types::DYDX);
		elemsum_2_y += solver->weights[n_ie + n_i + 1 + 3 * k] *
			kernel->basis_planar_planar(Parameter_Types::DYDY);
		elemsum_2_y += solver->weights[n_ie + n_i + 2 + 3 * k] *
			kernel->basis_planar_planar(Parameter_Types::DYDZ);
		elemsum_2_z += solver->weights[n_ie + n_i + 0 + 3 * k] *
			kernel->basis_planar_planar(Parameter_Types::DZDX);
		elemsum_2_z += solver->weights[n_ie + n_i + 1 + 3 * k] *
			kernel->basis_planar_planar(Parameter_Types::DZDY);
		elemsum_2_z += solver->weights[n_ie + n_i + 2 + 3 * k] *
			kernel->basis_planar_planar(Parameter_Types::DZDZ);
	}
	// tangent constraints
	for (int k = 0; k < n_t; k++) {
		kernel->set_points(p, constraints.tangent[k]);
		elemsum_3_x += solver->weights[n_ie + n_i + 3 * n_p + k] *
			kernel->basis_planar_tangent(Parameter_Types::DX);
		elemsum_3_y += solver->weights[n_ie + n_i + 3 * n_p + k] *
			kernel->basis_planar_tangent(Parameter_Types::DY);
		elemsum_3_z += solver->weights[n_ie + n_i + 3 * n_p + k] *
			kernel->basis_planar_tangent(Parameter_Types::DZ);
	}
	if (intern_params.poly_term) {
		Polynomial_Basis *p_basis_j = p_basis->clone();
		p_basis_j->set_point(p);
		VectorXd bx = p_basis_j->dx();
		VectorXd by = p_basis_j->dy();
		VectorXd bz = p_basis_j->dz();
		for (int k = 0; k < (int)bx.size(); k++) {
			poly_x += bx(k) * solver->weights[n_ie + n_i + 3 * n_p + n_t + k];
			poly_y += by(k) * solver->weights[n_ie + n_i + 3 * n_p + n_t + k];
			poly_z += bz(k) * solver->weights[n_ie + n_i + 3 * n_p + n_t + k];
		}
		delete p_basis_j;
	}
	double nx = elemsum_1_x + elemsum_2_x + elemsum_3_x + poly_x;
	double ny = elemsum_1_y + elemsum_2_y + elemsum_3_y + poly_y;
	double nz = elemsum_1_z + elemsum_2_z + elemsum_3_z + poly_z;
	p.set_vector_field(nx, ny, nz);
	delete kernel_j;
}

bool Single_Surface::get_equality_values(VectorXd &equality_values) {
	int j = 0;
	int k = 0;
	int l = 0;
	int m = 0;

	for (j = 0; j < (int)constraints.itrface.size(); j++)
		equality_values(j) = constraints.itrface[j].level();
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

bool Single_Surface::get_inequality_matrix(const MatrixXd &interpolation_matrix,
	MatrixXd &inequality_matrix) {
	if (inequality_matrix.rows() == 0 ||
		inequality_matrix.cols() != interpolation_matrix.cols())
		return false;

	// below inequality constraints have to be put in terms of s(x) >= level
	// so if s(x) <= level => -1.0*s(x) >= -1.0*level
	// for single surface level is always 0. so, -1.0*s(x) > 0
	int n_ie = intern_params.n_inequality;
	if (n_ie != 0) {
		for (int j = 0; j < inequality_matrix.rows(); j++) {
			for (int k = 0; k < inequality_matrix.cols(); k++) {
				if (constraints.inequality[j].level() > 0)
					inequality_matrix(j, k) = interpolation_matrix(j, k);
				else
					inequality_matrix(j, k) = -interpolation_matrix(j, k);
			}
		}
	}

	if (intern_params.restricted_range) {
		int n_i = intern_params.n_interface;
		int n_p = intern_params.n_planar;
		int n_t = intern_params.n_tangent;
		for (int j = 0; j < n_i; j++) {
			for (int k = 0; k < inequality_matrix.cols(); k++) {
				inequality_matrix(n_ie + 2 * j + 0, k) = interpolation_matrix(n_ie + j, k);  // lower bound constraint Ax >= lower_bound
				inequality_matrix(n_ie + 2 * j + 1, k) = -interpolation_matrix(n_ie + j, k);  // upper bound constraint -Ax >= -upper_bound
			}
		}
		for (int j = 0; j < n_p; j++) {
			for (int k = 0; k < inequality_matrix.cols(); k++) {
				inequality_matrix(n_ie + 2 * n_i + 6 * j + 0, k) =
					interpolation_matrix(n_ie + n_i + 3 * j + 0, k);  // lower bound constraint Ax >= lower_bound
				inequality_matrix(n_ie + 2 * n_i + 6 * j + 1, k) =
					-interpolation_matrix(n_ie + n_i + 3 * j + 0, k);  // upper bound constraint -Ax >= -upper_bound
				inequality_matrix(n_ie + 2 * n_i + 6 * j + 2, k) =
					interpolation_matrix(n_ie + n_i + 3 * j + 1, k);  // ...
				inequality_matrix(n_ie + 2 * n_i + 6 * j + 3, k) =
					-interpolation_matrix(n_ie + n_i + 3 * j + 1, k);
				inequality_matrix(n_ie + 2 * n_i + 6 * j + 4, k) =
					interpolation_matrix(n_ie + n_i + 3 * j + 2, k);
				inequality_matrix(n_ie + 2 * n_i + 6 * j + 5, k) =
					-interpolation_matrix(n_ie + n_i + 3 * j + 2, k);
			}
		}
		for (int j = 0; j < n_t; j++) {
			for (int k = 0; k < inequality_matrix.cols(); k++) {
				inequality_matrix(n_ie + 2 * n_i + 6 * n_p + 2 * j + 0, k) =
					interpolation_matrix(n_ie + n_i + 3 * n_p + j, k);  //  Ax >= angle_lower_bound (always positive)
				inequality_matrix(n_ie + 2 * n_i + 6 * n_p + 2 * j + 1, k) =
					-interpolation_matrix(n_ie + n_i + 3 * n_p + j, k);  // -Ax >= -angle_upper_bound (always positive)
			}
		}
	}

	return true;
}

bool Single_Surface::get_inequality_values(VectorXd &inequality_values) {
	int n_ie = intern_params.n_inequality;
	for (int j = 0; j < n_ie; j++) inequality_values(j) = 0.0;

	if (intern_params.restricted_range) {
		int n_i = intern_params.n_interface;
		int n_p = intern_params.n_planar;
		int n_t = intern_params.n_tangent;
		for (int j = 0; j < n_i; j++) {
			inequality_values(n_ie + 2 * j + 0) =
				constraints.itrface[j].level_lower_bound();  //  Ax >=  lower_bound
			inequality_values(n_ie + 2 * j + 1) =
				-constraints.itrface[j].level_upper_bound();  // -Ax >= -upper_bound
		}
		for (int j = 0; j < n_p; j++) {
			inequality_values(n_ie + 2 * n_i + 6 * j + 0) =
				constraints.planar[j].nx_lower_bound();  //  Ax >=  lower_bound
			inequality_values(n_ie + 2 * n_i + 6 * j + 1) =
				-constraints.planar[j].nx_upper_bound();  // -Ax >= -upper_bound
			inequality_values(n_ie + 2 * n_i + 6 * j + 2) =
				constraints.planar[j].ny_lower_bound();  // ...
			inequality_values(n_ie + 2 * n_i + 6 * j + 3) =
				-constraints.planar[j].ny_upper_bound();
			inequality_values(n_ie + 2 * n_i + 6 * j + 4) =
				constraints.planar[j].nz_lower_bound();
			inequality_values(n_ie + 2 * n_i + 6 * j + 5) =
				-constraints.planar[j].nz_upper_bound();
		}
		for (int j = 0; j < n_t; j++) {
			inequality_values(2 * n_i + 6 * n_p + 2 * j + 0) =
				constraints.tangent[j].angle_lower_bound();  //  Ax >=  angle_lower_bound
			inequality_values(2 * n_i + 6 * n_p + 2 * j + 1) =
				-constraints.tangent[j].angle_upper_bound();  // -Ax >= -angle_upper_bound
		}
	}
	return true;
}

bool Single_Surface::get_inequality_values(VectorXd &b, VectorXd &r) {
	int n_ie = intern_params.n_inequality;
	int n_i = intern_params.n_interface;
	int n_p = intern_params.n_planar;
	int n_t = intern_params.n_tangent;

	// REF:
	// minimize f = 1/2 xT H x
	// s.t. b <= Ax <= b + r

	// inequality/rock type data
	if (n_ie != 0) {
		// compute largest distance between points
		// this distance will represent the upper
		// bound for inequality constraints
		// NOTE: this notion depends on the norm of gradient of the scalar field
		// to be ~1. Which is not true in reality. Will affect results.
		// a possible direction for future work
		std::vector<Point> aggregated_pts = convert_constraints_to_points(constraints);
		double distance = get_largest_distance_between_points(aggregated_pts);
		for (int j = 0; j < n_ie; j++) {
			if (constraints.inequality[j].level() > 0)  // stratigraphically above
			{
				b(j) = 0.0;
				r(j) = distance;
			}
			else  // stratigraphically below
			{
				b(j) = -distance;
				r(j) = distance;
			}
		}
	}

	// interface data
	for (int j = 0; j < n_i; j++) {
		b(n_ie + j) = constraints.itrface[j].level_lower_bound();
		r(n_ie + j) = constraints.itrface[j].level_upper_bound() - constraints.itrface[j].level_lower_bound();
	}

	// planar data
	for (int j = 0; j < n_p; j++) {
		// x-component
		b(n_ie + n_i + 3 * j + 0) = constraints.planar[j].nx_lower_bound();
		r(n_ie + n_i + 3 * j + 0) = constraints.planar[j].nx_upper_bound() - constraints.planar[j].nx_lower_bound();
		// y-component
		b(n_ie + n_i + 3 * j + 1) = constraints.planar[j].ny_lower_bound();
		r(n_ie + n_i + 3 * j + 1) = constraints.planar[j].ny_upper_bound() - constraints.planar[j].ny_lower_bound();
		// z-component
		b(n_ie + n_i + 3 * j + 2) = constraints.planar[j].nz_lower_bound();
		r(n_ie + n_i + 3 * j + 2) = constraints.planar[j].nz_upper_bound() - constraints.planar[j].nz_lower_bound();
	}

	// tangent data
	for (int j = 0; j < n_t; j++) {
		b(n_ie + n_i + 3 * n_p + j) = constraints.tangent[j].angle_lower_bound();
		r(n_ie + n_i + 3 * n_p + j) =
			constraints.tangent[j].angle_upper_bound() - constraints.tangent[j].angle_lower_bound();
	}

	return true;
}

void Single_Surface::process_input_data() {
	if (!get_interface_data())
		throw GRBF_Exceptions::no_iterface_data;

	if (parameters.use_restricted_range) {
		for (auto &interface_pt : constraints.itrface) {
			interface_pt.setLevelBounds(parameters.interface_uncertainty);
			std::cout << " Oncontact Bounds: " << std::endl;
			std::cout << "	" << interface_pt.level_lower_bound()
				<< " <= 0 <= " << interface_pt.level_upper_bound()
				<< std::endl;
		}
		for (auto &planar_pt : constraints.planar) {
			planar_pt.setNormalBounds(parameters.angular_uncertainty, parameters.angular_uncertainty / 2);  // Need more ROBUST
			// METHOD. Try large statistical sampling from von Mises spherical distribution
			std::cout << " Planar[] Bounds: " << std::endl;
			std::cout << "	nx: " << planar_pt.nx_lower_bound()
				<< " <= " << planar_pt.nx()
				<< " <= " << planar_pt.nx_upper_bound() << std::endl;
			std::cout << "	ny: " << planar_pt.ny_lower_bound()
				<< " <= " << planar_pt.ny()
				<< " <= " << planar_pt.ny_upper_bound() << std::endl;
			std::cout << "	nz: " << planar_pt.nz_lower_bound()
				<< " <= " << planar_pt.nz()
				<< " <= " << planar_pt.nz_upper_bound() << std::endl;
		}
		for (auto &tangent_pt : constraints.tangent) {
			tangent_pt.setAngleBounds(parameters.angular_uncertainty);
			std::cout << " Tangent Bounds: " << std::endl;
			std::cout << "	" << tangent_pt.angle_lower_bound()
				<< " <= " << tangent_pt.inner_product_constraint()
				<< " <= " << tangent_pt.angle_upper_bound()
				<< std::endl;
		}
	}
}

bool Single_Surface::get_interpolation_matrix(MatrixXd &interpolation_matrix) {
	int n_ie = intern_params.n_inequality;
	int n_i = intern_params.n_interface;
	int n_p = intern_params.n_planar;
	int n_t = intern_params.n_tangent;

	// Row and Column constraint order : inequality (ine) -> interface (itr) ->
	// planar (p_x,p_y,p_z) -> tangent (t)

	// Base Matrix Structure
	// | ine/ine ine/itr ine/p_x ine/p_y ine/p_z ine/t |
	// | itr/ine itr/itr itr/p_y itr/p_y itr/p_z itr/t |
	// | p_x/ine p_x/itr p_x/p_x p_x/p_y p_x/p_z p_x/t |
	// | p_y/ine p_y/itr p_y/p_x p_y/p_y p_y/p_z p_y/t |
	// | p_z/ine p_z/itr p_z/p_x p_z/p_y p_z/p_z p_z/t |
	// |   t/ine   t/itr   t/p_x   t/p_y   t/p_z   t/t |

	// Inequality Constraints:
	for (int j = 0; j < n_ie; j++) {
		// Row:inequality/Column:inequality block
		for (int k = 0; k < n_ie; k++) {
			kernel->set_points(constraints.inequality[j], constraints.inequality[k]);
			interpolation_matrix(j, k) = kernel->basis_pt_pt();
		}
		// Row:inequality/Column:interface block
		for (int k = 0; k < n_i; k++) {
			kernel->set_points(constraints.inequality[j], constraints.itrface[k]);
			interpolation_matrix(j, k + n_ie) = kernel->basis_pt_pt();
		}
		// Row:inequality/Column:planar block
		for (int k = 0; k < n_p; k++) {
			kernel->set_points(constraints.inequality[j], constraints.planar[k]);
			interpolation_matrix(j, 3 * k + n_ie + n_i) = kernel->basis_pt_planar_x();
			interpolation_matrix(j, 3 * k + n_ie + n_i + 1) = kernel->basis_pt_planar_y();
			interpolation_matrix(j, 3 * k + n_ie + n_i + 2) = kernel->basis_pt_planar_z();
		}
		// Row:inequality/Column:tangent block
		for (int k = 0; k < n_t; k++) {
			kernel->set_points(constraints.inequality[j], constraints.tangent[k]);
			interpolation_matrix(j, n_ie + n_i + 3 * n_p + k) = kernel->basis_pt_tangent();
		}
	}
	// Interface Constraints:
	for (int j = 0; j < n_i; j++) {
		// Row:interface/Column:inequality block
		for (int k = 0; k < n_ie; k++) {
			kernel->set_points(constraints.itrface[j], constraints.inequality[k]);
			interpolation_matrix(j + n_ie, k) = kernel->basis_pt_pt();
		}
		// Row:interface/Column:interface block
		for (int k = 0; k < n_i; k++) {
			kernel->set_points(constraints.itrface[j], constraints.itrface[k]);
			interpolation_matrix(j + n_ie, k + n_ie) = kernel->basis_pt_pt();
		}
		// Row:interface/Column:planar block
		for (int k = 0; k < n_p; k++) {
			kernel->set_points(constraints.itrface[j], constraints.planar[k]);
			interpolation_matrix(j + n_ie, 3 * k + n_ie + n_i) = kernel->basis_pt_planar_x();
			interpolation_matrix(j + n_ie, 3 * k + n_ie + n_i + 1) = kernel->basis_pt_planar_y();
			interpolation_matrix(j + n_ie, 3 * k + n_ie + n_i + 2) = kernel->basis_pt_planar_z();
		}
		// Row:interface/Column:tangent block
		for (int k = 0; k < n_t; k++) {
			kernel->set_points(constraints.itrface[j], constraints.tangent[k]);
			interpolation_matrix(j + n_ie, n_ie + n_i + 3 * n_p + k) = kernel->basis_pt_tangent();
		}
	}
	// Planar Constraints
	for (int j = 0; j < n_p; j++) {
		// Row:planar/Column:inequality block
		for (int k = 0; k < n_ie; k++) {
			kernel->set_points(constraints.planar[j], constraints.inequality[k]);
			interpolation_matrix(3 * j + n_ie + n_i, k) = kernel->basis_planar_x_pt();
			interpolation_matrix(3 * j + n_ie + n_i + 1, k) = kernel->basis_planar_y_pt();
			interpolation_matrix(3 * j + n_ie + n_i + 2, k) = kernel->basis_planar_z_pt();
		}
		// Row:planar/Column:interface block
		for (int k = 0; k < n_i; k++) {
			kernel->set_points(constraints.planar[j], constraints.itrface[k]);
			interpolation_matrix(3 * j + n_ie + n_i, k + n_ie) = kernel->basis_planar_x_pt();
			interpolation_matrix(3 * j + n_ie + n_i + 1, k + n_ie) = kernel->basis_planar_y_pt();
			interpolation_matrix(3 * j + n_ie + n_i + 2, k + n_ie) = kernel->basis_planar_z_pt();
		}
		// Row:planar/Column:planar block
		for (int k = 0; k < n_p; k++) {
			kernel->set_points(constraints.planar[j], constraints.planar[k]);
			interpolation_matrix(3 * j + n_ie + n_i, 3 * k + n_ie + n_i) =
				kernel->basis_planar_planar(Parameter_Types::DXDX);
			interpolation_matrix(3 * j + n_ie + n_i, 3 * k + n_ie + n_i + 1) =
				kernel->basis_planar_planar(Parameter_Types::DXDY);
			interpolation_matrix(3 * j + n_ie + n_i, 3 * k + n_ie + n_i + 2) =
				kernel->basis_planar_planar(Parameter_Types::DXDZ);
			interpolation_matrix(3 * j + n_ie + n_i + 1, 3 * k + n_ie + n_i) =
				kernel->basis_planar_planar(Parameter_Types::DYDX);
			interpolation_matrix(3 * j + n_ie + n_i + 1, 3 * k + n_ie + n_i + 1) =
				kernel->basis_planar_planar(Parameter_Types::DYDY);
			interpolation_matrix(3 * j + n_ie + n_i + 1, 3 * k + n_ie + n_i + 2) =
				kernel->basis_planar_planar(Parameter_Types::DYDZ);
			interpolation_matrix(3 * j + n_ie + n_i + 2, 3 * k + n_ie + n_i) =
				kernel->basis_planar_planar(Parameter_Types::DZDX);
			interpolation_matrix(3 * j + n_ie + n_i + 2, 3 * k + n_ie + n_i + 1) =
				kernel->basis_planar_planar(Parameter_Types::DZDY);
			interpolation_matrix(3 * j + n_ie + n_i + 2, 3 * k + n_ie + n_i + 2) =
				kernel->basis_planar_planar(Parameter_Types::DZDZ);
		}
		// Row:planar/Column:tangent block
		for (int k = 0; k < n_t; k++) {
			kernel->set_points(constraints.planar[j], constraints.tangent[k]);
			interpolation_matrix(3 * j + n_ie + n_i, n_ie + n_i + 3 * n_p + k) =
				kernel->basis_planar_tangent(Parameter_Types::DX);
			interpolation_matrix(3 * j + n_ie + n_i + 1, n_ie + n_i + 3 * n_p + k) =
				kernel->basis_planar_tangent(Parameter_Types::DY);
			interpolation_matrix(3 * j + n_ie + n_i + 2, n_ie + n_i + 3 * n_p + k) =
				kernel->basis_planar_tangent(Parameter_Types::DZ);
		}
	}
	// Tangent Constraints
	for (int j = 0; j < n_t; j++) {
		// Row:tangent/Column:inequality block
		for (int k = 0; k < n_ie; k++) {
			kernel->set_points(constraints.tangent[j], constraints.inequality[k]);
			interpolation_matrix(j + n_ie + n_i + 3 * n_p, k) = kernel->basis_tangent_pt();
		}
		// Row:tangent/Column:interface block
		for (int k = 0; k < n_i; k++) {
			kernel->set_points(constraints.tangent[j], constraints.itrface[k]);
			interpolation_matrix(j + n_ie + n_i + 3 * n_p, k + n_ie) = kernel->basis_tangent_pt();
		}
		// Row:tangent/Column:planar block
		for (int k = 0; k < n_p; k++) {
			kernel->set_points(constraints.tangent[j], constraints.planar[k]);
			interpolation_matrix(j + n_ie + n_i + 3 * n_p, 3 * k + n_ie + n_i) =
				kernel->basis_tangent_planar(Parameter_Types::DX);
			interpolation_matrix(j + n_ie + n_i + 3 * n_p, 3 * k + n_ie + n_i + 1) =
				kernel->basis_tangent_planar(Parameter_Types::DY);
			interpolation_matrix(j + n_ie + n_i + 3 * n_p, 3 * k + n_ie + n_i + 2) =
				kernel->basis_tangent_planar(Parameter_Types::DZ);
		}
		// Row:tangent/Column:tangent block
		for (int k = 0; k < n_t; k++) {
			kernel->set_points(constraints.tangent[j], constraints.tangent[k]);
			interpolation_matrix(j + n_ie + n_i + 3 * n_p, n_ie + n_i + 3 * n_p + k) =
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