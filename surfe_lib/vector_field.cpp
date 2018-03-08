#include <vector_field.h>
#include <math_methods.h>
#include <matrix_solver.h>
#include <basis.h>
#include <algorithm>
#include <vector>

#include <iostream>
#include <iomanip>
#include <fstream>


bool Vector_Field::get_method_parameters()
{
	// # of constraints for each constraint type ...
	b_parameters.n_interface = 0;
	b_parameters.n_inequality = 0;
	b_parameters.n_planar = (int)b_input.planar->size();
	b_parameters.n_tangent = (int)b_input.tangent->size();
	// Total number of constraints ...
	b_parameters.n_constraints = b_parameters.n_interface +	b_parameters.n_inequality + 3*b_parameters.n_planar + b_parameters.n_tangent;
	// Total number of equality constraints
	//b_parameters.n_equality = b_parameters.n_interface + 3*b_parameters.n_planar + b_parameters.n_tangent;

// 	// polynomial parameters ...
// 
// 	b_parameters.poly_term = false;
// 	b_parameters.modified_basis = false;
// 	if (b_parameters.n_tangent != 0)
// 	b_parameters.problem_type = Parameter_Types::Linear;

//	b_parameters.n_poly_terms = 0;

	b_parameters.poly_term = false;
	if (m_parameters.use_restricted_range) b_parameters.restricted_range = true;
	else b_parameters.n_equality = b_parameters.n_interface + 3 * b_parameters.n_planar + b_parameters.n_tangent;

	// polynomial parameters ...
	if (b_parameters.n_inequality != 0 || m_parameters.use_restricted_range != 0)
	{
		b_parameters.modified_basis = true;
		b_parameters.problem_type = Parameter_Types::Quadratic;
	}
	else
	{
		b_parameters.modified_basis = false;
		b_parameters.problem_type = Parameter_Types::Linear;
	}

	int m = m_parameters.polynomial_order + 1;
	b_parameters.n_poly_terms = (int)(m*(m + 1)*(m + 2) / 6); // for 3D only...

	return true;
}

bool Vector_Field::process_input_data()
{
	if (m_parameters.use_restricted_range)
	{
		for (int j = 0; j < (int)b_input.tangent->size(); j++){
			b_input.tangent->at(j).setAngleBounds(m_parameters.angular_uncertainty);
			cout << " Tangent[" << j << "] Bounds: " << endl;
			cout << "	" << b_input.tangent->at(j).angle_lower_bound() << " <= " << b_input.tangent->at(j).inner_product_constraint() << " <= " << b_input.tangent->at(j).angle_upper_bound() << endl;
		}
	}

	return true;
}

bool Vector_Field::get_equality_values( VectorXd &equality_values )
{
	int j = 0;
	int k = 0;
	int l = 0;
	int m = 0;

	for (j = 0; j < (int)b_input.itrface->size(); j++) equality_values(j) = b_input.itrface->at(j).level();
	for (k = 0; k < (int)b_input.planar->size(); k++){
		equality_values(3 * k + j) = b_input.planar->at(k).nx();
		equality_values(3 * k + j + 1) = b_input.planar->at(k).ny();
		equality_values(3 * k + j + 2) = b_input.planar->at(k).nz();
	}
	for (l = 0; l < (int)b_input.tangent->size(); l++) equality_values(l + 3 * k + j) = b_input.tangent->at(l).inner_product_constraint();
	if (b_parameters.poly_term) for (m = 0; m < (int)b_parameters.n_poly_terms; m++) equality_values(m + l + 3 * k + j) = 0.0;

	return true;
}

bool Vector_Field::get_inequality_values(VectorXd &b, VectorXd &r)
{
	int n_t = b_parameters.n_tangent;
	int n_p = b_parameters.n_planar;

	// planar data
	for (int j = 0; j < n_p; j++){
		// x-component
		b(3 * j + 0) = b_input.planar->at(j).nx_lower_bound();
		r(3 * j + 0) = b_input.planar->at(j).nx_upper_bound() - b_input.planar->at(j).nx_lower_bound();
		// y-component
		b(3 * j + 1) = b_input.planar->at(j).ny_lower_bound();
		r(3 * j + 1) = b_input.planar->at(j).ny_upper_bound() - b_input.planar->at(j).ny_lower_bound();
		// z-component
		b(3 * j + 2) = b_input.planar->at(j).nz_lower_bound();
		r(3 * j + 2) = b_input.planar->at(j).nz_upper_bound() - b_input.planar->at(j).nz_lower_bound();
	}
	// tangent data
	for (int j = 0; j < n_t; j++){
		b(3 * n_p + j) = b_input.tangent->at(j).angle_lower_bound();
		r(3 * n_p + j) = b_input.tangent->at(j).angle_upper_bound() - b_input.tangent->at(j).angle_lower_bound();
	}

	return true;
}

bool Vector_Field::_get_polynomial_matrix_block(MatrixXd &poly_matrix)
{
	int n_ie = b_parameters.n_inequality;
	int n_i = b_parameters.n_interface;
	int n_p = b_parameters.n_planar;
	int n_t = b_parameters.n_tangent;

	int nl = n_i + n_ie; // n_ie should always be zero.

	p_basis = create_polynomial_basis(m_parameters.polynomial_order);

	if ((int)poly_matrix.rows() != b_parameters.n_poly_terms) return false;
	// for interface points ...
	for (int j = 0; j < nl; j++){
		p_basis->set_point(b_input.itrface->at(j));
		VectorXd b = p_basis->basis();
		if ((int)b.rows() != b_parameters.n_poly_terms) return false;
		for (int k = 0; k < (int)b.rows(); k++) poly_matrix(k, j) = b(k);
	}
	// for planar points ...
	for (int j = 0; j < n_p; j++){
		p_basis->set_point(b_input.planar->at(j));
		VectorXd bx = p_basis->dx();
		VectorXd by = p_basis->dy();
		VectorXd bz = p_basis->dz();
		if ((int)bx.rows() != b_parameters.n_poly_terms) return false;
		for (int k = 0; k < (int)bx.rows(); k++){
			poly_matrix(k, nl + 3 * j) = bx(k);
			poly_matrix(k, nl + 3 * j + 1) = by(k);
			poly_matrix(k, nl + 3 * j + 2) = bz(k);
		}
	}
	// for tangent points ...
	for (int j = 0; j < n_t; j++){
		p_basis->set_point(b_input.tangent->at(j));
		VectorXd bx = p_basis->dx();
		VectorXd by = p_basis->dy();
		VectorXd bz = p_basis->dz();
		if ((int)bx.rows() != b_parameters.n_poly_terms) return false;
		for (int k = 0; k < (int)bx.rows(); k++){
			poly_matrix(k, nl + 3 * n_p + j) = b_input.tangent->at(j).tx()*bx(k) + b_input.tangent->at(j).ty()*by(k) + b_input.tangent->at(j).tz()*bz(k);
		}
	}

	return true;
}

bool Vector_Field::_insert_polynomial_matrix_blocks_in_interpolation_matrix(const MatrixXd &poly_matrix, MatrixXd &interpolation_matrix)
{
	int n_ie = b_parameters.n_inequality;
	int n_i = b_parameters.n_interface;
	int n_p = b_parameters.n_planar;
	int n_t = b_parameters.n_tangent;

	// build polynomial blocks
	// | A PT |
	// | P 0  |
	// start with P
	for (int j = 0; j < (int)poly_matrix.rows(); j++){
		for (int k = 0; k < (int)poly_matrix.cols(); k++){
			interpolation_matrix(n_ie + n_i + 3 * n_p + n_t + j, k) = poly_matrix(j, k);
			interpolation_matrix(k, n_ie + n_i + 3 * n_p + n_t + j) = interpolation_matrix(n_ie + n_i + 3 * n_p + n_t + j, k);
		}
	}

	for (int j = 0; j < (int)poly_matrix.rows(); j++){
		for (int k = 0; k < (int)poly_matrix.rows(); k++){
			interpolation_matrix(n_ie + n_i + 3 * n_p + n_t + j, n_ie + n_i + 3 * n_p + n_t + k) = 0;
		}
	}

	return true;
}

Polynomial_Basis * Vector_Field::create_polynomial_basis(const int &poly_order)
{
	if (poly_order == 0) return new Poly_Zero;
	else if (poly_order == 1) return new Poly_First;
	else return new Poly_Second;
}

bool Vector_Field::get_interpolation_matrix( MatrixXd &interpolation_matrix )
{
	int n_p = b_parameters.n_planar;
	int n_t = b_parameters.n_tangent;

	// Base Matrix Structure
	// | p_x/p_x p_x/p_y p_x/p_z |
	// | p_y/p_x p_y/p_y p_y/p_z |
	// | p_z/p_x p_z/p_y p_z/p_z |

	// Planar Constraints
	for (int j = 0; j < n_p; j++ ){
		// Row:planar/Column:planar block
		for (int k = 0; k < n_p; k++ ){
			kernel->set_points(b_input.planar->at(j), b_input.planar->at(k));
			interpolation_matrix(3*j,3*k) = kernel->basis_planar_planar(Parameter_Types::DXDX);
			interpolation_matrix(3*j,3*k + 1) = kernel->basis_planar_planar(Parameter_Types::DXDY);
			interpolation_matrix(3*j,3*k + 2) = kernel->basis_planar_planar(Parameter_Types::DXDZ);
			interpolation_matrix(3*j + 1,3*k) = kernel->basis_planar_planar(Parameter_Types::DYDX);
			interpolation_matrix(3*j + 1,3*k + 1) = kernel->basis_planar_planar(Parameter_Types::DYDY);
			interpolation_matrix(3*j + 1,3*k + 2) = kernel->basis_planar_planar(Parameter_Types::DYDZ);
			interpolation_matrix(3*j + 2,3*k) = kernel->basis_planar_planar(Parameter_Types::DZDX);
			interpolation_matrix(3*j + 2,3*k + 1) = kernel->basis_planar_planar(Parameter_Types::DZDY);
			interpolation_matrix(3*j + 2,3*k + 2) = kernel->basis_planar_planar(Parameter_Types::DZDZ);
		}
		// Row:planar/Column:tangent block
		for (int k = 0; k < n_t; k++){
			kernel->set_points(b_input.planar->at(j), b_input.tangent->at(k));
			interpolation_matrix(3 * j, 3 * n_p + k) = kernel->basis_planar_tangent(Parameter_Types::DX);
			interpolation_matrix(3 * j + 1, 3 * n_p + k) = kernel->basis_planar_tangent(Parameter_Types::DY);
			interpolation_matrix(3 * j + 2, 3 * n_p + k) = kernel->basis_planar_tangent(Parameter_Types::DZ);
		}
	}
	// Tangent Constraints 
	for (int j = 0; j < n_t; j++){
		// Row:tangent/Column:planar block
		for (int k = 0; k < n_p; k++){
			kernel->set_points(b_input.tangent->at(j), b_input.planar->at(k));
			interpolation_matrix(j + 3 * n_p, 3 * k) = kernel->basis_tangent_planar(Parameter_Types::DX);
			interpolation_matrix(j + 3 * n_p, 3 * k + 1) = kernel->basis_tangent_planar(Parameter_Types::DY);
			interpolation_matrix(j + 3 * n_p, 3 * k + 2) = kernel->basis_tangent_planar(Parameter_Types::DZ);
		}
		// Row:tangent/Column:tangent block
		for (int k = 0; k < n_t; k++){
			kernel->set_points(b_input.tangent->at(j), b_input.tangent->at(k));
			interpolation_matrix(j + 3 * n_p, 3 * n_p + k) = kernel->basis_tangent_tangent();
		}
	}

	// build polynomial blocks if required
	// | A PT |
	// | P 0  |
	if (b_parameters.poly_term)
	{
		MatrixXd poly_matrix(b_parameters.n_poly_terms, b_parameters.n_constraints);
		if (!_get_polynomial_matrix_block(poly_matrix)) return false;
		// fill remaining matrix blocks (P, PT, 0)
		if (!_insert_polynomial_matrix_blocks_in_interpolation_matrix(poly_matrix, interpolation_matrix)) return false;
	}

	return true;
}

bool Vector_Field::setup_system_solver()
{

	int n = b_parameters.n_equality + b_parameters.n_poly_terms;

	if (b_parameters.restricted_range)
	{
		int n_c = b_parameters.n_constraints;
		VectorXd b(n_c);
		VectorXd r(n_c);
		get_inequality_values(b, r);

		MatrixXd interpolation_matrix(n_c, n_c);
		if (!get_interpolation_matrix(interpolation_matrix)) return false;

		MatrixXd inequality_matrix(n_c, n_c);
		inequality_matrix = interpolation_matrix;

		Quadratic_Predictor_Corrector_LOQO *qpc = new Quadratic_Predictor_Corrector_LOQO(interpolation_matrix, inequality_matrix, b, r);
		if (!qpc->solve())
		{
			error_msg.append(" LOQO Quadratic Solver failure.");
			return false;
		}
		solver = qpc;
	}
	else
	{
		VectorXd equality_values;
		get_equality_values(equality_values);

		MatrixXd interpolation_matrix(n, n);
		if (!get_interpolation_matrix(interpolation_matrix)) return false;

		Linear_LU_decomposition llu(interpolation_matrix, equality_values);
		if (!llu.solve()) return false;
		solver = &llu;
	}

	return true;
}

bool Vector_Field::convert_modified_kernel_to_rbf_kernel()
{
	if (rbf_kernel == NULL || kernel == NULL) return false;

	// prep for linear prob...
	// set the constraints...
	for (int j = 0; j < b_input.planar->size(); j++){
		eval_vector_interpolant_at_point(b_input.planar->at(j));
		// debug
		// 		cout<<" Planar["<<j<<"]: "<<endl;
		// 		cout<<"	Nx = "<<b_input.planar->at(j).nx()<<" Nx interpolated = "<<b_input.planar->at(j).nx_interp()<<endl;
		// 		cout<<"	Ny = "<<b_input.planar->at(j).ny()<<" Ny interpolated = "<<b_input.planar->at(j).ny_interp()<<endl;
		// 		cout<<"	Nz = "<<b_input.planar->at(j).nz()<<" Nz interpolated = "<<b_input.planar->at(j).nz_interp()<<endl;
		double normal[3] = { b_input.planar->at(j).nx_interp(), b_input.planar->at(j).ny_interp(), b_input.planar->at(j).nz_interp() };
		b_input.planar->at(j).setNormal(normal[0], normal[1], normal[2]);
	}
	for (int j = 0; j < b_input.tangent->size(); j++){
		eval_vector_interpolant_at_point(b_input.tangent->at(j));
		// 		cout<<" Tangent["<<j<<"]: "<<endl;
		// 		cout<<"	Tx = "<<b_input.tangent->at(j).tx()<<" Nx interpolated = "<<b_input.tangent->at(j).nx_interp()<<endl;
		// 		cout<<"	Ty = "<<b_input.tangent->at(j).ty()<<" Ny interpolated = "<<b_input.tangent->at(j).ny_interp()<<endl;
		// 		cout<<"	Tz = "<<b_input.tangent->at(j).tz()<<" Nz interpolated = "<<b_input.tangent->at(j).nz_interp()<<endl;
		// 		cout<<" Tx*nx + Ty*ny + Tz*nz = "<<b_input.tangent->at(j).tx()*b_input.tangent->at(j).nx_interp() + b_input.tangent->at(j).ty()*b_input.tangent->at(j).ny_interp() +
		// 			b_input.tangent->at(j).tz()*b_input.tangent->at(j).nz_interp()<<endl;
		double vf[3] = { b_input.tangent->at(j).nx_interp(), b_input.tangent->at(j).ny_interp(), b_input.tangent->at(j).nz_interp() };
		double inner_product = vf[0] * b_input.tangent->at(j).tx() + vf[1] * b_input.tangent->at(j).ty() + vf[2] * b_input.tangent->at(j).tz();
		b_input.tangent->at(j).setInnerProductConstraint(inner_product);
	}

	// switch from modified kernel to normal rbf kernel
	kernel = rbf_kernel;

	//int n_p = b_parameters.n_poly_terms;
	int n_p = 0;
	b_parameters.n_poly_terms = 0;
	if (m_parameters.use_restricted_range) b_parameters.restricted_range = false;
	b_parameters.n_interface = (int)b_input.itrface->size();
	b_parameters.n_inequality = (int)b_input.inequality->size();
	b_parameters.n_equality = b_parameters.n_interface + 3 * b_parameters.n_planar + b_parameters.n_tangent;
	//b_parameters.poly_term = true;
	b_parameters.poly_term = false;
	b_parameters.modified_basis = false;
	b_parameters.problem_type = Parameter_Types::Linear;
	int n_e = b_parameters.n_equality;
	VectorXd equality_values(n_e + n_p);
	get_equality_values(equality_values);

	MatrixXd interpolation_matrix(n_e + n_p, n_e + n_p);
	if (!get_interpolation_matrix(interpolation_matrix)) return false;

// 	std::cout << "Interpolation matrix:\n" << interpolation_matrix << std::endl;
// 	std::cout << "Solution vector:\n" << equality_values<< std::endl;

	Linear_LU_decomposition *llu = new Linear_LU_decomposition(interpolation_matrix, equality_values);
	if (!llu->solve()) return false;
	solver = llu;

	return true;
}

void Vector_Field::eval_scalar_interpolant_at_point( Point &p )
{
	int n_p = b_parameters.n_planar;
	int n_t = b_parameters.n_tangent;

	Kernel *kernel_j = kernel->clone();
	double elemsum_x = 0.0;
	double elemsum_y = 0.0;
	double elemsum_z = 0.0;
	double elemsum_3_x = 0.0;
	double elemsum_3_y = 0.0;
	double elemsum_3_z = 0.0;
	double poly_x = 0.0;
	double poly_y = 0.0;
	double poly_z = 0.0;
	for (int k = 0; k < n_p; k++ ){
		kernel_j->set_points(p, b_input.planar->at(k));
		elemsum_x += solver->weights[3*k] * kernel_j->basis_planar_planar(Parameter_Types::DXDX);
		elemsum_x += solver->weights[3*k + 1] * kernel_j->basis_planar_planar(Parameter_Types::DXDY);
		elemsum_x += solver->weights[3*k + 2] * kernel_j->basis_planar_planar(Parameter_Types::DXDZ);

		elemsum_y += solver->weights[3*k] * kernel_j->basis_planar_planar(Parameter_Types::DYDX);
		elemsum_y += solver->weights[3*k + 1] * kernel_j->basis_planar_planar(Parameter_Types::DYDY);
		elemsum_y += solver->weights[3*k + 2] * kernel_j->basis_planar_planar(Parameter_Types::DYDZ);

		elemsum_z += solver->weights[3*k] * kernel_j->basis_planar_planar(Parameter_Types::DZDX);
		elemsum_z += solver->weights[3*k + 1] * kernel_j->basis_planar_planar(Parameter_Types::DZDY);
		elemsum_z += solver->weights[3*k + 2] * kernel_j->basis_planar_planar(Parameter_Types::DZDZ);
	}
	// tangent constraints
	for (int k = 0; k < n_t; k++){
		kernel->set_points(p, b_input.tangent->at(k));
		elemsum_3_x += solver->weights[3 * n_p + k] * kernel->basis_planar_tangent(Parameter_Types::DX);
		elemsum_3_y += solver->weights[3 * n_p + k] * kernel->basis_planar_tangent(Parameter_Types::DY);
		elemsum_3_z += solver->weights[3 * n_p + k] * kernel->basis_planar_tangent(Parameter_Types::DZ);
	}
	if (b_parameters.poly_term)
	{
		Polynomial_Basis *p_basis_j = p_basis->clone();
		p_basis_j->set_point(p);
		VectorXd bx = p_basis_j->dx();
		VectorXd by = p_basis_j->dy();
		VectorXd bz = p_basis_j->dz();
		for (int k = 0; k < (int)bx.size(); k++){
			poly_x += bx(k) * solver->weights[3 * n_p + n_t + k];
			poly_y += by(k) * solver->weights[3 * n_p + n_t + k];
			poly_z += bz(k) * solver->weights[3 * n_p + n_t + k];
		}
		delete p_basis_j;
	}
	double nx = elemsum_x + elemsum_3_x + poly_x;
	double ny = elemsum_y + elemsum_3_y + poly_y;
	double nz = elemsum_z + elemsum_3_z + poly_z;
	p.set_vector_field(nx,ny,nz);
	delete kernel_j;
}

void Vector_Field::eval_vector_interpolant_at_point( Point &p )
{
	// this method is a copy of the eval_scalar_interpolant_at_points() method 
	// TO DO: Fix eval_scalar_interpolant_at_points() method to actually computer the scalar field and not the vector field.
	int n_p = b_parameters.n_planar;
	int n_t = b_parameters.n_tangent;

	Kernel *kernel_j = kernel->clone();
	double elemsum_x = 0.0;
	double elemsum_y = 0.0;
	double elemsum_z = 0.0;
	double elemsum_3_x = 0.0;
	double elemsum_3_y = 0.0;
	double elemsum_3_z = 0.0;
	double poly_x = 0.0;
	double poly_y = 0.0;
	double poly_z = 0.0;
	for (int k = 0; k < n_p; k++){
		kernel_j->set_points(p, b_input.planar->at(k));
		elemsum_x += solver->weights[3 * k] * kernel_j->basis_planar_planar(Parameter_Types::DXDX);
		elemsum_x += solver->weights[3 * k + 1] * kernel_j->basis_planar_planar(Parameter_Types::DXDY);
		elemsum_x += solver->weights[3 * k + 2] * kernel_j->basis_planar_planar(Parameter_Types::DXDZ);

		elemsum_y += solver->weights[3 * k] * kernel_j->basis_planar_planar(Parameter_Types::DYDX);
		elemsum_y += solver->weights[3 * k + 1] * kernel_j->basis_planar_planar(Parameter_Types::DYDY);
		elemsum_y += solver->weights[3 * k + 2] * kernel_j->basis_planar_planar(Parameter_Types::DYDZ);

		elemsum_z += solver->weights[3 * k] * kernel_j->basis_planar_planar(Parameter_Types::DZDX);
		elemsum_z += solver->weights[3 * k + 1] * kernel_j->basis_planar_planar(Parameter_Types::DZDY);
		elemsum_z += solver->weights[3 * k + 2] * kernel_j->basis_planar_planar(Parameter_Types::DZDZ);
	}
	// tangent constraints
	for (int k = 0; k < n_t; k++){
		kernel->set_points(p, b_input.tangent->at(k));
		elemsum_3_x += solver->weights[3 * n_p + k] * kernel->basis_planar_tangent(Parameter_Types::DX);
		elemsum_3_y += solver->weights[3 * n_p + k] * kernel->basis_planar_tangent(Parameter_Types::DY);
		elemsum_3_z += solver->weights[3 * n_p + k] * kernel->basis_planar_tangent(Parameter_Types::DZ);
	}
	if (b_parameters.poly_term)
	{
		Polynomial_Basis *p_basis_j = p_basis->clone();
		p_basis_j->set_point(p);
		VectorXd bx = p_basis_j->dx();
		VectorXd by = p_basis_j->dy();
		VectorXd bz = p_basis_j->dz();
		for (int k = 0; k < (int)bx.size(); k++){
			poly_x += bx(k) * solver->weights[3 * n_p + n_t + k];
			poly_y += by(k) * solver->weights[3 * n_p + n_t + k];
			poly_z += bz(k) * solver->weights[3 * n_p + n_t + k];
		}
		delete p_basis_j;
	}
	double nx = elemsum_x + elemsum_3_x + poly_x;
	double ny = elemsum_y + elemsum_3_y + poly_y;
	double nz = elemsum_z + elemsum_3_z + poly_z;
	p.set_vector_field(nx, ny, nz);
	delete kernel_j;
}
