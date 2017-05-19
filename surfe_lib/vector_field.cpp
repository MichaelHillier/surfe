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
	b_parameters.n_equality = b_parameters.n_interface + 3*b_parameters.n_planar + b_parameters.n_tangent;

	// polynomial parameters ...

	b_parameters.poly_term = false;
	b_parameters.modified_basis = false;
	b_parameters.problem_type = Parameter_Types::Linear;

	b_parameters.n_poly_terms = 0;

	return true;
}

bool Vector_Field::process_input_data()
{
	return true;
}

bool Vector_Field::get_equality_values( VectorXd &equality_values )
{
	int j = 0;
	int k = 0;
	for (j = 0; j < (int)b_input.planar->size(); j++){
		equality_values(3*j) = b_input.planar->at(j).nx();
		equality_values(3*j + 1) = b_input.planar->at(j).ny();
		equality_values(3*j + 2) = b_input.planar->at(j).nz();
	}
	for (k = 0; k < (int)b_input.tangent->size(); k++){
		equality_values(3 * j + k) = b_input.tangent->at(k).inner_product_constraint();
	}

	return true;
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

	return true;
}

bool Vector_Field::setup_system_solver()
{

	int n = b_parameters.n_equality + b_parameters.n_poly_terms;

	VectorXd equality_values(n);
	get_equality_values(equality_values);

	MatrixXd interpolation_matrix(n, n);
	if (!get_interpolation_matrix(interpolation_matrix)) return false;

	std::cout<<" Interpolation matrix:\n"<< interpolation_matrix << std::endl;
	std::cout<<" RHS:\n"<< equality_values << std::endl;

	Linear_LU_decomposition *llu = new Linear_LU_decomposition(interpolation_matrix,equality_values);
	if (!llu->solve()) return false;
	solver = llu;

	std::cout << " solution :\n" << solver->weights << std::endl;

	return true;
}

void Vector_Field::eval_scalar_interpolant_at_point( Point &p )
{
	int n_p = b_parameters.n_planar;
	int n_t = b_parameters.n_tangent;

	Kernel *kernel_j = kernel->clone();
	double elemsum_1 = 0.0;
	double elemsum_2 = 0.0;
	double elemsum_3 = 0.0;
	double elemsum_4 = 0.0;
	double poly = 0.0;
	for (int k = 0; k < n_p; k++ ){
		kernel_j->set_points(p, b_input.planar->at(k));
		elemsum_3 += solver->weights[3 * k] * kernel_j->basis_pt_planar_x();
		elemsum_3 += solver->weights[3 * k + 1] * kernel_j->basis_pt_planar_y();
		elemsum_3 += solver->weights[3 * k + 2] * kernel_j->basis_pt_planar_z();
	}
	for (int k = 0; k < n_t; k++){
		kernel_j->set_points(p, b_input.tangent->at(k));
		elemsum_4 += solver->weights[3 * n_p + k] * kernel_j->basis_pt_tangent();
	}
	p.set_scalar_field(elemsum_1 + elemsum_2 + elemsum_3 + elemsum_4 + poly);
	delete kernel_j;
}

void Vector_Field::eval_vector_interpolant_at_point( Point &p )
{
	// this method is a copy of the eval_scalar_interpolant_at_points() method 
	// TO DO: Fix eval_scalar_interpolant_at_points() method to actually computer the scalar field and not the vector field.
	int n_p = b_parameters.n_planar;

	Kernel *kernel_j = kernel->clone();
	double elemsum_x = 0.0;
	double elemsum_y = 0.0;
	double elemsum_z = 0.0;
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
	p.set_vector_field(elemsum_x,elemsum_y,elemsum_z);
	delete kernel_j;
}
