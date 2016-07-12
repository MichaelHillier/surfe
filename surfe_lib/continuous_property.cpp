#include <continuous_property.h>
#include <math_methods.h>
#include <matrix_solver.h>
#include <basis.h>
#include <algorithm>
#include <vector>

#include <iostream>
#include <iomanip>
#include <fstream>

bool Continuous_Property::_get_polynomial_matrix_block(std::vector< std::vector <double> > &poly_matrix)
{
	int n_ie = b_parameters.n_inequality;
	int n_i = b_parameters.n_interface;
	int n_p = b_parameters.n_planar;
	int n_t = b_parameters.n_tangent;

	int nl = n_i + n_ie; // n_ie should always be zero.

	p_basis = create_polynomial_basis(m_parameters.polynomial_order);

	if ((int)poly_matrix.size() != b_parameters.n_poly_terms ) return false;
	// for interface points ...
	for (int j = 0; j < nl; j++ ){
		p_basis->set_point(b_input.itrface->at(j));
		std::vector<double> b = p_basis->basis();
		if ((int)b.size() != b_parameters.n_poly_terms ) return false;
		for (int k = 0; k < (int)b.size(); k++ ) poly_matrix[k][j] = b[k];
	}
	// for planar points ...
	for (int j = 0; j < n_p; j++ ){
		p_basis->set_point(b_input.planar->at(j));
		std::vector<double> bx = p_basis->dx();
		std::vector<double> by = p_basis->dy();
		std::vector<double> bz = p_basis->dz();
		if ((int)bx.size() != b_parameters.n_poly_terms ) return false;
		for (int k = 0; k < (int)bx.size(); k++ ){
			poly_matrix[k][nl + 3*j] = bx[k];
			poly_matrix[k][nl + 3*j + 1] = by[k];
			poly_matrix[k][nl + 3*j + 2] = bz[k];
		}
	}
	// for tangent points ...
	for (int j = 0; j < n_t; j++ ){
		p_basis->set_point(b_input.tangent->at(j));
		std::vector<double> bx = p_basis->dx();
		std::vector<double> by = p_basis->dy();
		std::vector<double> bz = p_basis->dz();
		if ((int)bx.size() != b_parameters.n_poly_terms ) return false;
		for (int k = 0; k < (int)bx.size(); k++ ){
			poly_matrix[k][nl + 3 * n_p + j] = b_input.tangent->at(j).tx()*bx[k] + b_input.tangent->at(j).ty()*by[k] + b_input.tangent->at(j).tz()*bz[k];
		}
	}

	return true;
}

bool Continuous_Property::_insert_polynomial_matrix_blocks_in_interpolation_matrix( const std::vector< std::vector <double> > &poly_matrix, std::vector< std::vector <double> > &interpolation_matrix )
{
	int n_ie = b_parameters.n_inequality;
	int n_i = b_parameters.n_interface;
	int n_p = b_parameters.n_planar;
	int n_t = b_parameters.n_tangent;

	// build polynomial blocks
	// | A PT |
	// | P 0  |
	// start with P
	for (int j = 0; j < (int)poly_matrix.size(); j++ ){
		for (int k = 0; k < (int)poly_matrix[j].size(); k++ ){
			interpolation_matrix[n_ie + n_i + 3*n_p + n_t + j][k] = poly_matrix[j][k];
			interpolation_matrix[k][n_ie + n_i + 3*n_p + n_t + j] = interpolation_matrix[n_ie + n_i + 3*n_p + n_t + j][k];
		}
	}

	for (int j = 0; j < (int)poly_matrix.size(); j++ ){
		for (int k = 0; k < (int)poly_matrix.size(); k++ ){
			interpolation_matrix[n_ie + n_i + 3*n_p + n_t + j][n_ie + n_i + 3*n_p + n_t + k] = 0;
		}
	}

	return true;
}

Continuous_Property::Continuous_Property(const model_parameters& m_p, const Basic_input& basic_i)
{
	// set GUI parameters and basic input (inequality, interface, planar, tangent) data members to class
	m_parameters = m_p;
	b_input = basic_i;

	_iteration = 0;
}


Polynomial_Basis * Continuous_Property::create_polynomial_basis( const int &poly_order )
{
	if (poly_order == 0) return new Poly_Zero;
	else if (poly_order == 1) return new Poly_First;
	else return new Poly_Second;
}

bool Continuous_Property::get_method_parameters()
{
	// # of constraints for each constraint type ...
	b_parameters.n_interface = (int)b_input.itrface->size();
	b_parameters.n_inequality = 0;
	b_parameters.n_planar = 0;
	b_parameters.n_tangent = 0;
	// Total number of constraints ...
	b_parameters.n_constraints = b_parameters.n_interface +	b_parameters.n_inequality + 3*b_parameters.n_planar + b_parameters.n_tangent;
	// Total number of equality constraints
	b_parameters.n_equality = b_parameters.n_interface + 3*b_parameters.n_planar + b_parameters.n_tangent;

	// polynomial parameters ...
	if (b_parameters.n_inequality == 0)
	{
		b_parameters.poly_term = false; // NOTE: May want to have this as an option when using SPD functions
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
	b_parameters.n_poly_terms = 0; // for 3D only...

	return true;
}

bool Continuous_Property::setup_system_solver()
{

	int nrows = b_parameters.n_equality + b_parameters.n_poly_terms;
	std::vector<double> equality_values;
	get_equality_values(equality_values);
	if ((int)equality_values.size() != nrows) return false;
	std::vector< std::vector <double> > interpolation_matrix = Math_methods::make_std_matrix<double>(nrows,nrows);
	if (!get_interpolation_matrix(interpolation_matrix)) return false;
	Linear_LU_decomposition *llu = new Linear_LU_decomposition(interpolation_matrix,equality_values);
	if (!llu->solve()) return false;
	solver = llu;

	return true;
}

// bool Continuous_Property::get_minimial_greedy_input( Basic_input &greedy_input )
// {
// 	// get extremal interface points
// 	// Note:
// 	// if poly order = 1 find 4 interface points nicely sampling the volume 
// 	// if poly order = 2 find 6 interface points nicely sampling the volume
// 	std::vector < Interface > extremal_interace_pts;
// 	std::vector < Point > pts(b_input.interface.begin(),b_input.interface.end()); 
// 	std::vector < int > interface_indices = get_extremal_point_data_indices_from_points(pts);
// 	if ( (int)interface_indices.size() < b_parameters.n_poly_terms) return false;
// 	for (int j = 0; j < b_parameters.n_poly_terms; j++ ) extremal_interace_pts.push_back(b_input.interface[interface_indices[j]]);
// 
// 	greedy_input.interface = extremal_interace_pts;
// 
// 	return true;
// }

bool Continuous_Property::measure_residuals( Basic_input &input )
{
	if (solver == NULL) return false;

	// inequalities points
	for (int j = 0; j < (int)input.inequality->size(); j++ ){
		eval_scalar_interpolant_at_point(input.inequality->at(j));
		if (input.inequality->at(j).level() >= 0)
		{
			if (input.inequality->at(j).scalar_field() >= 0) input.inequality->at(j).setResidual(true);
			else input.inequality->at(j).setResidual(false);
		}
		else
		{
			if (input.inequality->at(j).scalar_field() < 0) input.inequality->at(j).setResidual(true);
			else input.inequality->at(j).setResidual(false);
		}
	}
	// interface points
	for (int j = 0; j < (int)input.itrface->size(); j++ ){
		eval_scalar_interpolant_at_point(input.itrface->at(j));
		input.itrface->at(j).setResidual(abs(input.itrface->at(j).scalar_field() - input.itrface->at(j).level()));
	}
	// planar points
	for (int j = 0; j < (int)input.planar->size(); j++ ){
		eval_vector_interpolant_at_point(input.planar->at(j));
		double angle = 0.0;
		std::vector<double> v1;
		std::vector<double> v2;
		v1.push_back(input.planar->at(j).nx());
		v1.push_back(input.planar->at(j).ny());
		v1.push_back(input.planar->at(j).nz());
		v2.push_back(input.planar->at(j).nx_interp());
		v2.push_back(input.planar->at(j).ny_interp());
		v2.push_back(input.planar->at(j).nz_interp());
		Math_methods::angle_btw_2_vectors(v1,v2,angle);
		input.planar->at(j).setResidual( angle);
	}
	// tangent points
	for (int j = 0; j < (int)input.tangent->size(); j++ ){
		eval_vector_interpolant_at_point(input.tangent->at(j));
		double angle = 0.0;
		std::vector<double> v1;
		std::vector<double> v2;
		v1.push_back(input.tangent->at(j).tx());
		v1.push_back(input.tangent->at(j).ty());
		v1.push_back(input.tangent->at(j).tz());
		v2.push_back(input.tangent->at(j).nx_interp());
		v2.push_back(input.tangent->at(j).ny_interp());
		v2.push_back(input.tangent->at(j).nz_interp());
		Math_methods::angle_btw_2_vectors<double>(v1,v2,angle);
		input.planar->at(j).setResidual( angle);
	}

	return true;
}

bool Continuous_Property::append_greedy_input(Basic_input &input )
{
	// planar > tangent > interface > inequalities

	double r2d = 57.29577951308232;

	// PLANAR Observations
	std::vector < double > large_planar_residuals;
	std::vector < int > large_planar_residuals_indices;
	for (int j = 0; j < (int)input.planar->size(); j++ ){
		double grad_err = input.planar->at(j).residual()*r2d;
		if ( grad_err > m_parameters.gradient_slack )
		{
			large_planar_residuals.push_back( grad_err );
			large_planar_residuals_indices.push_back(j);
		}
	}
	if ( large_planar_residuals.size() != 0)
	{
		Math_methods::sort_vector_w_index(large_planar_residuals,large_planar_residuals_indices);
		b_input.planar->push_back(input.planar->at(large_planar_residuals_indices[large_planar_residuals.size() - 1]));
		return true;
	}
	// TANGENT Observations
	for (int j = 0; j < (int)input.tangent->size(); j++){
		if ( input.tangent->at(j).residual()*r2d > m_parameters.gradient_slack )
		{
			this->b_input.tangent->push_back(input.tangent->at(j));
			return true;
		}
	}
	// INTERFACE Observations
	std::vector < double > large_interface_residuals;
	std::vector < int > large_interface_residuals_indices;
	for (int j = 0; j < (int)input.itrface->size(); j++){
		double interface_err = input.itrface->at(j).residual();
		if ( interface_err > m_parameters.interface_slack )
		{
			large_interface_residuals.push_back(interface_err);
			large_interface_residuals_indices.push_back(j);
		}
	}
	if ( large_interface_residuals.size() != 0)
	{
		Math_methods::sort_vector_w_index(large_interface_residuals,large_interface_residuals_indices);
		b_input.itrface->push_back(input.itrface->at(large_interface_residuals_indices[large_interface_residuals.size() - 1]));
		return true;
	}
	// INEQUALITIES Observations
	for (int j = 0; j < (int)input.inequality->size(); j++){
		if ( !input.inequality->at(j).residual() )
		{
			b_input.inequality->push_back(input.inequality->at(j));
			return true;
		}
	}

	return false;
}

void Continuous_Property::eval_scalar_interpolant_at_point( Point &p )
{
	int n_i = b_parameters.n_interface;
	int n_p = b_parameters.n_planar;
	int n_t = b_parameters.n_tangent;

	Kernel *kernel_j = kernel->clone();
	double elemsum_1 = 0.0;
	double elemsum_2 = 0.0;
	double elemsum_3 = 0.0;
	double elemsum_4 = 0.0;
	double poly = 0.0;
	for (int k = 0; k < n_i; k++ ){
		kernel_j->set_points(p, b_input.itrface->at(k));
		elemsum_2 += solver->weights[k] * kernel_j->basis_pt_pt();
	}
	for (int k = 0; k < n_p; k++ ){
		kernel_j->set_points(p, b_input.planar->at(k));
		elemsum_3 += solver->weights[n_i + 3*k] * kernel_j->basis_pt_planar_x();
		elemsum_3 += solver->weights[n_i + 3*k + 1] * kernel_j->basis_pt_planar_y();
		elemsum_3 += solver->weights[n_i + 3*k + 2] * kernel_j->basis_pt_planar_z();
	}
	for (int k = 0; k < n_t; k++ ){
		kernel_j->set_points(p, b_input.tangent->at(k));
		elemsum_4 += solver->weights[n_i + 3*n_p + k] * kernel_j->basis_pt_tangent();
	}
	if (b_parameters.poly_term)
	{
		Polynomial_Basis *p_basis_j = p_basis->clone();
		p_basis_j->set_point(p);
		std::vector<double> b = p_basis_j->basis();
		for (int k = 0; k < (int)b.size(); k++ ) poly += b[k] * solver->weights[n_i + 3*n_p + n_t + k];
		delete p_basis_j;
	}
	p.set_scalar_field(elemsum_1 + elemsum_2 + elemsum_3 +elemsum_4 + poly);
	delete kernel_j;
}

void Continuous_Property::eval_vector_interpolant_at_point( Point &p )
{
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
	for (int k = 0; k < n_i; k++ ){
		kernel->set_points(p, b_input.itrface->at(k));
		elemsum_1_x += solver->weights[k] * kernel->basis_planar_x_pt();
		elemsum_1_y += solver->weights[k] * kernel->basis_planar_y_pt();
		elemsum_1_z += solver->weights[k] * kernel->basis_planar_z_pt();
	}
	// normal constraints
	for (int k = 0; k < n_p; k++ ){
		kernel->set_points(p, b_input.planar->at(k));
		elemsum_2_x += solver->weights[n_i + 0 + 3*k] * kernel->basis_planar_planar(Parameter_Types::DXDX);
		elemsum_2_x += solver->weights[n_i + 1 + 3*k] * kernel->basis_planar_planar(Parameter_Types::DXDY);
		elemsum_2_x += solver->weights[n_i + 2 + 3*k] * kernel->basis_planar_planar(Parameter_Types::DXDZ);
		elemsum_2_y += solver->weights[n_i + 0 + 3*k] * kernel->basis_planar_planar(Parameter_Types::DYDX);
		elemsum_2_y += solver->weights[n_i + 1 + 3*k] * kernel->basis_planar_planar(Parameter_Types::DYDY);
		elemsum_2_y += solver->weights[n_i + 2 + 3*k] * kernel->basis_planar_planar(Parameter_Types::DYDZ);
		elemsum_2_z += solver->weights[n_i + 0 + 3*k] * kernel->basis_planar_planar(Parameter_Types::DZDX);
		elemsum_2_z += solver->weights[n_i + 1 + 3*k] * kernel->basis_planar_planar(Parameter_Types::DZDY);
		elemsum_2_z += solver->weights[n_i + 2 + 3*k] * kernel->basis_planar_planar(Parameter_Types::DZDZ);
	}
	// tangent constraints
	for (int k = 0; k < n_t; k++ ){
		kernel->set_points(p, b_input.tangent->at(k));
		elemsum_3_x += solver->weights[n_i + 3*n_p + k] * kernel->basis_planar_tangent(Parameter_Types::DX);
		elemsum_3_y += solver->weights[n_i + 3*n_p + k] * kernel->basis_planar_tangent(Parameter_Types::DY);
		elemsum_3_z += solver->weights[n_i + 3*n_p + k] * kernel->basis_planar_tangent(Parameter_Types::DZ);
	}
	if (b_parameters.poly_term)
	{
		Polynomial_Basis *p_basis_j = p_basis->clone();
		p_basis_j->set_point(p);
		std::vector<double> bx = p_basis_j->dx();
		std::vector<double> by = p_basis_j->dy();
		std::vector<double> bz = p_basis_j->dz();
		for (int k = 0; k < (int)bx.size(); k++ ){
			poly_x += bx[k] * solver->weights[n_i + 3*n_p + n_t + k];
			poly_y += by[k] * solver->weights[n_i + 3*n_p + n_t + k];
			poly_z += bz[k] * solver->weights[n_i + 3*n_p + n_t + k];
		}
		delete p_basis_j;
	}
	double nx = elemsum_1_x + elemsum_2_x + elemsum_3_x + poly_x;
	double ny = elemsum_1_y + elemsum_2_y + elemsum_3_y + poly_y;
	double nz = elemsum_1_z + elemsum_2_z + elemsum_3_z + poly_z;
	p.set_vector_field(nx,ny,nz);
	delete kernel_j;
}

bool Continuous_Property::get_equality_values( std::vector<double> &equality_values )
{
	for (int j = 0; j < (int)b_input.itrface->size(); j++) equality_values.push_back(b_input.itrface->at(j).level());
	for (int j = 0; j < (int)b_input.planar->size(); j++){
		equality_values.push_back(b_input.planar->at(j).nx());
		equality_values.push_back(b_input.planar->at(j).ny());
		equality_values.push_back(b_input.planar->at(j).nz());
	}
	for (int j = 0; j < (int)b_input.tangent->size(); j++) equality_values.push_back(0.0);
	if (b_parameters.poly_term) for (int j = 0; j < (int)b_parameters.n_poly_terms; j++ ) equality_values.push_back(0.0);

	return true;
}

bool Continuous_Property::process_input_data()
{
	if ((int)b_input.itrface->size() == 0) return false;

	//b_input.transform_data_to_4D_isotropic_space();
	return true;
}

bool Continuous_Property::get_interpolation_matrix( std::vector< std::vector <double> > &interpolation_matrix )
{
	int n_i = b_parameters.n_interface;
	int n_p = b_parameters.n_planar;
	int n_t = b_parameters.n_tangent;

	// Row and Column constraint order : interface (itr) -> planar (p_x,p_y,p_z) -> tangent (t)

	// Base Matrix Structure
	// | itr/itr itr/p_y itr/p_y itr/p_z itr/t |
	// | p_x/itr p_x/p_x p_x/p_y p_x/p_z p_x/t |
	// | p_y/itr p_y/p_x p_y/p_y p_y/p_z p_y/t |
	// | p_z/itr p_z/p_x p_z/p_y p_z/p_z p_z/t |
	// |   t/itr   t/p_x   t/p_y   t/p_z   t/t |

	// Interface Constraints:
	for (int j = 0; j < n_i; j++ ){
		// Row:interface/Column:interface block
		for (int k = 0; k < n_i; k++ ){
			kernel->set_points(b_input.itrface->at(j), b_input.itrface->at(k));
			interpolation_matrix[j][k] = kernel->basis_pt_pt();
		}
		// Row:interface/Column:planar block
		for (int k = 0; k < n_p; k++ ){
			kernel->set_points(b_input.itrface->at(j), b_input.planar->at(k));
			interpolation_matrix[j][3*k + n_i] = kernel->basis_pt_planar_x();
			interpolation_matrix[j][3*k + n_i + 1] = kernel->basis_pt_planar_y();
			interpolation_matrix[j][3*k + n_i + 2] = kernel->basis_pt_planar_z();
		}
		// Row:interface/Column:tangent block
		for (int k = 0; k < n_t; k++ ){
			kernel->set_points(b_input.itrface->at(j), b_input.tangent->at(k));
			interpolation_matrix[j][n_i + 3*n_p + k] = kernel->basis_pt_tangent();
		}
	}
	// Planar Constraints
	for (int j = 0; j < n_p; j++ ){
		// Row:planar/Column:interface block
		for (int k = 0; k < n_i; k++ ){
			kernel->set_points(b_input.planar->at(j), b_input.itrface->at(k));
			interpolation_matrix[3*j + n_i][k] = kernel->basis_planar_x_pt();
			interpolation_matrix[3*j + n_i + 1][k] = kernel->basis_planar_y_pt();
			interpolation_matrix[3*j + n_i + 2][k] = kernel->basis_planar_z_pt();
		}
		// Row:planar/Column:planar block
		for (int k = 0; k < n_p; k++ ){
			kernel->set_points(b_input.planar->at(j), b_input.planar->at(k));
			interpolation_matrix[3*j + n_i][3*k + n_i] = kernel->basis_planar_planar(Parameter_Types::DXDX);
			interpolation_matrix[3*j + n_i][3*k + n_i + 1] = kernel->basis_planar_planar(Parameter_Types::DXDY);
			interpolation_matrix[3*j + n_i][3*k + n_i + 2] = kernel->basis_planar_planar(Parameter_Types::DXDZ);
			interpolation_matrix[3*j + n_i + 1][3*k + n_i] = kernel->basis_planar_planar(Parameter_Types::DYDX);
			interpolation_matrix[3*j + n_i + 1][3*k + n_i + 1] = kernel->basis_planar_planar(Parameter_Types::DYDY);
			interpolation_matrix[3*j + n_i + 1][3*k + n_i + 2] = kernel->basis_planar_planar(Parameter_Types::DYDZ);
			interpolation_matrix[3*j + n_i + 2][3*k + n_i] = kernel->basis_planar_planar(Parameter_Types::DZDX);
			interpolation_matrix[3*j + n_i + 2][3*k + n_i + 1] = kernel->basis_planar_planar(Parameter_Types::DZDY);
			interpolation_matrix[3*j + n_i + 2][3*k + n_i + 2] = kernel->basis_planar_planar(Parameter_Types::DZDZ);
		}
		// Row:planar/Column:tangent block
		for (int k = 0; k < n_t; k++ ){
			kernel->set_points(b_input.planar->at(j), b_input.tangent->at(k));
			interpolation_matrix[3*j + n_i][n_i + 3*n_p + k] = kernel->basis_planar_tangent(Parameter_Types::DX);
			interpolation_matrix[3*j + n_i + 1][n_i + 3*n_p + k] = kernel->basis_planar_tangent(Parameter_Types::DY);
			interpolation_matrix[3*j + n_i + 2][n_i + 3*n_p + k] = kernel->basis_planar_tangent(Parameter_Types::DZ);
		}
	}
	// Tangent Constraints 
	for (int j = 0; j < n_t; j++ ){
		// Row:tangent/Column:interface block
		for (int k = 0; k < n_i; k++ ){
			kernel->set_points(b_input.tangent->at(j), b_input.itrface->at(k));
			interpolation_matrix[j + n_i + 3*n_p][k] = kernel->basis_tangent_pt();
		}
		// Row:tangent/Column:planar block
		for (int k = 0; k < n_p; k++ ){
			kernel->set_points(b_input.tangent->at(j), b_input.planar->at(k));
			interpolation_matrix[j + n_i + 3*n_p][3*k + n_i] = kernel->basis_tangent_planar(Parameter_Types::DX);
			interpolation_matrix[j + n_i + 3*n_p][3*k + n_i + 1] = kernel->basis_tangent_planar(Parameter_Types::DY);
			interpolation_matrix[j + n_i + 3*n_p][3*k + n_i + 2] = kernel->basis_tangent_planar(Parameter_Types::DZ);
		}
		// Row:tangent/Column:tangent block
		for (int k = 0; k < n_t; k++ ){
			kernel->set_points(b_input.tangent->at(j), b_input.tangent->at(k));
			interpolation_matrix[j + n_i + 3*n_p][n_i + 3*n_p + k] = kernel->basis_tangent_tangent();
		}
	}

	// build polynomial blocks if required
	// | A PT |
	// | P 0  |
	if (b_parameters.poly_term)
	{
		std::vector < std::vector <double> > poly_matrix = Math_methods::make_std_matrix<double>(b_parameters.n_poly_terms,b_parameters.n_constraints);
		if (!_get_polynomial_matrix_block(poly_matrix)) return false;
		// fill remaining matrix blocks (P, PT, 0)
		if (!_insert_polynomial_matrix_blocks_in_interpolation_matrix(poly_matrix,interpolation_matrix)) return false;
	}
	// for testing ...
	// 	if (m_parameters.use_smoothing)
	// 	{
	// 		double dist_err = sqrt(m_parameters.interface_slack * m_parameters.interface_slack / 3);
	// 		Point one(dist_err,dist_err,dist_err);
	// 		Point two(0,0,0);
	// 		kernel->set_points(one,two);
	// 		double err = kernel->basis_pt_pt();
	// 		for (int j = 0; j < (n_ie + n_i); j++ ) interpolation_matrix[j][j] += err;
	// 	}

	return true;
}