#include <lajaunie.h>
#include <math_methods.h>
#include <matrix_solver.h>
#include <basis.h>
#include <algorithm>
#include <vector>

#include <iostream>
#include <iomanip>
#include <fstream>

bool Lajaunie_Approach::_get_polynomial_matrix_block( std::vector< std::vector <double> > &poly_matrix )
{
	int n_ip = _n_increment_pair;
	int n_p = b_parameters.n_planar;
	int n_t = b_parameters.n_tangent;

	p_basis = create_polynomial_basis(m_parameters.polynomial_order);

	if ((int)poly_matrix.size() != b_parameters.n_poly_terms ) return false;

	// for Interface Increment Pair Constraints:
	for (int j = 0; j < (int)_increment_pairs->size();j++ ){
		p_basis->set_point(_increment_pairs->at(j)[0]);
		std::vector<double> b1 = p_basis->basis();
		p_basis->set_point(_increment_pairs->at(j)[1]);
		std::vector<double> b2 = p_basis->basis(); 
		if ((int)b1.size() != b_parameters.n_poly_terms || (int)b2.size() != b_parameters.n_poly_terms ) return false;
		for (int k = 0; k < (int)b1.size(); k++ ) poly_matrix[k][j] = b1[k] - b2[k];
	}
	// for planar points ...
	for (int j = 0; j < n_p; j++ ){
		p_basis->set_point(b_input.planar->at(j));
		std::vector<double> bx = p_basis->dx();
		std::vector<double> by = p_basis->dy();
		std::vector<double> bz = p_basis->dz();
		if ((int)bx.size() != b_parameters.n_poly_terms || (int)by.size() != b_parameters.n_poly_terms || (int)bz.size() != b_parameters.n_poly_terms) return false;
		for (int k = 0; k < (int)bx.size(); k++ ){
			poly_matrix[k][n_ip + 3*j] = bx[k];
			poly_matrix[k][n_ip + 3*j + 1] = by[k];
			poly_matrix[k][n_ip + 3*j + 2] = bz[k];
		}
	}
	// for tangent points ...
	for (int j = 0; j < n_t; j++ ){
		p_basis->set_point(b_input.tangent->at(j));
		std::vector<double> bx = p_basis->dx();
		std::vector<double> by = p_basis->dy();
		std::vector<double> bz = p_basis->dz();
		if ( (int)bx.size() != b_parameters.n_poly_terms || (int)by.size() != b_parameters.n_poly_terms || (int)bz.size() != b_parameters.n_poly_terms ) return false;
		for (int k = 0; k < (int)bx.size(); k++ ){
			poly_matrix[k][n_ip + 3 * n_p + j] = b_input.tangent->at(j).tx()*bx[k] + b_input.tangent->at(j).ty()*by[k] + b_input.tangent->at(j).tz()*bz[k];
		}
	}

	return true;
}

bool Lajaunie_Approach::_insert_polynomial_matrix_blocks_in_interpolation_matrix( const std::vector< std::vector <double> > &poly_matrix, std::vector< std::vector <double> > &interpolation_matrix )
{
	int n_ip = _n_increment_pair;
	int n_p = b_parameters.n_planar;
	int n_t = b_parameters.n_tangent;

	if ((n_ip +3*n_p + n_t + (int)poly_matrix.size()) > (int)interpolation_matrix.size() ||
		(n_ip +3*n_p + n_t + (int)poly_matrix.size()) > (int)interpolation_matrix[0].size()) return false;
	// build polynomial blocks
	// | A PT |
	// | P 0  |
	// start with P
	for (int j = 0; j < (int)poly_matrix.size(); j++ ){
		for (int k = 0; k < (int)poly_matrix[j].size(); k++ ){
			interpolation_matrix[n_ip + 3*n_p + n_t + j][k] = poly_matrix[j][k];
			interpolation_matrix[k][n_ip + 3*n_p + n_t + j] = interpolation_matrix[n_ip + 3*n_p + n_t + j][k];
		}
	}

	for (int j = 0; j < (int)poly_matrix.size(); j++ ){
		for (int k = 0; k < (int)poly_matrix.size(); k++ ){
			interpolation_matrix[n_ip + 3*n_p + n_t + j][n_ip + 3*n_p + n_t + k] = 0;
		}
	}

	return true;
}

bool Lajaunie_Approach::_get_increment_pairs()
{
	// the interface increment pairs
	int _n_interface_pairs = 0;
	for (int j = 0; j < (int)b_input.interface_point_lists->size(); j++) _n_interface_pairs += ((int)b_input.interface_point_lists->at(j).size() - 1);
	for (int j = 0; j < (int)b_input.interface_point_lists->size(); j++ ){
		for (int k = 0; k < ((int)b_input.interface_point_lists->at(j).size() - 1); k++ ){
			std::vector<Point> interface_incr_p;
			interface_incr_p.push_back(b_input.interface_point_lists->at(j)[0]);
			interface_incr_p.push_back(b_input.interface_point_lists->at(j)[k + 1]);
			_increment_pairs->push_back(interface_incr_p);
		}
	}
	_n_increment_pair = _n_interface_pairs;

	return true;
}

Lajaunie_Approach::Lajaunie_Approach(const model_parameters& m_p, const Basic_input& basic_i)
{
	// set GUI parameters and basic input (interface, planar, tangent) data members to class
	m_parameters = m_p;
	b_input = basic_i;

	_increment_pairs = new std::vector < std::vector < Point > >();
}

bool Lajaunie_Approach::get_method_parameters()
{
	// # of constraints for each constraint type ...
	b_parameters.n_interface = (int)b_input.interface->size();
	b_parameters.n_inequality = 0; // no support for inequality. if there is inequalities use the stratigraphic horizon method
	b_parameters.n_planar = (int)b_input.planar->size();
	b_parameters.n_tangent = (int)b_input.tangent->size();
	// Total number of constraints ...
	b_parameters.n_constraints = _n_increment_pair + 3*b_parameters.n_planar + b_parameters.n_tangent;
	// Total number of equality constraints
	b_parameters.n_equality = _n_increment_pair + 3*b_parameters.n_planar + b_parameters.n_tangent;

	// polynomial parameters ...
	b_parameters.poly_term = true; // NOTE: May want to have this as an option when using SPD functions
	b_parameters.modified_basis = false;
	b_parameters.problem_type = Parameter_Types::Linear;

	int m = m_parameters.polynomial_order + 1;
	b_parameters.n_poly_terms = (int)(m*(m + 1)*(m + 2) / 6) - 1; // minus 1 due to the nature of the indepentent pair constraints and the vanishing of the constant term in the polynomial

	return true;
}

bool Lajaunie_Approach::process_input_data()
{
	if ((int)b_input.interface->size() == 0) return false;
	else
	{
		if (!b_input.get_interface_data()) return false;
		if (!_get_increment_pairs()) return false;
	}
	return true;
}

bool Lajaunie_Approach::setup_system_solver()
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

	if (!_update_interface_iso_values()) return false;

	return true;
}

bool Lajaunie_Approach::get_equality_values( std::vector<double> &equality_values )
{
	for (int j = 0; j < _n_increment_pair; j++ ) equality_values.push_back(0.0);
	for (int j = 0; j < (int)b_input.planar->size(); j++){
		equality_values.push_back(b_input.planar->at(j).nx());
		equality_values.push_back(b_input.planar->at(j).ny());
		equality_values.push_back(b_input.planar->at(j).nz());
	}
	for (int j = 0; j < (int)b_input.tangent->size(); j++) equality_values.push_back(0.0);
	if (b_parameters.poly_term) for (int j = 0; j < (int)b_parameters.n_poly_terms; j++ ) equality_values.push_back(0.0);

	return true;
}

bool Lajaunie_Approach::get_interpolation_matrix( std::vector< std::vector <double> > &interpolation_matrix )
{
	int n_ip = _n_increment_pair;
	int n_p = b_parameters.n_planar;
	int n_t = b_parameters.n_tangent;

	// Row and Column constraint order : interface increment pair (iip) -> planar (p_x,p_y,p_z) -> tangent (t)

	// Base Matrix Structure
	// | iip/iip iip/p_y iip/p_y iip/p_z iip/t |
	// | p_x/iip p_x/p_x p_x/p_y p_x/p_z p_x/t |
	// | p_y/iip p_y/p_x p_y/p_y p_y/p_z p_y/t |
	// | p_z/iip p_z/p_x p_z/p_y p_z/p_z p_z/t |
	// |   t/iip  t/p_x   t/p_y   t/p_z   t/t |

	// Interface Increment Pair Constraints:
	for (int j = 0; j < (int)_increment_pairs->size(); j++ ){
		// Row:interface increment pair/Column:interface increment pair block
		for (int k = 0; k < (int)_increment_pairs->size(); k++ ){
			kernel->set_points(_increment_pairs->at(j)[0], _increment_pairs->at(k)[0]);
			double v1 = kernel->basis_pt_pt();
			kernel->set_points(_increment_pairs->at(j)[0], _increment_pairs->at(k)[1]);
			double v2 = kernel->basis_pt_pt();
			kernel->set_points(_increment_pairs->at(j)[1], _increment_pairs->at(k)[0]);
			double v3 = kernel->basis_pt_pt();
			kernel->set_points(_increment_pairs->at(j)[1], _increment_pairs->at(k)[1]);
			double v4 = kernel->basis_pt_pt();
			interpolation_matrix[j][k] = (v1 - v2) -  (v3 - v4);
		}
		// Row:interface increment pair/Column:planar block
		for (int k = 0; k < n_p; k++ ){
			kernel->set_points(_increment_pairs->at(j)[0], b_input.planar->at(k));
			double v1x = kernel->basis_pt_planar_x();
			double v1y = kernel->basis_pt_planar_y();
			double v1z = kernel->basis_pt_planar_z();
			kernel->set_points(_increment_pairs->at(j)[1], b_input.planar->at(k));
			double v2x = kernel->basis_pt_planar_x();
			double v2y = kernel->basis_pt_planar_y();
			double v2z = kernel->basis_pt_planar_z();
			interpolation_matrix[j][3*k + n_ip]     = v1x - v2x;
			interpolation_matrix[j][3*k + n_ip + 1] = v1y - v2y;
			interpolation_matrix[j][3*k + n_ip + 2] = v1z - v2z;
		}
		// Row:interface increment pair/Column:tangent block
		for (int k = 0; k < n_t; k++ ){
			kernel->set_points(_increment_pairs->at(j)[0], b_input.tangent->at(k));
			double v1 = kernel->basis_pt_tangent();
			kernel->set_points(_increment_pairs->at(j)[1], b_input.tangent->at(k));
			double v2 = kernel->basis_pt_tangent();
			interpolation_matrix[j][n_ip + 3*n_p + k] = v1 - v2;
		}
	}
	// Planar Constraints
	for (int j = 0; j < n_p; j++ ){
		// Row:planar/Column:interface increment pair
		for (int k = 0; k < (int)_increment_pairs->size();k++ ){
			kernel->set_points(b_input.planar->at(j), _increment_pairs->at(k)[0]);
			double v1x = kernel->basis_planar_x_pt();
			double v1y = kernel->basis_planar_y_pt();
			double v1z = kernel->basis_planar_z_pt();
			kernel->set_points(b_input.planar->at(j), _increment_pairs->at(k)[1]);
			double v2x = kernel->basis_planar_x_pt();
			double v2y = kernel->basis_planar_y_pt();
			double v2z = kernel->basis_planar_z_pt();
			interpolation_matrix[3*j + n_ip][k]     = v1x - v2x;
			interpolation_matrix[3*j + n_ip + 1][k] = v1y - v2y;
			interpolation_matrix[3*j + n_ip + 2][k] = v1z - v2z;
		}
		// Row:planar/Column:planar block
		for (int k = 0; k < n_p; k++ ){
			kernel->set_points(b_input.planar->at(j), b_input.planar->at(k));
			interpolation_matrix[3*j + n_ip][3*k + n_ip] = kernel->basis_planar_planar(Parameter_Types::DXDX);
			interpolation_matrix[3*j + n_ip][3*k + n_ip + 1] = kernel->basis_planar_planar(Parameter_Types::DXDY);
			interpolation_matrix[3*j + n_ip][3*k + n_ip + 2] = kernel->basis_planar_planar(Parameter_Types::DXDZ);
			interpolation_matrix[3*j + n_ip + 1][3*k + n_ip] = kernel->basis_planar_planar(Parameter_Types::DYDX);
			interpolation_matrix[3*j + n_ip + 1][3*k + n_ip + 1] = kernel->basis_planar_planar(Parameter_Types::DYDY);
			interpolation_matrix[3*j + n_ip + 1][3*k + n_ip + 2] = kernel->basis_planar_planar(Parameter_Types::DYDZ);
			interpolation_matrix[3*j + n_ip + 2][3*k + n_ip] = kernel->basis_planar_planar(Parameter_Types::DZDX);
			interpolation_matrix[3*j + n_ip + 2][3*k + n_ip + 1] = kernel->basis_planar_planar(Parameter_Types::DZDY);
			interpolation_matrix[3*j + n_ip + 2][3*k + n_ip + 2] = kernel->basis_planar_planar(Parameter_Types::DZDZ);
		}
		// Row:planar/Column:tangent block
		for (int k = 0; k < n_t; k++ ){
			kernel->set_points(b_input.planar->at(j), b_input.tangent->at(k));
			interpolation_matrix[3*j + n_ip][n_ip + 3*n_p + k] = kernel->basis_planar_tangent(Parameter_Types::DX);
			interpolation_matrix[3*j + n_ip + 1][n_ip + 3*n_p + k] = kernel->basis_planar_tangent(Parameter_Types::DY);
			interpolation_matrix[3*j + n_ip + 2][n_ip + 3*n_p + k] = kernel->basis_planar_tangent(Parameter_Types::DZ);
		}
	}
	// Tangent Constraints 
	for (int j = 0; j < n_t; j++ ){
		// Row:tangent/Column:interface increment pair block
		for (int k = 0; k < (int)_increment_pairs->size();k++ ){
			kernel->set_points(b_input.tangent->at(j), _increment_pairs->at(k)[0]);
			double v1 = kernel->basis_tangent_pt();
			kernel->set_points(b_input.tangent->at(j), _increment_pairs->at(k)[1]);
			double v2 = kernel->basis_tangent_pt();
			interpolation_matrix[j + n_ip + 3*n_p][k] = v1 - v2;
		}
		// Row:tangent/Column:planar block
		for (int k = 0; k < n_p; k++ ){
			kernel->set_points(b_input.tangent->at(j), b_input.planar->at(k));
			interpolation_matrix[j + n_ip + 3*n_p][3*k + n_ip] = kernel->basis_tangent_planar(Parameter_Types::DX);
			interpolation_matrix[j + n_ip + 3*n_p][3*k + n_ip + 1] = kernel->basis_tangent_planar(Parameter_Types::DY);
			interpolation_matrix[j + n_ip + 3*n_p][3*k + n_ip + 2] = kernel->basis_tangent_planar(Parameter_Types::DZ);
		}
		// Row:tangent/Column:tangent block
		for (int k = 0; k < n_t; k++ ){
			kernel->set_points(b_input.tangent->at(j), b_input.tangent->at(k));
			interpolation_matrix[j + n_ip + 3*n_p][n_ip + 3*n_p + k] = kernel->basis_tangent_tangent();
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

// 	if (m_parameters.use_smoothing)
// 	{
// 		double dist_err = sqrt(m_parameters.interface_slack * m_parameters.interface_slack / 3);
// 		Point one(dist_err,dist_err,dist_err);
// 		Point two(0,0,0);
// 		kernel->set_points(one,two);
// 		double err = kernel->basis_pt_pt();
// 		for (int j = 0; j < (int)interpolation_matrix.size(); j++ ) interpolation_matrix[j][j] += err;
// 	}

// 	std::ofstream out("E:\interpolation_mat.txt");
// 
// 	out.precision(15); // sets how many digits for each variable is outputted to the output file
// 
// 
// 	for (int j=0;j<(int)interpolation_matrix.size();j++){
// 		for (int k=0;k<(int)interpolation_matrix[j].size();k++){
// 			out<<"  "<<interpolation_matrix[j][k];
// 		}
// 		out<<std::endl;
// 	}


	return true;
}

Polynomial_Basis * Lajaunie_Approach::create_polynomial_basis( const int &poly_order )
{
	if (poly_order == 0) return new Poly_Zero(true);
	else if (poly_order == 1) return new Poly_First(true);
	else return new Poly_Second(true);
}


void Lajaunie_Approach::eval_scalar_interpolant_at_point( Point &p )
{
	int n_ip = _n_increment_pair;
	int n_p = b_parameters.n_planar;
	int n_t = b_parameters.n_tangent;

	Kernel *kernel_j = kernel->clone();
	double elemsum_1 = 0.0;
	double elemsum_2 = 0.0;
	double elemsum_3 = 0.0;
	double poly = 0.0;
	for (int k = 0; k < (int)_increment_pairs->size();k++ ){
		kernel_j->set_points(p, _increment_pairs->at(k)[0]);
		double v1 = kernel_j->basis_pt_pt();
		kernel_j->set_points(p, _increment_pairs->at(k)[1]);
		double v2 = kernel_j->basis_pt_pt();
		elemsum_1 += solver->weights[k] * (v1 - v2);
	}
	for (int k = 0; k < n_p; k++ ){
		kernel_j->set_points(p, b_input.planar->at(k));
		elemsum_2 += solver->weights[n_ip + 3*k] * kernel_j->basis_pt_planar_x();
		elemsum_2 += solver->weights[n_ip + 3*k + 1] * kernel_j->basis_pt_planar_y();
		elemsum_2 += solver->weights[n_ip + 3*k + 2] * kernel_j->basis_pt_planar_z();
	}
	for (int k = 0; k < n_t; k++ ){
		kernel_j->set_points(p, b_input.tangent->at(k));
		elemsum_3 += solver->weights[n_ip + 3*n_p + k] * kernel_j->basis_pt_tangent();
	}
	if (b_parameters.poly_term)
	{
		Polynomial_Basis *p_basis_j = p_basis->clone();
		p_basis_j->set_point(p);
		std::vector<double> b = p_basis_j->basis();
		for (int k = 0; k < (int)b.size(); k++ ) poly += b[k] * solver->weights[n_ip + 3*n_p + n_t + k];
		delete p_basis_j;
	}
	p.set_scalar_field(elemsum_1 + elemsum_2 + elemsum_3 + poly);
	delete kernel_j;
}

void Lajaunie_Approach::eval_vector_interpolant_at_point( Point &p )
{
	int n_ip = _n_increment_pair;
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
	for (int k = 0; k < n_ip;k++ ){
		kernel->set_points(p, _increment_pairs->at(k)[0]);
		double v1x = kernel->basis_planar_x_pt();
		double v1y = kernel->basis_planar_y_pt();
		double v1z = kernel->basis_planar_z_pt();
		kernel->set_points(p, _increment_pairs->at(k)[1]);
		double v2x = kernel->basis_planar_x_pt();
		double v2y = kernel->basis_planar_y_pt();
		double v2z = kernel->basis_planar_z_pt();
		elemsum_1_x += solver->weights[k]*(v1x - v2x);
		elemsum_1_y += solver->weights[k]*(v1y - v2y);
		elemsum_1_z += solver->weights[k]*(v1z - v2z);
	}
	// planar constraints
	for (int k = 0; k < n_p; k++ ){
		kernel->set_points(p, b_input.planar->at(k));
		elemsum_2_x += solver->weights[n_ip + 0 + 3*k]*kernel->basis_planar_planar(Parameter_Types::DXDX);
		elemsum_2_x += solver->weights[n_ip + 1 + 3*k]*kernel->basis_planar_planar(Parameter_Types::DXDY);
		elemsum_2_x += solver->weights[n_ip + 2 + 3*k]*kernel->basis_planar_planar(Parameter_Types::DXDZ);
		elemsum_2_y += solver->weights[n_ip + 0 + 3*k]*kernel->basis_planar_planar(Parameter_Types::DYDX);
		elemsum_2_y += solver->weights[n_ip + 1 + 3*k]*kernel->basis_planar_planar(Parameter_Types::DYDY);
		elemsum_2_y += solver->weights[n_ip + 2 + 3*k]*kernel->basis_planar_planar(Parameter_Types::DYDZ);
		elemsum_2_z += solver->weights[n_ip + 0 + 3*k]*kernel->basis_planar_planar(Parameter_Types::DZDX);
		elemsum_2_z += solver->weights[n_ip + 1 + 3*k]*kernel->basis_planar_planar(Parameter_Types::DZDY);
		elemsum_2_z += solver->weights[n_ip + 2 + 3*k]*kernel->basis_planar_planar(Parameter_Types::DZDZ);
	}
	// Tangent constraints
	for (int k = 0; k < n_t; k++ ){
		kernel->set_points(p, b_input.tangent->at(k));
		elemsum_3_x += solver->weights[n_ip + 3*n_p + k]*kernel->basis_planar_tangent(Parameter_Types::DX);
		elemsum_3_y += solver->weights[n_ip + 3*n_p + k]*kernel->basis_planar_tangent(Parameter_Types::DY);
		elemsum_3_z += solver->weights[n_ip + 3*n_p + k]*kernel->basis_planar_tangent(Parameter_Types::DZ);
	}
	if (b_parameters.poly_term)
	{
		Polynomial_Basis *p_basis_j = p_basis->clone();
		p_basis_j->set_point(p);
		std::vector<double> bx = p_basis_j->dx();
		std::vector<double> by = p_basis_j->dy();
		std::vector<double> bz = p_basis_j->dz();
		for (int k = 0; k < (int)bx.size(); k++ ){
			poly_x += bx[k] * solver->weights[n_ip + 3*n_p + n_t + k];
			poly_y += by[k] * solver->weights[n_ip + 3*n_p + n_t + k];
			poly_z += bz[k] * solver->weights[n_ip + 3*n_p + n_t + k];
		}
		delete p_basis_j;
	}
	double nx = elemsum_1_x + elemsum_2_x + elemsum_3_x + poly_x;
	double ny = elemsum_1_y + elemsum_2_y + elemsum_3_y + poly_y;
	double nz = elemsum_1_z + elemsum_2_z + elemsum_3_z + poly_z;
	p.set_vector_field(nx,ny,nz);
	delete kernel_j;
}
