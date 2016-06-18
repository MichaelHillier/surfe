#include <stratigraphic_surfaces.h>
#include <math_methods.h>
#include <matrix_solver.h>
#include <basis.h>
#include <algorithm>
#include <vector>

#include <iostream>
#include <iomanip>
#include <fstream>

bool Stratigraphic_Surfaces::_get_increment_pairs()
{
	// 1) sequenced contacts
	_n_sequenced_interface_pairs = (int)b_input.interface_test_points->size() - 1;
	_increment_pairs->resize(_n_sequenced_interface_pairs);
	for (int j = 0; j < _n_sequenced_interface_pairs; j++ ){
		_increment_pairs->at(j).push_back(b_input.interface_test_points->at(j));
		_increment_pairs->at(j).push_back(b_input.interface_test_points->at(j + 1));
	}
	// 2) lithostratigraphic inequalities 
	_n_sequenced_inequality_pairs = 0;
	for (int j = 0; j < (int)b_input.inequality->size(); j++ ){
		std::vector < std::vector <Point> > strati_seq_ine_pairs = _get_lithostratigraphic_increment_pairs_for_inequality_point(b_input.inequality->at(j));
		for (int k = 0; k < (int)strati_seq_ine_pairs.size(); k++ ) _increment_pairs->push_back(strati_seq_ine_pairs[k]);
		int j_size = (int)strati_seq_ine_pairs.size();
		_n_sequenced_inequality_pairs += j_size;
	}
	// 3) the interface increment pairs
	_n_interface_pairs = 0;
	for (int j = 0; j < (int)b_input.interface_point_lists->size(); j++) _n_interface_pairs += ((int)b_input.interface_point_lists->at(j).size() - 1);
	for (int j = 0; j < (int)b_input.interface_point_lists->size(); j++ ){
		for (int k = 0; k < ((int)b_input.interface_point_lists->at(j).size() - 1); k++){
			std::vector<Point> interface_incr_p;
			interface_incr_p.push_back(b_input.interface_point_lists->at(j)[0]);
			interface_incr_p.push_back(b_input.interface_point_lists->at(j)[k + 1]);
			_increment_pairs->push_back(interface_incr_p);
		}
	}
	_n_increment_pairs = _n_sequenced_interface_pairs + _n_sequenced_inequality_pairs + _n_interface_pairs;

	return true;
}

std::vector < std::vector < Point > > Stratigraphic_Surfaces::_get_lithostratigraphic_increment_pairs_for_inequality_point( const Inequality &ie_pt )
{
	std::vector < std::vector < Point > > strati_incr_p;
	// if horizon is above current pt (ie_pt)
	// strati_incr_p[0][0] = above_horizon_pt
	// strati_incr_p[0][1] = ie_pt
	// if horizon is below current pt
	// strati_incr_p[1][0] = ie_pt
	// strati_incr_p[1][1] = below_horizon_pt
	// the above structure is needed to ensure the inequality constraint is s(p1) > s (p2) e.g. p1 = strati_incr_p[1][0] p2 = strati_incr_p[1][1]
	int n_pair = 0;

	std::vector<double> horizon_levels;
	for (int j = 0; j < (int)b_input.interface_test_points->size(); j++) horizon_levels.push_back(b_input.interface_test_points->at(j).level());

	// are there horizons above ie_pt ?
	double level_above = _get_closest_horizon_level_above_given_level(ie_pt.level(),horizon_levels);
	if (level_above != NULL)
	{
		strati_incr_p.resize(n_pair + 1);
		for (int j = 0; j < (int)b_input.interface_test_points->size(); j++ ){
			if (b_input.interface_test_points->at(j).level() == level_above)
			{
				strati_incr_p[n_pair].push_back(b_input.interface_test_points->at(j));
				strati_incr_p[n_pair].push_back(ie_pt);
				break;
			}
		}
		n_pair++;
	}
	// are there horizons below ie_pt ?
	double level_below = _get_closest_horizon_level_below_given_level(ie_pt.level(),horizon_levels);
	if (level_below != NULL)
	{
		strati_incr_p.resize(n_pair + 1);
		for (int j = 0; j < (int)b_input.interface_test_points->size(); j++ ){
			if (b_input.interface_test_points->at(j).level() == level_below)
			{
				strati_incr_p[n_pair].push_back(ie_pt);
				strati_incr_p[n_pair].push_back(b_input.interface_test_points->at(j));
				break;
			}
		}
		n_pair++;
	}

	return strati_incr_p;
}

double Stratigraphic_Surfaces::_get_closest_horizon_level_above_given_level( const double &given_level,const std::vector<double> &horizon_levels )
{
	if ((int)horizon_levels.size() == 0) return NULL;
	std::vector<double> diff;
	std::vector<int> indx;

	for (int j = 0; j < (int)horizon_levels.size(); j++ ){
		diff.push_back(horizon_levels[j] - given_level);
		indx.push_back(j);
	}
	Math_methods::sort_vector_w_index(diff,indx);
	if (diff[0] > 0) return horizon_levels[indx[0]];
	else return NULL;
}

double Stratigraphic_Surfaces::_get_closest_horizon_level_below_given_level( const double &given_level,const std::vector<double> &horizon_levels )
{
	if ((int)horizon_levels.size() == 0) return NULL;
	std::vector<double> diff;
	std::vector<int> indx;

	for (int j = 0; j < (int)horizon_levels.size(); j++ ){
		diff.push_back(given_level - horizon_levels[j]);
		indx.push_back(j);
	}
	Math_methods::sort_vector_w_index(diff,indx);
	if (diff[0] > 0) return horizon_levels[indx[0]];
	else return NULL;
}

Stratigraphic_Surfaces::Stratigraphic_Surfaces(const model_parameters& m_p, const Basic_input& basic_i)
{
	// set GUI parameters and basic input (inequality, interface, planar, tangent) data members to class
	m_parameters = m_p;
	b_input = basic_i;

	_increment_pairs = new std::vector < std::vector < Point > >();
}

bool Stratigraphic_Surfaces::get_method_parameters()
{
	// # of constraints for each constraint type ...
	b_parameters.n_interface = (int)b_input.itrface->size();
	b_parameters.n_inequality = (int)b_input.inequality->size();
	b_parameters.n_planar = (int)b_input.planar->size();
	b_parameters.n_tangent = (int)b_input.tangent->size();
	// Total number of constraints ...
	b_parameters.n_constraints = _n_increment_pairs + 3*b_parameters.n_planar + b_parameters.n_tangent;
	// Total number of equality constraints
	b_parameters.n_equality = _n_interface_pairs + 3*b_parameters.n_planar + b_parameters.n_tangent;
	b_parameters.n_inequality = _n_sequenced_interface_pairs + _n_sequenced_inequality_pairs;

	// polynomial parameters ...
	b_parameters.poly_term = false;
	b_parameters.modified_basis = true;
	b_parameters.problem_type = Parameter_Types::Quadratic;

	int m = m_parameters.polynomial_order + 1;
	b_parameters.n_poly_terms = (int)(m*(m + 1)*(m + 2) / 6);

	return true;
}

bool Stratigraphic_Surfaces::process_input_data()
{
	if ((int)b_input.itrface->size() == 0) return false;
	else
	{
		if (!b_input.get_interface_data()) return false;
		if (!_get_increment_pairs()) return false;
	}
	return true;
}

bool Stratigraphic_Surfaces::get_equality_values( VectorXd &equality_values )
{
	int j = 0;
	int k = 0;
	int l = 0;
	for (j = 0; j < _n_interface_pairs; j++ ) equality_values(j) =  0.0;
	for (k = 0; k < (int)b_input.planar->size(); k++){
		equality_values(3 * k + j) = b_input.planar->at(j).nx();
		equality_values(3 * k + j + 1) = b_input.planar->at(j).ny();
		equality_values(3 * k + j + 2) = b_input.planar->at(j).nz();
	}
	for (l = 0; l < (int)b_input.tangent->size(); l++) equality_values(l + 3 * k + j) = 0.0;

	return true;
}

bool Stratigraphic_Surfaces::get_inequality_values( VectorXd &inequality_values )
{
	int j = 0;
	int k = 0;
	for (j = 0; j < _n_sequenced_interface_pairs; j++ ) inequality_values(j) = m_parameters.min_stratigraphic_thickness;
	for (k = 0; k < _n_sequenced_inequality_pairs; k++ ) inequality_values(k + j) = 0.0;

	return true;
}

bool Stratigraphic_Surfaces::get_interpolation_matrix( MatrixXd &interpolation_matrix )
{
	int n_ip = _n_increment_pairs;
	int n_p = b_parameters.n_planar;
	int n_t = b_parameters.n_tangent;

	// Row and Column constraint order : increment pair (ip) -> planar (p_x,p_y,p_z) -> tangent (t)
	// Note ip encapsulates 3 types of increment constraints :
	// 1) sequenced stratigraphic reference increments
	// 2) sequenced stratigraphic inequality increments
	// 3) interface increments
	// _increment_pair[] container must be built in the above order (e.g. 1),2),3))

	// Base Matrix Structure
	// |  ip/ip  ip/p_y  ip/p_y  ip/p_z  ip/t |
	// | p_x/ip p_x/p_x p_x/p_y p_x/p_z p_x/t |
	// | p_y/ip p_y/p_x p_y/p_y p_y/p_z p_y/t |
	// | p_z/ip p_z/p_x p_z/p_y p_z/p_z p_z/t |
	// |   t/ip   t/p_x   t/p_y   t/p_z   t/t |

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
			interpolation_matrix(j,k) = (v1 - v2) -  (v3 - v4);
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
			interpolation_matrix(j,3*k + n_ip)     = v1x - v2x;
			interpolation_matrix(j,3*k + n_ip + 1) = v1y - v2y;
			interpolation_matrix(j,3*k + n_ip + 2) = v1z - v2z;
		}
		// Row:interface increment pair/Column:tangent block
		for (int k = 0; k < n_t; k++ ){
			kernel->set_points(_increment_pairs->at(j)[0], b_input.tangent->at(k));
			double v1 = kernel->basis_pt_tangent();
			kernel->set_points(_increment_pairs->at(j)[1], b_input.tangent->at(k));
			double v2 = kernel->basis_pt_tangent();
			interpolation_matrix(j,n_ip + 3*n_p + k) = v1 - v2;
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
			interpolation_matrix(3*j + n_ip,k)     = v1x - v2x;
			interpolation_matrix(3*j + n_ip + 1,k) = v1y - v2y;
			interpolation_matrix(3*j + n_ip + 2,k) = v1z - v2z;
		}
		// Row:planar/Column:planar block
		for (int k = 0; k < n_p; k++ ){
			kernel->set_points(b_input.planar->at(j), b_input.planar->at(k));
			interpolation_matrix(3*j + n_ip,3*k + n_ip)         = kernel->basis_planar_planar(Parameter_Types::DXDX);
			interpolation_matrix(3*j + n_ip,3*k + n_ip + 1)     = kernel->basis_planar_planar(Parameter_Types::DXDY);
			interpolation_matrix(3*j + n_ip,3*k + n_ip + 2)     = kernel->basis_planar_planar(Parameter_Types::DXDZ);
			interpolation_matrix(3*j + n_ip + 1,3*k + n_ip)     = kernel->basis_planar_planar(Parameter_Types::DYDX);
			interpolation_matrix(3*j + n_ip + 1,3*k + n_ip + 1) = kernel->basis_planar_planar(Parameter_Types::DYDY);
			interpolation_matrix(3*j + n_ip + 1,3*k + n_ip + 2) = kernel->basis_planar_planar(Parameter_Types::DYDZ);
			interpolation_matrix(3*j + n_ip + 2,3*k + n_ip)     = kernel->basis_planar_planar(Parameter_Types::DZDX);
			interpolation_matrix(3*j + n_ip + 2,3*k + n_ip + 1) = kernel->basis_planar_planar(Parameter_Types::DZDY);
			interpolation_matrix(3*j + n_ip + 2,3*k + n_ip + 2) = kernel->basis_planar_planar(Parameter_Types::DZDZ);
		}
		// Row:planar/Column:tangent block
		for (int k = 0; k < n_t; k++ ){
			kernel->set_points(b_input.planar->at(j), b_input.tangent->at(k));
			interpolation_matrix(3*j + n_ip,n_ip + 3*n_p + k) = kernel->basis_planar_tangent(Parameter_Types::DX);
			interpolation_matrix(3*j + n_ip + 1,n_ip + 3*n_p + k) = kernel->basis_planar_tangent(Parameter_Types::DY);
			interpolation_matrix(3*j + n_ip + 2,n_ip + 3*n_p + k) = kernel->basis_planar_tangent(Parameter_Types::DZ);
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
			interpolation_matrix(j + n_ip + 3*n_p,k) = v1 - v2;
		}
		// Row:tangent/Column:planar block
		for (int k = 0; k < n_p; k++ ){
			kernel->set_points(b_input.tangent->at(j), b_input.planar->at(k));
			interpolation_matrix(j + n_ip + 3*n_p,3*k + n_ip)     = kernel->basis_tangent_planar(Parameter_Types::DX);
			interpolation_matrix(j + n_ip + 3*n_p,3*k + n_ip + 1) = kernel->basis_tangent_planar(Parameter_Types::DY);
			interpolation_matrix(j + n_ip + 3*n_p,3*k + n_ip + 2) = kernel->basis_tangent_planar(Parameter_Types::DZ);
		}
		// Row:tangent/Column:tangent block
		for (int k = 0; k < n_t; k++ ){
			kernel->set_points(b_input.tangent->at(j), b_input.tangent->at(k));
			interpolation_matrix(j + n_ip + 3*n_p,n_ip + 3*n_p + k) = kernel->basis_tangent_tangent();
		}
	}

	return true;
}

bool Stratigraphic_Surfaces::get_inequality_matrix( const MatrixXd &interpolation_matrix, MatrixXd &inequality_matrix )
{
	if (inequality_matrix.rows() == 0 || 
		inequality_matrix.rows()  > interpolation_matrix.rows() ||
		inequality_matrix.cols() != interpolation_matrix.cols() ) return false;

	for (int j = 0; j < (int)inequality_matrix.rows(); j++ ){
		for (int k = 0; k < (int)inequality_matrix.cols(); k++ ){
			inequality_matrix(j,k) = interpolation_matrix(j,k);
		}
	}

	return true;
}

bool Stratigraphic_Surfaces::setup_system_solver()
{
	// only way to solve this problem is via a quadratic optimization problem
	int n_ie = b_parameters.n_inequality;
	int n_e  = b_parameters.n_equality;
	int n_c  = b_parameters.n_constraints;

	VectorXd inequality_values(n_ie);
	get_inequality_values(inequality_values);

	VectorXd equality_values(n_e);
	get_equality_values(equality_values);

	MatrixXd interpolation_matrix(n_c,n_c);
	if (!get_interpolation_matrix(interpolation_matrix)) return false;

	MatrixXd inequality_matrix(n_ie,n_c);
	if (!get_inequality_matrix(interpolation_matrix,inequality_matrix)) return false;

	MatrixXd equality_matrix(n_e,n_c);
	if (!get_equality_matrix(interpolation_matrix,equality_matrix)) return false;

	Quadratic_Predictor_Corrector *qpc = new Quadratic_Predictor_Corrector(interpolation_matrix,equality_matrix,inequality_matrix,equality_values,inequality_values);
	if(!qpc->solve()) return false;
	solver = qpc;

	if (!_update_interface_iso_values()) return false;

	return true;
}

void Stratigraphic_Surfaces::eval_scalar_interpolant_at_point( Point &p )
{
	int n_ip = _n_increment_pairs;
	int n_p = b_parameters.n_planar;
	int n_t = b_parameters.n_tangent;

	Kernel *kernel_j = kernel->clone();
	double elemsum_1 = 0.0;
	double elemsum_2 = 0.0;
	double elemsum_3 = 0.0;
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
	p.set_scalar_field(elemsum_1 + elemsum_2 + elemsum_3);
	delete kernel_j;
}

void Stratigraphic_Surfaces::eval_vector_interpolant_at_point( Point &p )
{
	int n_ip = _n_increment_pairs;
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
	double nx = elemsum_1_x + elemsum_2_x + elemsum_3_x;
	double ny = elemsum_1_y + elemsum_2_y + elemsum_3_y;
	double nz = elemsum_1_z + elemsum_2_z + elemsum_3_z;
	p.set_vector_field(nx,ny,nz);
	delete kernel_j;
}
