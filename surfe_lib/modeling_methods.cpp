#include <math_methods.h>
#include <modeling_methods.h>
#include <matrix_solver.h>
#include <basis.h>
#include <algorithm>
#include <vector>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>

double round(double d)
{
	return floor(d + 0.5);
}

bool GRBF_Modelling_Methods::_update_interface_iso_values()
{
	// this is a messy method.
	// should really only be called if it is a Lajaunie method or Stratigraphic method
	// have put in safe guards to ensure no seg faults
	// function purpose:
	// When using increment approaches we do not know what the scalar field value will be
	// at interface points. Therefore, we wouldn't be able to do a iso-surface extraction.
	// To solve this issue, after the interpolant have been determined we evaluate the interpolant
	// at a interface point (in b_input.interface_test_points) for each interace. Then we
	// will know the right scalar field value for each interface to complete an iso-surface extraction.

	if ((int)b_input.interface_test_points->size() == 0) return false;

	// create evaluation point list from interface_test_points
	std::vector< Evaluation_Point > list;
	for (int j = 0; j < (int)b_input.interface_test_points->size(); j++){
		Evaluation_Point eval_pt(b_input.interface_test_points->at(j).x(), b_input.interface_test_points->at(j).y(), b_input.interface_test_points->at(j).z());
		list.push_back(eval_pt);
	}

	// evaluate the interpolant at these test points 
	if (solver != NULL) // check if we have a valid interpolant first
	{
		for (int j = 0; j < (int)b_input.interface_test_points->size(); j++) eval_scalar_interpolant_at_point(b_input.interface_test_points->at(j));
	}
	else return false;

	if ((int)b_input.interface_iso_values->size() != (int)b_input.interface_test_points->size()) return false;
	// update interface_iso_values to computed scalar field values
	for (int j = 0; j < (int)b_input.interface_iso_values->size(); j++) b_input.interface_iso_values->at(j) = b_input.interface_test_points->at(j).scalar_field();

	return true;
}

bool GRBF_Modelling_Methods::_output_greedy_debug_objects()
{
	if (solver == NULL) return false;
	else
	{
		if ((int)solver->weights.size() == 0) return false;
		else
		{
			if ((int)b_input.evaluation_pts->size()!= 0)	evaluate_scalar_interpolant();
			else return false;
		}
	}

	return true;
}

void GRBF_Modelling_Methods::_Progress( char message[], const int &step, const int &total )
{
	//progress width
	const int pwidth = 72;

	//minus label len
	int width = pwidth - strlen(message);
	int pos = (step * width) / total;
	int percent = (step * 100) / total;

	printf("%s[", message);

	//fill progress bar with =
	for (int i = 0; i < pos; i++)  printf("%c", '=');

	//fill progress bar with spaces
	printf("% *c", width - pos + 1, ']');
	printf(" %3d%%\r", percent);

}

void GRBF_Modelling_Methods::_Compute_Avg_NN_Dist_Data()
{
	int n_ie = b_parameters.n_inequality;
	int n_i = b_parameters.n_interface;
	int n_p = b_parameters.n_planar;
	int n_t = b_parameters.n_tangent;

	// for inequality points
	for (int j = 0; j < n_ie; j++ ){
		double min_dist = 10000000000000.0;
		for (int k = 0; k < n_ie; k++){
			if ( k != j )
			{
				double dist = distance_btw_pts(b_input.inequality->at(j),b_input.inequality->at(k));
				if ( dist < min_dist ) min_dist = dist;
			}
		}
		if ( n_ie == 1 ) min_dist = 0; // trap this edge case
		_avg_nn_dist_ie += min_dist;
	}
	if (n_ie != 0) _avg_nn_dist_ie /= n_ie;

	// for interface points
	for (int j = 0; j < n_i; j++ ){
		double min_dist = 10000000000000.0;
		for (int k = 0; k < n_i; k++){
			if ( k != j )
			{
				double dist = distance_btw_pts(b_input.itrface->at(j),b_input.itrface->at(k));
				if ( dist < min_dist ) min_dist = dist;
			}
		}
		if ( n_i == 1 ) min_dist = 0; // trap this edge case
		_avg_nn_dist_itr += min_dist;
	}
	if (n_i != 0) _avg_nn_dist_itr /= n_i;

	// for planar points
	for (int j = 0; j < n_p; j++ ){
		double min_dist = 10000000000000.0;
		for (int k = 0; k < n_p; k++){
			if ( k != j )
			{
				double dist = distance_btw_pts(b_input.planar->at(j),b_input.planar->at(k));
				if ( dist < min_dist ) min_dist = dist;
			}
		}
		if ( n_p == 1 ) min_dist = 0; // trap this edge case
		_avg_nn_dist_p += min_dist;
	}
	if (n_p != 0) _avg_nn_dist_p /= n_p;

	// for tangent points
	for (int j = 0; j < n_t; j++ ){
		double min_dist = 10000000000000.0;
		for (int k = 0; k < n_t; k++){
			if ( k != j )
			{
				double dist = distance_btw_pts(b_input.tangent->at(j),b_input.tangent->at(k));
				if ( dist < min_dist ) min_dist = dist;
			}
		}
		if ( n_t == 1 ) min_dist = 0; // trap this edge case
		_avg_nn_dist_t += min_dist;
	}
	if (n_t != 0) _avg_nn_dist_t /= n_t;
}

bool GRBF_Modelling_Methods::setup_basis_functions()
{
	rbf_kernel = this->create_rbf_kernel(m_parameters.basis_type,m_parameters.model_global_anisotropy);
	// check RBFKernel pointer
	if (rbf_kernel == NULL) return false;
	if (b_parameters.modified_basis)
	{
		if ((int)b_input.interface_point_lists->size() != 0) kernel = new Modified_Kernel(rbf_kernel, *b_input.interface_point_lists);
		else return false;
	}
	else kernel = rbf_kernel;

	return true;
}

bool GRBF_Modelling_Methods::evaluate_scalar_interpolant()
{
	if (solver == NULL) return false;
	else
	{
		if ((int)solver->weights.size() == 0) return false;
		else
		{
			int N = (int)b_input.evaluation_pts->size();
			int add = 0;
			int vv = round((double)N / 72.0); // 72.0 is the width of the progress bar 
			double factor = (100.0*(double)vv) / (double)N;

			#pragma omp parallel for schedule(dynamic)
			for (int j = 0; j < N; j++ ){
				eval_scalar_interpolant_at_point(b_input.evaluation_pts->at(j));
				if ( j % vv == 0 )
				{
#pragma omp atomic
					add++;
					int step = factor * add;
#pragma omp critical
					_Progress(" Computing Scalar field: ", step, 100 );
				}
			}
			cout<<endl;
		}
	}
	return true;
}

bool GRBF_Modelling_Methods::run_algorithm()
{
	clock_t tstart=clock();
	cout<<" Starting SURFE algorithm "<<endl;
	if ( !process_input_data()    ) return false;
	if ( !get_method_parameters() ) return false;
	if ( !setup_basis_functions() ) return false;
	if ( !setup_system_solver()   ) return false;
	if ( !evaluate_scalar_interpolant()  ) return false;
	cout<<" Total computation time = "<<((double)clock()-tstart)/CLOCKS_PER_SEC<<endl;

	return true;
}

bool GRBF_Modelling_Methods::run_greedy_algorithm()
{
	// check if there are non-zero errors permitted on the data
	if (m_parameters.interface_slack == 0 && m_parameters.gradient_slack == 0) return false;
	GRBF_Modelling_Methods *greedy_method = this->clone();

	greedy_method->_Compute_Avg_NN_Dist_Data();

	Basic_input greedy_input, excluded_input;
	// initialize starting data
	if (!get_minimial_and_excluded_input(greedy_input,excluded_input)) return false;
	greedy_method->b_input = greedy_input;
	greedy_method->b_input.evaluation_pts = b_input.evaluation_pts;

	bool converged = false;
	int iter = 0;
	while (!converged)
	{
		// run normal algorithm
		if ( !greedy_method->process_input_data()      ) return false;
		if ( !greedy_method->get_method_parameters()   ) return false;
		if ( !greedy_method->setup_basis_functions()   ) return false;
		if ( !greedy_method->setup_system_solver()     ) return false;

		// measure residuals
		if ( !greedy_method->measure_residuals(excluded_input)) return false;

		// debug: should output intermediate input constraints and modelled surface using those constraints
		// if ( !greedy_method->_output_greedy_debug_objects()) return false;

		// check topology : FUTURE

		// add appropriate data based on residuals
		if ( !greedy_method->append_greedy_input(excluded_input)) converged = true; // if no input is added convergence is assumed
		//iter++;
	}

	if ( !greedy_method->evaluate_scalar_interpolant() ) return false;

	b_input.evaluation_pts = greedy_method->b_input.evaluation_pts;

	return true;
}

bool GRBF_Modelling_Methods::get_equality_matrix( const std::vector< std::vector <double> > &interpolation_matrix, std::vector < std::vector < double > > &equality_matrix )
{
	if ((int)equality_matrix.size() == 0 || (int)equality_matrix.size() > (int)interpolation_matrix.size() || (int)equality_matrix[0].size() != (int)interpolation_matrix[0].size()) return false;
	int n_ie = (int)interpolation_matrix.size() - (int)equality_matrix.size();
	if (n_ie != b_parameters.n_inequality) return false;

	for (int j = 0; j < (int)equality_matrix.size(); j++ ){
		for (int k = 0; k < (int)equality_matrix[j].size(); k++ ){
			equality_matrix[j][k] = interpolation_matrix[j + n_ie][k];
		}
	}

	return true;
}

RBFKernel * GRBF_Modelling_Methods::create_rbf_kernel(const Parameter_Types::RBF &rbf_type, const bool &anisotropy)
{
	if (anisotropy)
	{
		if (rbf_type == Parameter_Types::Cubic) return new ACubic(*b_input.planar);
		else if (rbf_type == Parameter_Types::Gaussian) return new AGaussian(m_parameters.shape_parameter,*b_input.planar);
		else if (rbf_type == Parameter_Types::IMQ) return new AIMQ(m_parameters.shape_parameter, *b_input.planar);
		else if (rbf_type == Parameter_Types::MQ) return new AMQ(m_parameters.shape_parameter, *b_input.planar);
		else if (rbf_type == Parameter_Types::R) return new AR(*b_input.planar);
		else return new ATPS(*b_input.planar);
	}
	else
	{
		//if (b_input._weights.size() != 0) return new Scaled_Cubic(b_input._weights,b_input._points);
		if (rbf_type == Parameter_Types::Cubic) return new Cubic;
		else if (rbf_type == Parameter_Types::Gaussian) return new Gaussian(m_parameters.shape_parameter);
		else if (rbf_type == Parameter_Types::IMQ) return new IMQ(m_parameters.shape_parameter);
		else if (rbf_type == Parameter_Types::MQ) return new MQ(m_parameters.shape_parameter);
		else if (rbf_type == Parameter_Types::R) return new R;
		else return new TPS;
	}
}
