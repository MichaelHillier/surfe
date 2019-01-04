#include <math_methods.h>
#include <modeling_methods.h>
#include <matrix_solver.h>
#include <basis.h>
#include <single_surface.h>
#include <lajaunie.h>
#include <stratigraphic_surfaces.h>
#include <continuous_property.h>

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

	// evaluate the interpolant at these interface_test_points
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

bool GRBF_Modelling_Methods::setup_basis_functions()
{
	rbf_kernel = this->create_rbf_kernel(m_parameters.basis_type,m_parameters.model_global_anisotropy);
	// check RBFKernel pointer
	if (rbf_kernel == NULL) return false;
	if (b_parameters.modified_basis)
	{
		if ((int)b_input.interface_point_lists->size() != 0) kernel = new Modified_Kernel(rbf_kernel, *b_input.interface_point_lists);
		else if ((int)b_input.tangent->size() != 0)
		{
			std::vector<Point> pts;
			for (int j = 0; j < b_input.tangent->size(); j++){
				Point apt(b_input.tangent->at(j).x(), b_input.tangent->at(j).y(), b_input.tangent->at(j).z());
				pts.push_back(apt);
			}
			kernel = new Modified_Kernel(rbf_kernel, pts);
		}
		else return false;
	}
	else kernel = rbf_kernel;

	return true;
}

bool GRBF_Modelling_Methods::check_interpolant()
{
	for (int j = 0; j < b_input.itrface->size(); j++ ){
		cout<<" Interface["<<j<<"]: "<<endl;
		eval_scalar_interpolant_at_point(b_input.itrface->at(j));
		cout<<"	Scalar field = "<<b_input.itrface->at(j).scalar_field()<<endl;
	}

	for (int j = 0; j < b_input.planar->size(); j++ ){
		eval_vector_interpolant_at_point(b_input.planar->at(j));
		double vf[3] = {b_input.planar->at(j).nx_interp(),b_input.planar->at(j).ny_interp(),b_input.planar->at(j).nz_interp()};
		cout<<" Planar["<<j<<"]: "<<endl;
		cout<<"	Nx = "<<b_input.planar->at(j).nx()<<" Nx interpolated = "<<b_input.planar->at(j).nx_interp()<<endl;
		cout<<"	Ny = "<<b_input.planar->at(j).ny()<<" Ny interpolated = "<<b_input.planar->at(j).ny_interp()<<endl;
		cout<<"	Nz = "<<b_input.planar->at(j).nz()<<" Nz interpolated = "<<b_input.planar->at(j).nz_interp()<<endl;
	}
	for (int j = 0; j < b_input.tangent->size(); j++ ){
		eval_vector_interpolant_at_point(b_input.tangent->at(j));
		double vf[3] = {b_input.tangent->at(j).nx_interp(),b_input.tangent->at(j).ny_interp(),b_input.tangent->at(j).nz_interp()};
		cout<<" Tangent["<<j<<"]: "<<endl;
		cout<<"	Tx = "<<b_input.tangent->at(j).tx()<<" Nx interpolated = "<<b_input.tangent->at(j).nx_interp()<<endl;
		cout<<"	Ty = "<<b_input.tangent->at(j).ty()<<" Ny interpolated = "<<b_input.tangent->at(j).ny_interp()<<endl;
		cout<<"	Tz = "<<b_input.tangent->at(j).tz()<<" Nz interpolated = "<<b_input.tangent->at(j).nz_interp()<<endl;
		cout<<" Tx*nx + Ty*ny + Tz*nz = "<<b_input.tangent->at(j).tx()*b_input.tangent->at(j).nx_interp() + b_input.tangent->at(j).ty()*b_input.tangent->at(j).ny_interp() +
			b_input.tangent->at(j).tz()*b_input.tangent->at(j).nz_interp()<<endl;
	}

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
			if (b_parameters.modified_basis)
			{
				if (!convert_modified_kernel_to_rbf_kernel())
				{
					error_msg.append(" QPP solution conversion to Linear Failure.");
					return false;
				}
			}

			int N = (int)b_input.evaluation_pts->size();
			int add = 0;
			int vv = round((double)N / 72.0); // 72.0 is the width of the progress bar 
			double factor = (100.0*(double)vv) / (double)N;
			
			#pragma omp parallel for schedule(dynamic)
			for (int j = 0; j < N; j++ ){
				eval_scalar_interpolant_at_point(b_input.evaluation_pts->at(j));
				//eval_vector_interpolant_at_point(b_input.evaluation_pts->at(j));
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

	// set OpenMP parameters
	const int nthreads = omp_get_max_threads();
	omp_set_dynamic(false);
	if (nthreads >= 8) omp_set_num_threads(nthreads - 2);
	else omp_set_num_threads(nthreads - 1);
	///////////////////////////////////////

	cout<<" Starting SURFE algorithm "<<endl;
	cout<<" Processing input data...";
	if (!process_input_data())
	{
		error_msg = "Error processing input data";
		return false;
	}
	cout<<"done!"<<endl;
	cout<<" Get method parameters...";
	if (!get_method_parameters())
	{
		error_msg.append(" Error getting method parameters.");
		return false;
	}
	cout<<"done!"<<endl;
	cout<<" Setup basis functions...";
	if (!setup_basis_functions())
	{
		error_msg.append(" Error setting up basis functions.");
		return false;
	}
	cout<<"done!"<<endl;
	cout<<" Solve mathematical problem...";
	if (!setup_system_solver())
	{
		error_msg.append(" Error solving mathematical equations.");
		return false;
	}
	cout<<"done!"<<endl;
	cout<<" Evaluate scalar interpolant at grid nodes...";
	if (!evaluate_scalar_interpolant())
	{
		error_msg.append(" Error evaluating interpolant in grid of points.");
		return false;
	}
	cout<<"done!"<<endl;
	cout<<" Total computation time = "<<((double)clock()-tstart)/CLOCKS_PER_SEC<<endl;

	return true;
}

bool GRBF_Modelling_Methods::get_equality_matrix( const MatrixXd &interpolation_matrix, MatrixXd &equality_matrix )
{
	if (equality_matrix.rows() == 0 ||
		equality_matrix.rows() > interpolation_matrix.rows() ||
		equality_matrix.cols() != interpolation_matrix.cols()) return false;
	int n_ie = (int)interpolation_matrix.rows() - (int)equality_matrix.rows();
	if (n_ie != b_parameters.n_inequality) return false;

	for (int j = 0; j < equality_matrix.rows(); j++ ){
		for (int k = 0; k < equality_matrix.cols(); k++ ){
			equality_matrix(j,k) = interpolation_matrix(j + n_ie,k);
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

bool GRBF_Modelling_Methods::_output_greedy_debug_objects()
{
	if (solver == NULL) return false;
	else
	{
		if ((int)solver->weights.size() == 0) return false;
		else
		{
			if ((int)b_input.evaluation_pts->size() != 0)	evaluate_scalar_interpolant();
			else return false;
		}
	}

	return true;
}

GRBF_Modelling_Methods * GRBF_Modelling_Methods::get_method( const model_parameters& m_parameters,const Basic_input& input )
{
	if (m_parameters.model_type == Parameter_Types::Single_surface) return new Single_Surface(m_parameters,input);
	else if (m_parameters.model_type == Parameter_Types::Lajaunie_approach) return new Lajaunie_Approach(m_parameters,input);
	else if (m_parameters.model_type == Parameter_Types::Stratigraphic_horizons) return new Stratigraphic_Surfaces(m_parameters,input);
	else return new Continuous_Property(m_parameters,input);
}

bool GRBF_Modelling_Methods::run_greedy_algorithm()
{
	// check if there are non-zero errors permitted on the data
	if (m_parameters.interface_uncertainty == 0 && m_parameters.angular_uncertainty == 0) return false;

	GRBF_Modelling_Methods *greedy_method = get_method(m_parameters,b_input);

	greedy_method->b_input.compute_avg_nn_distances();

	Basic_input greedy_input, excluded_input;
	// initialize starting data
	if (!get_minimial_and_excluded_input(greedy_input, excluded_input)) return false;
	greedy_input.SetInequalityAvgNNDist(greedy_method->b_input.GetInequalityAvgNNDist());
	greedy_input.SetInterfaceAvgNNDist(greedy_method->b_input.GetInterfaceAvgNNDist());
	greedy_input.SetPlanarAvgNNDist(greedy_method->b_input.GetPlanarAvgNNDist());
	greedy_input.SetTangentAvgNNDist(greedy_method->b_input.GetTangentAvgNNDist());
	greedy_method->b_input = greedy_input;

	bool converged = false;
	int iter = 0;
	while (!converged)
	{
		// run normal algorithm
		if (!greedy_method->process_input_data()) return false;
		if (!greedy_method->get_method_parameters()) return false;
		if (!greedy_method->setup_basis_functions()) return false;
		if (!greedy_method->setup_system_solver()) return false;

		// measure residuals
		if (!greedy_method->measure_residuals(excluded_input)) return false;

		// debug: should output intermediate input constraints and modelled surface using those constraints
		// if ( !greedy_method->_output_greedy_debug_objects()) return false;

		// add appropriate data based on residuals
		if (!greedy_method->append_greedy_input(excluded_input)) converged = true; // if no input is added convergence is assumed
		iter++;
		greedy_method->_SetIteration(iter);
	}

	greedy_method->b_input.evaluation_pts = b_input.evaluation_pts;
	if (!greedy_method->evaluate_scalar_interpolant()) return false;

	b_input.evaluation_pts = greedy_method->b_input.evaluation_pts;
	b_input.interface_iso_values->clear();
	b_input.interface_iso_values = new std::vector<double>(*greedy_method->b_input.interface_iso_values);

	return true;
}

unsigned int Tensor_Methods::_get_nearest_orientation_pt_index(const TensorEvaluationPoints &pt)
{
	double nn_dist = 0;
	unsigned int index = 0;
	input.orientation->at(0).x();
	for (unsigned int j = 0; j < (unsigned int)input.orientation->size(); j++){
		double dx = pt.x() - input.orientation->at(j).x();
		double dy = pt.y() - input.orientation->at(j).y();
		double dz = pt.z() - input.orientation->at(j).z();
		double d = sqrt(dx*dx + dy*dy + dz*dz);
		if (j == 0)
		{
			nn_dist = d;
			index = j;
		}
		else
		{
			if (d < nn_dist)
			{
				nn_dist = d;
				index = j;
			}
		}
	}
	return index;
}

void Tensor_Methods::_sort_eigensystem(Matrix3d &evectors, Vector3d &evalues)
{
	// sort the eigenvectors  & eigenvalues according to ascending eigenvalue order
	// create temp storage ...
	double EigenValues[3];
	double EigenVectors[3][3];
	std::vector<double> eigenvalues;
	for (int j = 0; j < 3; j++){
		EigenValues[j] = evalues[j];
		eigenvalues.push_back(EigenValues[j]);
		for (int k = 0; k < 3; k++){
			EigenVectors[j][k] = evectors(j, k);
		}
	}

	sort(eigenvalues.begin(), eigenvalues.end());

	int Ind[3] = { 0 };
	if (eigenvalues.at(0) == EigenValues[0])Ind[2] = 0;
	if (eigenvalues.at(0) == EigenValues[1])Ind[2] = 1;
	if (eigenvalues.at(0) == EigenValues[2])Ind[2] = 2;
	if (eigenvalues.at(1) == EigenValues[0])Ind[1] = 0;
	if (eigenvalues.at(1) == EigenValues[1])Ind[1] = 1;
	if (eigenvalues.at(1) == EigenValues[2])Ind[1] = 2;
	if (eigenvalues.at(2) == EigenValues[0])Ind[0] = 0;
	if (eigenvalues.at(2) == EigenValues[1])Ind[0] = 1;
	if (eigenvalues.at(2) == EigenValues[2])Ind[0] = 2;

	double tEigenVectors[3][3] = { 0 };
	double tEigenValues[3];
	for (int l = 0; l < 3; l++){
		for (int k = 0; k < 3; k++){
			tEigenVectors[l][k] = EigenVectors[l][Ind[k]];
		}
		tEigenValues[l] = EigenValues[Ind[l]];
	}
	for (int l = 0; l < 3; l++){
		for (int k = 0; k < 3; k++){
			EigenVectors[l][k] = tEigenVectors[l][k];
		}
		EigenValues[l] = abs(tEigenValues[l]);
	}
	evectors(0, 0) = tEigenVectors[0][2];
	evectors(1, 0) = tEigenVectors[1][2];
	evectors(2, 0) = tEigenVectors[2][2];
	evectors(0, 1) = tEigenVectors[0][1];
	evectors(1, 1) = tEigenVectors[1][1];
	evectors(2, 1) = tEigenVectors[2][1];
	evectors(0, 2) = tEigenVectors[0][0];
	evectors(1, 2) = tEigenVectors[1][0];
	evectors(2, 2) = tEigenVectors[2][0];

	evalues[0] = EigenValues[2];
	evalues[1] = EigenValues[1];
	evalues[2] = EigenValues[0];
}

double Tensor_Methods::_GetLargestInterPointDistance(const std::vector<Orientation> &pts)
{
	double lipd = 1.0;
	std::vector<double> LD;
	for (int j = 0; j < pts.size(); j++){
		std::vector<double> d;
		for (int k = 0; k < pts.size(); k++){
			double dx = pts.at(j).x() - pts.at(k).x();
			double dy = pts.at(j).y() - pts.at(k).y();
			double dz = pts.at(j).z() - pts.at(k).z();
			double r2 = dx*dx + dy*dy + dz*dz;
			d.push_back(r2);
		}
		std::sort(d.begin(), d.end());
		LD.push_back(d[d.size() - 1]);
	}
	std::sort(LD.begin(), LD.end());

	return sqrt(LD[LD.size() - 1]);
}

bool Tensor_Methods::run_algorithm()
{
	// get the local anisotropy @ every point (input.orientation[])
	//if (!input.GetLocalAnisotropy(parameters)) return false;


// 	for (int j = 0; j < (int)input.orientation->size(); j++){
// 		std::cout << " Orientation[" << j << "]" << endl;
// 		std::cout << " Eigenvectors:\n" << input.orientation->at(j).eigenvectors << std::endl;
// 		std::cout << " Eigenvalues:\n" << input.orientation->at(j).eigenvalues << std::endl;
// 	}
	distance = _GetLargestInterPointDistance(*input.orientation);

	#pragma omp parallel for 
	for (int j = 0; j < (int)input.evalpts->size(); j++) interpolate_tensor_field_at_pt(input.evalpts->at(j));
	return true;
}

bool Tensor_Methods::interpolate_tensor_field_at_pt(TensorEvaluationPoints& eval_pt)
{
	// get closest input pt
	unsigned int index = _get_nearest_orientation_pt_index(eval_pt);
	//cout<<" Nearest orientation point= ["<<index<<"]"<<endl;

	// initialize this object 
	Orientation cur_tensor_est = input.orientation->at(index);
	// compute sigma_t^(-1/2) & sigma_t^(1/2)
	// for first iteration t, sigma_t will be defined from the closest input pt ( Centers[index] )
	Matrix3d sigma_t_inv_sqrt;
	Matrix3d sigma_t_sqrt;

	// spectral decomposition of tensor @ Centers[index] -> U * Diag(d_i) * U^T
	// U = | ev1_x ev2_x ev3_x |  Diag(d_i) = | ev1  0   0  |
	//     | ev1_y ev2_y ev3_y |			  | 0   ev2  0  |
	//     | ev1_z ev2_z ev3_z |              | 0    0  ev3 | 

	// sigma_t^(-1/2) = U * Diag{ (d_i)^(-1/2) } * U^T
	// sigma t^(1/2) = U * Diag { (d_i)^(1/2) } * U^T

	// declare temporary containers
	Matrix3d  sigma_t_inv_sqrt_temp;
	Matrix3d  sigma_t_sqrt_temp;
	Matrix3d  U;
	Matrix3d  UT;
	Vector3d  Diag_inv_sqrt;
	Vector3d  Diag_sqrt;
	Matrix3d  Uj;
	Matrix3d  UTj;
	Vector3d  Diagj;
	Matrix3d  temp_mat;
	Matrix3d  temp_mat2;

	bool converged = false;
	int niter = 0;
	double last_avg_err = 0;
	double last_avg = 0;
	double last_fin_mat_add = 0;
	while (!converged)
	{
		//cout<<" niter = "<<niter<<endl;
		niter++;
		U = cur_tensor_est.eigenvectors;
		UT = U.transpose();
		for (int j = 0; j < 3; j++){
			Diag_sqrt(j) = sqrt(cur_tensor_est.eigenvalues(j));
			Diag_inv_sqrt(j) = 1.0 / sqrt(cur_tensor_est.eigenvalues(j));
		}
// 		cout << " U:\n" << U << endl;
// 		cout << " UT:\n" << UT << endl;
// 		cout << " Diag_sqrt:\n" << Diag_sqrt << endl;
// 		cout << " Diag_inv_sqrt:\n" << Diag_inv_sqrt << endl;

		// compute the sigma matrices
		// sigma_t^(1/2)
		//cout<<" Multiplying Eigenvector Matrix U by Diag(sqrt(ev[1,2,3])) * UT of current estimated tensor "<<endl;
		sigma_t_sqrt = U * Diag_sqrt.asDiagonal() * UT;
		//cout << " sigma_t_sqrt:\n" << sigma_t_sqrt << endl;

		// sigma_t^(-1/2)
		//cout<<" Multiplying Eigenvector Matrix U by Diag(1/sqrt(ev[1,2,3])) * UT of current estimated tensor "<<endl;
		sigma_t_inv_sqrt = U * Diag_inv_sqrt.asDiagonal() * UT;

		//cout << " sigma_t_inv_sqrt:\n" << sigma_t_inv_sqrt << endl;

		Matrix3d cur_matrix;
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				cur_matrix(j,k) = 0.0;
			}
		}
		double total_weight = 0;

		for (int j = 0; j < (int)input.orientation->size(); j++){
			// cout<<"	WORKING ON Centers["<<j<<"]"<<endl;
			// compute sigma_t^(-1/2) * sigma_j * sigma_t^(-1/2) 
			// sigma[j]  = U[j] * Diag (eigenvalues[1,2,3]) * UT[j]
			// get U[j] & UT[j]
			//cout<<" spectral decomposition @ ["<<j<<"]"<<endl;
			Uj = input.orientation->at(j).eigenvectors;
			UTj = Uj.transpose();
			Diagj = input.orientation->at(j).eigenvalues;

// 			cout << " Uj:\n" << Uj << endl;
// 			cout << " UTj:\n" << UTj << endl;
// 			cout << " Diagj:\n" << Diagj << endl;

			//cout<<" Getting sigma^(-1/2)*Uj*Diagj*UTj by sigma^(-1/2) "<<endl;
			temp_mat2 = sigma_t_inv_sqrt*Uj*Diagj.asDiagonal()*UTj*sigma_t_inv_sqrt;
			//cout << " temp_mat2:\n" << temp_mat2 << endl;

			// complete a spectral decomposition on the inversion product
			SelfAdjointEigenSolver<Matrix3d> ea;
			ea.compute(temp_mat2);
			Matrix3d U_j;
			Vector3d d_j;
			U_j = ea.eigenvectors().real();
			d_j = ea.eigenvalues().real();
			Matrix3d U_j_T;
			Vector3d D_j;
			Vector3d D_j_log;
			U_j_T = U_j.transpose();
			D_j = d_j;
			//cout << " d_j:\n" <<d_j << endl;
			for (int k = 0; k < 3; k++) D_j_log[k] = log(d_j[k]);

// 			cout << " U_j:\n" << U_j << endl;
// 			cout << " U_j_T:\n" << U_j_T << endl;
// 			cout << " D_j_log:\n" << D_j_log << endl;

			//cout<<" Eigenvalues = {"<<d_j[0]<<","<<d_j[1]<<","<<d_j[2]<<"}"<<endl;
			//cout<<" Multiplying U_j*D_j_log * U_j_T"<<endl;
			temp_mat2 = U_j*D_j_log.asDiagonal()*U_j_T;

			//std::cout << " temp_mat2:\n" << temp_mat2 << std::endl;

			double dx = input.orientation->at(j).x() - eval_pt.x();
			double dy = input.orientation->at(j).y() - eval_pt.y();
			double dz = input.orientation->at(j).z() - eval_pt.z();
			double r = sqrt(dx*dx + dy*dy + dz*dz);
			//cout<<" r = "<<r<<endl;
			double j_weight = 1 / (pow(r,parameters.idw_power));//exp(-1.0*value);//1 - 3.0 * sr*sr + 2.0 * sr*sr*sr;//1.0 / s;//exp(-1.0*s);// 1 - 3 * sr*sr + 2 * sr*sr*sr; //1 / (sr*sr);
			//cout<<" basis weight = "<<basis<<endl;
			//cout<<" multiplying U_j*D_j_log*U_j_T by basis weight"<<endl;
			temp_mat = temp_mat2 * j_weight;
			//std::cout << " temp_mat:\n" << temp_mat << std::endl;
			//cout<<" Adding current matrix to previous iteration "<<endl;
			cur_matrix += temp_mat;
			total_weight += j_weight;
		}
		//std::cout << " cur_matrix:\n" << cur_matrix<< std::endl;
		//cout<<" total weight = "<<total_weight<<endl;
		Matrix3d fin_mat;
		//cout<<" Muliplying cur_matrix by 1/total_weight"<<endl;
		fin_mat = cur_matrix*(1.0 / total_weight);
		double fin_mat_add = 0;
		//cout<<" fin_matrix = "<<endl;
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				//cout<< setw(15) <<fin_mat[j][k];
				fin_mat_add += abs(fin_mat(j,k));
			}
			//cout<<endl;
		}
		//cout<<"	FINAL MATRIX ADDED = "<<fin_mat_add<<endl;
		// complete last spectral decomposition on fin_mat
		//cout<<" completing last spectral decomposition !!!!!!!!!!!!!"<<endl;
		SelfAdjointEigenSolver<Matrix3d> eaf;
		eaf.compute(fin_mat);
		Matrix3d U_f;
		Vector3d d_f;
		U_f = eaf.eigenvectors().real();
		d_f = eaf.eigenvalues().real();

		//Math_methods::sort_eigensystem(U_f,d_f);
		Matrix3d U_f_T;
		Vector3d D_f;
		Vector3d D_f_exp;
		U_f_T = U_f.transpose();
		D_f = d_f;
		for (int j = 0; j < 3; j++) D_f_exp(j) = exp(d_f(j));
		//cout<<" eigenvalues= {"<<d_f[0]<<","<<d_f[1]<<","<<d_f[2]<<"}"<<endl;

		//cout<<" Multiplying U_f*D_f_exp by U_f_T"<<endl;
		temp_mat = U_f * D_f_exp.asDiagonal() * U_f_T;
		//cout<<" Multiplying sigma^(1/2)*U_f*D_f_exp*U_f_T by sigma^(1/2)"<<endl;
		temp_mat2 = sigma_t_sqrt * temp_mat * sigma_t_sqrt;

		//cout<<" current estimated tensor "<<endl;
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				//cout<< setw(15) <<temp_mat2[j][k];
			}
			//cout<<endl;
		}


		//cout<<" completing FINAL spectral decomposition !!!!!!!!!!! "<<endl;
		SelfAdjointEigenSolver<Matrix3d> eafinal;
		eafinal.compute(temp_mat2);
		Matrix3d u_tp1;
		Vector3d d_tp1;
		u_tp1 = eafinal.eigenvectors().real();
		d_tp1 = eafinal.eigenvalues().real();

		// look at diff between current estimation of the tensor field and the last estimation 
		// determined if convergence have been achieved

		// check if convergence has been achieved
		// take an average of the all the elements in the difference matrix
		if (fin_mat_add < 1e-6 || niter > 500) converged = true;
		if (niter > 1)
		{
			if (fin_mat_add > last_fin_mat_add) converged = true;
		}
		if (!converged)
		{
			cur_tensor_est.eigenvectors = u_tp1;
			cur_tensor_est.eigenvalues = d_tp1;
			for (int j = 0; j < 3; j++){
				for (int k = 0; k < 3; k++){
					//cout<< setw(15) << u_tp1[j][k];
				}
				//cout<<endl;
			}
			//cout<<" eigenvalues = {"<<d_tp1[0]<<","<<d_tp1[1]<<","<<d_tp1[2]<<"}"<<endl;
		}
		last_fin_mat_add = fin_mat_add;
	}


	Matrix3d A = cur_tensor_est.eigenvectors;
	Vector3d b = cur_tensor_est.eigenvalues;
	_sort_eigensystem(A, b);
	cur_tensor_est.eigenvectors = A;
	cur_tensor_est.eigenvalues = b;

	// set interpolation properties to node 
	eval_pt.eigenvalues = cur_tensor_est.eigenvalues;
	eval_pt.eigenvectors = cur_tensor_est.eigenvectors;

	//cout<<" eval_pt -> eigenvector matrix: "<<endl;
	for (int j = 0; j < 3; j++){
		for (int k = 0; k < 3; k++){
			//cout<< setw(15) << eval_pt.eigenvectors(j,k);
		}
		//cout<<endl;
	}

	//cout<<" eval_pt -> eigenvalues {"<<eval_pt.Eigenvalue[0]<<","<<eval_pt.Eigenvalue[1]<<","<<eval_pt.Eigenvalue[2]<<"}"<<endl;

	//cout<<" convergence achieved@ iter["<<niter<<"]"<<endl;


	return true;
}
