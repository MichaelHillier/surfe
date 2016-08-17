#ifndef math_methods_h
#define math_methods_h

#include <math_lib_module.h>
#include <gmpxx.h>

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

using namespace Eigen;

class MATH_LIB_EXPORT Math_methods{
private:
	template <class T> static T _find_step_length(const Matrix <T, Dynamic, 1> &a, const Matrix <T, Dynamic, 1> &da,
		const Matrix <T, Dynamic, 1> &b, const Matrix <T, Dynamic, 1> &db);
	template <class T> static void _rot(std::vector< std::vector < T > > &a, const T &s, const T &tau, const int &i, const int &j, const int &k, const int &l);
	static double _get_double(const double &d) { return d;}
	static double _get_double(const mpf_class &d) { return d.get_d(); }
	static double _find_step(const VectorXd &da, const VectorXd &a);
	static double _find_positivity_step( const VectorXd &da, const VectorXd &a,
		                          const VectorXd &db, const VectorXd &b,
								  const VectorXd &dc, const VectorXd &c,
								  const VectorXd &dd, const VectorXd &d);
public:
	template <class T> static bool sort_vector_w_index(std::vector<T> &arr, std::vector<int> &brr);
	template <class T> static T max_element_wrt_zero(const T &a, const T &b);
	template <class T> static void SWAP(T &a, T &b){ T x = a; a = b; b = x; }
	template <class T> static bool angle_btw_2_vectors( const std::vector < T > &v1, const std::vector < T > &v2, T &angle);
	static double RandomDouble(const double &min, const double&max);
	template <class T> static bool quadratic_solver(const Matrix <T, Dynamic, Dynamic> &H,
		const Matrix <T, Dynamic, Dynamic> &A,
		const Matrix <T, Dynamic, Dynamic> &C,
		const Matrix <T, Dynamic, 1> &b,
		const Matrix <T, Dynamic, 1> &d,
		Matrix <T, Dynamic, 1> &fvalues);
	static bool Math_methods::quadratic_solver_loqo( const MatrixXd &H, const MatrixXd &A, const VectorXd &b, const VectorXd &r, VectorXd &fvalues );
// 	static bool quadratic_solver_loqo(const Matrix <double, Dynamic, Dynamic> &H,
// 		const Matrix <double, Dynamic, Dynamic> &A,
// 		const Matrix <double, Dynamic, 1> &b,
// 		const Matrix <double, Dynamic, 1> &r,
// 		Matrix <double, Dynamic, 1> &fvalues);
};

template <class T>
bool Math_methods::angle_btw_2_vectors( const std::vector < T > &v1, const std::vector < T > &v2, T &angle )
{
	if (v1.size() != v2.size() ) return false;

	T dotp = 0.0;
	T v1norm = 0.0;
	T v2norm = 0.0;
	for (int j = 0; j < (int)v1.size(); j++ ){
		dotp += v1[j] * v2[j];
		v1norm += v1[j] * v1[j];
		v2norm += v2[j] * v2[j];
	}

	angle = acos(dotp / (sqrt(v1norm)*sqrt(v2norm)));
	return true;
}

template <class T>
T Math_methods::max_element_wrt_zero( const T &a, const T &b)
{
	T max_value = a;
	T c = 0.0;
	if (b > max_value) max_value = b;
	if (c > max_value) max_value = c;
	
	return max_value;
}

template <class T>
T Math_methods::_find_step_length(const Matrix <T, Dynamic, 1> &a, const Matrix <T, Dynamic, 1> &da,
	const Matrix <T, Dynamic, 1> &b, const Matrix <T, Dynamic, 1> &db)
{
	int n = (int)a.rows();
	//if (a.rows() != b.rows() || a.rows() != da.rows() || a.rows() != db.rows() || da.rows() != db.rows()) throw -1;
	// find step length (a) ... this is the most important step
	// has to satisfy all the non-negativity conditions ...
	T min_alpha_a = 100.0; // improper initialization. proper initialization will occur in first iteration of below for loop 
	T min_alpha_b = 100.0; // ?????		why do i do this? seems dumb
	T max_alpha_a = 0.0; // improper initialization. proper initialization will occur in first iteration of below for loop 
	T max_alpha_b = 0.0;
	T alpha = 0.0;
	Matrix <T, Dynamic, 1> alpha_a(n);
	Matrix <T, Dynamic, 1> alpha_b(n);
	for (int j = 0; j < n; j++){
		cout<<" z["<<j<<"]= "<<b[j].get_d()<<endl;
		cout<<" s["<<j<<"]= "<<a[j].get_d()<<endl;
		cout<<" ["<<j<<"] dz_aff = "<<db[j].get_d()<<" ds_aff= "<<da[j].get_d()<<endl;
		// do alpha_b first ...
		if (b(j) > 0.0) alpha_b(j) =  b(j) / db(j) ;
		else alpha_b(j) = -1.0*b(j) / db(j);
		// alpha_a
		if (a(j) > 0.0) alpha_a(j) =  a(j) / da(j);
		else alpha_a(j) = -1.0*a(j) / da(j);

		if (alpha_b(j) < min_alpha_b && alpha_b(j) > 0.00000000000001) min_alpha_b = alpha_b(j);
		if (alpha_a(j) < min_alpha_a && alpha_a(j) > 0.00000000000001) min_alpha_a = alpha_a(j);
		cout<<" alpha_s["<<j<<"]= "<<alpha_a[j].get_d()<<endl;
		cout<<" alpha_z["<<j<<"]= "<<alpha_b[j].get_d()<<endl;
		//cout<<" max_alpha_s = "<<max_alpha_a.get_d()<<endl;
		//cout<<" max_alpha_z = "<<max_alpha_b.get_d()<<endl;
		cout<<" min_alpha_s = "<<min_alpha_a.get_d()<<endl;
		cout<<" min_alpha_z = "<<min_alpha_b.get_d()<<endl;
	}

	if (min_alpha_b < min_alpha_a) alpha = min_alpha_b;
	else alpha = min_alpha_a;

	if (alpha > 1.0 || alpha == 0) alpha = 1.0;
	return alpha;
}

template <class T>
bool Math_methods::quadratic_solver(const Matrix <T, Dynamic, Dynamic> &H,
	const Matrix <T, Dynamic, Dynamic> &A,
	const Matrix <T, Dynamic, Dynamic> &C,
	const Matrix <T, Dynamic, 1> &b,
	const Matrix <T, Dynamic, 1> &d,
	Matrix <T, Dynamic, 1> &fvalues)
{
	// Describe the Quadratic Optimization problem
	// min (w.r.t. x) f(x) = 1/2 * xT * H * x => our objective function
	// s.t.            T * x = b
	//                C * x >= d
	// NOTE: H must be a positive definite matrix !!! There is a check for this in the constructor for the matrix solver
	// The Interior Point method using a Predictor-Corrector method will be used to solve the QP
	// Lagrangian (Primal)
	// L(x,y,z) = 1/2 * xT * H * x - yT * (T * x - b) - zT * (C * x - d)
	// KKT Conditions
	// dL/dx = H * x - AT * y - CT * z = 0
	// dL/dy = T * x - b               = 0
	// dL/dz = C * x - s - d           = 0  s (slack variable)
	//         z >= 0 , s >= 0
	//         z_i * s_i = 0
	// Solve Block Matrix
	// | H  -AT  -CT   0 || dx |     | rh |  -> r1      
	// | T    0    0   0 || dy | = - | ra |  -> r2   EQN (1)
	// | C    0    0  -I || dz |     | rc |  -> r3 
	// | 0    0    S   Z || ds |     | rsz|  -> r4 
	// where S = diag (s0,s1,...,sn)
	//       Z = diag (z0,z1,...,zn)
	//       I  = Identity matrix
	//       rh = H * x - AT * y - CT * z
	//       ra = T * x - b
	//       rc = C * x - s - d
	//       rsz = Z * S
	//       n := # of rows in H
	//       na := # of rows in T
	//       nc := # of rows in C
	//       H[n x n]
	//       T[na x n]
	//       C[nc x n]
	//       S[nc x nc]
	//       Z[nc x nc]
	//       I[nc x nc]
	// EQN (1) can be reduced in size by block elimination since S and Z are strictly positive
	// New EQN becomes
	// | H  AT      CT   || dx|    |     -rh     | 
	// | T   0	     0   ||-dy| =  |     -ra     |
	// | C   0 -Z^(-1)*S ||-dz|    |-rc-Z^(-1)rsz|

	int n = (int)H.rows();
	int na = (int)A.rows();
	int nc = (int)C.rows();

	///////////////////////////////////////////////////////////////////////
	/////////////////////////// Initialization ////////////////////////////
	///////////////////////////////////////////////////////////////////////
	// Calculation Helper Variables
	T add, temp, temp2, temp_p, minus_one, elemsum, datanorm, zTs, mu, zero, prev_mu;
	zero = 0.0;
	elemsum = 0.0;
	// KKT System
	Matrix <T, Dynamic, 1> x(n);
	x.setZero();
	Matrix <T, Dynamic, 1> y(na);
	y.setZero();
	Matrix <T, Dynamic, 1> z(nc);
	z.setZero();
	Matrix <T, Dynamic, 1> s(nc);
	s.setZero();
	Matrix <T, Dynamic, 1> rh(n);
	rh.setZero();
	Matrix <T, Dynamic, 1> ra(na);
	ra.setZero();
	Matrix <T, Dynamic, 1> rc(nc);
	rc.setZero();
	Matrix <T, Dynamic, 1> rsz(nc);
	rsz.setZero();

	Matrix <T, Dynamic, 1> t1(n);
	t1.setZero();
	Matrix <T, Dynamic, 1> t2(n);
	t2.setZero();
	Matrix <T, Dynamic, 1> t3(n);
	t3.setZero();
	Matrix <T, Dynamic, 1> t4(na);
	t4.setZero();
	Matrix <T, Dynamic, 1> t5(nc);
	t5.setZero();

	Matrix <T, Dynamic, Dynamic> AT = A.transpose(); 
	Matrix <T, Dynamic, Dynamic> CT = C.transpose();
	Matrix <T, Dynamic, Dynamic> KKT(n + na + 2 * nc, n + na + 2 * nc);
	Matrix <T, Dynamic, Dynamic> KKT_predictor(n + na + 2 * nc, n + na + 2 * nc);
	Matrix <T, Dynamic, Dynamic> KKT_corrector(n + na + 2 * nc, n + na + 2 * nc);

	Matrix <T, Dynamic, 1> solution_vector(n + na + 2 * nc);
	solution_vector.setZero();

	// step containers and related helper variables
	T alpha, min_alpha_z, min_alpha_s, max_alpha_z, max_alpha_s, mu_aff, sigma;
	Matrix <T, Dynamic, 1> alpha_z(nc);
	alpha_z.setZero();
	Matrix <T, Dynamic, 1> alpha_s(nc);
	alpha_s.setZero();
	Matrix <T, Dynamic, 1> dx_aff(n);
	dx_aff.setZero();
	Matrix <T, Dynamic, 1> dy_aff(na);
	dy_aff.setZero();
	Matrix <T, Dynamic, 1> dz_aff(nc);
	dz_aff.setZero();
	Matrix <T, Dynamic, 1> ds_aff(nc);
	ds_aff.setZero();
	Matrix <T, Dynamic, 1> rsz_aff(nc);
	rsz_aff.setZero();
	Matrix <T, Dynamic, 1> dvector(n + na + 2 * nc);
	dvector.setZero();

	Matrix <T, Dynamic, 1> dx(n);
	dx.setZero();
	Matrix <T, Dynamic, 1> dy(na);
	dy.setZero();
	Matrix <T, Dynamic, 1> dz(nc);
	dz.setZero();
	Matrix <T, Dynamic, 1> ds(nc);
	ds.setZero();
	Matrix <T, Dynamic, 1> dvector_corr(n + na + 2 * nc);
	dvector_corr.setZero();

	///////////////////////////////////////////////////////////////////////
	////////////////////// End of Initialization //////////////////////////
	///////////////////////////////////////////////////////////////////////

	// get norm of matrix 
	datanorm = sqrt(H.maxCoeff());

	// setup initial iterate for z and s...
	for (int j = 0; j < nc; j++){
		z(j) = datanorm;
		s(j) = datanorm;
	}

	// Build KKT matrix blocks [1][1], [1][2], [1][3], [1][4]
	for (int j = 0; j < n; j++){ // Matrix rows...
		// Matrix columns...
		// [1][1] Block
		for (int k = 0; k < n; k++) KKT(j,k) = H(j,k);
		// [1][2] Block
		for (int k = 0; k < na; k++) KKT(j,n + k) = -1.0*AT(j,k);
		// [1][3] Block
		for (int k = 0; k < nc; k++) KKT(j,n + na + k) = -1.0*CT(j,k);
		// [1][4] Block
		for (int k = 0; k < nc; k++) KKT(j,n + na + nc + k) = 0.0;
	}
	// Build KKT matrix blocks [2][1], [2][2], [2][3], [2][4]
	for (int j = 0; j < na; j++){
		// [2][1] Block
		for (int k = 0; k < n; k++)  KKT(n + j,k) = A(j,k);
		// [2][2] Block
		for (int k = 0; k < na; k++) KKT(n + j,n + k) = 0.0;
		// [2][3] Block
		for (int k = 0; k < nc; k++) KKT(n + j,n + na + k) = 0.0;
		// [2][4] Block
		for (int k = 0; k < nc; k++) KKT(n + j,n + na + nc + k) = 0.0;
	}
	// Build KKT matrix blocks [3][1], [3][2], [3][3], [3][4]
	for (int j = 0; j < nc; j++){
		// [3][1] Block
		for (int k = 0; k < n; k++)  KKT(n + na + j, k) = C(j, k);
		// [3][2] Block
		for (int k = 0; k < na; k++) KKT(n + na + j,n + k) = 0.0;
		// [3][3] Block
		for (int k = 0; k < nc; k++) KKT(n + na + j,n + na + k) = 0.0;
		// [3][4] Block
		for (int k = 0; k < nc; k++){
			if (j == k) KKT(n + na + j,n + na + nc + k) = -1.0;
			else KKT(n + na + j,n + na + nc + k) = 0.0;
		}
	}
	// Build KKT matrix blocks [4][1], [4][2]
	for (int j = 0; j < nc; j++){
		// [4][1] Block
		for (int k = 0; k < n; k++)  KKT(n + na + nc + j,k) = 0.0;
		// [4][2] Block
		for (int k = 0; k < na; k++) KKT(n + na + nc + j,n + k) = 0.0;
	}

	// Build KKT matrix blocks: [4][3], [4][4]
	for (int j = 0; j < nc; j++){
		// [4][3] Block
		for (int k = 0; k < nc; k++){
			if (j == k) KKT(n + na + nc + j,n + na + k) = s(j);
			else KKT(n + na + nc + j,n + na + k) = 0.0;
		}
		// [4][4] Block
		for (int k = 0; k < nc; k++){
			if (j == k) KKT(n + na + nc + j,n + na + nc + k) = z(j);
			else KKT(n + na + nc + j,n + na + nc + k) = 0.0;
		}
	}

// 	std::ofstream out("E:\interpolation_mat.txt");
// 	out.precision(15);
// 	for (int j = 0; j < n + na + 2*nc; j++ ){
// 		for (int k = 0; k < n + na + 2*nc; k++ ){
// 			out<<"  "<<KKT(j,k).get_d();
// 		}
// 		out<<endl;
// 	}

	//////////////////////////////////////////////
	// Determine Starting Point for the iterate //
	//////////////////////////////////////////////
	// 1st compute residuals for the affine system
	t1 = H*x;
	t2 = AT*y;
	t3 = CT*z;
	t4 = A*x;
	t5 = C*x;
	// 	_output_vector(t3,"t3.txt"); // debug
	// 	_output_vector(z,"z.txt");
	// 	_output_matrix(CT,"CT.txt");
	for (int j = 0; j < n; j++){ // rh residual = H*x - AT*y - CT*z
		rh(j) = t1(j) - t2(j) - t3(j);
		solution_vector(j) = -1.0*rh(j);
	}
	for (int j = 0; j < na; j++){ // ra residual = A*x - b
		ra(j) = t4(j) - b(j);
		solution_vector(n + j) = -1.0*ra(j);
	}
	zTs = 0.0; // re-initialize this adder
	for (int j = 0; j < nc; j++){ // rc residual = C*x - s - d
		rc(j) = t5(j) - s(j) - d(j);
		rsz(j) = s(j) * z(j); // rsz residual = Z[]*S[]*e
		solution_vector(n + na + j) = -1.0*rc(j);
		zTs += rsz(j);
		solution_vector(n + na + nc + j) = -1.0*rsz(j);
	}
	// Solve Affine system ...
//	for (int j = 0; j < n + na + 2 * nc; j++) dvector(j) = solution_vector(j);
	for (int j = 0; j < n + na + 2 * nc; j++){
		for (int k = 0; k < n + na + 2 * nc; k++) KKT_predictor(j,k) = KKT(j,k);
	}
	dvector = KKT_predictor.partialPivLu().solve(solution_vector);

	for (int j = 0; j < n + na + 2 * nc; j++ ){
		cout<<" dvector["<<j<<"]= "<<dvector[j].get_d()<<endl;
	}

	for (int j = 0; j < n;  j++) dx(j) = dvector(j);
	for (int j = 0; j < na; j++) dy(j) = dvector(j + n);
	for (int j = 0; j < nc; j++){
		dz(j) = dvector(j + n + na);
		ds(j) = dvector(j + n + na + nc);
	}

	// Update iterate using full affine scaling
	// (x,y,z,s)->(x,y,z,s) + (dx_aff,dy_aff,dz_aff,ds_aff)
	for (int j = 0; j < n; j++){
		x(j) += dx(j);
		cout<<" x["<<j<<"]= "<<dx(j).get_d()<<endl;
	}
	for (int j = 0; j < na; j++){ 
		y(j) += dy(j);
		cout<<" y["<<j<<"]= "<<dy(j).get_d()<<endl;
	}
	for (int j = 0; j < nc; j++) {
		z(j) += dz(j);
		cout<<" z["<<j<<"]= "<<dz(j).get_d()<<endl;
	}
	for (int j = 0; j < nc; j++) {
		s(j) += ds(j);
		cout<<" s["<<j<<"]= "<<ds(j).get_d()<<endl;
	}

	// above iterate likely infeasible
	// measure violation
	std::vector < T > max_violation_list;
	for (int j = 0; j < nc; j++){
		T v1 = -1.0*z(j);
		T v2 = -1.0*s(j);
		T maximum = max_element_wrt_zero(v1, v2);
		if (maximum < 0) maximum = 0.0;
		max_violation_list.push_back(maximum);
	}
	// sort max_violation_list
	std::sort(max_violation_list.begin(), max_violation_list.end());
	T max_violation = max_violation_list[nc - 1];
	// compute the iterate shift for complementary variables z & s
	T shift = 1000.0 + 2.0*max_violation;
	// update complementary variables ...
	for (int j = 0; j < nc; j++){
		z(j) += shift;
		s(j) += shift;
	}

	bool found_soln = false;
	int iter = 0;
	while (!found_soln)
	{
		// Update KKT matrix blocks: [4][3], [4][4]
		for (int j = 0; j < nc; j++){
			// [4][3] Block
			for (int k = 0; k < nc; k++){
				if (j == k) KKT(n + na + nc + j,n + na + k) = s(j);
				else KKT(n + na + nc + j,n + na + k) = 0.0;
			}
			// [4][4] Block
			for (int k = 0; k < nc; k++){
				if (j == k) KKT(n + na + nc + j,n + na + nc + k) = z(j);
				else KKT(n + na + nc + j,n + na + nc + k) = 0.0;
			}
		}

		//_output_matrix(KKT,"KKT_matrix.txt"); // debug

		// COMPUTE RESIDUALS ...
		t1 = H*x;
		t2 = AT*y;
		t3 = CT*z;
		t4 = A*x;
		t5 = C*x;
		for (int j = 0; j < n; j++){ // rh residual = H*x - AT*y - CT*z
			rh(j) = t1(j) - t2(j) - t3(j);
			solution_vector(j) = -1.0*rh(j);
		}
		for (int j = 0; j < na; j++){ // ra residual = A*x - b
			ra(j) = t4(j) - b(j);
			solution_vector(n + j) = -1.0*ra(j);
		}
		zTs = 0.0; // re-initialize this adder
		for (int j = 0; j < nc; j++){ // rc residual = C*x - s - d
			rc(j) = t5(j) - s(j) - d(j);
			rsz(j) = s(j) * z(j); // rsz residual = Z[]*S[]*e
			solution_vector(n + na + j) = -1.0*rc(j);
			zTs += rsz(j);
			solution_vector(n + na + nc + j) = -1.0*rsz(j);
		}
		//calculate mu
		mu = zTs / nc;
		cout<<" mu = "<<mu.get_d() <<endl;
		//double mu_d = _get_double(mu); // debug
		if (iter > 5 && mu > prev_mu) return false;
		prev_mu = mu;
		if (mu < 0.00000001)
		{
			found_soln = true;
			for (int j = 0; j < n; j++) fvalues(j) = x(j); // get solution
			continue;
		}
		///////////////////////////////////////////////////////
		/////////////////// Predictor Step ////////////////////
		///////////////////////////////////////////////////////
		// Solve Affine system ...
		//cout<<" solution vector for affine system "<<endl;
//		for (int j = 0; j < n + na + 2 * nc; j++) dvector(j) = solution_vector(j);
		for (int j = 0; j < n + na + 2 * nc; j++){
			for (int k = 0; k < n + na + 2 * nc; k++) KKT_predictor(j,k) = KKT(j,k);
		}
		dvector = KKT_predictor.partialPivLu().solve(solution_vector);

		for (int j = 0; j < n + na + 2 * nc; j++ ){
			cout<<" dvector_pred["<<j<<"]= "<<dvector[j].get_d()<<endl;
		}


		// get step vectors...
		for (int j = 0; j < nc; j++){
			dz_aff(j) = dvector(j + n + na);
			ds_aff(j) = dvector(j + n + na + nc);
		}
		for (int j = 0; j < nc; j++) {
			cout<<" dz_aff["<<j<<"] = "<<dz_aff(j).get_d()<<endl;
		}
		for (int j = 0; j < nc; j++) {
			cout<<" ds_aff["<<j<<"] = "<<ds_aff(j).get_d()<<endl;
		}

		alpha = _find_step_length(s, ds_aff, z, dz_aff);
		cout<< " alpha = "<<alpha.get_d()<<endl;

		// calculate mu_aff
		elemsum = 0.0;
		for (int j = 0; j < nc; j++) elemsum += (z(j) + alpha*dz_aff(j)) * (s(j) + alpha*ds_aff(j));
		mu_aff = elemsum / nc;
		sigma = (mu_aff / mu)*(mu_aff / mu)*(mu_aff / mu);

		///////////////////////////////////////////////////////
		/////////////////// Corrector Step ////////////////////
		///////////////////////////////////////////////////////
		// modify rsz vector using the calculated correctors ...
		// recompute rsz residuals...
		// Z*S*e - sigma*mu*e + dz_aff*ds_aff*e
		for (int j = 0; j < nc; j++){
			rsz_aff(j) = rsz(j) - sigma*mu + dz_aff(j) * ds_aff(j);
			solution_vector(n + na + nc + j) = -1.0*rsz_aff(j);
		}
//		for (int j = 0; j < n + na + 2 * nc; j++) dvector_corr(j) = solution_vector(j);

		// Solve the corrected linear system ...
		for (int j = 0; j < n + na + 2 * nc; j++){
			for (int k = 0; k < n + na + 2 * nc; k++) KKT_corrector(j,k) = KKT(j,k);
		}
		dvector_corr = KKT_corrector.partialPivLu().solve(solution_vector);

		for (int j = 0; j < n + na + 2 * nc; j++ ){
			cout<<" dvector_corr["<<j<<"]= "<<dvector_corr[j].get_d()<<endl;
		}


		// get step vectors from corrector step...
		for (int j = 0; j < n; j++){
			dx(j) = dvector_corr(j);
			if (j < na) dy(j) = dvector_corr(j + n);
			if (j < nc)
			{
				dz(j) = dvector_corr(j + n + na);
				ds(j) = dvector_corr(j + n + na + nc);
			}
		}
		for (int j = 0; j < nc; j++) {
			cout<<" dz["<<j<<"] = "<<dz(j).get_d()<<endl;
		}
		for (int j = 0; j < nc; j++) {
			cout<<" ds["<<j<<"] = "<<ds(j).get_d()<<endl;
		}

		alpha = _find_step_length(s, ds, z, dz);

		cout<< " alpha = "<<alpha.get_d()<<endl;

		// update x,y,z,s vectors using alpha step
		for (int j = 0; j < n; j++) x(j) += alpha*dx(j);
		for (int j = 0; j < na; j++) y(j) += alpha*dy(j);
		for (int j = 0; j < nc; j++) z(j) += alpha*dz(j);
		for (int j = 0; j < nc; j++) s(j) += alpha*ds(j);
		iter++;
	}
	return true;
}

template <class T>
bool Math_methods::sort_vector_w_index( std::vector<T> &arr, std::vector<int> &brr )
{
	if (arr.size() != brr.size()) return false;

	const int M = 7, NSTACK = 50;
	int i, ir, j, k, jstack = -1, l = 0;
	T a;
	int b;
	std::vector<int> istack;
	istack.resize(NSTACK);

	int n = (int)arr.size();
	ir = n - 1;
	for (;;) {
		if (ir - l < M) {
			for (j = l + 1;j <= ir;j++) {
				a = arr[j];
				b = brr[j];
				for (i = j - 1;i >= l;i--) {
					if (arr[i] <= a) break;
					arr[i + 1] = arr[i];
					brr[i + 1] = brr[i];
				}
				arr[i + 1] = a;
				brr[i + 1] = b;
			}
			if (jstack < 0) break;
			ir = istack[jstack--];
			l = istack[jstack--];
		}
		else {
			k = (l + ir) >> 1;
			SWAP(arr[k], arr[l + 1]);
			SWAP(brr[k], brr[l + 1]);
			if (arr[l] > arr[ir]) {
				SWAP(arr[l], arr[ir]);
				SWAP(brr[l], brr[ir]);
			}
			if (arr[l + 1] > arr[ir]) {
				SWAP(arr[l + 1], arr[ir]);
				SWAP(brr[l + 1], brr[ir]);
			}
			if (arr[l] > arr[l + 1]) {
				SWAP(arr[l], arr[l + 1]);
				SWAP(brr[l], brr[l + 1]);
			}
			i = l + 1;
			j = ir;
			a = arr[l + 1];
			b = brr[l + 1];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i], arr[j]);
				SWAP(brr[i], brr[j]);
			}
			arr[l + 1] = arr[j];
			arr[j] = a;
			brr[l + 1] = brr[j];
			brr[j] = b;
			jstack += 2;
			if (jstack >= NSTACK) return false;
			if (ir - i + 1 >= j - l) {
				istack[jstack] = ir;
				istack[jstack - 1] = i;
				ir = j - 1;
			}
			else {
				istack[jstack] = j - 1;
				istack[jstack - 1] = l;
				l = i;
			}
		}
	}

	return true;
}

#endif