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

#include <cstdlib>
#include <math_methods.h>

double Math_methods::_find_step(const VectorXd &da, const VectorXd &a) {
	int n = (int)da.rows();
	double max = DBL_MIN;
	for (int j = 0; j < n; j++) {
		double alpha = std::max(abs(da[j] / a[j]) / 0.95, 1.0);
		if (alpha > max)
			max = alpha;
		// std::cout<<" alpha = "<<alpha<<" max= "<<max<<std::endl;
	}
	return max;
}

double Math_methods::_find_positivity_step(
	const VectorXd &da, const VectorXd &a,
	const VectorXd &db, const VectorXd &b,
	const VectorXd &dc, const VectorXd &c,
	const VectorXd &dd, const VectorXd &d)
{
	double max = DBL_MIN;

	double max_a = _find_step(da, a);
	double max_b = _find_step(db, b);
	double max_c = _find_step(dc, c);
	double max_d = _find_step(dd, d);

	if (max_a > max)
		max = max_a;
	if (max_b > max)
		max = max_b;
	if (max_c > max)
		max = max_c;
	if (max_d > max)
		max = max_d;

	return max;
}

double Math_methods::RandomDouble(const double &min, const double &max)
{
	double f = (double)rand() / RAND_MAX;
	return min + f * (max - min);
}

bool Math_methods::quadratic_solver_loqo(
	const MatrixXd &H, const MatrixXd &A,
	const VectorXd &b, const VectorXd &r,
	VectorXd &fvalues)
{
	// minimize f = 1/2 xT H x
	// s.t. b <= Ax <= b + r
	int n = (int)H.rows();

	// double max_ele = H.maxCoeff();

	MatrixXd KKT(2 * n, 2 * n);
	VectorXd c(n);
	c.setZero();
	VectorXd rhs(2 * n);
	rhs << c, b;
	// std::cout<<" rhs:\n"<< rhs << std::endl;

	KKT << -(H + MatrixXd::Identity(n, n)), A.transpose(), A, MatrixXd::Identity(n, n);
	// std::cout<<" Initial KKT matrix:\n"<< KKT << std::endl;

	VectorXd soln(2 * n);
	soln = KKT.partialPivLu().solve(rhs);
	// std::cout<<" Soln from Initial KKT matrix:\n"<< soln << std::endl;
	VectorXd x(n);
	x = soln.head(n);
	VectorXd y(n);
	y = soln.tail(n);
	VectorXd g(n);
	VectorXd z(n);
	VectorXd t(n);
	VectorXd s(n);
	VectorXd v(n);
	VectorXd w(n);
	VectorXd p(n);
	VectorXd q(n);
	for (int j = 0; j < n; j++) {
		g(j) = std::max(abs(x(j)), 100.0);
		z(j) = g(j);
		t(j) = g(j);
		s(j) = g(j);
		v(j) = std::max(abs(y(j)), 100.0);
		w(j) = v(j);
		p(j) = std::max(abs(r(j) - w(j)), 100.0);
		q(j) = v(j);
	}
	MatrixXd DebugMatrixV(n, 10);
	DebugMatrixV.col(0) = x;
	DebugMatrixV.col(1) = y;
	DebugMatrixV.col(2) = g;
	DebugMatrixV.col(3) = z;
	DebugMatrixV.col(4) = t;
	DebugMatrixV.col(5) = s;
	DebugMatrixV.col(6) = v;
	DebugMatrixV.col(7) = w;
	DebugMatrixV.col(8) = p;
	DebugMatrixV.col(9) = q;
	// std::cout<<" Current Variable matrix:\n"<< DebugMatrixV << std::endl;

	double mu = (z.dot(g) + v.dot(w) + s.dot(t) + p.dot(q)) / 4 * n;

	// construct G Z, V W, S T, and P Q diagonal matrices
	MatrixXd G(n, n);
	MatrixXd Z(n, n);
	MatrixXd V(n, n);
	MatrixXd W(n, n);
	MatrixXd S(n, n);
	MatrixXd T(n, n);
	MatrixXd P(n, n);
	MatrixXd Q(n, n);
	// below matrices are computed after predictor step
	MatrixXd dG(n, n);
	MatrixXd dV(n, n);
	MatrixXd dT(n, n);
	MatrixXd dP(n, n);

	// helper variables
	VectorXd rho(n);
	VectorXd nu(n);
	VectorXd alpha(n);
	VectorXd sigma(n);
	VectorXd tau(n);
	VectorXd beta(n);
	VectorXd gamma_z(n);
	VectorXd gamma_w(n);
	VectorXd gamma_s(n);
	VectorXd gamma_q(n);
	VectorXd tauh(n);
	VectorXd betah(n);
	VectorXd alphah(n);
	VectorXd nuh(n);
	VectorXd ev(n);
	ev.setOnes();

	// step directions
	VectorXd dx(n);
	VectorXd dy(n);
	VectorXd dg(n);
	VectorXd dz(n);
	VectorXd dt(n);
	VectorXd ds(n);
	VectorXd dv(n);
	VectorXd dw(n);
	VectorXd dp(n);
	VectorXd dq(n);
	// MatrixXd DebugMatrixStepV(n,10);

	bool converged = false;
	double last_iterations_sig_fig = 0.0;
	double primal_infeasibility = sqrt(rho.dot(rho) + tau.dot(tau) + alpha.dot(alpha) + nu.dot(nu)) / (sqrt(b.dot(b)) + 1.0);
	double dual_infeasibilitiy = sqrt(sigma.dot(sigma) + beta.dot(beta));
	double last_primal_infeasibility = primal_infeasibility;
	double last_dual_infeasibility = dual_infeasibilitiy;
	double last_primal = 0;
	double last_dual = 0;
	VectorXd iter_minus_one_x = x;
	int iter = 0;
	while (!converged) {
		iter++;
		// 		MatrixXd DebugV(n,10);
		// 		DebugV.col(0) = x;
		// 		DebugV.col(1) = y;
		// 		DebugV.col(2) = g;
		// 		DebugV.col(3) = z;
		// 		DebugV.col(4) = t;
		// 		DebugV.col(5) = s;
		// 		DebugV.col(6) = v;
		// 		DebugV.col(7) = w;
		// 		DebugV.col(8) = p;
		// 		DebugV.col(9) = q;
		// 		std::cout<<" Current Variable matrix (x,y,g,z,t,s,v,w,p,q):\n"<< DebugV
		// << std::endl;

		double primal_obj = 0.5 * x.transpose() * H * x;
		double dual_obj = b.dot(y) - 0.5 * x.transpose() * H * x - r.dot(q);

		double sigfig = std::max(-std::log10(abs(primal_obj - dual_obj) / (abs(primal_obj) + 1.0)), 0.0);

		primal_infeasibility = sqrt(rho.dot(rho) + tau.dot(tau) + alpha.dot(alpha) + nu.dot(nu)) / (sqrt(b.dot(b)) + 1.0);
		dual_infeasibilitiy = sqrt(sigma.dot(sigma) + beta.dot(beta));

		std::cout << " Iteration[" << iter << "]" << std::endl;
		std::cout << "	Primal_obj = " << primal_obj << std::endl;
		std::cout << "	Dual_obj = " << dual_obj << std::endl;
		std::cout << "	Significant figures = " << sigfig << std::endl;
		std::cout << "	Primal Infeasibility = " << primal_infeasibility << std::endl;
		std::cout << "	Dual Infeasibility = " << dual_infeasibilitiy << std::endl;

		if (sigfig > 6)
		{
			converged = true;
			fvalues = x; // strong duality gap
			break;
		}
		if (dual_obj > primal_obj || sigfig < last_iterations_sig_fig)
			return false;

		last_iterations_sig_fig = sigfig;
		last_primal_infeasibility = primal_infeasibility;
		last_dual_infeasibility = dual_infeasibilitiy;
		last_primal = primal_obj;
		last_dual = dual_obj;
		iter_minus_one_x = x;

		G = g.asDiagonal();
		Z = z.asDiagonal();
		V = v.asDiagonal();
		W = w.asDiagonal();
		S = s.asDiagonal();
		T = t.asDiagonal();
		P = p.asDiagonal();
		Q = q.asDiagonal();

		MatrixXd D(n, n);
		D = (S.inverse() * T + G * Z.inverse()).inverse();
		MatrixXd E(n, n);
		E = (V * W.inverse() + P.inverse() * Q).inverse();

		rho = b - A * x + w;
		nu = -x + g - t;
		alpha = r - w - p;
		sigma = -A.transpose() * y - z + H * x;
		tau = -z - s;
		beta = y + q - v;
		// predictor nonlinearities
		gamma_z = -z;
		gamma_w = -w;
		gamma_s = -s;
		gamma_q = -q;

		tauh = tau - gamma_s;
		betah = beta - V * W.inverse() * gamma_w;
		alphah = alpha - P * Q.inverse() * gamma_q;
		nuh = nu + G * Z.inverse() * gamma_z;

		// 		MatrixXd RHSV(n,14);
		// 		RHSV.col(0) = rho;
		// 		RHSV.col(1) = nu;
		// 		RHSV.col(2) = alpha;
		// 		RHSV.col(3) = sigma;
		// 		RHSV.col(4) = tau;
		// 		RHSV.col(5) = beta;
		// 		RHSV.col(6) = gamma_z;
		// 		RHSV.col(7) = gamma_w;
		// 		RHSV.col(8) = gamma_s;
		// 		RHSV.col(9) = gamma_q;
		// 		RHSV.col(10) = tauh;
		// 		RHSV.col(11) = betah;
		// 		RHSV.col(12) = alphah;
		// 		RHSV.col(13) = nuh;
		// 		std::cout<<" Current RHS matrix
		// (rho,nu,alpha,sigma,tau,beta,gamma_z,gamma_w,gamma_s,gamma_q,tauh,betah,alphah,nuh):\n"<<
		// RHSV << std::endl;

		rhs << (sigma - D * (nuh + S.inverse() * T * tauh)), (rho - E * (betah - P.inverse() * Q * alphah));

		KKT << -(H + D), A.transpose(), A, E;
		// 		std::cout<<" Predictor KKT matrix:\n"<< KKT << std::endl;
		// 		std::cout<<" Predictor rhs:\n"<< rhs << std::endl;

		soln = KKT.partialPivLu().solve(rhs);

		if (!soln.allFinite())
		{
			std::cout << " Numerical issue with solving linear system..." << std::endl;
			return false;
		}

		// get "delta" variables for predictor system ...
		dx = soln.head(n);
		dy = soln.tail(n);

		// 		std::cout<<" Soln from Predictor KKT matrix:\n"<< soln <<
		// std::endl; 		std::cout<<" dx:\n"<< dx << std::endl; 		std::cout<<" dy:\n"<<
		// dy << std::endl;

		dw = -E * (betah - P.inverse() * Q * alphah + dy);
		dt = -D * S.inverse() * T * (G * Z.inverse() * tauh - nuh + dx);
		dz = G.inverse() * Z * (nuh - dx - dt);
		dq = P.inverse() * Q * (dw - alphah);
		dv = V * W.inverse() * (gamma_w - dw);
		ds = gamma_s - S * T.inverse() * dt;
		dp = P * Q.inverse() * (gamma_q - dq);
		dg = G * Z.inverse() * (gamma_z - dz);

		// Debug
		//  		DebugMatrixStepV.col(0) = dx;
		// 		DebugMatrixStepV.col(1) = dy;
		//  		DebugMatrixStepV.col(2) = dg;
		//  		DebugMatrixStepV.col(3) = dz;
		//  		DebugMatrixStepV.col(4) = dt;
		//  		DebugMatrixStepV.col(5) = ds;
		//  		DebugMatrixStepV.col(6) = dv;
		//  		DebugMatrixStepV.col(7) = dw;
		//  		DebugMatrixStepV.col(8) = dp;
		//  		DebugMatrixStepV.col(9) = dq;
		//  		std::cout<<" Current Step Variable matrix (after
		//  predictor step) (dx,dy,dg,dz,dt,ds,dv,dw,dp,dq):\n"<< DebugMatrixStepV
		//  << std::endl;

		// compute step lengths for primal and dual systems
		double alpha_p = _find_positivity_step(dg, g, dw, w, dt, t, dp, p);
		double alpha_d = _find_positivity_step(dz, z, dv, v, ds, s, dq, q);
		// std::cout<<" predictor alpha_p = "<<alpha_p <<" alpha_d= "<<alpha_d <<
		// std::endl;
		// alpha_pd will be the maximum of the two above values
		double alpha_pd = std::max(alpha_p, alpha_d);
		// std::cout<<" alpha_pd = "<<alpha_pd << std::endl;
		double fraction = pow(((alpha_pd - 1.0) / (alpha_pd + 10.0)), 2);

		// update mu
		mu = (z.dot(g) + v.dot(w) + s.dot(t) + p.dot(q)) * (fraction) / (4 * n);
		std::cout << "	mu (predictor) = " << mu << " fraction = " << fraction << std::endl;

		// update rhs variables rho,nu,alpha,sigma,tau,beta,gamma's
		// first compute dG,dV,dT,and dP matrices from dg,dv,dt, and dp vectors
		// above
		dG = dg.asDiagonal();
		dV = dv.asDiagonal();
		dT = dt.asDiagonal();
		dP = dp.asDiagonal();
		gamma_z = mu * G.inverse() * ev - z - G.inverse() * dG * dz;
		gamma_w = mu * V.inverse() * ev - w - V.inverse() * dV * dw;
		gamma_s = mu * T.inverse() * ev - s - T.inverse() * dT * ds;
		gamma_q = mu * P.inverse() * ev - q - P.inverse() * dP * dq;
		tauh = tau - gamma_s;
		betah = beta - V * W.inverse() * gamma_w;
		alphah = alpha - P * Q.inverse() * gamma_q;
		nuh = nu + G * Z.inverse() * gamma_z;

		rhs << (sigma - D * (nuh + S.inverse() * T * tauh)),
			(rho - E * (betah - P.inverse() * Q * alphah));

		// 		std::cout<<" Corrector KKT matrix:\n"<< KKT << std::endl;
		// 		std::cout<<" Corrector rhs:\n"<< rhs << std::endl;

		soln = KKT.partialPivLu().solve(rhs);

		if (!soln.allFinite()) {
			std::cout << " Numerical issue with solving linear system..." << std::endl;
			return false;
		}

		// get "delta" variables for corrector system ...
		dx = soln.head(n);
		dy = soln.tail(n);

		// 		std::cout<<" Soln from Corrector KKT matrix:\n"<< soln <<
		// std::endl; 		std::cout<<" dx:\n"<< dx << std::endl; 		std::cout<<" dy:\n"<<
		// dy << std::endl;

		dw = -E * (betah - P.inverse() * Q * alphah + dy);
		dt = -D * S.inverse() * T * (G * Z.inverse() * tauh - nuh + dx);
		dz = G.inverse() * Z * (nuh - dx - dt);
		dq = P.inverse() * Q * (dw - alphah);
		dv = V * W.inverse() * (gamma_w - dw);
		ds = gamma_s - S * T.inverse() * dt;
		dp = P * Q.inverse() * (gamma_q - dq);
		dg = G * Z.inverse() * (gamma_z - dz);

		// Debug
		// 		DebugMatrixStepV.col(0) = dx;
		// 		DebugMatrixStepV.col(1) = dy;
		// 		DebugMatrixStepV.col(2) = dg;
		// 		DebugMatrixStepV.col(3) = dz;
		// 		DebugMatrixStepV.col(4) = dt;
		// 		DebugMatrixStepV.col(5) = ds;
		// 		DebugMatrixStepV.col(6) = dv;
		// 		DebugMatrixStepV.col(7) = dw;
		// 		DebugMatrixStepV.col(8) = dp;
		// 		DebugMatrixStepV.col(9) = dq;
		// 		std::cout<<" Current Step Variable matrix (after corrector step)
		// (dx,dy,dg,dz,dt,ds,dv,dw,dp,dq):\n"<< DebugMatrixStepV << std::endl;

		// compute step lengths for primal and dual systems
		alpha_p = _find_positivity_step(dg, g, dw, w, dt, t, dp, p);
		alpha_d = _find_positivity_step(dz, z, dv, v, ds, s, dq, q);

		// std::cout<<" alpha_p= "<<alpha_p<<" alpha_d= "<<alpha_d<<std::endl;

		// update solution
		x += (1 / alpha_p) * dx;
		g += (1 / alpha_p) * dg;
		w += (1 / alpha_p) * dw;
		t += (1 / alpha_p) * dt;
		p += (1 / alpha_p) * dp;

		y += (1 / alpha_d) * dy;
		z += (1 / alpha_d) * dz;
		v += (1 / alpha_d) * dv;
		s += (1 / alpha_d) * ds;
		q += (1 / alpha_d) * dq;
	}

	return true;
}

bool Math_methods::angle_btw_2_vectors(
	const std::vector<double> &v1, const std::vector<double> &v2, double &angle)
{
	if (v1.size() != v2.size())
		return false;

	double dotp = 0.0;
	double v1norm = 0.0;
	double v2norm = 0.0;
	for (int j = 0; j < (int)v1.size(); j++) {
		dotp += v1[j] * v2[j];
		v1norm += v1[j] * v1[j];
		v2norm += v2[j] * v2[j];
	}

	angle = acos(dotp / (sqrt(v1norm) * sqrt(v2norm)));
	return true;
}

double Math_methods::max_element_wrt_zero(const double &a, const double &b) {
	double max_value = a;
	double c = 0.0;
	if (b > max_value)
		max_value = b;
	if (c > max_value)
		max_value = c;

	return max_value;
}

double Math_methods::_find_step_length(
	const VectorXd &a,
	const VectorXd &da,
	const VectorXd &b,
	const VectorXd &db)
{
	int n = (int)a.rows();
	// if (a.rows() != b.rows() || a.rows() != da.rows() || a.rows() != db.rows()
	// || da.rows() != db.rows()) throw -1;
	// find step length (a) ... this is the most important step
	// has to satisfy all the non-negativity conditions ...
	double min_alpha_a = 100.0; // improper initialization. proper initialization will
							// occur in first iteration of below for loop
	double min_alpha_b = 100.0; // ?????		why do i do this? seems dumb
	double max_alpha_a = 0.0;   // improper initialization. proper initialization will
							// occur in first iteration of below for loop
	double max_alpha_b = 0.0;
	double alpha = 0.0;
	VectorXd alpha_a(n);
	VectorXd alpha_b(n);
	for (int j = 0; j < n; j++) {
		// cout<<" z["<<j<<"]= "<<mpf_get_d(z[j])<<endl;
		// cout<<" s["<<j<<"]= "<<mpf_get_d(s[j])<<endl;
		// cout<<" ["<<j<<"] dz_aff = "<<mpf_get_d(dz_aff[j])<<" ds_aff=
		// "<<mpf_get_d(ds_aff[j])<<endl;
		// do alpha_b first ...
		if (b(j) > 0.0)
			alpha_b(j) = b(j) / db(j);
		else
			alpha_b(j) = -1.0 * b(j) / db(j);
		// alpha_a
		if (a(j) > 0.0)
			alpha_a(j) = a(j) / da(j);
		else
			alpha_a(j) = -1.0 * a(j) / da(j);

		if (alpha_b(j) < min_alpha_b && alpha_b(j) > 0.00000000000001)
			min_alpha_b = alpha_b(j);
		if (alpha_a(j) < min_alpha_a && alpha_a(j) > 0.00000000000001)
			min_alpha_a = alpha_a(j);
		// cout<<" alpha_s["<<j<<"]= "<<mpf_get_d(alpha_s[j])<<endl;
		// cout<<" alpha_z["<<j<<"]= "<<mpf_get_d(alpha_z[j])<<endl;
		// cout<<" max_alpha_s = "<<mpf_get_d(max_alpha_s)<<endl;
		// cout<<" max_alpha_z = "<<mpf_get_d(max_alpha_z)<<endl;
		// cout<<" min_alpha_s = "<<mpf_get_d(min_alpha_s)<<endl;
		// cout<<" min_alpha_z = "<<mpf_get_d(min_alpha_z)<<endl;
	}

	if (min_alpha_b < min_alpha_a)
		alpha = min_alpha_b;
	else
		alpha = min_alpha_a;

	if (alpha > 1.0 || alpha == 0)
		alpha = 1.0;
	return alpha;
}

bool Math_methods::quadratic_solver(
	const MatrixXd &H,
	const MatrixXd &A,
	const MatrixXd &C,
	const VectorXd &b,
	const VectorXd &d,
	VectorXd &fvalues)
{
	// Describe the Quadratic Optimization problem
	// min (w.r.t. x) f(x) = 1/2 * xT * H * x => our objective function
	// s.t.            T * x = b
	//                C * x >= d
	// NOTE: H must be a positive definite matrix !!! There is a check for this in
	// the constructor for the matrix solver The Interior Point method using a
	// Predictor-Corrector method will be used to solve the QP Lagrangian (Primal)
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
	// EQN (1) can be reduced in size by block elimination since S and Z are
	// strictly positive New EQN becomes | H  AT      CT   || dx|    |     -rh |
	// | T   0	     0   ||-dy| =  |     -ra     |
	// | C   0 -Z^(-1)*S ||-dz|    |-rc-Z^(-1)rsz|

	int n = (int)H.rows();
	int na = (int)A.rows();
	int nc = (int)C.rows();

	///////////////////////////////////////////////////////////////////////
	/////////////////////////// Initialization ////////////////////////////
	///////////////////////////////////////////////////////////////////////
	// Calculation Helper Variables
	double elemsum, datanorm, zTs, mu, zero, prev_mu;
	zero = 0.0;
	elemsum = 0.0;
	// KKT System
	VectorXd x(n);
	x.setZero();
	VectorXd y(na);
	y.setZero();
	VectorXd z(nc);
	z.setZero();
	VectorXd s(nc);
	s.setZero();
	VectorXd rh(n);
	rh.setZero();
	VectorXd ra(na);
	ra.setZero();
	VectorXd rc(nc);
	rc.setZero();
	VectorXd rsz(nc);
	rsz.setZero();

	VectorXd t1(n);
	t1.setZero();
	VectorXd t2(n);
	t2.setZero();
	VectorXd t3(n);
	t3.setZero();
	VectorXd t4(na);
	t4.setZero();
	VectorXd t5(nc);
	t5.setZero();

	MatrixXd AT = A.transpose();
	MatrixXd CT = C.transpose();
	MatrixXd KKT(n + na + 2 * nc, n + na + 2 * nc);
	MatrixXd KKT_predictor(n + na + 2 * nc, n + na + 2 * nc);
	MatrixXd KKT_corrector(n + na + 2 * nc, n + na + 2 * nc);

	VectorXd solution_vector(n + na + 2 * nc);
	solution_vector.setZero();

	// step containers and related helper variables
	double alpha, mu_aff, sigma;
	VectorXd alpha_z(nc);
	alpha_z.setZero();
	VectorXd alpha_s(nc);
	alpha_s.setZero();
	VectorXd dx_aff(n);
	dx_aff.setZero();
	VectorXd dy_aff(na);
	dy_aff.setZero();
	VectorXd dz_aff(nc);
	dz_aff.setZero();
	VectorXd ds_aff(nc);
	ds_aff.setZero();
	VectorXd rsz_aff(nc);
	rsz_aff.setZero();
	VectorXd dvector(n + na + 2 * nc);
	dvector.setZero();

	VectorXd dx(n);
	dx.setZero();
	VectorXd dy(na);
	dy.setZero();
	VectorXd dz(nc);
	dz.setZero();
	VectorXd ds(nc);
	ds.setZero();
	VectorXd dvector_corr(n + na + 2 * nc);
	dvector_corr.setZero();

	///////////////////////////////////////////////////////////////////////
	////////////////////// End of Initialization //////////////////////////
	///////////////////////////////////////////////////////////////////////

	// get norm of matrix
	datanorm = sqrt(H.maxCoeff());

	// setup initial iterate for z and s...
	for (int j = 0; j < nc; j++) {
		z(j) = datanorm;
		s(j) = datanorm;
	}

	// Build KKT matrix blocks [1][1], [1][2], [1][3], [1][4]
	for (int j = 0; j < n; j++) { // Matrix rows...
		// Matrix columns...
		// [1][1] Block
		for (int k = 0; k < n; k++)
			KKT(j, k) = H(j, k);
		// [1][2] Block
		for (int k = 0; k < na; k++)
			KKT(j, n + k) = -1.0 * AT(j, k);
		// [1][3] Block
		for (int k = 0; k < nc; k++)
			KKT(j, n + na + k) = -1.0 * CT(j, k);
		// [1][4] Block
		for (int k = 0; k < nc; k++)
			KKT(j, n + na + nc + k) = 0.0;
	}
	// Build KKT matrix blocks [2][1], [2][2], [2][3], [2][4]
	for (int j = 0; j < na; j++) {
		// [2][1] Block
		for (int k = 0; k < n; k++)
			KKT(n + j, k) = A(j, k);
		// [2][2] Block
		for (int k = 0; k < na; k++)
			KKT(n + j, n + k) = 0.0;
		// [2][3] Block
		for (int k = 0; k < nc; k++)
			KKT(n + j, n + na + k) = 0.0;
		// [2][4] Block
		for (int k = 0; k < nc; k++)
			KKT(n + j, n + na + nc + k) = 0.0;
	}
	// Build KKT matrix blocks [3][1], [3][2], [3][3], [3][4]
	for (int j = 0; j < nc; j++) {
		// [3][1] Block
		for (int k = 0; k < n; k++)
			KKT(n + na + j, k) = C(j, k);
		// [3][2] Block
		for (int k = 0; k < na; k++)
			KKT(n + na + j, n + k) = 0.0;
		// [3][3] Block
		for (int k = 0; k < nc; k++)
			KKT(n + na + j, n + na + k) = 0.0;
		// [3][4] Block
		for (int k = 0; k < nc; k++) {
			if (j == k)
				KKT(n + na + j, n + na + nc + k) = -1.0;
			else
				KKT(n + na + j, n + na + nc + k) = 0.0;
		}
	}
	// Build KKT matrix blocks [4][1], [4][2]
	for (int j = 0; j < nc; j++) {
		// [4][1] Block
		for (int k = 0; k < n; k++)
			KKT(n + na + nc + j, k) = 0.0;
		// [4][2] Block
		for (int k = 0; k < na; k++)
			KKT(n + na + nc + j, n + k) = 0.0;
	}

	// Build KKT matrix blocks: [4][3], [4][4]
	for (int j = 0; j < nc; j++) {
		// [4][3] Block
		for (int k = 0; k < nc; k++) {
			if (j == k)
				KKT(n + na + nc + j, n + na + k) = s(j);
			else
				KKT(n + na + nc + j, n + na + k) = 0.0;
		}
		// [4][4] Block
		for (int k = 0; k < nc; k++) {
			if (j == k)
				KKT(n + na + nc + j, n + na + nc + k) = z(j);
			else
				KKT(n + na + nc + j, n + na + nc + k) = 0.0;
		}
	}

	//////////////////////////////////////////////
	// Determine Starting Point for the iterate //
	//////////////////////////////////////////////
	// 1st compute residuals for the affine system
	t1 = H * x;
	t2 = AT * y;
	t3 = CT * z;
	t4 = A * x;
	t5 = C * x;
	// 	_output_vector(t3,"t3.txt"); // debug
	// 	_output_vector(z,"z.txt");
	// 	_output_matrix(CT,"CT.txt");
	for (int j = 0; j < n; j++) { // rh residual = H*x - AT*y - CT*z
		rh(j) = t1(j) - t2(j) - t3(j);
		solution_vector(j) = -1.0 * rh(j);
	}
	for (int j = 0; j < na; j++) { // ra residual = A*x - b
		ra(j) = t4(j) - b(j);
		solution_vector(n + j) = -1.0 * ra(j);
	}
	zTs = 0.0;                     // re-initialize this adder
	for (int j = 0; j < nc; j++) { // rc residual = C*x - s - d
		rc(j) = t5(j) - s(j) - d(j);
		rsz(j) = s(j) * z(j); // rsz residual = Z[]*S[]*e
		solution_vector(n + na + j) = -1.0 * rc(j);
		zTs += rsz(j);
		solution_vector(n + na + nc + j) = -1.0 * rsz(j);
	}
	// Solve Affine system ...
	//	for (int j = 0; j < n + na + 2 * nc; j++) dvector(j) =
	//solution_vector(j);
	for (int j = 0; j < n + na + 2 * nc; j++) {
		for (int k = 0; k < n + na + 2 * nc; k++)
			KKT_predictor(j, k) = KKT(j, k);
	}
	dvector = KKT_predictor.partialPivLu().solve(solution_vector);

	// 	for (int j = 0; j < n + na + 2 * nc; j++ ){
	// 		cout<<" dvector["<<j<<"]= "<<dvector[j].get_d()<<endl;
	// 	}

	for (int j = 0; j < n; j++)
		dx(j) = dvector(j);
	for (int j = 0; j < na; j++)
		dy(j) = dvector(j + n);
	for (int j = 0; j < nc; j++) {
		dz(j) = dvector(j + n + na);
		ds(j) = dvector(j + n + na + nc);
	}

	// Update iterate using full affine scaling
	// (x,y,z,s)->(x,y,z,s) + (dx_aff,dy_aff,dz_aff,ds_aff)
	for (int j = 0; j < n; j++)
		x(j) += dx(j);
	for (int j = 0; j < na; j++)
		y(j) += dy(j);
	for (int j = 0; j < nc; j++)
		z(j) += dz(j);
	for (int j = 0; j < nc; j++)
		s(j) += ds(j);

	// above iterate likely infeasible
	// measure violation
	std::vector<double> max_violation_list;
	for (int j = 0; j < nc; j++) {
		double v1 = -1.0 * z(j);
		double v2 = -1.0 * s(j);
		double maximum = max_element_wrt_zero(v1, v2);
		if (maximum < 0)
			maximum = 0.0;
		max_violation_list.push_back(maximum);
	}
	// sort max_violation_list
	std::sort(max_violation_list.begin(), max_violation_list.end());
	double max_violation = max_violation_list[nc - 1];
	// compute the iterate shift for complementary variables z & s
	double shift = 1000.0 + 2.0 * max_violation;
	// update complementary variables ...
	for (int j = 0; j < nc; j++) {
		z(j) += shift;
		s(j) += shift;
	}

	bool found_soln = false;
	int iter = 0;
	while (!found_soln) {
		// Update KKT matrix blocks: [4][3], [4][4]
		for (int j = 0; j < nc; j++) {
			// [4][3] Block
			for (int k = 0; k < nc; k++) {
				if (j == k)
					KKT(n + na + nc + j, n + na + k) = s(j);
				else
					KKT(n + na + nc + j, n + na + k) = 0.0;
			}
			// [4][4] Block
			for (int k = 0; k < nc; k++) {
				if (j == k)
					KKT(n + na + nc + j, n + na + nc + k) = z(j);
				else
					KKT(n + na + nc + j, n + na + nc + k) = 0.0;
			}
		}

		//_output_matrix(KKT,"KKT_matrix.txt"); // debug

		// COMPUTE RESIDUALS ...
		t1 = H * x;
		t2 = AT * y;
		t3 = CT * z;
		t4 = A * x;
		t5 = C * x;
		for (int j = 0; j < n; j++) { // rh residual = H*x - AT*y - CT*z
			rh(j) = t1(j) - t2(j) - t3(j);
			solution_vector(j) = -1.0 * rh(j);
		}
		for (int j = 0; j < na; j++) { // ra residual = A*x - b
			ra(j) = t4(j) - b(j);
			solution_vector(n + j) = -1.0 * ra(j);
		}
		zTs = 0.0;                     // re-initialize this adder
		for (int j = 0; j < nc; j++) { // rc residual = C*x - s - d
			rc(j) = t5(j) - s(j) - d(j);
			rsz(j) = s(j) * z(j); // rsz residual = Z[]*S[]*e
			solution_vector(n + na + j) = -1.0 * rc(j);
			zTs += rsz(j);
			solution_vector(n + na + nc + j) = -1.0 * rsz(j);
		}
		// calculate mu
		mu = zTs / nc;
		double mu_d = _get_double(mu); // debug
		std::cout << " mu[" << iter << "]= " << mu_d << std::endl;
		if (iter > 5 && mu > prev_mu) {
			found_soln = true;
			for (int j = 0; j < n; j++)
				fvalues(j) = x(j); // get solution
			continue;
			// 			std::cout << " QPP did not converge. Due to precision problem..." << std::endl;
			// 			return false;
		}
		prev_mu = mu;
		if (mu < 0.00000001) {
			found_soln = true;
			for (int j = 0; j < n; j++)
				fvalues(j) = x(j); // get solution
			continue;
		}
		///////////////////////////////////////////////////////
		/////////////////// Predictor Step ////////////////////
		///////////////////////////////////////////////////////
		// Solve Affine system ...
		// cout<<" solution vector for affine system "<<endl;
		//		for (int j = 0; j < n + na + 2 * nc; j++) dvector(j) =
		//solution_vector(j);
		for (int j = 0; j < n + na + 2 * nc; j++) {
			for (int k = 0; k < n + na + 2 * nc; k++)
				KKT_predictor(j, k) = KKT(j, k);
		}
		dvector = KKT_predictor.partialPivLu().solve(solution_vector);

		// get step vectors...
		for (int j = 0; j < nc; j++) {
			dz_aff(j) = dvector(j + n + na);
			ds_aff(j) = dvector(j + n + na + nc);
		}

		alpha = _find_step_length(s, ds_aff, z, dz_aff);

		// calculate mu_aff
		elemsum = 0.0;
		for (int j = 0; j < nc; j++)
			elemsum += (z(j) + alpha * dz_aff(j)) * (s(j) + alpha * ds_aff(j));
		mu_aff = elemsum / nc;
		sigma = (mu_aff / mu) * (mu_aff / mu) * (mu_aff / mu);

		///////////////////////////////////////////////////////
		/////////////////// Corrector Step ////////////////////
		///////////////////////////////////////////////////////
		// modify rsz vector using the calculated correctors ...
		// recompute rsz residuals...
		// Z*S*e - sigma*mu*e + dz_aff*ds_aff*e
		for (int j = 0; j < nc; j++) {
			rsz_aff(j) = rsz(j) - sigma * mu + dz_aff(j) * ds_aff(j);
			solution_vector(n + na + nc + j) = -1.0 * rsz_aff(j);
		}
		//		for (int j = 0; j < n + na + 2 * nc; j++) dvector_corr(j) =
		//solution_vector(j);

		// Solve the corrected linear system ...
		for (int j = 0; j < n + na + 2 * nc; j++) {
			for (int k = 0; k < n + na + 2 * nc; k++)
				KKT_corrector(j, k) = KKT(j, k);
		}
		dvector_corr = KKT_corrector.partialPivLu().solve(solution_vector);

		// get step vectors from corrector step...
		for (int j = 0; j < n; j++) {
			dx(j) = dvector_corr(j);
			if (j < na)
				dy(j) = dvector_corr(j + n);
			if (j < nc) {
				dz(j) = dvector_corr(j + n + na);
				ds(j) = dvector_corr(j + n + na + nc);
			}
		}

		alpha = _find_step_length(s, ds, z, dz);

		// update x,y,z,s vectors using alpha step
		for (int j = 0; j < n; j++)
			x(j) += alpha * dx(j);
		for (int j = 0; j < na; j++)
			y(j) += alpha * dy(j);
		for (int j = 0; j < nc; j++)
			z(j) += alpha * dz(j);
		for (int j = 0; j < nc; j++)
			s(j) += alpha * ds(j);
		iter++;
	}
	return true;
}

bool Math_methods::sort_vector_w_index(std::vector<double> &arr,
	std::vector<int> &brr) {
	if (arr.size() != brr.size())
		return false;

	const int M = 7, NSTACK = 50;
	int i, ir, j, k, jstack = -1, l = 0;
	double a;
	int b;
	std::vector<int> istack;
	istack.resize(NSTACK);

	int n = (int)arr.size();
	ir = n - 1;
	for (;;) {
		if (ir - l < M) {
			for (j = l + 1; j <= ir; j++) {
				a = arr[j];
				b = brr[j];
				for (i = j - 1; i >= l; i--) {
					if (arr[i] <= a)
						break;
					arr[i + 1] = arr[i];
					brr[i + 1] = brr[i];
				}
				arr[i + 1] = a;
				brr[i + 1] = b;
			}
			if (jstack < 0)
				break;
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
				do
					i++;
				while (arr[i] < a);
				do
					j--;
				while (arr[j] > a);
				if (j < i)
					break;
				SWAP(arr[i], arr[j]);
				SWAP(brr[i], brr[j]);
			}
			arr[l + 1] = arr[j];
			arr[j] = a;
			brr[l + 1] = brr[j];
			brr[j] = b;
			jstack += 2;
			if (jstack >= NSTACK)
				return false;
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