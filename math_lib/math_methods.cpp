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
  int n = da.rows();
  double max = DBL_MIN;
  for (int j = 0; j < n; j++) {
    double alpha = std::max(abs(da[j] / a[j]) / 0.95, 1.0);
    if (alpha > max)
      max = alpha;
    // std::cout<<" alpha = "<<alpha<<" max= "<<max<<std::endl;
  }
  return max;
}

double Math_methods::_find_positivity_step(const VectorXd &da, const VectorXd &a,
                                    const VectorXd &db, const VectorXd &b,
                                    const VectorXd &dc, const VectorXd &c,
                                    const VectorXd &dd, const VectorXd &d) {
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

bool Math_methods::quadratic_solver_loqo(const MatrixXd &H, const MatrixXd &A,
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

	KKT << -(H + MatrixXd::Identity(n, n)), A.transpose(), A,
		MatrixXd::Identity(n, n);
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
	double primal_infeasibility = sqrt(rho.dot(rho) + tau.dot(tau) + alpha.dot(alpha) + nu.dot(nu)) /(sqrt(b.dot(b)) + 1.0);
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

		double sigfig = std::max(
			-std::log10(abs(primal_obj - dual_obj) / (abs(primal_obj) + 1.0)), 0.0);

		primal_infeasibility = sqrt(rho.dot(rho) + tau.dot(tau) + alpha.dot(alpha) + nu.dot(nu)) /(sqrt(b.dot(b)) + 1.0);
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

		rhs << (sigma - D * (nuh + S.inverse() * T * tauh)),
			(rho - E * (betah - P.inverse() * Q * alphah));

		KKT << -(H + D), A.transpose(), A, E;
		// 		std::cout<<" Predictor KKT matrix:\n"<< KKT << std::endl;
		// 		std::cout<<" Predictor rhs:\n"<< rhs << std::endl;

		soln = KKT.partialPivLu().solve(rhs);

		if (!soln.allFinite())
		{
			std::cout << " Numerical issue with solving linear system..."
					<< std::endl;
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
		std::cout << "	mu (predictor) = " << mu << " fraction = " << fraction
					<< std::endl;

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
			std::cout << " Numerical issue with solving linear system..."
					<< std::endl;
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
