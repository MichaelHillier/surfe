#include <math_methods.h>
#include <cstdlib>

double Math_methods::_find_step( const VectorXd &da, const VectorXd &a )
{
	int n = da.rows();
	double max = DBL_MIN;
	for (int j = 0; j < n; j++ ){
		double alpha = -da[j]/(a[j]*0.95);
		if ( alpha > max && alpha < 1) max = alpha;
	}
	return max;
}

double Math_methods::_find_positivity_step( const VectorXd &da, const VectorXd &a, const VectorXd &db, const VectorXd &b, const VectorXd &dc, const VectorXd &c, const VectorXd &dd, const VectorXd &d )
{
	double max = DBL_MIN;
	
	double max_a = _find_step(da,a);
	double max_b = _find_step(db,b);
	double max_c = _find_step(dc,c);
	double max_d = _find_step(dd,d);

	if ( max_a > max ) max = max_a;
	if ( max_b > max ) max = max_b;
	if ( max_c > max ) max = max_c;
	if ( max_d > max ) max = max_d;

	return max;
}

double Math_methods::RandomDouble( const double &min, const double&max )
{
	double f = (double)rand()/RAND_MAX;
	return min + f*(max - min);
}

bool Math_methods::quadratic_solver_loqo( const MatrixXd &H, const MatrixXd &A, const VectorXd &b, const VectorXd &r, VectorXd &fvalues )
{
	int n = (int)H.rows();

	MatrixXd KKT(2*n,2*n);
	VectorXd c(n);
	c.setZero();
	VectorXd rhs(2*n);
	rhs << c, b;
	std::cout<<" rhs:\n"<< rhs << std::endl;
	
	KKT << -(H + MatrixXd::Identity(n,n)), A.transpose(), A, MatrixXd::Identity(n,n);
	std::cout<<" Initial KKT matrix:\n"<< KKT << std::endl;

	VectorXd soln(2*n);
	soln = KKT.partialPivLu().solve(rhs);
	std::cout<<" Soln from Initial KKT matrix:\n"<< soln << std::endl;
	VectorXd x(n);
	x = soln.segment(0,n-1);
	VectorXd y(n);
	y = soln.segment(n,2*n - 1); 
    VectorXd g(n);
	VectorXd z(n);
	VectorXd t(n);
	VectorXd s(n);
	VectorXd v(n);
	VectorXd w(n);
	VectorXd p(n);
	VectorXd q(n);
	for (int j = 0; j < n; j++ ){
		g(j) = std::max(abs(x(j)),100.0);
		z(j) = g(j);
		t(j) = g(j);
		s(j) = g(j);
		v(j) = std::max(abs(y(j)),100.0);
		w(j) = v(j);
		p(j) = std::max(abs(r(j) - w(j)),100.0);
		q(j) = v(j);
	}
	MatrixXd DebugMatrixV(n,10);
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
	std::cout<<" Current Variable matrix:\n"<< DebugMatrixV << std::endl;

	double mu = (z.dot(g) + v.dot(w) + s.dot(t) + p.dot(q))/4*n; 

 	// construct G Z, V W, S T, and P Q diagonal matrices
    MatrixXd G(n,n);
	MatrixXd Z(n,n);
 	MatrixXd V(n,n);
 	MatrixXd W(n,n);
 	MatrixXd S(n,n);
 	MatrixXd T(n,n);
 	MatrixXd P(n,n);
 	MatrixXd Q(n,n);
	Q = q.asDiagonal();
	// below matrices are computed after predictor step
	MatrixXd dG(n,n);
	MatrixXd dV(n,n);
	MatrixXd dT(n,n);
	MatrixXd dP(n,n);


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
	MatrixXd DebugMatrixStepV(n,10);

	bool converged = false;
	while (!converged)
	{
		double primal_obj = 0.5*x.transpose()*H*x;
		double   dual_obj = b.dot(y) - 0.5*x.transpose()*H*x - r.dot(q);

		double sigfig = std::max(-std::log10(abs(primal_obj - dual_obj)/(abs(primal_obj) + 1)),0.0);

		if ( sigfig > 6 )
		{
			converged = true;
			break;
		}
 
		G = g.asDiagonal();
		Z = z.asDiagonal();
		V = v.asDiagonal();
		W = w.asDiagonal();
		S = s.asDiagonal();
		T = t.asDiagonal();
		P = p.asDiagonal();
		Q = q.asDiagonal();

		MatrixXd D(n,n);
		D = (S.inverse()*T + G*Z.inverse()).inverse();
		MatrixXd E(n,n);
		E = (V*W.inverse() + P.inverse()*Q).inverse();

		rho     = b - A*x + w;
		nu      = -x + g - t;
		alpha   = r - w - p;
		sigma   = -A.transpose()*y - z + H*x;
		tau     = -z -s;
		beta    = y + q - v;
		// predictor nonlinearities
		gamma_z = -z;
		gamma_w = -w;
		gamma_s = -s;
		gamma_q = -q;

		tauh = tau - gamma_s;
		betah = beta - V*W.inverse()*gamma_w;
		alphah = alpha - P*Q.inverse()*gamma_q;
		nuh = nu + G*Z.inverse()*gamma_z;

		rhs << (sigma - D*(nuh + S.inverse()*T*tauh)),(rho - E*(betah - P.inverse()*Q*alphah));

		KKT << -(H + D), A.transpose(), A, E;

		soln = KKT.partialPivLu().solve(rhs);

		// get "delta" variables for predictor system ...
		dx = soln.segment(0,n-1);
		dy = soln.segment(n,2*n - 1); 
		dw = -E*(betah -P.inverse()*Q*alphah +dy);
		dt = -D*S.inverse()*T*(G*Z.inverse()*tauh - nuh + dx);
		dz = G.inverse()*Z*(nuh - dx - dt);
		dq = P.inverse()*Q*(dw - alphah);
		dv = V*W.inverse()*(gamma_w - dw);
		ds = gamma_s - S*T.inverse()*dt;
		dp = P*Q.inverse()*(gamma_q - dq);
		dg = G*Z.inverse()*(gamma_z - dz);

		// Debug
		DebugMatrixStepV.col(0) = dx;
		DebugMatrixStepV.col(1) = dy;
		DebugMatrixStepV.col(2) = dg;
		DebugMatrixStepV.col(3) = dz;
		DebugMatrixStepV.col(4) = dt;
		DebugMatrixStepV.col(5) = ds;
		DebugMatrixStepV.col(6) = dv;
		DebugMatrixStepV.col(7) = dw;
		DebugMatrixStepV.col(8) = dp;
		DebugMatrixStepV.col(9) = dq;
		std::cout<<" Current Step Variable matrix (after predictor step):\n"<< DebugMatrixStepV << std::endl;

		// compute step lengths for primal and dual systems
		double alpha_p = _find_positivity_step(dg,g,dw,w,dt,t,dp,p);
		double alpha_d = _find_positivity_step(dz,z,dv,v,ds,s,dq,q);
		// alpha_pd will be the maximum of the two above values
		double alpha_pd = std::max(alpha_p,alpha_d);
		double fraction = pow(((alpha_pd - 1.0)/(alpha_pd + 10.0)),2);

		// update mu
		mu = (z.dot(g) + v.dot(w) + s.dot(t) + p.dot(q))*(fraction)/4*n;

		// update rhs variables rho,nu,alpha,sigma,tau,beta,gamma's
		// first compute dG,dV,dT,and dP matrices from dg,dv,dt, and dp vectors above
		dG = dg.asDiagonal();
		dV = dv.asDiagonal();
		dT = dt.asDiagonal();
		dP = dp.asDiagonal();
		gamma_z = mu*G.inverse()*ev - z - G.inverse()*dG*dz;
		gamma_w = mu*V.inverse()*ev - w - V.inverse()*dV*dw;
		gamma_s = mu*T.inverse()*ev - s - T.inverse()*dT*ds;
		gamma_q = mu*P.inverse()*ev - q - P.inverse()*dP*dq;
		tauh = tau - gamma_s;
		betah = beta - V*W.inverse()*gamma_w;
		alphah = alpha - P*Q.inverse()*gamma_q;
		nuh = nu + G*Z.inverse()*gamma_z;

		rhs << (sigma - D*(nuh + S.inverse()*T*tauh)),(rho - E*(betah - P.inverse()*Q*alphah));

		soln = KKT.partialPivLu().solve(rhs);

		// get "delta" variables for corrector system ...
		dx = soln.segment(0,n-1);
		dy = soln.segment(n,2*n - 1); 
		dw = -E*(betah -P.inverse()*Q*alphah +dy);
		dt = -D*S.inverse()*T*(G*Z.inverse()*tauh - nuh + dx);
		dz = G.inverse()*Z*(nuh - dx - dt);
		dq = P.inverse()*Q*(dw - alphah);
		dv = V*W.inverse()*(gamma_w - dw);
		ds = gamma_s - S*T.inverse()*dt;
		dp = P*Q.inverse()*(gamma_q - dq);
		dg = G*Z.inverse()*(gamma_z - dz);

		// Debug
		DebugMatrixStepV.col(0) = dx;
		DebugMatrixStepV.col(1) = dy;
		DebugMatrixStepV.col(2) = dg;
		DebugMatrixStepV.col(3) = dz;
		DebugMatrixStepV.col(4) = dt;
		DebugMatrixStepV.col(5) = ds;
		DebugMatrixStepV.col(6) = dv;
		DebugMatrixStepV.col(7) = dw;
		DebugMatrixStepV.col(8) = dp;
		DebugMatrixStepV.col(9) = dq;
		std::cout<<" Current Step Variable matrix (after corrector step):\n"<< DebugMatrixStepV << std::endl;

		// compute step lengths for primal and dual systems
		alpha_p = _find_positivity_step(dg,g,dw,w,dt,t,dp,p);
		alpha_d = _find_positivity_step(dz,z,dv,v,ds,s,dq,q);

		// update solution
		x += (1/alpha_p)*dx;
		g += (1/alpha_p)*dg;
		w += (1/alpha_p)*dw;
		t += (1/alpha_p)*dt;
		p += (1/alpha_p)*dp;

		y += (1/alpha_d)*dy;
		z += (1/alpha_d)*dz;
		v += (1/alpha_d)*dv;
		s += (1/alpha_d)*ds;
		q += (1/alpha_d)*dq;
		
	}

	return true;
}
