#include <math_methods.h>
#include <cstdlib>

double Math_methods::RandomDouble( const double &min, const double&max )
{
	double f = (double)rand()/RAND_MAX;
	return min + f*(max - min);
}

bool Math_methods::quadratic_solver_loqo( const Matrix <double, Dynamic, Dynamic> &H, const Matrix <double, Dynamic, Dynamic> &A, const Matrix <double, Dynamic, 1> &b, const Matrix <double, Dynamic, 1> &r, Matrix <double, Dynamic, 1> &fvalues )
{
	int n = (int)H.rows();

	Matrix<double,Dynamic,Dynamic> KKT(2*n,2*n);
	VectorXd c(n);
	c.setZero();
	VectorXd rhs(2*n);
	rhs << c, b;
	std::cout<<" rhs:\n"<< rhs << std::endl;
	
	KKT << -(H + MatrixXd::Identity(n,n)), A.transpose(), A, MatrixXd::Identity(n,n);
	std::cout<<" Initial KKT matrix:\n"<< KKT << std::endl;

	VectorXd soln(2*n);
	soln = KKT.partialPivLu().solve(rhs);
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

	double mu = (z.dot(g) + v.dot(w) + s.dot(t) + p.dot(q))/4*n; 

 	// construct G Z, V W, S T, and P Q diagonal matrices
 	MatrixXd G(n);
	G = g.asDiagonal();
	MatrixXd Z(n);
	Z = z.asDiagonal();
 	MatrixXd V(n);
	V = v.asDiagonal();
 	MatrixXd W(n);
	W = w.asDiagonal();
 	MatrixXd S(n);
	S = s.asDiagonal();
 	MatrixXd T(n);
	T = t.asDiagonal();
 	MatrixXd P(n);
	P = p.asDiagonal();
 	MatrixXd Q(n);
	Q = q.asDiagonal();


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

	bool converged = false;
	while (!converged)
	{
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

		soln.partialPivLu().solve(rhs);

		dx = soln.segment(0,n-1);
		dy = soln.segment(n,2*n - 1); 

	}

	return true;
}
