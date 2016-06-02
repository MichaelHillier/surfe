#ifndef math_methods_h
#define math_methods_h

#include <mpirxx.h>

#include <vector>
#include <algorithm>

#define SQR(x)   ((x)*(x)) // x^2

class Math_methods{
private:
	template <class T> static T _find_step_length(const std::vector < T > &a,const std::vector < T > &da,
													const std::vector < T > &b,const std::vector < T > &db);
	template <class T> static void _rot(std::vector< std::vector < T > > &a, const T &s, const T &tau, const int &i, const int &j, const int &k, const int &l);
	// convenience functions for debugging ...
	static void _output_matrix(const std::vector < std::vector < double > > & matrix, const std::string &filename);
	static void _output_matrix(const std::vector < std::vector < mpf_class > > & matrix, const std::string &filename);
	static void _output_vector(const std::vector<mpf_class > &vector, const std::string &filename);
	static void _output_vector(const std::vector<double > &vector, const std::string &filename);
	static double _get_double(const double &d) { return d;}
	static double _get_double(const mpf_class &d) { return d.get_d(); }
public:
	template <class T> static bool sort_vector_w_index(std::vector<T> &arr, std::vector<int> &brr);
	template <class T> static bool eigenanalysis( std::vector < std::vector <T> > &A, std::vector < std::vector < T > > &Q, std::vector< T > &w);
	template <class T> static bool ludcmp(std::vector< std::vector < T > > &a, std::vector<int> &indx, T &d); 
	template <class T> static void lubksb(std::vector< std::vector < T > > &a, std::vector<int> &indx, std::vector<T> &b);
	template <class T> static bool solve_sym_linear_system_via_LU_dcmp(std::vector< std::vector <T> > &a, std::vector <T> &b);
	template <class T> static std::vector< std::vector < T > > make_std_matrix(const int &rows, const int&cols);
	template <class T> static T find_max_element_in_matrix(const std::vector < std::vector < T > > &mat);
	template <class T> static bool quadratic_solver(const std::vector < std:: vector < T > > &H,
													const std::vector < std:: vector < T > > &A,
													const std::vector < std:: vector < T > > &C,
													const std::vector < T > &b,
													const std::vector < T > &d,
													std::vector <T> &fvalues);
	template <class T> static std::vector< std::vector < T > > matrix_transpose(const std::vector < std::vector < T > > &mat);
	template <class T> static bool matrix_vector_multiply(const std::vector < std::vector < T > > &A, const std::vector< T > &b, std::vector < T > &c);
	template <class T> static T max_element_wrt_zero(const T &a, const T &b);
	template <class T> static void SWAP(T &a, T &b){ T x = a; a = b; b = x; }
	template <class T> static bool angle_btw_2_vectors( const std::vector < T > &v1, const std::vector < T > &v2, T &angle);
	static double RandomDouble(const double &min, const double&max);
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
void Math_methods::_rot( std::vector< std::vector < T > > &a, const T &s, const T &tau, const int &i, const int &j, const int &k, const int &l )
{
	T g,h;

	g=a[i][j];
	h=a[k][l];
	a[i][j]=g-s*(h+g*tau); // if T is mpf_class type then explosion could happen here. Maybe do explicit cast from int to mpf_class ?
	a[k][l]=h+s*(g-h*tau);
}

template <class T>
bool Math_methods::eigenanalysis( std::vector < std::vector <T> > &A, std::vector < std::vector < T > > &Q, std::vector< T > &w )
{
	// computation via Jacobi method.

 	// check is matrix a is square
 	const int n = (int)A.size();
 	for (int j = 0; j < n; j++ ){
 		if ((int)A[j].size() != n ) return false;
 	}

	T tresh,theta,tau,t,sm,s,h,g,c;

	std::vector<T> b;
	b.resize(n);
	std::vector<T> z;
	z.resize(n);
	for (int ip=0;ip<n;ip++) {
		for (int iq=0;iq<n;iq++) Q[ip][iq]=0.0;
		Q[ip][ip]=1.0;
	}
	for (int ip=0;ip<n;ip++) {
		b[ip]=w[ip]=A[ip][ip];
		z[ip]=0.0;
	}
	int nrot=0;
	for (int i=1;i<=50;i++) {
		sm=0.0;
		for (int ip=0;ip<n-1;ip++) {
			for (int iq=ip+1;iq<n;iq++)
				sm += abs(A[ip][iq]);
		}
		if (sm == 0.0)
			return 0;
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (int ip=0; ip<n-1; ip++) {
			for (int iq=ip+1; iq<n; iq++) {
				g=100.0*abs(A[ip][iq]);
				if (i > 4 && (abs(w[ip])+g) == abs(w[ip])
					&& (abs(w[iq])+g) == abs(w[iq]))
					A[ip][iq]=0.0;
				else if (abs(A[ip][iq]) > tresh) {
					h=w[iq]-w[ip];
					if ((abs(h)+g) == abs(h))
						t=(A[ip][iq])/h;
					else {
						theta=0.5*h/(A[ip][iq]);
						t=1.0/(abs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*A[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					w[ip] -= h;
					w[iq] += h;
					A[ip][iq]=0.0;
					for (int j=0; j<ip; j++)
						_rot(A,s,tau,j,ip,j,iq);
					for (int j=ip+1; j<iq; j++)
						_rot(A,s,tau,ip,j,j,iq);
					for (int j=iq+1; j<n; j++)
						_rot(A,s,tau,ip,j,iq,j);
					for (int j=0; j<n; j++)
						_rot(Q,s,tau,j,ip,j,iq);
					++nrot;
				}
			}
		}
		for (int ip=0; ip<n; ip++) {
			b[ip] += z[ip];
			w[ip]=b[ip];
			z[ip]=0.0;
		}
	}

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
T Math_methods::find_max_element_in_matrix( const std::vector < std::vector < T > > &mat )
{
	if ((int)mat[0].size() == 0) throw -1;
	T max_element = mat[0][0];
	for (int j = 0; j < (int)mat.size(); j++ ){
		for (int k = 0; k < (int)mat[j].size(); k++ ){
			if (mat[j][k] > max_element) max_element = mat[j][k];
		}
	}
	return max_element;
}

template <class T>
T Math_methods::_find_step_length(const std::vector < T > &a,const std::vector < T > &da,
											const std::vector < T > &b,const std::vector < T > &db)
{
	int n = (int)a.size();
	if (a.size() != b.size() || a.size() != da.size() || a.size() != db.size() || da.size() != db.size()) throw -1;
	// find step length (a) ... this is the most important step
	// has to satisfy all the non-negativity conditions ...
	T min_alpha_a = 100.0; // improper initialization. proper initialization will occur in first iteration of below for loop 
	T min_alpha_b = 100.0; // ?????		why do i do this? seems dumb
	T max_alpha_a = 0.0; // improper initialization. proper initialization will occur in first iteration of below for loop 
	T max_alpha_b = 0.0;
	T alpha = 0.0;
	std::vector < T > alpha_a;
	std::vector < T > alpha_b;
	for (int j = 0; j < n; j++){
		//cout<<" z["<<j<<"]= "<<mpf_get_d(z[j])<<endl;
		//cout<<" s["<<j<<"]= "<<mpf_get_d(s[j])<<endl;
		//cout<<" ["<<j<<"] dz_aff = "<<mpf_get_d(dz_aff[j])<<" ds_aff= "<<mpf_get_d(ds_aff[j])<<endl;
		// do alpha_b first ...
		if (b[j] > 0.0) alpha_b.push_back( b[j] / db[j] );
		else alpha_b.push_back(-1.0*b[j] / db[j]);
		// alpha_a
		if (a[j] > 0.0) alpha_a.push_back( a[j] / da[j] );
		else alpha_a.push_back( -1.0*a[j] / da[j] );

		//if (mpf_cmp(alpha_s[j],max_alpha_s) > 0 && mpf_cmp_d(alpha_s[j],1.0) <= 0) mpf_set(max_alpha_s,alpha_s[j]);
		//if (mpf_cmp(alpha_z[j],max_alpha_z) > 0 && mpf_cmp_d(alpha_z[j],1.0) <= 0) mpf_set(max_alpha_z,alpha_z[j]);
		if (alpha_b[j] < min_alpha_b && alpha_b[j] > 0.00000000000001) min_alpha_b = alpha_b[j];
		if (alpha_a[j] < min_alpha_a && alpha_a[j] > 0.00000000000001) min_alpha_a = alpha_a[j];
		// 			if (mpf_cmp(alpha_s[j],min_alpha_s) < 0) mpf_set(min_alpha_s,alpha_s[j]);
		// 			if (mpf_cmp(alpha_z[j],min_alpha_z) < 0) mpf_set(min_alpha_z,alpha_z[j]);
		//cout<<" alpha_s["<<j<<"]= "<<mpf_get_d(alpha_s[j])<<endl;
		//cout<<" alpha_z["<<j<<"]= "<<mpf_get_d(alpha_z[j])<<endl;
		// 			cout<<" max_alpha_s = "<<mpf_get_d(max_alpha_s)<<endl;
		// 			cout<<" max_alpha_z = "<<mpf_get_d(max_alpha_z)<<endl;
		//cout<<" min_alpha_s = "<<mpf_get_d(min_alpha_s)<<endl;
		//cout<<" min_alpha_z = "<<mpf_get_d(min_alpha_z)<<endl;
	}

	if (min_alpha_b < min_alpha_a) alpha = min_alpha_b;
	else alpha = min_alpha_a;

	if (alpha > 1.0 || alpha == 0) alpha = 1.0;
	return alpha;
}

template <class T>
bool Math_methods::quadratic_solver( const std::vector < std:: vector < T > > &H, const std::vector < std:: vector < T > > &A, const std::vector < std:: vector < T > > &C, const std::vector < T > &b, const std::vector < T > &d, std::vector <T> &fvalues )
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

int n = (int)H.size();
int na = (int)A.size();
int nc = (int)C.size();

///////////////////////////////////////////////////////////////////////
/////////////////////////// Initialization ////////////////////////////
///////////////////////////////////////////////////////////////////////
// Calculation Helper Variables
T add,temp,temp2,temp_p,minus_one,elemsum,datanorm,zTs,mu,zero,prev_mu;
zero = 0.0;
elemsum = 0.0;
// KKT System
std::vector<T> x,y,z,s,rh,ra,rc,rsz;
x.resize(n,0);
y.resize(na,0);
z.resize(nc,0);
s.resize(nc,0);
rh.resize(n,0);
ra.resize(na,0);
rc.resize(nc,0);
rsz.resize(nc,0);
std::vector<T> t1,t2,t3,t4,t5;
t1.resize(n,0);
t2.resize(n,0);
t3.resize(n,0);
t4.resize(na,0);
t5.resize(nc,0);
std::vector< std::vector < T > > AT = Math_methods::matrix_transpose(A);
std::vector< std::vector < T > > CT = Math_methods::matrix_transpose(C);
std::vector< std::vector < T > > KKT = Math_methods::make_std_matrix<T>(n + na + 2*nc,n + na + 2*nc);
std::vector< std::vector < T > > KKT_predictor = Math_methods::make_std_matrix<T>(n + na + 2*nc,n + na + 2*nc);
std::vector< std::vector < T > > KKT_corrector = Math_methods::make_std_matrix<T>(n + na + 2*nc,n + na + 2*nc);
std::vector< T > solution_vector;
solution_vector.resize(n + na + 2*nc,0);

// step containers and related helper variables
T alpha,min_alpha_z,min_alpha_s,max_alpha_z,max_alpha_s,mu_aff,sigma;
std::vector < T > alpha_z,alpha_s,dx_aff,dy_aff,dz_aff,ds_aff,rsz_aff,dvector;
alpha_z.resize(nc,0);
alpha_s.resize(nc,0);
dx_aff.resize(n,0);
dy_aff.resize(na,0);
dz_aff.resize(nc,0);
ds_aff.resize(nc,0);
rsz_aff.resize(nc,0);
dvector.resize(n + na + 2*nc,0);
std::vector < T > dx,dy,dz,ds,dvector_corr;
dx.resize(n,0);
dy.resize(na,0);
dz.resize(nc,0);
ds.resize(nc,0);
dvector_corr.resize(n + na + 2*nc,0);

///////////////////////////////////////////////////////////////////////
////////////////////// End of Initialization //////////////////////////
///////////////////////////////////////////////////////////////////////

// get norm of matrix 
datanorm = sqrt(Math_methods::find_max_element_in_matrix(H));

// setup initial iterate for z and s...
for (int j = 0; j < nc; j++ ){
	z[j] = datanorm;
	s[j] = datanorm;
}

// Build KKT matrix blocks [1][1], [1][2], [1][3], [1][4]
for (int j = 0; j < n; j++){ // Matrix rows...
	// Matrix columns...
	// [1][1] Block
	for (int k = 0; k < n; k++) KKT[j][k] = H[j][k];
	// [1][2] Block
	for (int k = 0; k < na; k++) KKT[j][n + k] = -1.0*AT[j][k];
	// [1][3] Block
	for (int k = 0; k < nc; k++) KKT[j][n + na + k] = -1.0*CT[j][k]; 
	// [1][4] Block
	for (int k = 0; k < nc; k++) KKT[j][n + na + nc + k] = 0.0;
}
// Build KKT matrix blocks [2][1], [2][2], [2][3], [2][4]
for (int j = 0; j < na; j++){
	// [2][1] Block
	for (int k = 0; k < n; k++)  KKT[n + j][k] = A[j][k];
	// [2][2] Block
	for (int k = 0; k < na; k++) KKT[n + j][n + k] = 0.0;
	// [2][3] Block
	for (int k = 0; k < nc; k++) KKT[n + j][n + na + k] = 0.0;
	// [2][4] Block
	for (int k = 0; k < nc; k++) KKT[n + j][n + na + nc + k] = 0.0;
}
// Build KKT matrix blocks [3][1], [3][2], [3][3], [3][4]
for (int j = 0; j < nc; j++){
	// [3][1] Block
	for (int k = 0; k < n; k++)  KKT[n + na + j][k] = C[j][k];
	// [3][2] Block
	for (int k = 0; k < na; k++) KKT[n + na + j][n + k] = 0.0;
	// [3][3] Block
	for (int k = 0; k < nc; k++) KKT[n + na + j][n + na + k] = 0.0;
	// [3][4] Block
	for (int k = 0; k < nc; k++){
		if (j == k) KKT[n + na + j][n + na + nc + k] = -1.0;
		else KKT[n + na + j][n + na + nc + k] = 0.0;
	}
}
// Build KKT matrix blocks [4][1], [4][2]
for (int j = 0; j < nc; j++){
	// [4][1] Block
	for (int k = 0; k < n; k++)  KKT[n + na + nc + j][k] = 0.0;
	// [4][2] Block
	for (int k = 0; k < na; k++) KKT[n + na + nc + j][n + k] = 0.0;
}

// Build KKT matrix blocks: [4][3], [4][4]
for (int j = 0; j < nc; j++){
	// [4][3] Block
	for (int k = 0; k < nc; k++){
		if (j == k) KKT[n + na + nc + j][n + na + k] = s[j];
		else KKT[n + na + nc + j][n + na + k] = 0.0;
	}
	// [4][4] Block
	for (int k = 0; k < nc; k++){
		if (j == k) KKT[n + na + nc + j][n + na + nc + k] = z[j];
		else KKT[n + na + nc + j][n + na + nc + k] = 0.0;
	}
}

//////////////////////////////////////////////
// Determine Starting Point for the iterate //
//////////////////////////////////////////////
// 1st compute residuals for the affine system
Math_methods::matrix_vector_multiply(H,x,t1);
Math_methods::matrix_vector_multiply(AT,y,t2);
Math_methods::matrix_vector_multiply(CT,z,t3);
// 	_output_vector(t3,"t3.txt"); // debug
// 	_output_vector(z,"z.txt");
// 	_output_matrix(CT,"CT.txt");
for (int j = 0; j < n; j ++ ){ // rh residual = H*x - AT*y - CT*z
	rh[j] = t1[j] - t2[j] - t3[j];
	solution_vector[j] = -1.0*rh[j];
}
Math_methods::matrix_vector_multiply(A,x,t4);
for (int j = 0; j < na; j++ ){ // ra residual = A*x - b
	ra[j] = t4[j] - b[j];
	solution_vector[n + j] = -1.0*ra[j];
}
zTs = 0.0; // re-initialize this adder
Math_methods::matrix_vector_multiply(C,x,t5);
for (int j = 0; j < nc; j++ ){ // rc residual = C*x - s - d
	rc[j] = t5[j] - s[j] - d[j];
	rsz[j] = s[j]*z[j]; // rsz residual = Z[]*S[]*e
	solution_vector[n + na + j] = -1.0*rc[j];
	zTs += rsz[j];
	solution_vector[n + na + nc + j] = -1.0*rsz[j];
}
// Solve Affine system ...
for (int j = 0; j < n + na + 2*nc; j++ ) dvector[j] = solution_vector[j];
for (int j = 0; j < n + na + 2*nc; j++ ){
	for (int k = 0; k < n + na + 2*nc; k++ ) KKT_predictor[j][k] = KKT[j][k];
}
// 	_output_matrix(KKT_predictor,"KKT_predictor.txt"); // debug
// 	_output_vector(dvector,"dvector_b4.txt");
Math_methods::solve_sym_linear_system_via_LU_dcmp(KKT_predictor,dvector);
//	_output_vector(dvector,"dvector_af.txt"); // debug

for (int j = 0; j < n; j++ ){
	dx[j] = dvector[j];
	if (j < na) dy[j] = dvector[j + n];
	if (j < nc)
	{
		dz[j] = dvector[j + n + na];
		ds[j] = dvector[j + n + na + nc];
	}
}
// Update iterate using full affine scaling
// (x,y,z,s)->(x,y,z,s) + (dx_aff,dy_aff,dz_aff,ds_aff)
for (int j = 0; j < n;  j++ ) x[j] += dx[j];
for (int j = 0; j < na; j++ ) y[j] += dy[j];
for (int j = 0; j < nc; j++ ) z[j] += dz[j];
for (int j = 0; j < nc; j++ ) s[j] += ds[j];

// above iterate likely infeasible
// measure violation
std::vector < T > max_violation_list;
for (int j = 0; j < nc; j++ ){
	T v1 = -1.0*z[j];
	T v2 = -1.0*s[j];
	T maximum = max_element_wrt_zero(v1,v2);
	if ( maximum < 0) maximum = 0.0;
	max_violation_list.push_back( maximum );
}
// sort max_violation_list
std::sort(max_violation_list.begin(),max_violation_list.end());
T max_violation = max_violation_list[nc - 1];
// compute the iterate shift for complementary variables z & s
T shift = 1000.0 + 2.0*max_violation;
// update complementary variables ...
for (int j = 0; j < nc; j++ ){
	z[j] += shift;
	s[j] += shift;
}

bool found_soln = false;
int iter = 0;
while (!found_soln)
{
	// Update KKT matrix blocks: [4][3], [4][4]
	for (int j = 0; j < nc; j++){
		// [4][3] Block
		for (int k = 0; k < nc; k++){
			if (j == k) KKT[n + na + nc + j][n + na + k] = s[j];
			else KKT[n + na + nc + j][n + na + k] = 0.0;
		}
		// [4][4] Block
		for (int k = 0; k < nc; k++){
			if (j == k) KKT[n + na + nc + j][n + na + nc + k] = z[j];
			else KKT[n + na + nc + j][n + na + nc + k] = 0.0;
		}
	}

	//_output_matrix(KKT,"KKT_matrix.txt"); // debug

	// COMPUTE RESIDUALS ...
	Math_methods::matrix_vector_multiply(H,x,t1);
	Math_methods::matrix_vector_multiply(AT,y,t2);
	Math_methods::matrix_vector_multiply(CT,z,t3);
	for (int j = 0; j < n; j ++ ){ // rh residual = H*x - AT*y - CT*z
		rh[j] = t1[j] - t2[j] - t3[j];
		solution_vector[j] = -1.0*rh[j];
	}
	Math_methods::matrix_vector_multiply(A,x,t4);
	for (int j = 0; j < na; j++ ){ // ra residual = A*x - b
		ra[j] = t4[j] - b[j];
		solution_vector[n + j] = -1.0*ra[j];
	}
	zTs = 0.0; // re-initialize this adder
	Math_methods::matrix_vector_multiply(C,x,t5);
	for (int j = 0; j < nc; j++ ){ // rc residual = C*x - s - d
		rc[j] = t5[j] - s[j] - d[j];
		rsz[j] = s[j]*z[j]; // rsz residual = Z[]*S[]*e
		solution_vector[n + na + j] = -1.0*rc[j];
		zTs += rsz[j];
		solution_vector[n + na + nc + j] = -1.0*rsz[j];
	}
	//calculate mu
	mu = zTs / nc;
	//double mu_d = _get_double(mu); // debug
	if (iter > 5 && mu > prev_mu) return false;
	prev_mu = mu;
	if (mu < 0.00000001)
	{
		found_soln = true;
		for (int j = 0; j < n; j++) fvalues[j] = x[j]; // get solution
		continue;
	}
	///////////////////////////////////////////////////////
	/////////////////// Predictor Step ////////////////////
	///////////////////////////////////////////////////////
	// Solve Affine system ...
	//cout<<" solution vector for affine system "<<endl;
	for (int j = 0; j < n + na + 2*nc; j++ ) dvector[j] = solution_vector[j];
	for (int j = 0; j < n + na + 2*nc; j++ ){
		for (int k = 0; k < n + na + 2*nc; k++ ) KKT_predictor[j][k] = KKT[j][k];
	}
	Math_methods::solve_sym_linear_system_via_LU_dcmp(KKT_predictor,dvector);

	// get step vectors...
	for (int j = 0; j < nc; j++ ){
		dz_aff[j] = dvector[j + n + na];
		ds_aff[j] = dvector[j + n + na + nc];
	}
	
	alpha = _find_step_length(s,ds_aff,z,dz_aff);
	
	// calculate mu_aff
	elemsum = 0.0;
	for (int j = 0; j < nc; j++ ) elemsum += (z[j] + alpha*dz_aff[j]) * (s[j] + alpha*ds_aff[j]); 
	mu_aff = elemsum / nc;
	sigma = (mu_aff/mu)*(mu_aff/mu)*(mu_aff/mu);

	///////////////////////////////////////////////////////
	/////////////////// Corrector Step ////////////////////
	///////////////////////////////////////////////////////
	// modify rsz vector using the calculated correctors ...
	// recompute rsz residuals...
	// Z*S*e - sigma*mu*e + dz_aff*ds_aff*e
	for (int j = 0; j < nc; j ++){
		rsz_aff[j] = rsz[j] - sigma*mu + dz_aff[j]*ds_aff[j];
		solution_vector[n + na + nc + j] = -1.0*rsz_aff[j];
	}
	for (int j = 0; j < n + na + 2*nc; j++ ) dvector_corr[j] = solution_vector[j];

	// Solve the corrected linear system ...
	for (int j = 0; j < n + na + 2*nc; j++ ){
		for (int k = 0; k < n + na + 2*nc; k++ ) KKT_corrector[j][k] = KKT[j][k];
	}
	Math_methods::solve_sym_linear_system_via_LU_dcmp(KKT_corrector,dvector_corr);

	// get step vectors from corrector step...
	for (int j = 0; j < n; j++ ){
		dx[j] = dvector_corr[j];
		if (j < na) dy[j] = dvector_corr[j + n];
		if (j < nc)
		{
			dz[j] = dvector_corr[j + n + na];
			ds[j] = dvector_corr[j + n + na + nc];
		}
	}

	alpha = _find_step_length(s,ds,z,dz);

	// update x,y,z,s vectors using alpha step
	for (int j = 0; j < n;  j++ ) x[j] += alpha*dx[j];
	for (int j = 0; j < na; j++ ) y[j] += alpha*dy[j];
	for (int j = 0; j < nc; j++ ) z[j] += alpha*dz[j];
	for (int j = 0; j < nc; j++ ) s[j] += alpha*ds[j];
	iter++;
}
return true;
}

template <class T>
bool Math_methods::matrix_vector_multiply( const std::vector < std::vector < T > > &A, const std::vector< T > &b, std::vector < T > &c )
{
	int nrows = (int)A.size();
	int ncols = (int)A[0].size();
	int nvrows= (int)b.size();
	if (nvrows != ncols) return false;

	for (int j=0;j<nrows;j++){
		T elemsum = 0.0;
		for (int k=0;k<ncols;k++){
			elemsum+=A[j][k]*b[k];
		}
		c[j]=elemsum;
	}
	return true;
}

template <class T> std::vector< std::vector < T > >
Math_methods::matrix_transpose( const std::vector < std::vector < T > > &mat )
{
	int nrows = (int)mat[0].size();
	int ncols = (int)mat.size();
	std::vector < std::vector < T > > trans_mat = make_std_matrix<T>(nrows,ncols);
	for (int j = 0; j < nrows; j++ ){
		for (int k = 0; k < ncols; k++ ){
			trans_mat[j][k] = mat[k][j];
		}
	}
	return trans_mat;
}

template <class T>
std::vector< std::vector < T > > Math_methods::make_std_matrix( const int &nrows, const int&ncols )
{
	std::vector< std::vector < T > > matrix;
	matrix.resize(nrows);
	for (int j = 0; j < nrows; j++ ) matrix[j].resize(ncols);
	return matrix;
}

template <class T>
bool Math_methods::ludcmp( std::vector< std::vector < T > > &a, std::vector<int> &indx, T &d )
{
	const T TINY = 1.0e-20;
	int i, imax = 0, j, k;
	T big, dum, sum, temp;

	int n = (int)a.size();
	std::vector< T > vv;
	vv.resize(n);
	d = 1.0;
	for (i = 0;i<n;i++) {
		big = 0.0;
		for (j = 0;j<n;j++)
			if ((temp = abs(a[i][j])) > big) big = temp;
		if (big == 0.0) return false;
		vv[i] = 1.0 / big;
	}
	for (j = 0;j<n;j++) {
		for (i = 0;i<j;i++) {
			sum = a[i][j];
			for (k = 0;k<i;k++) sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
		}
		big = 0.0;
		for (i = j;i<n;i++) {
			sum = a[i][j];
			for (k = 0;k<j;k++) sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			if ((dum = vv[i] * abs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 0;k<n;k++) {
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			d = -d;
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a[j][j] == 0.0) a[j][j] = TINY;
		if (j != n - 1) {
			dum = 1.0 / (a[j][j]);
			for (i = j + 1;i<n;i++) a[i][j] *= dum;
		}
	}

	return true;
}

template <class T>
void Math_methods::lubksb( std::vector< std::vector < T > > &a, std::vector<int> &indx, std::vector<T> &b )
{
	int i, ii = 0, ip, j;
	T sum;

	int n = (int)a.size();
	for (i = 0;i<n;i++) {
		ip = indx[i];
		sum = b[ip];
		//cout<<" ii= "<<ii<<" ip= "<<ip<<" index["<<i<<"]= "<<indx[i]<<" sum = "<<sum<<endl; // debug
		b[ip] = b[i];
		if (ii != 0)
			for (j = ii - 1;j<i;j++) sum -= a[i][j] * b[j];
		else if (sum != 0.0)
			ii = i + 1;
		b[i] = sum;
	}
	for (i = n - 1;i >= 0;i--) {
		sum = b[i];
		//cout<<" Sum= "<<sum<<endl;
		for (j = i + 1;j<n;j++) sum -= a[i][j] * b[j];
		//cout<<" a[i][i]= "<<a[i][i]<<endl;  // debug
		b[i] = sum / a[i][i];
		//cout<<" b["<<i<<"]= "<<b[i]<<" sum= "<<sum<<" a[i][i]= "<<a[i][i]<<endl;  // debug
	}
}

template <class T>
bool Math_methods::solve_sym_linear_system_via_LU_dcmp( std::vector< std::vector <T> > &a, std::vector <T> &b )
{
	std::vector<int> indx;
	indx.resize((int)a.size(),0);
	T d;
	// perform LU decomposition
	if (!ludcmp<T>(a, indx, d)) return false;
	// perform LU back substitution to obtain solution
	lubksb<T>(a, indx, b);

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