#include <anisotropy_input.h>

bool Orientation::_compute_strike_dip_polarity_from_normal()
{
	// get dip first 
	_dip = acos(_normal[2])*R2D; // could do better. e.g. for overturn cases puts _dip > 90. but there formula's for getting normals works so sticking with it for now
	// get dip_direction
	double dip_direction = atan2(_normal[1], _normal[0])*R2D;


	// if negative azimuth get positive angle
	if (dip_direction < 0) dip_direction += 360;

	// get strike
	_strike = 360 - dip_direction;

	return true; // check this computation
}

bool Orientation::_compute_normal_from_strike_dip_polarity()
{
	// Get down dip vector - v
	double vx = cos(-1.0*(_strike * D2R)) * cos(-1.0*(_dip * D2R));
	double vy = sin(-1.0*(_strike * D2R)) * cos(-1.0*(_dip * D2R));
	double vz = sin(-1.0*(_dip * D2R));
	double l_dv = sqrt(vx*vx + vy*vy + vz*vz);
	vx /= l_dv;
	vy /= l_dv;
	vz /= l_dv;
	_dipvector[0] = vx;
	_dipvector[1] = vy;
	_dipvector[2] = vz;
	
	// Get strike vector - vp
	double vpx = -1.0 * vy;
	double vpy = vx;
	double vpz = 0;
	double l_sv = sqrt(vpx*vpx + vpy*vpy + vpz*vpz);
	vpx /= l_sv;
	vpy /= l_sv;
	vpz /= l_sv;
	_strikevector[0] = vpx;
	_strikevector[1] = vpy;
	_strikevector[2] = vpz;

	// Compute cross product of v x vp
	double Nx = (vy * vpz) - (vz * vpy);
	double Ny = (vz * vpx) - (vx * vpz);
	double Nz = (vx * vpy) - (vy * vpx);

	// Compute Length
	double Length = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);

	// Normalize
	Nx = Nx / Length;
	Ny = Ny / Length;
	Nz = Nz / Length;

	// Check polarity of computed normal is consistent with the property "polarity"
	// polarity == 1 -> Is overturned: Must point downward
	// polarity == 0 -> Is upright: Must point upward
	// polarity != 0/1 -> Is unknown: Can point upward or downward, it doesn't matter.
	//                    Retains its initialized polarity set above.
	int polarity = 0; // assume upright .. in the end it doesn't matter for this type of data
	if ((polarity == 1 && Nz > 0) || (polarity == 0 && Nz < 0))
	{
		// Flip vector
		Nx = -Nx;
		Ny = -Ny;
		Nz = -Nz;
	}

	// assign normal
	_normal[0] = Nx;
	_normal[1] = Ny;
	_normal[2] = Nz;

	return true;
}


void Orientation::_set_eigensystem()
{
	eigenvectors[0][0] = _strikevector[0];
	for (int j = 0; j < 3; j++){
		eigenvectors[j][0] = _strikevector[j];
		eigenvectors[j][1] = _dipvector[j];
		eigenvectors[j][2] = _normal[j];
	}
	eigenvalues[0] = 1;
	eigenvalues[1] = 1;
	eigenvalues[2] = 10;
}

void local_anisotropy_input::setPCA()
{
	unsigned int n = (unsigned int)orientation->size();
	for (int j = 0; j < n; j++){

	}
}

void local_anisotropy_input::getMeanTensor(Matrix3d &u_mean, Matrix3d &d_mean)
{
	Orientation cur_pt = orientation->at(0);

	Matrix3d sigma_t_inv_sqrt;
	Matrix3d sigma_t_sqrt;

	Matrix3d sigma_t_inv_sqrt_temp;
	Matrix3d sigma_t_sqrt_temp;
	Matrix3d U;
	Matrix3d UT;
	Matrix3d Diag_inv_sqrt;
	Matrix3d Diag_sqrt;
	Matrix3d Uj;
	Matrix3d UTj;
	Matrix3d Diagj;
	Matrix3d temp_mat;
	Matrix3d temp_mat2;

	bool converged = false;
	int niter = 0;
	double last_avg_err = 0;
	double last_avg = 0;
	double last_fin_mat_add = 0;
	while (!converged)
	{
		niter++;
		U = cur_pt.eigenvectors;
		UT = cur_pt.eigenvalues.transpose();
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				if (j == k)
				{
					Diag_sqrt(j,k) = sqrt(cur_pt.eigenvalues[j]);
					Diag_inv_sqrt(j,k) = 1 / sqrt(cur_pt.eigenvalues[j]);
				}
				else
				{
					Diag_sqrt(j, k) = 0.0;
					Diag_inv_sqrt(j,k) = 0.0;
				}
			}
		}

		sigma_t_sqrt = U*Diag_sqrt*UT;
		sigma_t_inv_sqrt = U*Diag_inv_sqrt*UT;

		Matrix3d cur_matrix;
		for (int j = 0; j < orientation->size(); j++){
			Uj = orientation->at(j).eigenvectors;
			UTj = orientation->at(j).eigenvectors.transpose();
			for (int k = 0; k < 3; k++){
				for (int l = 0; l < 3; l++){
					if (k == l) Diagj(k,l) = orientation->at(j).eigenvalues[k];
					else Diagj(k,l) = 0;
				}
			}
			temp_mat2 = sigma_t_inv_sqrt*Uj*Diagj*UTj*sigma_t_inv_sqrt;
			SelfAdjointEigenSolver<Matrix3d> es;
			es.compute(temp_mat2);
			Matrix3d U_j;
			Vector3d d_j;
			U_j = es.eigenvectors().real();
			d_j = es.eigenvalues().real();
			Matrix3d U_j_T;
			Matrix3d D_j;
			Matrix3d D_j_log;
			U_j_T = U_j.transpose();
			for (int k = 0; k < 3; k++){
				for (int l = 0; l < 3; l++){
					if (k == l)
					{
						D_j(k,l) = d_j[k];
						D_j_log(k, l) = log(d_j[k]);
					}
					else
					{
						D_j(k, l) = 0;
						D_j_log(k, l) = 0;
					}
				}
			}
			temp_mat = U_j*D_j_log*U_j_T;
			cur_matrix += temp_mat;
		}
		Matrix3d fin_mat;
		fin_mat = (1.0 / (double)orientation->size())*cur_matrix;
		double fin_mat_add = 0.0;
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				fin_mat_add += abs(fin_mat(j,k));
			}
		}

}
