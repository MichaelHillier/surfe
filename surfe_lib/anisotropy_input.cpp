#include <anisotropy_input.h>
#include <math_methods.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>

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

bool Orientation::_compute_a_normal_from_strike_dip()
{
	// Please note that the polarity of the normal here DOES NOT matter
	// Both the +ve/-ve are valid since these type of tensor only respect the axial nature of this 'direction'

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

	// assign normal
	_normal[0] = Nx;
	_normal[1] = Ny;
	_normal[2] = Nz;

	return true;
}


void Orientation::_set_eigensystem()
{

	for (int j = 0; j < 3; j++){
		eigenvectors(j, 0) = _strikevector[j];
		eigenvectors(j, 1) = _dipvector[j];
		eigenvectors(j, 2) = _normal[j];
	}
// 	std::cout << " strikevector:\n" << _strikevector << std::endl;
// 	std::cout << " dipvector:\n" << _dipvector << std::endl;
// 	std::cout << " normal vector:\n" << _normal << std::endl;
// 	std::cout << " Eigenvector matrix:\n" << eigenvectors << std::endl;
	eigenvalues[0] = 0.1;
	eigenvalues[1] = 0.5;
	eigenvalues[2] = 1.4;
	double det = eigenvectors.determinant();
	if (det < 0) // flip plunge vector
	{
		//std::cout << " Had to flip" << std::endl;
		eigenvectors(0, 0) *= -1.0;
		eigenvectors(1, 0) *= -1.0;
		eigenvectors(2, 0) *= -1.0;
	}
	U = eigenvectors;
	for (int j = 0; j < 3; j++){
		for (int k = 0; k < 3; k++){
			if (j == k)
			{
				D(j, k) = eigenvalues[j];
				S(j, k) = sqrt(eigenvalues[j] / eigenvalues[0]);
			}
			else
			{
				D(j, k) = 0.0;
				S(j, k) = 0.0;
			}
		}
	}
	Tensor = U * D * U.transpose();
	Transform = U * S * U.transpose();
}

bool Orientation::getDipVector(double(&vector)[3])
{
	if (_dipvector[0] >= -1.0 && _dipvector[0] <= 1.0) vector[0] = _dipvector[0];
	else return false;
	if (_dipvector[1] >= -1.0 && _dipvector[1] <= 1.0) vector[1] = _dipvector[1];
	else return false;
	if (_dipvector[2] >= -1.0 && _dipvector[2] <= 1.0) vector[2] = _dipvector[2];
	else return false;

	return true;
}

bool Orientation::getStrikeVector(double(&vector)[3])
{
	if (_strikevector[0] >= -1.0 && _strikevector[0] <= 1.0) vector[0] = _strikevector[0];
	else return false;
	if (_strikevector[1] >= -1.0 && _strikevector[1] <= 1.0) vector[1] = _strikevector[1];
	else return false;
	if (_strikevector[2] >= -1.0 && _strikevector[2] <= 1.0) vector[2] = _strikevector[2];
	else return false;

	return true;
}

void TensorInput::_sort_eigensystem(Matrix3d &evectors, Vector3d &evalues)
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

bool TensorInput::_get_neighbourhoods(const int &n_neighbors)
{
	if (orientation->size() == 0) return false;

	neighbourhoods.resize(orientation->size());

	for (int j = 0; j < orientation->size(); j++){
		std::vector<double> d;
		std::vector<int> idx;
		for (int k = 0; k < orientation->size(); k++){
			double dx = orientation->at(j).x() - orientation->at(k).x();
			double dy = orientation->at(j).y() - orientation->at(k).y();
			double dz = orientation->at(j).z() - orientation->at(k).z();
			double r2 = dx*dx + dy*dy + dz*dz;
			d.push_back(r2);
			idx.push_back(k);
		}
		if (d.size() < n_neighbors) return false;
		Math_methods::sort_vector_w_index(d, idx);
		for (int k = 0; k < n_neighbors; k++){
			neighbourhoods[j].push_back(orientation->at(idx[k]));
		}
	}

	return true;
}

void TensorInput::_get_local_anisotropy_from_neighbourhoods(std::vector < std::vector < Orientation > > &nh)
{
	for (int j = 0; j < nh.size(); j++){
		Matrix3d OrientationMatrix;
		double SumXX = 0; // Sum(x_i *  x_i)
		double SumXY = 0; // Sum(x_i *  y_i)
		double SumXZ = 0; // Sum(x_i *  z_i)
		double SumYY = 0; // Sum(y_i *  y_i)
		double SumYZ = 0; // Sum(y_i *  z_i)
		double SumZZ = 0; // Sum(z_i *  z_i)
		int nsize = (int)nh[j].size();
		for (int k = 0; k < nh[j].size(); k++){
			//std::cout << "Normal[" << k << "] = {" << nh[j].at(k).nx() << "," << nh[j].at(k).ny() << "," << nh[j].at(k).nz() << "}" << std::endl;
			SumXX += nh[j].at(k).nx() * nh[j].at(k).nx();
			SumYY += nh[j].at(k).ny() * nh[j].at(k).ny();
			SumZZ += nh[j].at(k).nz() * nh[j].at(k).nz();
			SumXY += nh[j].at(k).nx() * nh[j].at(k).ny();
			SumXZ += nh[j].at(k).nx() * nh[j].at(k).nz();
			SumYZ += nh[j].at(k).ny() * nh[j].at(k).nz();
		}
		OrientationMatrix(0,0) = SumXX;
		OrientationMatrix(0,1) = SumXY;
		OrientationMatrix(0,2) = SumXZ;
		OrientationMatrix(1,0) = SumXY;
		OrientationMatrix(1,1) = SumYY;
		OrientationMatrix(1,2) = SumYZ;
		OrientationMatrix(2,0) = SumXZ;
		OrientationMatrix(2,1) = SumYZ;
		OrientationMatrix(2,2) = SumZZ;
		//std::cout << " OrientationMatrix[" << j << "]:\n" << OrientationMatrix << std::endl;
		SelfAdjointEigenSolver<Matrix3d> ea;
		ea.compute(OrientationMatrix);
		Matrix3d evectors;
		Vector3d evalues;
		evectors = ea.eigenvectors().real();
		evalues = ea.eigenvalues().real();
// 		std::cout << " Eigenvectors:\n" << evectors << std::endl;
// 		std::cout << " Eigenvalues:\n" << evalues << std::endl;
		// The only negative eigenvalues that should come 
		// out of the eigen analysis is values within machine
		// precision to 0 (e.g. -5e-16 etc) from tensor matrices
		for (int k = 0; k < 3; k++){
			if (evalues[k] < 0) evalues[k] = abs(evalues[k]);
		}
		double zero_limit = 0.001;
		if (evalues[0] < zero_limit || evalues[1] < zero_limit || evalues[2] < zero_limit)
		{
			for (int k = 0; k < 3; k++){
				if (evalues[k] < zero_limit) evalues[k] = zero_limit;
			}
		}
		_sort_eigensystem(evectors, evalues);
		double det = evectors.determinant();
		if (det < 0) // flip plunge vector
		{
			//std::cout << " Had to flip" << std::endl;
			evectors(0,0) *= -1;
			evectors(1,0) *= -1;
			evectors(2,0) *= -1;
		}
		// the following code is to fix numerical errors that arise from scenarios where
		// there is no discernible plunge detected: occurs when the two smallest eigenvalues ~ 0 
		if (evalues[1] < 0.1 && evalues[2] < 0.1)
		{
			if (abs(evalues[1] - evalues[2]) < 0.1)
			{
				evalues[1] += 0.1;
				evalues[2] += 0.3;
				// ensure eigenvalues add up to the required value: total# of members within the local neighbourhood
				double cur_eig_add = evalues[0] + evalues[1] + evalues[2];
				double ratio = nsize / cur_eig_add;
				for (int k = 0; k < 3; k++) evalues[k] *= ratio;
			}
		}
		_sort_eigensystem(evectors, evalues);
		orientation->at(j).eigenvectors = evectors;
		orientation->at(j).eigenvalues = evalues;
		orientation->at(j).U = evectors;
		for (int k = 0; k < 3; k++){
			for (int l = 0; l < 3; l++){
				if (k == l)
				{
					orientation->at(j).D(k, l) = evalues[k];
					orientation->at(j).S(k, l) = sqrt(evalues[k]/evalues[0]);
				}
				else
				{
					orientation->at(j).D(k, l) = 0.0;
					orientation->at(j).S(k, l) = 0.0;
				}
			}
		}
		orientation->at(j).Tensor = orientation->at(j).U * orientation->at(j).D * orientation->at(j).U.transpose();
		orientation->at(j).Transform = orientation->at(j).U * orientation->at(j).S * orientation->at(j).U.transpose();
	}
}

bool TensorInput::GetLocalAnisotropy(const model_parameters &parameters)
{
	// 1st get neighbourhoods
	if (!_get_neighbourhoods(parameters.nearest_neighbours)) return false;

	_get_local_anisotropy_from_neighbourhoods(neighbourhoods);

	return true;
}

Vector3d TensorEvaluationPoints::GetDipVector()
{
	Vector3d dipvector;
	dipvector[0] = eigenvectors(0, 1);
	dipvector[1] = eigenvectors(1, 1);
	dipvector[2] = eigenvectors(2, 1);
	return dipvector;
}

Vector3d TensorEvaluationPoints::GetStrikeVector()
{
	Vector3d strikevector;
	strikevector[0] = eigenvectors(0, 0);
	strikevector[1] = eigenvectors(1, 0);
	strikevector[2] = eigenvectors(2, 0);
	return strikevector;
}

Vector3d TensorEvaluationPoints::GetNormalVector()
{
	Vector3d normalvector;
	normalvector[0] = eigenvectors(0, 2);
	normalvector[1] = eigenvectors(1, 2);
	normalvector[2] = eigenvectors(2, 2);
	return normalvector;
}
