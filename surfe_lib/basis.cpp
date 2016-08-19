#include <basis.h>
#include <math_methods.h>

#include <algorithm>

using namespace std;

inline void RBFKernel::radius()
{
	// compute euclidean distance
	// compute component differences
	_x_delta = _p1->x() - _p2->x();
	_y_delta = _p1->y() - _p2->y();
	_z_delta = _p1->z() - _p2->z();
	_c_delta = _p1->c() - _p2->c();

	_radius = sqrt( _x_delta*_x_delta + _y_delta*_y_delta + _z_delta*_z_delta + _c_delta*_c_delta );
}

inline void RBFKernel::scaled_radius()
{
	// compute scaled euclidean distance
	// compute component differences
	double dx = _p1->x() - _p2->x();
	double dy = _p1->y() - _p2->y();
	double dz = _p1->z() - _p2->z();
	// compute scaled component differences
	_x_delta = _Transform(0,0)*dx + _Transform(0,1)*dy + _Transform(0,2)*dz;
	_y_delta = _Transform(1,0)*dx + _Transform(1,1)*dy + _Transform(1,2)*dz;
	_z_delta = _Transform(2,0)*dx + _Transform(2,1)*dy + _Transform(2,2)*dz;

	_radius = sqrt( _x_delta*_x_delta + _y_delta*_y_delta + _z_delta*_z_delta );
}

bool RBFKernel::get_global_anisotropy(const std::vector<Planar> &planar)
{
	if ((int)planar.size() < 2) return false;

	double SumXX = 0; // Sum(x_i *  x_i)
	double SumXY = 0; // Sum(x_i *  y_i)
	double SumXZ = 0; // Sum(x_i *  z_i)
	double SumYY = 0; // Sum(y_i *  y_i)
	double SumYZ = 0; // Sum(y_i *  z_i)
	double SumZZ = 0; // Sum(z_i *  z_i)
	Matrix3f covMat;
	for (int j = 0; j < (int)planar.size(); j++) {
		double VNormal[3] = { planar[j].nx(),planar[j].ny(),planar[j].nz() };
		// Normals...
		SumXX += VNormal[0] * VNormal[0];
		SumYY += VNormal[1] * VNormal[1];
		SumZZ += VNormal[2] * VNormal[2];
		SumXY += VNormal[0] * VNormal[1];
		SumXZ += VNormal[0] * VNormal[2];
		SumYZ += VNormal[1] * VNormal[2];
	}
	covMat(0,0) = SumXX;
	covMat(0,1) = SumXY;
	covMat(0,2) = SumXZ;
	covMat(1,0) = SumXY;
	covMat(1,1) = SumYY;
	covMat(1,2) = SumYZ;
	covMat(2,0) = SumXZ;
	covMat(2,1) = SumYZ;
	covMat(2,2) = SumZZ;

	SelfAdjointEigenSolver<Matrix3f> es;
	es.compute(covMat);
	Vector3f eVals = es.eigenvalues().real();
	Matrix3f eVectors = es.eigenvectors().real();
	//cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
	//cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;

	_Global_Plunge[0] = eVectors(0,0);
	_Global_Plunge[1] = eVectors(1,0);
	_Global_Plunge[2] = eVectors(2,0);

	//double plunge_vector[3]       = { eVectors(0,0),eVectors(1,0),eVectors(2,0) };
	//double facing_vector[3]       = { eVectors(0,1),eVectors(1,1),eVectors(2,1) };
	//double cross_plunge_vector[3] = { eVectors(0,2),eVectors(1,2),eVectors(2,2) };

    // might have to check eval(0) smallest eigenvalue
	// if there are normals sampled from a perfect cylinderical fold
	// eval(0) = ~1e-16 machine precision.
	// perhaps set to eval(0) = 0.0001 ...

	Matrix3f ISQR;
	ISQR(0,0) = 1;
	ISQR(1,0) = 0;
	ISQR(2,0) = 0;
	ISQR(0,1) = 0;
	ISQR(1,1) = sqrt(eVals(1)/eVals(0));
	ISQR(2,1) = 0;
	ISQR(0,2) = 0;
	ISQR(1,2) = 0;
	ISQR(2,2) = sqrt(eVals(2)/eVals(0));

	_Transform = eVectors*ISQR*eVectors.transpose();
	//cout << "The _Transform is:" << endl << _Transform << endl << endl;

	return true; // should really validate the computation here
}

double RBFKernel::basis_pt_pt()
{
	return this->basis();
}

double RBFKernel::basis_pt_planar_x()
{
	return this->dx_p2();
}

double RBFKernel::basis_planar_x_pt()
{
	return this->dx_p1();
}

double RBFKernel::basis_pt_planar_y()
{
	return this->dy_p2();
}

double RBFKernel::basis_planar_y_pt()
{
	return this->dy_p1();
}

double RBFKernel::basis_pt_planar_z()
{
	return this->dz_p2();
}

double RBFKernel::basis_planar_z_pt()
{
	return this->dz_p1();
}

double RBFKernel::basis_pt_tangent()
{
	// It is your responsibility to supply the right type
	// It is not worth the cost of dynamic_cast + if
	Tangent *t = static_cast<Tangent*>(this->p2());
	return this->dx_p2()*t->tx() + this->dy_p2()*t->ty() + this->dz_p2()*t->tz();
}

double RBFKernel::basis_tangent_pt()
{
	Tangent *t = static_cast<Tangent*>(this->p1());
	return this->dx_p1()*t->tx() + this->dy_p1()*t->ty() + this->dz_p1()*t->tz();
}

double RBFKernel::basis_planar_planar( const Parameter_Types::SecondDerivatives& sd )
{
	if (sd == Parameter_Types::DXDX) return this->dxx();
	else if (sd == Parameter_Types::DXDY) return this->dxy();
	else if (sd == Parameter_Types::DXDZ) return this->dxz();
	else if (sd == Parameter_Types::DYDX) return this->dyx();
	else if (sd == Parameter_Types::DYDY) return this->dyy();
	else if (sd == Parameter_Types::DYDZ) return this->dyz();
	else if (sd == Parameter_Types::DZDX) return this->dzx();
	else if (sd == Parameter_Types::DZDY) return this->dzy();
	else return this->dzz();
}

double RBFKernel::basis_tangent_tangent()
{
	Tangent *t1 = static_cast<Tangent*>(this->p1());
	Tangent *t2 = static_cast<Tangent*>(this->p2());

	return t1->tx()*t2->tx() * this->dxx() +
		   t1->tx()*t2->ty() * this->dxy() +
		   t1->tx()*t2->tz() * this->dxz() +
		   t1->ty()*t2->tx() * this->dyx() +
		   t1->ty()*t2->ty() * this->dyy() +
		   t1->ty()*t2->tz() * this->dyz() +
		   t1->tz()*t2->tx() * this->dzx() +
		   t1->tz()*t2->ty() * this->dzy() +
		   t1->tz()*t2->tz() * this->dzz();
}

double RBFKernel::basis_planar_tangent(const Parameter_Types::FirstDerivatives& fd)
{
	Tangent *t = static_cast<Tangent*>(this->p2());
	double tdx = 0;
	double tdy = 0;
	double tdz = 0;
	if (fd == Parameter_Types::DX)
	{
		tdx = this->dxx();
		tdy = this->dxy();
		tdz = this->dxz();
	}
	if (fd == Parameter_Types::DY)
	{
		tdx = this->dyx();
		tdy = this->dyy();
		tdz = this->dyz();
	}
	if (fd == Parameter_Types::DZ)
	{
		tdx = this->dzx();
		tdy = this->dzy();
		tdz = this->dzz();
	}
	return t->tx()*tdx + t->ty()*tdy + t->tz()*tdz;
}

double RBFKernel::basis_tangent_planar(const Parameter_Types::FirstDerivatives& fd)
{
	Tangent *t = static_cast<Tangent*>(this->p1());
	double tdx = 0;
	double tdy = 0;
	double tdz = 0;
	if (fd == Parameter_Types::DX)
	{
		tdx = this->dxx();
		tdy = this->dyx();
		tdz = this->dzx();
	}
	if (fd == Parameter_Types::DY)
	{
		tdx = this->dxy();
		tdy = this->dyy();
		tdz = this->dzy();
	}
	if (fd == Parameter_Types::DZ)
	{
		tdx = this->dxz();
		tdy = this->dyz();
		tdz = this->dzz();
	}
	return t->tx()*tdx + t->ty()*tdy + t->tz()*tdz;
}



double Cubic::basis()
{
	radius();
	return _radius*_radius*_radius;
}

double Cubic::dx_p1()
{
	radius();
	return 3*_radius*_x_delta;
}

double Cubic::dx_p2()
{
	radius();
	return -3*_radius*_x_delta;
}

double Cubic::dy_p1()
{
	radius();
	return 3*_radius*_y_delta;
}

double Cubic::dy_p2()
{
	radius();
	return -3*_radius*_y_delta;
}

double Cubic::dz_p1()
{
	radius();
	return 3*_radius*_z_delta;
}

double Cubic::dz_p2()
{
	radius();
	return -3*_radius*_z_delta;
}

double Cubic::dxx()
{
	radius();
	if (_radius == 0) return 0.0;
	else return -3.0* ( ((_x_delta*_x_delta)/_radius) + _radius );
}

double Cubic::dxy()
{
	radius();
	if (_radius == 0) return 0.0;
	else return -3.0*((_x_delta*_y_delta)/_radius);
}

double Cubic::dxz()
{
	radius();
	if (_radius == 0) return 0.0;
	else return -3.0*((_x_delta*_z_delta) / _radius);
}

double Cubic::dyx()
{
	return dxy();
}

double Cubic::dyy()
{
	radius();
	if (_radius == 0) return 0.0;
	else return -3.0* (((_y_delta*_y_delta) / _radius) + _radius);
}

double Cubic::dyz()
{
	radius();
	if (_radius == 0) return 0.0;
	else return -3.0*((_y_delta*_z_delta) / _radius);
}

double Cubic::dzx()
{
	return dxz();
}

double Cubic::dzy()
{
	return dyz();
}

double Cubic::dzz()
{
	radius();
	if (_radius == 0) return 0.0;
	else return -3.0* (((_z_delta*_z_delta) / _radius) + _radius);
}

double ACubic::basis()
{
	scaled_radius();
	return _radius*_radius*_radius;
}

double ACubic::dx_p1()
{
	scaled_radius();
	return 3.0*_radius*(_Transform(0,0)*_x_delta + _Transform(1,0)*_y_delta + _Transform(2,0)*_z_delta);
}

double ACubic::dx_p2()
{
	scaled_radius();
	return -3.0 * _radius*(_Transform(0,0) * _x_delta + _Transform(1,0) * _y_delta + _Transform(2,0) * _z_delta);
}

double ACubic::dy_p1()
{
	scaled_radius();
	return 3.0 * _radius*(_Transform(0,1) * _x_delta + _Transform(1,1) * _y_delta + _Transform(2,1) * _z_delta);
}

double ACubic::dy_p2()
{
	scaled_radius();
	return -3.0 * _radius*(_Transform(0,1) * _x_delta + _Transform(1,1) * _y_delta + _Transform(2,1) * _z_delta);
}

double ACubic::dz_p1()
{
	scaled_radius();
	return 3.0 * _radius*(_Transform(0,2) * _x_delta + _Transform(1,2) * _y_delta + _Transform(2,2) * _z_delta);
}

double ACubic::dz_p2()
{
	scaled_radius();
	return -3.0 * _radius*(_Transform(0,2) * _x_delta + _Transform(1,2) * _y_delta + _Transform(2,2) * _z_delta);
}

double ACubic::dxx()
{
	scaled_radius();
	if (_radius == 0) return 0.0;
	else
	{
		double a = _Transform(0,0) * _x_delta + _Transform(1,0) * _y_delta + _Transform(2,0) * _z_delta;
		double b = _Transform(0,0) * _Transform(0,0) + _Transform(1,0) * _Transform(1,0) + _Transform(2,0) * _Transform(2,0);
		return -3.0* (((a*a) / _radius) + b*_radius);
	}
}

double ACubic::dxy()
{
	scaled_radius();
	if (_radius == 0) return 0.0;
	else
	{
		double a = _Transform(0,0) * _x_delta + _Transform(1,0) * _y_delta + _Transform(2,0) * _z_delta;
		double b = _Transform(0,1) * _x_delta + _Transform(1,1) * _y_delta + _Transform(2,1) * _z_delta;
		double c = _Transform(0,0) * _Transform(0,1) + _Transform(1,0) * _Transform(1,1) + _Transform(2,0) * _Transform(2,1);
		return -3.0* (((a*b) / _radius) + c*_radius);
	}
}

double ACubic::dxz()
{
	scaled_radius();
	if (_radius == 0) return 0.0;
	else
	{
		double a = _Transform(0,0) * _x_delta + _Transform(1,0) * _y_delta + _Transform(2,0) * _z_delta;
		double b = _Transform(0,2) * _x_delta + _Transform(1,2) * _y_delta + _Transform(2,2) * _z_delta;
		double c = _Transform(0,0) * _Transform(0,2) + _Transform(1,0) * _Transform(1,2) + _Transform(2,0) * _Transform(2,2);
		return -3.0* (((a*b) / _radius) + c*_radius);
	}
}

double ACubic::dyx()
{
	return dxy();
}

double ACubic::dyy()
{
	scaled_radius();
	if (_radius == 0) return 0.0;
	else
	{
		double a = _Transform(0,1) * _x_delta + _Transform(1,1) * _y_delta + _Transform(2,1) * _z_delta;
		double b = _Transform(0,1) * _Transform(0,1) + _Transform(1,1) * _Transform(1,1) + _Transform(2,1) * _Transform(2,1);
		return -3.0* (((a*a) / _radius) + b*_radius);
	}
}

double ACubic::dyz()
{
	scaled_radius();
	if (_radius == 0) return 0.0;
	else
	{
		double a = _Transform(0,1) * _x_delta + _Transform(1,1) * _y_delta + _Transform(2,1) * _z_delta;
		double b = _Transform(0,2) * _x_delta + _Transform(1,2) * _y_delta + _Transform(2,2) * _z_delta;
		double c = _Transform(0,1) * _Transform(0,2) + _Transform(1,1) * _Transform(1,2) + _Transform(2,1) * _Transform(2,2);
		return -3.0* (((a*b) / _radius) + c*_radius);
	}
}

double ACubic::dzx()
{
	return dxz();
}

double ACubic::dzy()
{
	return dyz();
}

double ACubic::dzz()
{
	scaled_radius();
	if (_radius == 0) return 0.0;
	else
	{
		double a = _Transform(0,2) * _x_delta + _Transform(1,2) * _y_delta + _Transform(2,2) * _z_delta;
		double b = _Transform(0,2) * _Transform(0,2) + _Transform(1,2) * _Transform(1,2) + _Transform(2,2) * _Transform(2,2);
		return -3.0* (((a*a) / _radius) + b*_radius);
	}
}

double Gaussian::basis()
{
	radius();
	return exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
}

double Gaussian::dx_p1()
{
	radius();
	return -2.0 * _shape_parameter*_shape_parameter*_x_delta*exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
}

double Gaussian::dx_p2()
{
	radius();
	return 2.0 * _shape_parameter*_shape_parameter*_x_delta*exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
}

double Gaussian::dy_p1()
{
	radius();
	return -2.0 * _shape_parameter*_shape_parameter*_y_delta*exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
}

double Gaussian::dy_p2()
{
	radius();
	return 2.0 * _shape_parameter*_shape_parameter*_y_delta*exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
}

double Gaussian::dz_p1()
{
	radius();
	return -2.0 * _shape_parameter*_shape_parameter*_z_delta*exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
}

double Gaussian::dz_p2()
{
	radius();
	return 2.0 * _shape_parameter*_shape_parameter*_z_delta*exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
}

double Gaussian::dxx()
{
	radius();
	return (2.0*_shape_parameter*_shape_parameter - 4.0*pow(_shape_parameter,4)*_x_delta*_x_delta)*exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
}

double Gaussian::dxy()
{
	radius();
	return (-4.0*pow(_shape_parameter, 4)*_x_delta*_y_delta)*exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
}
double Gaussian::dxz()
{
	radius();
	return (-4.0*pow(_shape_parameter, 4)*_x_delta*_z_delta)*exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
}

double Gaussian::dyx()
{
	return dxy();
}

double Gaussian::dyy()
{
	radius();
	return (2.0*_shape_parameter*_shape_parameter - 4.0*pow(_shape_parameter, 4)*_y_delta*_y_delta)*exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
}

double Gaussian::dyz()
{
	radius();
	return (-4.0*pow(_shape_parameter, 4)*_y_delta*_z_delta)*exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
}

double Gaussian::dzx()
{
	return dxz();
}

double Gaussian::dzy()
{
	return dyz();
}

double Gaussian::dzz()
{
	radius();
	return (2.0*_shape_parameter*_shape_parameter - 4.0*pow(_shape_parameter, 4)*_z_delta*_z_delta)*exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
}

double AGaussian::basis()
{
	scaled_radius();
	return exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
}

double AGaussian::dx_p1()
{
	scaled_radius();
	return -2.0*_shape_parameter*_shape_parameter*(_Transform(0,0) * _x_delta + _Transform(1,0) * _y_delta + 
		_Transform(2,0) * _z_delta) * exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
}

double AGaussian::dx_p2()
{
	scaled_radius();
	return 2.0*_shape_parameter*_shape_parameter*(_Transform(0,0) * _x_delta + _Transform(1,0) * _y_delta +
		_Transform(2,0) * _z_delta) * exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
}

double AGaussian::dy_p1()
{
	scaled_radius();
	return -2.0*_shape_parameter*_shape_parameter*(_Transform(0,1) * _x_delta + _Transform(1,1) * _y_delta +
		_Transform(2,1) * _z_delta) * exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
}

double AGaussian::dy_p2()
{
	scaled_radius();
	return 2.0*_shape_parameter*_shape_parameter*(_Transform(0,1) * _x_delta + _Transform(1,1) * _y_delta +
		_Transform(2,1) * _z_delta) * exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
}

double AGaussian::dz_p1()
{
	scaled_radius();
	return -2.0*_shape_parameter*_shape_parameter*(_Transform(0,2) * _x_delta + _Transform(1,2) * _y_delta +
		_Transform(2,2) * _z_delta) * exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
}

double AGaussian::dz_p2()
{
	scaled_radius();
	return 2.0*_shape_parameter*_shape_parameter*(_Transform(0,2) * _x_delta + _Transform(1,2) * _y_delta +
		_Transform(2,2) * _z_delta) * exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
}

double AGaussian::dxx()
{
	scaled_radius();
	double a = _Transform(0,0) * _Transform(0,0) + _Transform(1,0) * _Transform(1,0) + _Transform(2,0) * _Transform(2,0);
	double b = _Transform(0,0) * _x_delta + _Transform(1,0) * _y_delta + _Transform(2,0) * _z_delta;
	double c = _shape_parameter*_shape_parameter* exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
	return 2.0*c*(a - 2.0*_shape_parameter*_shape_parameter*b*b);
}

double AGaussian::dxy()
{
	scaled_radius();
	double a = _Transform(0,0) * _Transform(0,1) + _Transform(1,0) * _Transform(1,1) + _Transform(2,0) * _Transform(2,1);
	double b = _Transform(0,0) * _x_delta + _Transform(1,0) * _y_delta + _Transform(2,0) * _z_delta;
	double c = _Transform(0,1) * _x_delta + _Transform(1,1) * _y_delta + _Transform(2,1) * _z_delta;
	double d = _shape_parameter*_shape_parameter* exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
	return 2.0*d*(a - 2.0*_shape_parameter*_shape_parameter*b*c);
}

double AGaussian::dxz()
{
	scaled_radius();
	double a = _Transform(0,0) * _Transform(0,2) + _Transform(1,0) * _Transform(1,2) + _Transform(2,0) * _Transform(2,2);
	double b = _Transform(0,0) * _x_delta + _Transform(1,0) * _y_delta + _Transform(2,0) * _z_delta;
	double c = _Transform(0,2) * _x_delta + _Transform(1,2) * _y_delta + _Transform(2,2) * _z_delta;
	double d = _shape_parameter*_shape_parameter* exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
	return 2.0*d*(a - 2.0*_shape_parameter*_shape_parameter*b*c);
}

double AGaussian::dyx()
{
	return dxy();
}

double AGaussian::dyy()
{
	scaled_radius();
	double a = _Transform(0,1) * _Transform(0,1) + _Transform(1,1) * _Transform(1,1) + _Transform(2,1) * _Transform(2,1);
	double b = _Transform(0,1) * _x_delta + _Transform(1,1) * _y_delta + _Transform(2,1) * _z_delta;
	double c = _shape_parameter*_shape_parameter* exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
	return 2.0*c*(a - 2.0*_shape_parameter*_shape_parameter*b*b);
}

double AGaussian::dyz()
{
	scaled_radius();
	double a = _Transform(0,1) * _Transform(0,2) + _Transform(1,1) * _Transform(1,2) + _Transform(2,1) * _Transform(2,2);
	double b = _Transform(0,1) * _x_delta + _Transform(1,1) * _y_delta + _Transform(2,1) * _z_delta;
	double c = _Transform(0,2) * _x_delta + _Transform(1,2) * _y_delta + _Transform(2,2) * _z_delta;
	double d = _shape_parameter*_shape_parameter* exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
	return 2.0*d*(a - 2.0*_shape_parameter*_shape_parameter*b*c);
}

double AGaussian::dzx()
{
	return dxz();
}

double AGaussian::dzy()
{
	return dyz();
}

double AGaussian::dzz()
{
	scaled_radius();
	double a = _Transform(0,2) * _Transform(0,2) + _Transform(1,2) * _Transform(1,2) + _Transform(2,2) * _Transform(2,2);
	double b = _Transform(0,2) * _x_delta + _Transform(1,2) * _y_delta + _Transform(2,2) * _z_delta;
	double c = _shape_parameter*_shape_parameter* exp(-(_shape_parameter*_shape_parameter*_radius*_radius));
	return 2.0*c*(a - 2.0*_shape_parameter*_shape_parameter*b*b);
}

double MQ::basis()
{
	radius();
	return pow(_shape_parameter + _radius*_radius, 0.5);
}

double MQ::dx_p1()
{
	radius();
	return _x_delta/pow(_shape_parameter + _radius*_radius, 0.5);
}

double MQ::dx_p2()
{
	radius();
	return (-1.0*_x_delta) / pow(_shape_parameter + _radius*_radius, 0.5);
}

double MQ::dy_p1()
{
	radius();
	return _y_delta / pow(_shape_parameter + _radius*_radius, 0.5);
}

double MQ::dy_p2()
{
	radius();
	return (-1.0*_y_delta) / pow(_shape_parameter + _radius*_radius, 0.5);
}

double MQ::dz_p1()
{
	radius();
	return _z_delta / pow(_shape_parameter + _radius*_radius, 0.5);
}

double MQ::dz_p2()
{
	radius();
	return (-1.0*_z_delta) / pow(_shape_parameter + _radius*_radius, 0.5);
}

double MQ::dxx()
{
	radius();
	return (_x_delta*_x_delta) / pow(_shape_parameter + _radius*_radius, 1.5) - (1.0 / pow(_shape_parameter + _radius*_radius, 0.5));
}
double MQ::dxy()
{
	radius();
	return (_x_delta*_y_delta) / pow(_shape_parameter + _radius*_radius, 1.5);
}

double MQ::dxz()
{
	radius();
	return (_x_delta*_z_delta) / pow(_shape_parameter + _radius*_radius, 1.5);
}

double MQ::dyx()
{
	return dxy();
}

double MQ::dyy()
{
	radius();
	return (_y_delta*_y_delta) / pow(_shape_parameter + _radius*_radius, 1.5) - (1.0 / pow(_shape_parameter + _radius*_radius, 0.5));
}

double MQ::dyz()
{
	radius();
	return (_y_delta*_z_delta) / pow(_shape_parameter + _radius*_radius, 1.5);
}

double MQ::dzx()
{
	return dxz();
}

double MQ::dzy()
{
	return dyz();
}

double MQ::dzz()
{
	radius();
	return (_z_delta*_z_delta) / pow(_shape_parameter + _radius*_radius, 1.5) - (1.0 / pow(_shape_parameter + _radius*_radius, 0.5));
}

double MQ3::basis()
{
	_p1->set_c(_c);
	_p2->set_c(0.0);
	radius();
	return _radius*_radius*_radius;
}

double MQ3::dx_p1()
{
	_p1->set_c(_c);
	_p2->set_c(0.0);
	radius();
	return 3 * _radius*_x_delta;
}

double MQ3::dx_p2()
{
	_p1->set_c(_c);
	_p2->set_c(0.0);
	radius();
	return -3 * _radius*_x_delta;
}

double MQ3::dy_p1()
{
	_p1->set_c(_c);
	_p2->set_c(0.0);
	radius();
	return 3 * _radius*_y_delta;
}

double MQ3::dy_p2()
{
	_p1->set_c(_c);
	_p2->set_c(0.0);
	radius();
	return -3 * _radius*_y_delta;
}

double MQ3::dz_p1()
{
	_p1->set_c(_c);
	_p2->set_c(0.0);
	radius();
	return 3 * _radius*_z_delta;
}

double MQ3::dz_p2()
{
	_p1->set_c(_c);
	_p2->set_c(0.0);
	radius();
	return -3 * _radius*_z_delta;
}

double MQ3::dxx()
{
	_p1->set_c(_c);
	_p2->set_c(0.0);
	radius();
	if (_radius == 0) return 0.0;
	else return -3.0* (((_x_delta*_x_delta) / _radius) + _radius);
}

double MQ3::dxy()
{
	_p1->set_c(_c);
	_p2->set_c(0.0);
	radius();
	if (_radius == 0) return 0.0;
	else return -3.0*((_x_delta*_y_delta) / _radius);
}

double MQ3::dxz()
{
	_p1->set_c(_c);
	_p2->set_c(0.0);
	radius();
	if (_radius == 0) return 0.0;
	else return -3.0*((_x_delta*_z_delta) / _radius);
}

double MQ3::dyx()
{
	_p1->set_c(_c);
	_p2->set_c(0.0);
	return dxy();
}

double MQ3::dyy()
{
	_p1->set_c(_c);
	_p2->set_c(0.0);
	radius();
	if (_radius == 0) return 0.0;
	else return -3.0* (((_y_delta*_y_delta) / _radius) + _radius);
}

double MQ3::dyz()
{
	_p1->set_c(_c);
	_p2->set_c(0.0);
	radius();
	if (_radius == 0) return 0.0;
	else return -3.0*((_y_delta*_z_delta) / _radius);
}

double MQ3::dzx()
{
	_p1->set_c(_c);
	_p2->set_c(0.0);
	return dxz();
}

double MQ3::dzy()
{
	_p1->set_c(_c);
	_p2->set_c(0.0);
	return dyz();
}

double MQ3::dzz()
{
	_p1->set_c(_c);
	_p2->set_c(0.0);
	radius();
	if (_radius == 0) return 0.0;
	else return -3.0* (((_z_delta*_z_delta) / _radius) + _radius);
}


double AMQ::basis()
{
	scaled_radius();
	return pow(_shape_parameter + _radius*_radius, 0.5);
}

double AMQ::dx_p1()
{
	scaled_radius();
	return (_Transform(0,0) * _x_delta + _Transform(1,0) * _y_delta + _Transform(2,0) * _z_delta) / pow(_shape_parameter + _radius*_radius, 0.5);
}

double AMQ::dx_p2()
{
	scaled_radius();
	return -1.0*(_Transform(0,0) * _x_delta + _Transform(1,0) * _y_delta + _Transform(2,0) * _z_delta) / pow(_shape_parameter + _radius*_radius, 0.5);
}

double AMQ::dy_p1()
{
	scaled_radius();
	return (_Transform(0,1) * _x_delta + _Transform(1,1) * _y_delta + _Transform(2,1) * _z_delta) / pow(_shape_parameter + _radius*_radius, 0.5);
}

double AMQ::dy_p2()
{
	scaled_radius();
	return -1.0*(_Transform(0,1) * _x_delta + _Transform(1,1) * _y_delta + _Transform(2,1) * _z_delta) / pow(_shape_parameter + _radius*_radius, 0.5);
}

double AMQ::dz_p1()
{
	scaled_radius();
	return (_Transform(0,2) * _x_delta + _Transform(1,2) * _y_delta + _Transform(2,2) * _z_delta) / pow(_shape_parameter + _radius*_radius, 0.5);
}

double AMQ::dz_p2()
{
	scaled_radius();
	return -1.0*(_Transform(0,2) * _x_delta + _Transform(1,2) * _y_delta + _Transform(2,2) * _z_delta) / pow(_shape_parameter + _radius*_radius, 0.5);
}

double AMQ::dxx()
{
	scaled_radius();
	double a = (-2.0*_Transform(0,0) * _x_delta - 2.0*_Transform(1,0) * _y_delta - 2.0*_Transform(2,0) * _z_delta);
	double b = (2.0*_Transform(0,0) * _x_delta + 2.0*_Transform(1,0) * _y_delta + 2.0*_Transform(2,0) * _z_delta);
	double c = 4.0*pow(_shape_parameter + _radius*_radius, 1.5);
	double d = -1.0*(a * b) / c;
	double f = (-_Transform(0,0) * _Transform(0,0) - _Transform(1,0) * _Transform(1,0) - _Transform(2,0) * _Transform(2,0)) / pow(_shape_parameter + _radius*_radius, 0.5);
	return d + f;
}

double AMQ::dxy()
{
	scaled_radius();
	double a = (2.0*_Transform(0,0) * _x_delta + 2.0*_Transform(1,0) * _y_delta + 2.0*_Transform(2,0) * _z_delta);
	double b = (-2.0*_Transform(0,1) * _x_delta - 2.0*_Transform(1,1) * _y_delta - 2.0*_Transform(2,1) * _z_delta);
	double c = 4.0*pow(_shape_parameter + _radius*_radius, 1.5);
	double d = -1.0*(a * b) / c;
	double f = (-_Transform(0,0) * _Transform(0,1) - _Transform(1,0) * _Transform(1,1) - _Transform(2,0) * _Transform(2,1)) / pow(_shape_parameter + _radius*_radius, 0.5);
	return d + f;
}

double AMQ::dxz()
{
	scaled_radius();
	double a = (2.0*_Transform(0,0) * _x_delta + 2.0*_Transform(1,0) * _y_delta + 2.0*_Transform(2,0) * _z_delta);
	double b = (-2.0*_Transform(0,2) * _x_delta - 2.0*_Transform(1,2) * _y_delta - 2.0*_Transform(2,2) * _z_delta);
	double c = 4.0*pow(_shape_parameter + _radius*_radius, 1.5);
	double d = -1.0*(a * b) / c;
	double f = (-_Transform(0,0) * _Transform(0,2) - _Transform(1,0) * _Transform(1,2) - _Transform(2,0) * _Transform(2,2)) / pow(_shape_parameter + _radius*_radius, 0.5);
	return d + f;
}

double AMQ::dyx()
{
	return dxy();
}

double AMQ::dyy()
{
	scaled_radius();
	double a = (-2.0*_Transform(0,1) * _x_delta - 2.0*_Transform(1,1) * _y_delta - 2.0*_Transform(2,1) * _z_delta);
	double b = (2.0*_Transform(0,1) * _x_delta + 2.0*_Transform(1,1) * _y_delta + 2.0*_Transform(2,1) * _z_delta);
	double c = 4.0*pow(_shape_parameter + _radius*_radius, 1.5);
	double d = -1.0*(a * b) / c;
	double f = (-_Transform(0,1) * _Transform(0,1) - _Transform(1,1) * _Transform(1,1) - _Transform(2,1) * _Transform(2,1)) / pow(_shape_parameter + _radius*_radius, 0.5);
	return d + f;
}

double AMQ::dyz()
{
	scaled_radius();
	double a = (2.0*_Transform(0,1) * _x_delta + 2.0*_Transform(1,1) * _y_delta + 2.0*_Transform(2,1) * _z_delta);
	double b = (-2.0*_Transform(0,2) * _x_delta - 2.0*_Transform(1,2) * _y_delta - 2.0*_Transform(2,2) * _z_delta);
	double c = 4.0*pow(_shape_parameter + _radius*_radius, 1.5);
	double d = -1.0*(a * b) / c;
	double f = (-_Transform(0,1) * _Transform(0,2) - _Transform(1,1) * _Transform(1,2) - _Transform(2,1) * _Transform(2,2)) / pow(_shape_parameter + _radius*_radius, 0.5);
	return d + f;
}

double AMQ::dzx()
{
	return dxz();
}

double AMQ::dzy()
{
	return dyz();
}

double AMQ::dzz()
{
	scaled_radius();
	double a = (-2.0*_Transform(0,2) * _x_delta - 2.0*_Transform(1,2) * _y_delta - 2.0*_Transform(2,2) * _z_delta);
	double b = (2.0*_Transform(0,2) * _x_delta + 2.0*_Transform(1,2) * _y_delta + 2.0*_Transform(2,2) * _z_delta);
	double c = 4.0*pow(_shape_parameter + _radius*_radius, 1.5);
	double d = -1.0*(a * b) / c;
	double f = (-_Transform(0,2) * _Transform(0,2) - _Transform(1,2) * _Transform(1,2) - _Transform(2,2) * _Transform(2,2)) / pow(_shape_parameter + _radius*_radius, 0.5);
	return d + f;
}

double TPS::basis()
{
	radius();
	if (_radius != 0) return pow(_radius,4)*log(_radius);
	else return 0;
}

double TPS::dx_p1()
{
	radius();
	if (_radius != 0) return _x_delta*_radius*_radius + 4.0*_x_delta*_radius*_radius*log(_radius);
	else return 0;
}

double TPS::dx_p2()
{
	radius();
	if (_radius != 0) return -_x_delta*_radius*_radius - 4.0*_x_delta*_radius*_radius*log(_radius);
	else return 0;
}

double TPS::dy_p1()
{
	radius();
	if (_radius != 0) return _y_delta*_radius*_radius + 4.0*_y_delta*_radius*_radius*log(_radius);
	else return 0;
}

double TPS::dy_p2()
{
	radius();
	if (_radius != 0) return -_y_delta*_radius*_radius - 4.0*_y_delta*_radius*_radius*log(_radius);
	else return 0;
}

double TPS::dz_p1()
{
	radius();
	if (_radius != 0) return _z_delta*_radius*_radius + 4.0*_z_delta*_radius*_radius*log(_radius);
	else return 0;
}

double TPS::dz_p2()
{
	radius();
	if (_radius != 0) return -_z_delta*_radius*_radius - 4.0*_z_delta*_radius*_radius*log(_radius);
	else return 0;
}

double TPS::dxx()
{
	radius();
	if (_radius != 0) return -7.0*_x_delta*_x_delta - _y_delta*_y_delta - _z_delta*_z_delta - 8.0*_x_delta*_x_delta*log(_radius) - 4.0*_radius*_radius*log(_radius);
	else return 0;
}

double TPS::dxy()
{
	radius();
	if (_radius != 0) return -6.0*_x_delta*_y_delta - 8.0*_x_delta*_y_delta*log(_radius);
	else return 0;
}

double TPS::dxz()
{
	radius();
	if (_radius != 0) return -6.0*_x_delta*_z_delta - 8.0*_x_delta*_z_delta*log(_radius);
	else return 0;
}

double TPS::dyx()
{
	return dxy();
}

double TPS::dyy()
{
	radius();
	if (_radius != 0) return -7.0*_y_delta*_y_delta - _x_delta*_x_delta - _z_delta*_z_delta - 8.0*_y_delta*_y_delta*log(_radius) - 4.0*_radius*_radius*log(_radius);
	else return 0;
}

double TPS::dyz()
{
	radius();
	if (_radius != 0) return -6.0*_y_delta*_z_delta - 8.0*_y_delta*_z_delta*log(_radius);
	else return 0;
}

double TPS::dzx()
{
	return dxz();
}

double TPS::dzy()
{
	return dyz();
}

double TPS::dzz()
{
	radius();
	if (_radius != 0) return -7.0*_z_delta*_z_delta - _y_delta*_y_delta - _x_delta*_x_delta - 8.0*_z_delta*_z_delta*log(_radius) - 4.0*_radius*_radius*log(_radius);
	else return 0;
}

double ATPS::basis()
{
	scaled_radius();
	if (_radius != 0) return pow(_radius, 4)*log(_radius);
	else return 0;
}

double ATPS::dx_p1()
{
	scaled_radius();
	if (_radius != 0)
	{
		double a = (_Transform(0,0) * _x_delta + _Transform(1,0) * _y_delta + _Transform(2,0) * _z_delta)*_radius*_radius;
		double b = 4.0*a*log(_radius);
		return a + b;
	}
	else return 0;
}

double ATPS::dx_p2()
{
	scaled_radius();
	if (_radius != 0)
	{
		double a = -1.0*(_Transform(0,0) * _x_delta + _Transform(1,0) * _y_delta + _Transform(2,0) * _z_delta)*_radius*_radius;
		double b = 4.0*a*log(_radius);
		return a + b;
	}
	else return 0;
}

double ATPS::dy_p1()
{
	scaled_radius();
	if (_radius != 0)
	{
		double a = (_Transform(0,1) * _x_delta + _Transform(1,1) * _y_delta + _Transform(2,1) * _z_delta)*_radius*_radius;
		double b = 4.0*a*log(_radius);
		return a + b;
	}
	else return 0;
}

double ATPS::dy_p2()
{
	scaled_radius();
	if (_radius != 0)
	{
		double a = -1.0*(_Transform(0,1) * _x_delta + _Transform(1,1) * _y_delta + _Transform(2,1) * _z_delta)*_radius*_radius;
		double b = 4.0*a*log(_radius);
		return a + b;
	}
	else return 0;
}

double ATPS::dz_p1()
{
	scaled_radius();
	if (_radius != 0)
	{
		double a = (_Transform(0,2) * _x_delta + _Transform(1,2) * _y_delta + _Transform(2,2) * _z_delta)*_radius*_radius;
		double b = 4.0*a*log(_radius);
		return a + b;
	}
	else return 0;
}

double ATPS::dz_p2()
{
	scaled_radius();
	if (_radius != 0)
	{
		double a = -1.0*(_Transform(0,2) * _x_delta + _Transform(1,2) * _y_delta + _Transform(2,2) * _z_delta)*_radius*_radius;
		double b = 4.0*a*log(_radius);
		return a + b;
	}
	else return 0;
}

double ATPS::dxx()
{
	scaled_radius();
	if (_radius != 0)
	{
		double a = _Transform(0,0) * _x_delta + _Transform(1,0) * _y_delta + _Transform(2,0) * _z_delta;
		double b = -6.0*a*a;
		double c = -1.0*(_Transform(0,0) * _Transform(0,0) + _Transform(1,0) * _Transform(1,0) + _Transform(2,0) * _Transform(2,0))*_radius*_radius;
		double d = -8.0*a*a*log(_radius);
		double f = -4.0*(_Transform(0,0) * _Transform(0,0) + _Transform(1,0) * _Transform(1,0) + _Transform(2,0) * _Transform(2,0))*_radius*_radius*log(_radius);
		return b + c + d + f;
	}
	else return 0;
}

double ATPS::dxy()
{
	scaled_radius();
	if (_radius != 0)
	{
		double a = _Transform(0,0) * _x_delta + _Transform(1,0) * _y_delta + _Transform(2,0) * _z_delta;
		double b = _Transform(0,1) * _x_delta + _Transform(1,1) * _y_delta + _Transform(2,1) * _z_delta;
		double c = -6.0*a*b;
		double d = -1.0*(_Transform(0,0) * _Transform(0,1) + _Transform(1,0) * _Transform(1,1) + _Transform(2,0) * _Transform(2,1))*_radius*_radius;
		double f = -8.0*a*b*log(_radius);
		double g = -4.0*(_Transform(0,0) * _Transform(0,1) + _Transform(1,0) * _Transform(1,1) + _Transform(2,0) * _Transform(2,1))*_radius*_radius*log(_radius);
		return c + d + f + g;
	}
	else return 0;
}

double ATPS::dxz()
{
	scaled_radius();
	if (_radius != 0)
	{
		double a = _Transform(0,0) * _x_delta + _Transform(1,0) * _y_delta + _Transform(2,0) * _z_delta;
		double b = _Transform(0,2) * _x_delta + _Transform(1,2) * _y_delta + _Transform(2,2) * _z_delta;
		double c = -6.0*a*b;
		double d = -1.0*(_Transform(0,0) * _Transform(0,2) + _Transform(1,0) * _Transform(1,2) + _Transform(2,0) * _Transform(2,2))*_radius*_radius;
		double f = -8.0*a*b*log(_radius);
		double g = -4.0*(_Transform(0,0) * _Transform(0,2) + _Transform(1,0) * _Transform(1,2) + _Transform(2,0) * _Transform(2,2))*_radius*_radius*log(_radius);
		return c + d + f + g;
	}
	else return 0;
}

double ATPS::dyx()
{
	return dxy();
}

double ATPS::dyy()
{
	scaled_radius();
	if (_radius != 0)
	{
		double a = _Transform(0,1) * _x_delta + _Transform(1,1) * _y_delta + _Transform(2,1) * _z_delta;
		double b = -6.0*a*a;
		double c = -1.0*(_Transform(0,1) * _Transform(0,1) + _Transform(1,1) * _Transform(1,1) + _Transform(2,1) * _Transform(2,1))*_radius*_radius;
		double d = -8.0*a*a*log(_radius);
		double f = -4.0*(_Transform(0,1) * _Transform(0,1) + _Transform(1,1) * _Transform(1,1) + _Transform(2,1) * _Transform(2,1))*_radius*_radius*log(_radius);
		return b + c + d + f;
	}
	else return 0;
}

double ATPS::dyz()
{
	scaled_radius();
	if (_radius != 0)
	{
		double a = _Transform(0,1) * _x_delta + _Transform(1,1) * _y_delta + _Transform(2,1) * _z_delta;
		double b = _Transform(0,2) * _x_delta + _Transform(1,2) * _y_delta + _Transform(2,2) * _z_delta;
		double c = -6.0*a*b;
		double d = -1.0*(_Transform(0,1) * _Transform(0,2) + _Transform(1,1) * _Transform(1,2) + _Transform(2,1) * _Transform(2,2))*_radius*_radius;
		double f = -8.0*a*b*log(_radius);
		double g = -4.0*(_Transform(0,1) * _Transform(0,2) + _Transform(1,1) * _Transform(1,2) + _Transform(2,1) * _Transform(2,2))*_radius*_radius*log(_radius);
		return c + d + f + g;
	}
	else return 0;
}

double ATPS::dzx()
{
	return dxz();
}

double ATPS::dzy()
{
	return dyz();
}

double ATPS::dzz()
{
	scaled_radius();
	if (_radius != 0)
	{
		double a = _Transform(0,2) * _x_delta + _Transform(1,2) * _y_delta + _Transform(2,2) * _z_delta;
		double b = -6.0*a*a;
		double c = -1.0*(_Transform(0,2) * _Transform(0,2) + _Transform(1,2) * _Transform(1,2) + _Transform(2,2) * _Transform(2,2))*_radius*_radius;
		double d = -8.0*a*a*log(_radius);
		double f = -4.0*(_Transform(0,2) * _Transform(0,2) + _Transform(1,2) * _Transform(1,2) + _Transform(2,2) * _Transform(2,2))*_radius*_radius*log(_radius);
		return b + c + d + f;
	}
	else return 0;
}

double IMQ::basis()
{
	radius();
	return 1.0 / pow(_shape_parameter + _radius*_radius, 0.5);
}

double IMQ::dx_p1()
{
	radius();
	return -1.0*_x_delta/ pow(_shape_parameter + _radius*_radius, 1.5);
}

double IMQ::dx_p2()
{
	radius();
	return _x_delta / pow(_shape_parameter + _radius*_radius, 1.5);
}

double IMQ::dy_p1()
{
	radius();
	return -1.0*_y_delta / pow(_shape_parameter + _radius*_radius, 1.5);
}

double IMQ::dy_p2()
{
	radius();
	return _y_delta / pow(_shape_parameter + _radius*_radius, 1.5);
}

double IMQ::dz_p1()
{
	radius();
	return -1.0*_z_delta / pow(_shape_parameter + _radius*_radius, 1.5);
}

double IMQ::dz_p2()
{
	radius();
	return _z_delta / pow(_shape_parameter + _radius*_radius, 1.5);
}

double IMQ::dxx()
{
	radius();
	return (-3.0*_x_delta*_x_delta / pow(_shape_parameter + _radius*_radius, 2.5)) + 1.0 / pow(_shape_parameter + _radius*_radius, 1.5);
}

double IMQ::dxy()
{
	radius();
	return -3.0*_x_delta*_y_delta/ pow(_shape_parameter + _radius*_radius, 2.5);
}

double IMQ::dxz()
{
	radius();
	return -3.0*_x_delta*_z_delta / pow(_shape_parameter + _radius*_radius, 2.5);
}

double IMQ::dyx()
{
	return dxy();
}

double IMQ::dyy()
{
	radius();
	return (-3.0*_y_delta*_y_delta / pow(_shape_parameter + _radius*_radius, 2.5)) + 1.0 / pow(_shape_parameter + _radius*_radius, 1.5);
}

double IMQ::dyz()
{
	radius();
	return -3.0*_y_delta*_z_delta / pow(_shape_parameter + _radius*_radius, 2.5);
}

double IMQ::dzx()
{
	return dxz();
}

double IMQ::dzy()
{
	return dyz();
}

double IMQ::dzz()
{
	radius();
	return (-3.0*_z_delta*_z_delta / pow(_shape_parameter + _radius*_radius, 2.5)) + 1.0 / pow(_shape_parameter + _radius*_radius, 1.5);
}

double AIMQ::basis()
{
	scaled_radius();
	return 1 / pow(_shape_parameter + _radius*_radius, 0.5);
}

double AIMQ::dx_p1()
{
	scaled_radius();
	return -1.0*(_Transform(0,0) * _x_delta + _Transform(1,0) * _y_delta + _Transform(2,0) * _z_delta) / pow(_shape_parameter + _radius*_radius, 1.5);
}

double AIMQ::dx_p2()
{
	scaled_radius();
	return (_Transform(0,0) * _x_delta + _Transform(1,0) * _y_delta + _Transform(2,0) * _z_delta) / pow(_shape_parameter + _radius*_radius, 1.5);
}

double AIMQ::dy_p1()
{
	scaled_radius();
	return -1.0*(_Transform(0,1) * _x_delta + _Transform(1,1) * _y_delta + _Transform(2,1) * _z_delta) / pow(_shape_parameter + _radius*_radius, 1.5);
}

double AIMQ::dy_p2()
{
	scaled_radius();
	return (_Transform(0,1) * _x_delta + _Transform(1,1) * _y_delta + _Transform(2,1) * _z_delta) / pow(_shape_parameter + _radius*_radius, 1.5);
}

double AIMQ::dz_p1()
{
	scaled_radius();
	return -1.0*(_Transform(0,2) * _x_delta + _Transform(1,2) * _y_delta + _Transform(2,2) * _z_delta) / pow(_shape_parameter + _radius*_radius, 1.5);
}
double AIMQ::dz_p2()
{
	scaled_radius();
	return (_Transform(0,2) * _x_delta + _Transform(1,2) * _y_delta + _Transform(2,2) * _z_delta) / pow(_shape_parameter + _radius*_radius, 1.5);
}

double AIMQ::dxx()
{
	scaled_radius();
	double a = _Transform(0,0) * _x_delta + _Transform(1,0) * _y_delta + _Transform(2,0) * _z_delta;
	double b = -3.0*a*a / pow(_shape_parameter + _radius*_radius, 2.5);
	double c = -1.0*(_Transform(0,0) * _Transform(0,0) + _Transform(1,0) * _Transform(1,0) + _Transform(2,0) * _Transform(2,0)) / pow(_shape_parameter + _radius*_radius, 1.5);
	return b - c;
}

double AIMQ::dxy()
{
	scaled_radius();
	double a = _Transform(0,0) * _x_delta + _Transform(1,0) * _y_delta + _Transform(2,0) * _z_delta;
	double b = _Transform(0,1) * _x_delta + _Transform(1,1) * _y_delta + _Transform(2,1) * _z_delta;
	double c = -3.0*a*b / pow(_shape_parameter + _radius*_radius, 2.5);
	double d = -1.0*(_Transform(0,0) * _Transform(0,1) + _Transform(1,0) * _Transform(1,1) + _Transform(2,0) * _Transform(2,1)) / pow(_shape_parameter + _radius*_radius, 1.5);
	return c - d;
}

double AIMQ::dxz()
{
	scaled_radius();
	double a = _Transform(0,0) * _x_delta + _Transform(1,0) * _y_delta + _Transform(2,0) * _z_delta;
	double b = _Transform(0,2) * _x_delta + _Transform(1,2) * _y_delta + _Transform(2,2) * _z_delta;
	double c = -3.0*a*b / pow(_shape_parameter + _radius*_radius, 2.5);
	double d = -1.0*(_Transform(0,0) * _Transform(0,2) + _Transform(1,0) * _Transform(1,2) + _Transform(2,0) * _Transform(2,2)) / pow(_shape_parameter + _radius*_radius, 1.5);
	return c - d;
}

double AIMQ::dyx()
{
	return dxy();
}

double AIMQ::dyy()
{
	scaled_radius();
	double a = _Transform(0,1) * _x_delta + _Transform(1,1) * _y_delta + _Transform(2,1) * _z_delta;
	double b = -3.0*a*a / pow(_shape_parameter + _radius*_radius, 2.5);
	double c = -1.0*(_Transform(0,1) * _Transform(0,1) + _Transform(1,1) * _Transform(1,1) + _Transform(2,1) * _Transform(2,1)) / pow(_shape_parameter + _radius*_radius, 1.5);
	return b - c;
}

double AIMQ::dyz()
{
	scaled_radius();
	double a = _Transform(0,1) * _x_delta + _Transform(1,1) * _y_delta + _Transform(2,1) * _z_delta;
	double b = _Transform(0,2) * _x_delta + _Transform(1,2) * _y_delta + _Transform(2,2) * _z_delta;
	double c = -3.0*a*b / pow(_shape_parameter + _radius*_radius, 2.5);
	double d = -1.0*(_Transform(0,1) * _Transform(0,2) + _Transform(1,1) * _Transform(1,2) + _Transform(2,1) * _Transform(2,2)) / pow(_shape_parameter + _radius*_radius, 1.5);
	return c - d;
}

double AIMQ::dzx()
{
	return dxz();
}

double AIMQ::dzy()
{
	return dyz();
}

double AIMQ::dzz()
{
	scaled_radius();
	double a = _Transform(0,2) * _x_delta + _Transform(1,2) * _y_delta + _Transform(2,2) * _z_delta;
	double b = -3.0*a*a / pow(_shape_parameter + _radius*_radius, 2.5);
	double c = -1.0*(_Transform(0,2) * _Transform(0,2) + _Transform(1,2) * _Transform(1,2) + _Transform(2,2) * _Transform(2,2)) / pow(_shape_parameter + _radius*_radius, 1.5);
	return b - c;
}

double R::basis()
{
	radius();
	return _radius;
}

double R::dx_p1()
{
	throw - 666;
	return 0;
}

double R::dx_p2()
{
	throw - 666;
	return 0;
}

double R::dy_p1()
{
	throw - 666;
	return 0;
}

double R::dy_p2()
{
	throw - 666;
	return 0;
}

double R::dz_p1()
{
	throw - 666;
	return 0;
}

double R::dz_p2()
{
	throw - 666;
	return 0;
}

double R::dxx()
{
	throw - 666;
	return 0;
}

double R::dxy()
{
	throw - 666;
	return 0;
}

double R::dxz()
{
	throw - 666;
	return 0;
}

double R::dyx()
{
	return dxy();
}

double R::dyy()
{
	throw - 666;
	return 0;
}

double R::dyz()
{
	throw - 666;
	return 0;
}

double R::dzx()
{
	return dxz();
}

double R::dzy()
{
	return dyz();
}

double R::dzz()
{
	throw - 666;
	return 0;
}

double AR::basis()
{
	scaled_radius();
	return _radius;
}

double AR::dx_p1()
{
	throw - 666;
	return 0;
}

double AR::dx_p2()
{
	throw - 666;
	return 0;
}

double AR::dy_p1()
{
	throw - 666;
	return 0;
}

double AR::dy_p2()
{
	throw - 666;
	return 0;
}

double AR::dz_p1()
{
	throw - 666;
	return 0;
}

double AR::dz_p2()
{
	throw - 666;
	return 0;
}

double AR::dxx()
{
	throw - 666;
	return 0;
}

double AR::dxy()
{
	throw - 666;
	return 0;
}

double AR::dxz()
{
	throw - 666;
	return 0;
}

double AR::dyx()
{
	return dxy();
}

double AR::dyy()
{
	throw - 666;
	return 0;
}

double AR::dyz()
{
	throw - 666;
	return 0;
}

double AR::dzx()
{
	return dxz();
}

double AR::dzy()
{
	return dyz();
}

double AR::dzz()
{
	throw - 666;
	return 0;
}

bool Lagrangian_Polynomial_Basis::_get_unisolvent_subset(const std::vector < std::vector < Interface > > &interface_point_lists)
{
	// NOTE : Currently only supporting 1st order polynomials ( have p(x) = A*x + B*y + C*z + D )
	// Tried implementing 2nd order however it is not practical
	// e.g. finding the algrebraic solution in mathematical (via Thomas algo) resulted in ~ 100 page equation.  

	// use interface_point_lists data
	if ((int)interface_point_lists.size() == 0) return false;

	// find horizon with largest number of points
	int index = 0;
	for (int j = 1; j < (int)interface_point_lists.size(); j++) {
		if ((int)interface_point_lists[j].size() >(int)interface_point_lists[index].size()) index = j;
	}
	if ((int)interface_point_lists[index].size() < 4) return false; // not enough points to create the 1st order Lagrangian Polynomial Basis

	int n = (int)interface_point_lists[index].size();

	std::vector<double> Xcoord_array;
	std::vector<int> Index_Xcoord_array;
	std::vector<double> Ycoord_array;
	std::vector<int> Index_Ycoord_array; 
	std::vector<double> Zcoord_array;
	std::vector<int> Index_Zcoord_array;

	for (int j = 0; j < n; j++) {
		Index_Xcoord_array.push_back(j);
		Index_Ycoord_array.push_back(j);
		Index_Zcoord_array.push_back(j);
		Xcoord_array.push_back(interface_point_lists[index][j].x());
		Ycoord_array.push_back(interface_point_lists[index][j].y());
		Zcoord_array.push_back(interface_point_lists[index][j].z());
	}
	// Sort those arrays
	Math_methods::sort_vector_w_index<double>(Xcoord_array, Index_Xcoord_array);
	Math_methods::sort_vector_w_index<double>(Ycoord_array, Index_Ycoord_array);
	Math_methods::sort_vector_w_index<double>(Zcoord_array, Index_Zcoord_array);

	// 	cout<<" sorted x array"<<endl;
	// 	for (int j = 0; j < n; j++ ) cout<<" ["<<j<<"]= "<<Xcoord_array[j]<<" index= "<<Index_Xcoord_array[j]<<endl;
	// 	cout<<" sorted y array"<<endl;
	// 	for (int j = 0; j < n; j++ ) cout<<" ["<<j<<"]= "<<Ycoord_array[j]<<" index= "<<Index_Ycoord_array[j]<<endl;
	// 	cout<<" sorted z array"<<endl;
	// 	for (int j = 0; j < n; j++ ) cout<<" ["<<j<<"]= "<<Zcoord_array[j]<<" index= "<<Index_Zcoord_array[j]<<endl;

	double dx = Xcoord_array[n - 1] - Xcoord_array[0];
	double dy = Ycoord_array[n - 1] - Ycoord_array[0];
	double dz = Zcoord_array[n - 1] - Zcoord_array[0];

	/*	cout<<" dx = "<<dx<<" dy= "<<dy<<" dz= "<<dz<<endl;*/

	// Find the axis that has the largest sampling
	if (dx >= dy && dx >= dz)
	{
		unisolvent_subset_points.push_back(interface_point_lists[index][Index_Xcoord_array[0]]);
		unisolvent_subset_points.push_back(interface_point_lists[index][Index_Xcoord_array[n - 1]]);
		if (dy >= dz)
		{
			for (int j = 0; j < n; j++) {
				if (Index_Ycoord_array[j] != Index_Xcoord_array[0] && Index_Ycoord_array[j] != Index_Xcoord_array[n - 1])
				{
					unisolvent_subset_points.push_back(interface_point_lists[index][Index_Ycoord_array[j]]);
					break;
				}
			}
			for (int j = 0; j < n; j++) {
				if (Index_Ycoord_array[n - 1 - j] != Index_Xcoord_array[0] && Index_Ycoord_array[n - 1 - j] != Index_Xcoord_array[n - 1])
				{
					unisolvent_subset_points.push_back(interface_point_lists[index][Index_Ycoord_array[n - 1 - j]]);
					break;
				}
			}
		}
		else
		{
			for (int j = 0; j < n; j++) {
				if (Index_Zcoord_array[j] != Index_Xcoord_array[0] && Index_Zcoord_array[j] != Index_Xcoord_array[n - 1])
				{
					unisolvent_subset_points.push_back(interface_point_lists[index][Index_Zcoord_array[j]]);
					break;
				}
			}
			for (int j = 0; j < n; j++) {
				if (Index_Zcoord_array[n - 1 - j] != Index_Xcoord_array[0] && Index_Zcoord_array[n - 1 - j] != Index_Xcoord_array[n - 1])
				{
					unisolvent_subset_points.push_back(interface_point_lists[index][Index_Zcoord_array[n - 1 - j]]);
					break;
				}
			}
		}
	}

	if (dy >= dx && dy >= dz)
	{
		unisolvent_subset_points.push_back(interface_point_lists[index][Index_Ycoord_array[0]]);
		unisolvent_subset_points.push_back(interface_point_lists[index][Index_Ycoord_array[n - 1]]);
		if (dx >= dz)
		{
			for (int j = 0; j < n; j++) {
				if (Index_Xcoord_array[j] != Index_Ycoord_array[0] && Index_Xcoord_array[j] != Index_Ycoord_array[n - 1])
				{
					unisolvent_subset_points.push_back(interface_point_lists[index][Index_Xcoord_array[j]]);
					break;
				}
			}
			for (int j = 0; j < n; j++) {
				if (Index_Xcoord_array[n - 1 - j] != Index_Ycoord_array[0] && Index_Xcoord_array[n - 1 - j] != Index_Ycoord_array[n - 1])
				{
					unisolvent_subset_points.push_back(interface_point_lists[index][Index_Xcoord_array[n - 1 - j]]);
					break;
				}
			}
		}
		else
		{
			for (int j = 0; j < n; j++) {
				if (Index_Zcoord_array[j] != Index_Ycoord_array[0] && Index_Zcoord_array[j] != Index_Ycoord_array[n - 1])
				{
					unisolvent_subset_points.push_back(interface_point_lists[index][Index_Zcoord_array[j]]);
					break;
				}
			}
			for (int j = 0; j < n; j++) {
				if (Index_Zcoord_array[n - 1 - j] != Index_Ycoord_array[0] && Index_Zcoord_array[n - 1 - j] != Index_Ycoord_array[n - 1])
				{
					unisolvent_subset_points.push_back(interface_point_lists[index][Index_Zcoord_array[n - 1 - j]]);
					break;
				}
			}
		}
	}

	if (dz >= dx && dz >= dy)
	{
		unisolvent_subset_points.push_back(interface_point_lists[index][Index_Zcoord_array[0]]);
		unisolvent_subset_points.push_back(interface_point_lists[index][Index_Zcoord_array[n - 1]]);
		if (dx >= dy)
		{
			for (int j = 0; j < n; j++) {
				if (Index_Xcoord_array[j] != Index_Zcoord_array[0] && Index_Xcoord_array[j] != Index_Zcoord_array[n - 1])
				{
					unisolvent_subset_points.push_back(interface_point_lists[index][Index_Xcoord_array[j]]);
					break;
				}
			}
			for (int j = 0; j < n; j++) {
				if (Index_Xcoord_array[n - 1 - j] != Index_Zcoord_array[0] && Index_Xcoord_array[n - 1 - j] != Index_Zcoord_array[n - 1])
				{
					unisolvent_subset_points.push_back(interface_point_lists[index][Index_Xcoord_array[n - 1 - j]]);
					break;
				}
			}
		}
		else
		{
			for (int j = 0; j < n; j++) {
				if (Index_Ycoord_array[j] != Index_Zcoord_array[0] && Index_Ycoord_array[j] != Index_Zcoord_array[n - 1])
				{
					unisolvent_subset_points.push_back(interface_point_lists[index][Index_Ycoord_array[j]]);
					break;
				}
			}
			for (int j = 0; j < n; j++) {
				if (Index_Ycoord_array[n - 1 - j] != Index_Zcoord_array[0] && Index_Ycoord_array[n - 1 - j] != Index_Zcoord_array[n - 1])
				{
					unisolvent_subset_points.push_back(interface_point_lists[index][Index_Ycoord_array[n - 1 - j]]);
					break;
				}
			}
		}
	}

	// check if Subset list exists on a plane
	bool exists_on_plane = false;
	int add_x = 0;
	int add_y = 0;
	int add_z = 0;
	int Q = (int)unisolvent_subset_points.size();
	for (int j = 1; j < Q; j++) {
		// compare each coord of [0] with the same coord of the other pts []
		// if coord is same as all other pts, pts lie on a plane
		if (unisolvent_subset_points[j].x() == unisolvent_subset_points[0].x()) add_x++;
		if (unisolvent_subset_points[j].y() == unisolvent_subset_points[0].y()) add_y++;
		if (unisolvent_subset_points[j].z() == unisolvent_subset_points[0].z()) add_z++;
	}

	if (add_x == (Q - 1) || add_y == (Q - 1) || add_z == (Q - 1)) exists_on_plane = true;
	if (exists_on_plane)
	{
		bool found_unique_point = false;
		if (add_x == (Q - 1)) // search for a pt in pointset with a different x
		{
			for (int j = 0; j < (int)interface_point_lists[index].size(); j++) {
				if (interface_point_lists[index][j].x() != unisolvent_subset_points[0].x())
				{
					unisolvent_subset_points.pop_back();
					unisolvent_subset_points.push_back(interface_point_lists[index][j]);
					found_unique_point = true;
					break;
				}
			}
		}
		if (add_y == (Q - 1)) // search for a pt in pointset with a different y
		{
			for (int j = 0; j < (int)interface_point_lists[index].size(); j++) {
				if (interface_point_lists[index][j].y() != unisolvent_subset_points[0].y())
				{
					unisolvent_subset_points.pop_back();
					unisolvent_subset_points.push_back(interface_point_lists[index][j]);
					found_unique_point = true;
					break;
				}
			}
		}
		if (add_z == (Q - 1)) // search for a pt in pointset with a different z
		{
			for (int j = 0; j < (int)interface_point_lists[index].size(); j++) {
				if (interface_point_lists[index][j].z() != unisolvent_subset_points[0].z())
				{
					unisolvent_subset_points.pop_back();
					unisolvent_subset_points.push_back(interface_point_lists[index][j]);
					found_unique_point = true;
					break;
				}
			}
		}
		if (!found_unique_point)
		{
			if (add_x == (Q - 1)) unisolvent_subset_points[0].set_x(unisolvent_subset_points[0].x() + 0.0001);
			if (add_y == (Q - 1)) unisolvent_subset_points[0].set_y(unisolvent_subset_points[0].y() + 0.0001);
			if (add_z == (Q - 1)) unisolvent_subset_points[0].set_z(unisolvent_subset_points[0].z() + 0.0001);
		}
	}

	if ((int)unisolvent_subset_points.size() == 4) return true;
	else return false;
}

void Lagrangian_Polynomial_Basis::_initialize_basis()
{
	mpf_class x1 = unisolvent_subset_points[0].x();
	mpf_class y1 = unisolvent_subset_points[0].y();
	mpf_class z1 = unisolvent_subset_points[0].z();

	mpf_class x2 = unisolvent_subset_points[1].x();
	mpf_class y2 = unisolvent_subset_points[1].y();
	mpf_class z2 = unisolvent_subset_points[1].z();

	mpf_class x3 = unisolvent_subset_points[2].x();
	mpf_class y3 = unisolvent_subset_points[2].y();
	mpf_class z3 = unisolvent_subset_points[2].z();

	mpf_class x4 = unisolvent_subset_points[3].x();
	mpf_class y4 = unisolvent_subset_points[3].y();
	mpf_class z4 = unisolvent_subset_points[3].z();

	mpf_class d = (x2*y3*z1 - x2*y4*z1 - x1*y3*z2 + x1*y4*z2 - x2*y1*z3 + x1*y2*z3 - x1*y4*z3 + x2*y4*z3 + 
		x4*(y2*z1 - y3*z1 - y1*z2 + y3*z2 + y1*z3 - y2*z3) + x2*y1*z4 - x1*y2*z4 + x1*y3*z4 - x2*y3*z4 + 
		x3*(y4*z1 + y1*z2 - y4*z2 - y1*z4 + y2*(-z1 + z4)));

	_polynomial_constants.resize(16);

	// Coefficents for lagrange polynomial basis 1 
	_polynomial_constants(0) = ( (x4*y3*z2 - x3*y4*z2 - x4*y2*z3 + x2*y4*z3 + x3*y2*z4 - x2*y3*z4)/d ); // constant coef
	_polynomial_constants(1) = ( (-(y3*z2) + y4*z2 + y2*z3 - y4*z3 - y2*z4 + y3*z4)/d ); // x coef
	_polynomial_constants(2) = ( (x3*z2 - x4*z2 - x2*z3 + x4*z3 + x2*z4 - x3*z4)/d ); // y coef
	_polynomial_constants(3) = ( (-(x3*y2) + x4*y2 + x2*y3 - x4*y3 - x2*y4 + x3*y4)/d ); // z coef

	// Coefficents for lagrange polynomial basis 2  
	_polynomial_constants(4) = ( (-(x4*y3*z1) + x3*y4*z1 + x4*y1*z3 - x1*y4*z3 - x3*y1*z4 + x1*y3*z4)/d ); // constant coef
	_polynomial_constants(5) = ( (y3*z1 - y4*z1 - y1*z3 + y4*z3 + y1*z4 - y3*z4)/d ); // x coef
	_polynomial_constants(6) = ( (-(x3*z1) + x4*z1 + x1*z3 - x4*z3 - x1*z4 + x3*z4)/d ); // y coef
	_polynomial_constants(7) = ( (x3*y1 - x4*y1 - x1*y3 + x4*y3 + x1*y4 - x3*y4)/d ); // z coef

	// Coefficents for lagrange polynomial basis 3 
	_polynomial_constants(8) = ( (x4*y2*z1 - x2*y4*z1 - x4*y1*z2 + x1*y4*z2 + x2*y1*z4 - x1*y2*z4)/d ); // constant coef
	_polynomial_constants(9) = ( (-(y2*z1) + y4*z1 + y1*z2 - y4*z2 - y1*z4 + y2*z4)/d ); // x coef
	_polynomial_constants(10) = ( (x2*z1 - x4*z1 - x1*z2 + x4*z2 + x1*z4 - x2*z4)/d ); // y coef
	_polynomial_constants(11) = ( (-(x2*y1) + x4*y1 + x1*y2 - x4*y2 - x1*y4 + x2*y4)/d ); // z coef

	// Coefficents for lagrange polynomial basis 4
	_polynomial_constants(12) = ( (-(x3*y2*z1) + x2*y3*z1 + x3*y1*z2 - x1*y3*z2 - x2*y1*z3 + x1*y2*z3)/d ); // constant coef
	_polynomial_constants(13) = ( (y2*z1 - y3*z1 - y1*z2 + y3*z2 + y1*z3 - y2*z3)/d ); // x coef
	_polynomial_constants(14) = ( (-(x2*z1) + x3*z1 + x1*z2 - x3*z2 - x1*z3 + x2*z3)/d ); // y coef
	_polynomial_constants(15) = ( (x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3)/d ); // z coef

	//std::vector<double> poly_const_double;
	//for (int j = 0; j< 16;j++) poly_const_double.push_back(_polynomial_constants[j].get_d());

	_derivative_polynomial_constants.resize(3,4);

	_derivative_polynomial_constants(0,0) = _polynomial_constants[1]; // basis 1
	_derivative_polynomial_constants(0,1) = _polynomial_constants[5]; // basis 2
	_derivative_polynomial_constants(0,2) = _polynomial_constants[9]; // basis 3
	_derivative_polynomial_constants(0,3) = _polynomial_constants[13]; // basis 4

	_derivative_polynomial_constants(1,0) = _polynomial_constants[2]; // basis 1 
	_derivative_polynomial_constants(1,1) = _polynomial_constants[6]; // basis 2 
	_derivative_polynomial_constants(1,2) = _polynomial_constants[10]; // basis 3 
	_derivative_polynomial_constants(1,3) = _polynomial_constants[14]; // basis 4

	_derivative_polynomial_constants(2,0) = _polynomial_constants[3]; // basis 1
	_derivative_polynomial_constants(2,1) = _polynomial_constants[7]; // basis 2 
	_derivative_polynomial_constants(2,2) = _polynomial_constants[11]; // basis 3
	_derivative_polynomial_constants(2,3) = _polynomial_constants[15]; // basis 4
}

Matrix <mpf_class, Dynamic, 1> Lagrangian_Polynomial_Basis::poly(const Point *p)
{
	Matrix <mpf_class, Dynamic, 1> basis;
	basis.resize(4);

	basis(0) = _polynomial_constants(0)  + _polynomial_constants(1)*p->x()  + _polynomial_constants(2)*p->y()  + _polynomial_constants(3)*p->z();
	basis(1) = _polynomial_constants(4)  + _polynomial_constants(5)*p->x()  + _polynomial_constants(6)*p->y()  + _polynomial_constants(7)*p->z();
	basis(2) = _polynomial_constants(8)  + _polynomial_constants(9)*p->x()  + _polynomial_constants(10)*p->y() + _polynomial_constants(11)*p->z();
	basis(3) = _polynomial_constants(12) + _polynomial_constants(13)*p->x() + _polynomial_constants(14)*p->y() + _polynomial_constants(15)*p->z();

	return basis;
}

Matrix <mpf_class, Dynamic, 1> Lagrangian_Polynomial_Basis::poly_dx( const Point *p )
{
	Matrix <mpf_class, Dynamic, 1> basis;
	basis.resize(4);

	basis(0) = _derivative_polynomial_constants(0,0);
	basis(1) = _derivative_polynomial_constants(0,1);
	basis(2) = _derivative_polynomial_constants(0,2);
	basis(3) = _derivative_polynomial_constants(0,3);

	return basis;
}

Matrix <mpf_class, Dynamic, 1> Lagrangian_Polynomial_Basis::poly_dy( const Point *p )
{
	Matrix <mpf_class, Dynamic, 1> basis;
	basis.resize(4);

	basis(0) = _derivative_polynomial_constants(1,0);
	basis(1) = _derivative_polynomial_constants(1,1);
	basis(2) = _derivative_polynomial_constants(1,2);
	basis(3) = _derivative_polynomial_constants(1,3);

	return basis;
}

Matrix <mpf_class, Dynamic, 1> Lagrangian_Polynomial_Basis::poly_dz( const Point *p )
{
	Matrix <mpf_class, Dynamic, 1> basis;
	basis.resize(4);

	basis(0) = _derivative_polynomial_constants(2,0);
	basis(1) = _derivative_polynomial_constants(2,1);
	basis(2) = _derivative_polynomial_constants(2,2);
	basis(3) = _derivative_polynomial_constants(2,3);

	return basis;
}

double Modified_Kernel::basis_pt_pt()
{
	mpf_class t1,t2,t3,t4;

	Matrix <mpf_class, Dynamic, 1> p1 = this->_aLPB->poly(this->p1());
	Matrix <mpf_class, Dynamic, 1> p2 = this->_aLPB->poly(this->p2());
	
	for (int j = 0; j < 4; j++){
		this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],*this->p2());
		mpf_class b1(this->_aRBFKernel->basis());
		this->_aRBFKernel->set_points(*this->p1(),this->_aLPB->unisolvent_subset_points[j]);
		mpf_class b2(this->_aRBFKernel->basis());
		t1 += p1(j) * b1;
		t2 += p2(j) * b2;
		t3 += p1(j) * p2(j);
		for (int k = 0; k < 4; k++ ){
 			if (k != j)
 			{
				this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],this->_aLPB->unisolvent_subset_points[k]);
				mpf_class b3(this->_aRBFKernel->basis());
				t4 += p1(j) * p2(k) * b3;
			}
		}
	}
	this->_aRBFKernel->set_points(*this->_p1,*this->_p2);
	mpf_class sum = this->_aRBFKernel->basis() - t1 - t2 + t3 + t4;
	return sum.get_d();
}

double Modified_Kernel::basis_pt_planar_x()
{
	mpf_class t1x,t2x,t3x,t4x;

	Matrix <mpf_class, Dynamic, 1> p1 = this->_aLPB->poly(this->p1());
	Matrix <mpf_class, Dynamic, 1> p2x = this->_aLPB->poly_dx(this->p2());

	for (int j = 0; j < 4; j++){
		this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],*this->p2());
		mpf_class b1x(this->_aRBFKernel->dx_p2());
		this->_aRBFKernel->set_points(*this->p1(),this->_aLPB->unisolvent_subset_points[j]);
		mpf_class b2(this->_aRBFKernel->basis());
		t1x += p1(j)*b1x;  
		t2x += p2x(j)*b2;
		t3x += p1(j)*p2x(j);
		for (int k = 0; k < 4; k++ ){
			if ( k!=j )
			{
				this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],this->_aLPB->unisolvent_subset_points[k]);
				mpf_class b3(this->_aRBFKernel->basis());
				t4x += p1(j) * p2x(k) * b3;
			}
		}
	}

	this->_aRBFKernel->set_points(*this->_p1,*this->_p2);
	mpf_class sum = this->_aRBFKernel->dx_p2() - t1x - t2x + t3x + t4x;
	return sum.get_d();
}

double Modified_Kernel::basis_planar_x_pt()
{
	mpf_class t1x,t2x,t3x,t4x;

	Matrix <mpf_class, Dynamic, 1> p1x = this->_aLPB->poly_dx(this->p1());
	Matrix <mpf_class, Dynamic, 1> p2 = this->_aLPB->poly(this->p2());

	for (int j = 0; j < 4; j++){
		this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],*this->p2());
		mpf_class b1(this->_aRBFKernel->basis());
		this->_aRBFKernel->set_points(*this->p1(),this->_aLPB->unisolvent_subset_points[j]);
		mpf_class b2x(this->_aRBFKernel->dx_p1());
		t1x += p1x(j)*b1;  
		t2x += p2(j)*b2x;
		t3x += p1x(j)*p2(j);
		for (int k = 0; k < 4; k++ ){
			if (k != j)
			{
				this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],this->_aLPB->unisolvent_subset_points[k]);
				mpf_class b3(this->_aRBFKernel->basis());
				t4x += p1x(j) * p2(k) * b3;
			}
		}
	}

	this->_aRBFKernel->set_points(*this->_p1,*this->_p2);
	mpf_class sum = this->_aRBFKernel->dx_p1() - t1x - t2x + t3x + t4x;
	return sum.get_d();
}

double Modified_Kernel::basis_pt_planar_y()
{
	mpf_class t1y,t2y,t3y,t4y;

	Matrix <mpf_class, Dynamic, 1> p1 = this->_aLPB->poly(this->p1());
	Matrix <mpf_class, Dynamic, 1> p2y = this->_aLPB->poly_dy(this->p2());

	for (int j = 0; j < 4; j++){
		this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],*this->p2());
		mpf_class b1y(this->_aRBFKernel->dy_p2());
		this->_aRBFKernel->set_points(*this->p1(),this->_aLPB->unisolvent_subset_points[j]);
		mpf_class b2(this->_aRBFKernel->basis());
		t1y += p1(j)*b1y;  
		t2y += p2y(j)*b2;
		t3y += p1(j)*p2y(j);
		for (int k = 0; k < 4; k++ ){
			if (k != j)
			{
				this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],this->_aLPB->unisolvent_subset_points[k]);
				mpf_class b3(this->_aRBFKernel->basis());
				t4y += p1(j) * p2y(k) * b3;
			}
		}
	}

	this->_aRBFKernel->set_points(*this->_p1,*this->_p2);
	mpf_class sum = this->_aRBFKernel->dy_p2() - t1y - t2y + t3y + t4y;
	return sum.get_d();
}

double Modified_Kernel::basis_planar_y_pt()
{
	mpf_class t1y,t2y,t3y,t4y;

	Matrix <mpf_class, Dynamic, 1> p1y = this->_aLPB->poly_dy(this->p1());
	Matrix <mpf_class, Dynamic, 1> p2 = this->_aLPB->poly(this->p2());

	for (int j = 0; j < 4; j++){
		this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],*this->p2());
		mpf_class b1(this->_aRBFKernel->basis());
		this->_aRBFKernel->set_points(*this->p1(),this->_aLPB->unisolvent_subset_points[j]);
		mpf_class b2y(this->_aRBFKernel->dy_p1());
		t1y += p1y(j)*b1;  
		t2y += p2(j)*b2y;
		t3y += p1y(j)*p2(j);
		for (int k = 0; k < 4; k++ ){
			if (k != j)
			{
				this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],this->_aLPB->unisolvent_subset_points[k]);
				mpf_class b3(this->_aRBFKernel->basis());
				t4y += p1y(j) * p2(k) * b3;
			}
		}
	}
	
	this->_aRBFKernel->set_points(*this->_p1,*this->_p2);
	mpf_class sum = this->_aRBFKernel->dy_p1() - t1y - t2y + t3y + t4y;
	return sum.get_d();
}

double Modified_Kernel::basis_pt_planar_z()
{
	mpf_class t1z,t2z,t3z,t4z;

	Matrix <mpf_class, Dynamic, 1> p1 = this->_aLPB->poly(this->p1());
	Matrix <mpf_class, Dynamic, 1> p2z = this->_aLPB->poly_dz(this->p2());

	for (int j = 0; j < 4; j++){
		this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],*this->p2());
		mpf_class b1z(this->_aRBFKernel->dz_p2());
		this->_aRBFKernel->set_points(*this->p1(),this->_aLPB->unisolvent_subset_points[j]);
		mpf_class b2(this->_aRBFKernel->basis());
		t1z += p1(j)*b1z;  
		t2z += p2z(j)*b2;
		t3z += p1(j)*p2z(j);
		for (int k = 0; k < 4; k++ ){
			if ( k!= j)
			{
				this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],this->_aLPB->unisolvent_subset_points[k]);
				mpf_class b3(this->_aRBFKernel->basis());
				t4z += p1(j) * p2z(k) * b3;
			}
		}
	}
	
	this->_aRBFKernel->set_points(*this->_p1,*this->_p2);
	mpf_class sum = this->_aRBFKernel->dz_p2() - t1z - t2z + t3z + t4z;
	return sum.get_d();
}

double Modified_Kernel::basis_planar_z_pt()
{
	mpf_class t1z,t2z,t3z,t4z;

	Matrix <mpf_class, Dynamic, 1> p1z = this->_aLPB->poly_dz(this->p1());
	Matrix <mpf_class, Dynamic, 1> p2 = this->_aLPB->poly(this->p2());

	for (int j = 0; j < 4; j++){
		this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],*this->p2());
		mpf_class b1(this->_aRBFKernel->basis());
		this->_aRBFKernel->set_points(*this->p1(),this->_aLPB->unisolvent_subset_points[j]);
		mpf_class b2z(this->_aRBFKernel->dz_p1());
		t1z += p1z(j)*b1;  
		t2z += p2(j)*b2z;
		t3z += p1z(j)*p2(j);
		for (int k = 0; k < 4; k++ ){
			if (k != j)
			{
				this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],this->_aLPB->unisolvent_subset_points[k]);
				mpf_class b3(this->_aRBFKernel->basis());
				t4z += p1z(j) * p2(k) * b3;
			}
		}
	}
	
	this->_aRBFKernel->set_points(*this->_p1,*this->_p2);
	mpf_class sum = this->_aRBFKernel->dz_p1() - t1z - t2z + t3z + t4z;
	return sum.get_d();
}

double Modified_Kernel::basis_pt_tangent()
{
	mpf_class t1x,t2x,t3x,t4x,t1y,t2y,t3y,t4y,t1z,t2z,t3z,t4z;

	Matrix <mpf_class, Dynamic, 1> p1 = this->_aLPB->poly(this->p1());

	Matrix <mpf_class, Dynamic, 1> p2x = this->_aLPB->poly_dx(this->p2());
	Matrix <mpf_class, Dynamic, 1> p2y = this->_aLPB->poly_dy(this->p2());
	Matrix <mpf_class, Dynamic, 1> p2z = this->_aLPB->poly_dz(this->p2());

	for (int j = 0; j < 4; j++){
		this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],*this->p2());
		mpf_class b1x(this->_aRBFKernel->dx_p2());
		mpf_class b1y(this->_aRBFKernel->dy_p2());
		mpf_class b1z(this->_aRBFKernel->dz_p2());
		this->_aRBFKernel->set_points(*this->p1(),this->_aLPB->unisolvent_subset_points[j]);
		mpf_class b2(this->_aRBFKernel->basis());
		// for dx component
		t1x += p1(j) * b1x;  
		t2x += p2x(j) * b2;
		t3x += p1(j) * p2x(j);
		// for dy component
		t1y += p1(j) * b1y;  
		t2y += p2y(j) * b2;
		t3y += p1(j) * p2y(j);
		// for dz component
		t1z += p1(j) * b1z;  
		t2z += p2z(j) * b2;
		t3z += p1(j) * p2z(j);
		for (int k = 0; k < 4; k++ ){
			if (k != j)
			{
				this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],this->_aLPB->unisolvent_subset_points[k]);
				mpf_class b3(this->_aRBFKernel->basis());
				t4x += p1(j) * p2x(k) * b3;
				t4y += p1(j) * p2y(k) * b3;
				t4z += p1(j) * p2z(k) * b3;
			}
		}
	}

	this->_aRBFKernel->set_points(*this->_p1,*this->_p2);
	mpf_class dx = this->_aRBFKernel->dx_p2() - t1x - t2x + t3x + t4x;
	mpf_class dy = this->_aRBFKernel->dy_p2() - t1y - t2y + t3y + t4y;
	mpf_class dz = this->_aRBFKernel->dz_p2() - t1z - t2z + t3z + t4z;

	// It is your responsibility to supply the right type
	// It is not worth the cost of dynamic_cast + if
	Tangent *t = static_cast<Tangent*>(this->p2());
	mpf_class sum = dx * t->tx() + dy * t->ty() + dz * t->tz();

	delete t;
	return sum.get_d();
}

double Modified_Kernel::basis_tangent_pt()
{
	mpf_class t1x,t2x,t3x,t4x,t1y,t2y,t3y,t4y,t1z,t2z,t3z,t4z;

	Matrix <mpf_class, Dynamic, 1> p2 = this->_aLPB->poly(this->p2());

	Matrix <mpf_class, Dynamic, 1> p1x = this->_aLPB->poly_dx(this->p1());
	Matrix <mpf_class, Dynamic, 1> p1y = this->_aLPB->poly_dy(this->p1());
	Matrix <mpf_class, Dynamic, 1> p1z = this->_aLPB->poly_dz(this->p1());

	for (int j = 0; j < 4; j++){
		this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],*this->p2());
		mpf_class b2x(this->_aRBFKernel->dx_p1());
		mpf_class b2y(this->_aRBFKernel->dy_p1());
		mpf_class b2z(this->_aRBFKernel->dz_p1());
		this->_aRBFKernel->set_points(*this->p1(),this->_aLPB->unisolvent_subset_points[j]);
		mpf_class b1(this->_aRBFKernel->basis());
		// for dx component
		t1x += p1x(j) * b1;  
		t2x += p2(j) * b2x;
		t3x += p1x(j) * p2(j);
		// for dy component
		t1y += p1y(j) * b1;  
		t2y += p2(j) * b2y;
		t3y += p1y(j) * p2(j);
		// for dz component
		t1z += p1z(j) * b1;  
		t2z += p2(j) * b2z;
		t3z += p1z(j) * p2(j);
		for (int k = 0; k < 4; k++ ){
			if ( k != j)
			{
				this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],this->_aLPB->unisolvent_subset_points[k]);
				mpf_class b3(this->_aRBFKernel->basis());
				t4x += p1x(j) * p2(k) * b3;
				t4y += p1y(j) * p2(k) * b3;
				t4z += p1z(j) * p2(k) * b3;
			}
		}
	}
	
	this->_aRBFKernel->set_points(*this->_p1,*this->_p2);
	mpf_class dx = this->_aRBFKernel->dx_p1() - t1x - t2x + t3x + t4x;
	mpf_class dy = this->_aRBFKernel->dy_p1() - t1y - t2y + t3y + t4y;
	mpf_class dz = this->_aRBFKernel->dz_p1() - t1z - t2z + t3z + t4z;

	// It is your responsibility to supply the right type
	// It is not worth the cost of dynamic_cast + if
	Tangent *t = static_cast<Tangent*>(this->p1());
	mpf_class sum = dx * t->tx() + dy * t->ty() + dz * t->tz();

	delete t;
	return sum.get_d();
}

double Modified_Kernel::basis_planar_planar( const Parameter_Types::SecondDerivatives& sd )
{
	mpf_class t1,t2,t3,t4;

	Matrix <mpf_class, Dynamic, 1> p1x = this->_aLPB->poly_dx(this->p1());
	Matrix <mpf_class, Dynamic, 1> p1y = this->_aLPB->poly_dy(this->p1());
	Matrix <mpf_class, Dynamic, 1> p1z = this->_aLPB->poly_dz(this->p1());

	Matrix <mpf_class, Dynamic, 1> p2x = this->_aLPB->poly_dx(this->p2());
	Matrix <mpf_class, Dynamic, 1> p2y = this->_aLPB->poly_dy(this->p2());
	Matrix <mpf_class, Dynamic, 1> p2z = this->_aLPB->poly_dz(this->p2());

	for (int j = 0; j < 4; j++){
		this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],*this->p2());
		mpf_class b1x(this->_aRBFKernel->dx_p2());
		mpf_class b1y(this->_aRBFKernel->dy_p2());
		mpf_class b1z(this->_aRBFKernel->dz_p2());
		this->_aRBFKernel->set_points(*this->p1(),this->_aLPB->unisolvent_subset_points[j]);
		mpf_class b2x(this->_aRBFKernel->dx_p1());
		mpf_class b2y(this->_aRBFKernel->dy_p1());
		mpf_class b2z(this->_aRBFKernel->dz_p1());
		if ( sd == Parameter_Types::DXDX)
		{
			t1 += p1x(j) * b1x;  
			t2 += p2x(j) * b2x;
			t3 += p1x(j) * p2x(j);
		}
		if ( sd == Parameter_Types::DYDY)
		{
			t1 += p1y(j) * b1y;  
			t2 += p2y(j) * b2y;
			t3 += p1y(j) * p2y(j);
		}
		if ( sd == Parameter_Types::DZDZ)
		{
			t1 += p1z(j) * b1z;  
			t2 += p2z(j) * b2z;
			t3 += p1z(j) * p2z(j);
		}
		if ( sd == Parameter_Types::DXDY)
		{
			t1 += p1x(j) * b1y;  
			t2 += p2y(j) * b2x;
			t3 += p1x(j) * p2y(j);
		}
		if ( sd == Parameter_Types::DXDZ)
		{
			t1 += p1x(j) * b1z;  
			t2 += p2z(j) * b2x;
			t3 += p1x(j) * p2z(j);
		}
		if ( sd == Parameter_Types::DYDZ)
		{
			t1 += p1y(j) * b1z;  
			t2 += p2z(j) * b2y;
			t3 += p1y(j) * p2z(j);
		}
		if ( sd == Parameter_Types::DYDX )
		{
			t1 += p1y(j) * b1x;  
			t2 += p2x(j) * b2y;
			t3 += p1y(j) * p2x(j);
		}
		if ( sd == Parameter_Types::DZDX )
		{
			t1 += p1z(j) * b1x;  
			t2 += p2x(j) * b2z;
			t3 += p1z(j) * p2x(j);
		}
		if ( sd == Parameter_Types::DZDY )
		{
			t1 += p1z(j) * b1y;  
			t2 += p2y(j) * b2z;
			t3 += p1z(j) * p2y(j);
		}
		for (int k = 0; k < 4; k++ ){
			if ( k != j)
			{
				this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],this->_aLPB->unisolvent_subset_points[k]);
				mpf_class b3(this->_aRBFKernel->basis());
				if ( sd == Parameter_Types::DXDX ) t4 += p2x(k) * b3 * p1x(j);
				if ( sd == Parameter_Types::DYDY ) t4 += p2y(k) * b3 * p1y(j);
				if ( sd == Parameter_Types::DZDZ ) t4 += p2z(k) * b3 * p1z(j);
				if ( sd == Parameter_Types::DXDY ) t4 += p2y(k) * b3 * p1x(j);
				if ( sd == Parameter_Types::DXDZ ) t4 += p2z(k) * b3 * p1x(j);
				if ( sd == Parameter_Types::DYDZ ) t4 += p2z(k) * b3 * p1y(j);
				if ( sd == Parameter_Types::DYDX ) t4 += p2x(k) * b3 * p1y(j);
				if ( sd == Parameter_Types::DZDX ) t4 += p2x(k) * b3 * p1z(j);
				if ( sd == Parameter_Types::DZDY ) t4 += p2y(k) * b3 * p1z(j);
			}
		}
	}
	
	this->_aRBFKernel->set_points(*this->_p1,*this->_p2);
	mpf_class base;
	if ( sd == Parameter_Types::DXDX ) base = this->_aRBFKernel->dxx();
	if ( sd == Parameter_Types::DYDY ) base = this->_aRBFKernel->dyy();
	if ( sd == Parameter_Types::DZDZ ) base = this->_aRBFKernel->dzz();
	if ( sd == Parameter_Types::DXDY ) base = this->_aRBFKernel->dxy();
	if ( sd == Parameter_Types::DXDZ ) base = this->_aRBFKernel->dxz();
	if ( sd == Parameter_Types::DYDZ ) base = this->_aRBFKernel->dyz();
	if ( sd == Parameter_Types::DYDX ) base = this->_aRBFKernel->dyx();
	if ( sd == Parameter_Types::DZDX ) base = this->_aRBFKernel->dzx();
	if ( sd == Parameter_Types::DZDY ) base = this->_aRBFKernel->dzy();

	mpf_class sum = base - t1 - t2 + t3 + t4;
	return sum.get_d();

}

double Modified_Kernel::basis_tangent_tangent()
{
	Matrix <mpf_class, Dynamic, 1> p1x = this->_aLPB->poly_dx(this->p1());
	Matrix <mpf_class, Dynamic, 1> p1y = this->_aLPB->poly_dy(this->p1());
	Matrix <mpf_class, Dynamic, 1> p1z = this->_aLPB->poly_dz(this->p1());

	Matrix <mpf_class, Dynamic, 1> p2x = this->_aLPB->poly_dx(this->p2());
	Matrix <mpf_class, Dynamic, 1> p2y = this->_aLPB->poly_dy(this->p2());
	Matrix <mpf_class, Dynamic, 1> p2z = this->_aLPB->poly_dz(this->p2());

	mpf_class t1xx,t2xx,t3xx,t4xx;
	mpf_class t1yy,t2yy,t3yy,t4yy;
	mpf_class t1zz,t2zz,t3zz,t4zz;
	mpf_class t1xy,t2xy,t3xy,t4xy;
	mpf_class t1xz,t2xz,t3xz,t4xz;
	mpf_class t1yz,t2yz,t3yz,t4yz;
	mpf_class t1yx,t2yx,t3yx,t4yx;
	mpf_class t1zx,t2zx,t3zx,t4zx;
	mpf_class t1zy,t2zy,t3zy,t4zy;

	for (int j = 0; j < 4; j++){
		this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],*this->p2());
		mpf_class b1x(this->_aRBFKernel->dx_p2());
		mpf_class b1y(this->_aRBFKernel->dy_p2());
		mpf_class b1z(this->_aRBFKernel->dz_p2());
		this->_aRBFKernel->set_points(*this->p1(),this->_aLPB->unisolvent_subset_points[j]);
		mpf_class b2x(this->_aRBFKernel->dx_p1());
		mpf_class b2y(this->_aRBFKernel->dy_p1());
		mpf_class b2z(this->_aRBFKernel->dz_p1());
		// dxx
		t1xx += p1x(j) * b1x;  
		t2xx += p2x(j) * b2x;
		t3xx += p1x(j) * p2x(j);
		// dyy
		t1yy += p1y(j) * b1y;  
		t2yy += p2y(j) * b2y;
		t3yy += p1y(j) * p2y(j);
		// dzz
		t1zz += p1z(j) * b1z;  
		t2zz += p2z(j) * b2z;
		t3zz += p1z(j) * p2z(j);
		// dxy
		t1xy += p1x(j) * b1y;  
		t2xy += p2y(j) * b2x;
		t3xy += p1x(j) * p2y(j);
		// dxz
		t1xz += p1x(j) * b1z;  
		t2xz += p2z(j) * b2x;
		t3xz += p1x(j) * p2z(j);
		// dyz
		t1yz += p1y(j) * b1z;  
		t2yz += p2z(j) * b2y;
		t3yz += p1y(j) * p2z(j);
		// dyx
		t1yx += p1y(j) * b1x;  
		t2yx += p2x(j) * b2y;
		t3yx += p1y(j) * p2x(j);
		// dzx
		t1zx += p1z(j) * b1x;  
		t2zx += p2x(j) * b2z;
		t3zx += p1z(j) * p2x(j);
		// dzy
		t1zy += p1z(j) * b1y;  
		t2zy += p2y(j) * b2z;
		t3zy += p1z(j) * p2y(j);
		for (int k = 0; k < 4; k++ ){
			if ( k != j)
			{
				this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],this->_aLPB->unisolvent_subset_points[k]);
				mpf_class b3(this->_aRBFKernel->basis());
				t4xx += p2x(k) * b3 * p1x(j);
				t4yy += p2y(k) * b3 * p1y(j);
				t4zz += p2z(k) * b3 * p1z(j);
				t4xy += p2y(k) * b3 * p1x(j);
				t4xz += p2z(k) * b3 * p1x(j);
				t4yz += p2z(k) * b3 * p1y(j);
				t4yx += p2x(k) * b3 * p1y(j);
				t4zx += p2x(k) * b3 * p1z(j);
				t4zy += p2y(k) * b3 * p1z(j);
			}
		}
	}
	
	this->_aRBFKernel->set_points(*this->_p1,*this->_p2);
	mpf_class xx = this->_aRBFKernel->dxx() - t1xx - t2xx + t3xx + t4xx;
	mpf_class yy = this->_aRBFKernel->dyy() - t1yy - t2yy + t3yy + t4yy;
	mpf_class zz = this->_aRBFKernel->dzz() - t1zz - t2zz + t3zz + t4zz;
	mpf_class xy = this->_aRBFKernel->dxy() - t1xy - t2xy + t3xy + t4xy;
	mpf_class xz = this->_aRBFKernel->dxz() - t1xz - t2xz + t3xz + t4xz;
	mpf_class yx = this->_aRBFKernel->dyx() - t1yx - t2yx + t3yx + t4yx;
	mpf_class yz = this->_aRBFKernel->dyz() - t1yz - t2yz + t3yz + t4yz;
	mpf_class zx = this->_aRBFKernel->dzx() - t1zx - t2zx + t3zx + t4zx;
	mpf_class zy = this->_aRBFKernel->dzy() - t1zy - t2zy + t3zy + t4zy;

	Tangent *t1 = static_cast<Tangent*>(this->p1());
	Tangent *t2 = static_cast<Tangent*>(this->p2());

	mpf_class value = t1->tx()*t2->tx()*xx + t1->tx()*t2->ty()*xy + t1->tx()*t2->tz()*xz + 
					  t1->ty()*t2->tx()*yx + t1->ty()*t2->ty()*yy + t1->ty()*t2->tz()*yz + 
					  t1->tz()*t2->tx()*zx + t1->tz()*t2->ty()*zy + t1->tz()*t2->tz()*zz;

	delete t1,t2;
	return value.get_d();
}

double Modified_Kernel::basis_planar_tangent( const Parameter_Types::FirstDerivatives& fd )
{
	if (fd == Parameter_Types::DX)
	{
		Matrix <mpf_class, Dynamic, 1> p1x = this->_aLPB->poly_dx(this->p1());
		Matrix <mpf_class, Dynamic, 1> p2x = this->_aLPB->poly_dx(this->p2());
		Matrix <mpf_class, Dynamic, 1> p2y = this->_aLPB->poly_dy(this->p2());
		Matrix <mpf_class, Dynamic, 1> p2z = this->_aLPB->poly_dz(this->p2());

		mpf_class t1xx,t2xx,t3xx,t4xx;
		mpf_class t1xy,t2xy,t3xy,t4xy;
		mpf_class t1xz,t2xz,t3xz,t4xz;

		for (int j = 0; j < 4; j++){
			this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],*this->p2());
			mpf_class b1x(this->_aRBFKernel->dx_p2());
			mpf_class b1y(this->_aRBFKernel->dy_p2());
			mpf_class b1z(this->_aRBFKernel->dz_p2());
			this->_aRBFKernel->set_points(*this->p1(),this->_aLPB->unisolvent_subset_points[j]);
			mpf_class b2x(this->_aRBFKernel->dx_p1());
			// dxx
			t1xx += p1x(j) * b1x;  
			t2xx += p2x(j) * b2x;
			t3xx += p1x(j) * p2x(j);
			// dxy
			t1xy += p1x(j) * b1y;  
			t2xy += p2y(j) * b2x;
			t3xy += p1x(j) * p2y(j);
			// dxz
			t1xz += p1x(j) * b1z;  
			t2xz += p2z(j) * b2x;
			t3xz += p1x(j) * p2z(j);
			for (int k = 0; k < 4; k++ ){
				if ( k != j)
				{
					this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],this->_aLPB->unisolvent_subset_points[k]);
					mpf_class b3(this->_aRBFKernel->basis());
					t4xx += p2x(k) * b3 * p1x(j);
					t4xy += p2y(k) * b3 * p1x(j);
					t4xz += p2z(k) * b3 * p1x(j);
				}
			}
		}
		
		this->_aRBFKernel->set_points(*this->_p1,*this->_p2);
		mpf_class xx = this->_aRBFKernel->dxx() - t1xx - t2xx + t3xx + t4xx;
		mpf_class xy = this->_aRBFKernel->dxy() - t1xy - t2xy + t3xy + t4xy;
		mpf_class xz = this->_aRBFKernel->dxz() - t1xz - t2xz + t3xz + t4xz;

		Tangent *t = static_cast<Tangent*>(this->p2());

		mpf_class value = t->tx()*xx + t->ty()*xy + t->tz()*xz;

		delete t;
		return value.get_d();
	}
	else if (fd == Parameter_Types::DY)
	{
		Matrix <mpf_class, Dynamic, 1> p1y = this->_aLPB->poly_dy(this->p1());
		Matrix <mpf_class, Dynamic, 1> p2x = this->_aLPB->poly_dx(this->p2());
		Matrix <mpf_class, Dynamic, 1> p2y = this->_aLPB->poly_dy(this->p2());
		Matrix <mpf_class, Dynamic, 1> p2z = this->_aLPB->poly_dz(this->p2());

		mpf_class t1yy,t2yy,t3yy,t4yy;
		mpf_class t1yz,t2yz,t3yz,t4yz;
		mpf_class t1yx,t2yx,t3yx,t4yx;

		for (int j = 0; j < 4; j++){
			this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],*this->p2());
			mpf_class b1x(this->_aRBFKernel->dx_p2());
			mpf_class b1y(this->_aRBFKernel->dy_p2());
			mpf_class b1z(this->_aRBFKernel->dz_p2());
			this->_aRBFKernel->set_points(*this->p1(),this->_aLPB->unisolvent_subset_points[j]);
			mpf_class b2y(this->_aRBFKernel->dy_p1());
			// dyy
			t1yy += p1y(j) * b1y;  
			t2yy += p2y(j) * b2y;
			t3yy += p1y(j) * p2y(j);
			// dyz
			t1yz += p1y(j) * b1z;  
			t2yz += p2z(j) * b2y;
			t3yz += p1y(j) * p2z(j);
			// dyx
			t1yx += p1y(j) * b1x;  
			t2yx += p2x(j) * b2y;
			t3yx += p1y(j) * p2x(j);
			for (int k = 0; k < 4; k++ ){
				if ( k != j)
				{
					this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],this->_aLPB->unisolvent_subset_points[k]);
					mpf_class b3(this->_aRBFKernel->basis());
					t4yy += p2y(k) * b3 * p1y(j);
					t4yz += p2z(k) * b3 * p1y(j);
					t4yx += p2x(k) * b3 * p1y(j);
				}
			}
		}
		
		this->_aRBFKernel->set_points(*this->_p1,*this->_p2);
		mpf_class yy = this->_aRBFKernel->dyy() - t1yy - t2yy + t3yy + t4yy;
		mpf_class yx = this->_aRBFKernel->dyx() - t1yx - t2yx + t3yx + t4yx;
		mpf_class yz = this->_aRBFKernel->dyz() - t1yz - t2yz + t3yz + t4yz;

		Tangent *t = static_cast<Tangent*>(this->p2());

		mpf_class value = t->tx()*yx + t->ty()*yy + t->tz()*yz;

		delete t;
		return value.get_d();
	}
	else // fd == DZ
	{
		Matrix <mpf_class, Dynamic, 1> p1z = this->_aLPB->poly_dz(this->p1());
		Matrix <mpf_class, Dynamic, 1> p2x = this->_aLPB->poly_dx(this->p2());
		Matrix <mpf_class, Dynamic, 1> p2y = this->_aLPB->poly_dy(this->p2());
		Matrix <mpf_class, Dynamic, 1> p2z = this->_aLPB->poly_dz(this->p2());

		mpf_class t1zz,t2zz,t3zz,t4zz;
		mpf_class t1zx,t2zx,t3zx,t4zx;
		mpf_class t1zy,t2zy,t3zy,t4zy;

		for (int j = 0; j < 4; j++){
			this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],*this->p2());
			mpf_class b1x(this->_aRBFKernel->dx_p2());
			mpf_class b1y(this->_aRBFKernel->dy_p2());
			mpf_class b1z(this->_aRBFKernel->dz_p2());
			this->_aRBFKernel->set_points(*this->p1(),this->_aLPB->unisolvent_subset_points[j]);
			mpf_class b2z(this->_aRBFKernel->dz_p1());
			// dzz
			t1zz += p1z(j) * b1z;  
			t2zz += p2z(j) * b2z;
			t3zz += p1z(j) * p2z(j);
			// dzx
			t1zx += p1z(j) * b1x;  
			t2zx += p2x(j) * b2z;
			t3zx += p1z(j) * p2x(j);
			// dzy
			t1zy += p1z(j) * b1y;  
			t2zy += p2y(j) * b2z;
			t3zy += p1z(j) * p2y(j);
			for (int k = 0; k < 4; k++ ){
				if ( k != j)
				{
					this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],this->_aLPB->unisolvent_subset_points[k]);
					mpf_class b3(this->_aRBFKernel->basis());
					t4zz += p2z(k) * b3 * p1z(j);
					t4zx += p2x(k) * b3 * p1z(j);
					t4zy += p2y(k) * b3 * p1z(j);
				}
			}
		}
		
		this->_aRBFKernel->set_points(*this->_p1,*this->_p2);
		mpf_class zz = this->_aRBFKernel->dzz() - t1zz - t2zz + t3zz + t4zz;
		mpf_class zx = this->_aRBFKernel->dzx() - t1zx - t2zx + t3zx + t4zx;
		mpf_class zy = this->_aRBFKernel->dzy() - t1zy - t2zy + t3zy + t4zy;

		Tangent *t = static_cast<Tangent*>(this->p2());

		mpf_class value = t->tx()*zx + t->ty()*zy + t->tz()*zz;

		delete t;
		return value.get_d();
	}
}

double Modified_Kernel::basis_tangent_planar( const Parameter_Types::FirstDerivatives& fd )
{
	if (fd == Parameter_Types::DX)
	{

		Matrix <mpf_class, Dynamic, 1> p1x = this->_aLPB->poly_dx(this->p1());
		Matrix <mpf_class, Dynamic, 1> p1y = this->_aLPB->poly_dy(this->p1());
		Matrix <mpf_class, Dynamic, 1> p1z = this->_aLPB->poly_dz(this->p1());
		Matrix <mpf_class, Dynamic, 1> p2x = this->_aLPB->poly_dx(this->p2());

		mpf_class t1xx,t2xx,t3xx,t4xx;
		mpf_class t1yx,t2yx,t3yx,t4yx;
		mpf_class t1zx,t2zx,t3zx,t4zx;

		for (int j = 0; j < 4; j++){
			this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],*this->p2());
			mpf_class b1x(this->_aRBFKernel->dx_p2());
			this->_aRBFKernel->set_points(*this->p1(),this->_aLPB->unisolvent_subset_points[j]);
			mpf_class b2x(this->_aRBFKernel->dx_p1());
			mpf_class b2y(this->_aRBFKernel->dy_p1());
			mpf_class b2z(this->_aRBFKernel->dz_p1());
			// dxx
			t1xx += p1x(j) * b1x;  
			t2xx += p2x(j) * b2x;
			t3xx += p1x(j) * p2x(j);
			// dyx
			t1yx += p1y(j) * b1x;  
			t2yx += p2x(j) * b2y;
			t3yx += p1y(j) * p2x(j);
			// dzx
			t1zx += p1z(j) * b1x;  
			t2zx += p2x(j) * b2z;
			t3zx += p1z(j) * p2x(j);
			for (int k = 0; k < 4; k++ ){
				if ( k != j)
				{
					this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],this->_aLPB->unisolvent_subset_points[k]);
					mpf_class b3(this->_aRBFKernel->basis());
					t4xx += p2x(k) * b3 * p1x(j);
					t4yx += p2x(k) * b3 * p1y(j);
					t4zx += p2x(k) * b3 * p1z(j);
				}
			}
		}
		
		this->_aRBFKernel->set_points(*this->_p1,*this->_p2);
		mpf_class xx = this->_aRBFKernel->dxx() - t1xx - t2xx + t3xx + t4xx;
		mpf_class yx = this->_aRBFKernel->dyx() - t1yx - t2yx + t3yx + t4yx;
		mpf_class zx = this->_aRBFKernel->dzx() - t1zx - t2zx + t3zx + t4zx;

		Tangent *t = static_cast<Tangent*>(this->p1());

		mpf_class value = t->tx()*xx + t->ty()*yx + t->tz()*zx;

		delete t;
		return value.get_d();
	}
	else if (fd == Parameter_Types::DY)
	{
		Matrix <mpf_class, Dynamic, 1> p1x = this->_aLPB->poly_dx(this->p1());
		Matrix <mpf_class, Dynamic, 1> p1y = this->_aLPB->poly_dy(this->p1());
		Matrix <mpf_class, Dynamic, 1> p1z = this->_aLPB->poly_dz(this->p1());
		Matrix <mpf_class, Dynamic, 1> p2y = this->_aLPB->poly_dy(this->p2());

		mpf_class t1yy,t2yy,t3yy,t4yy;
		mpf_class t1zy,t2zy,t3zy,t4zy;
		mpf_class t1xy,t2xy,t3xy,t4xy;

		for (int j = 0; j < 4; j++){
			this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],*this->p2());
			mpf_class b1y(this->_aRBFKernel->dy_p2());
			this->_aRBFKernel->set_points(*this->p1(),this->_aLPB->unisolvent_subset_points[j]);
			mpf_class b2x(this->_aRBFKernel->dx_p1());
			mpf_class b2y(this->_aRBFKernel->dy_p1());
			mpf_class b2z(this->_aRBFKernel->dz_p1());
			// dyy
			t1yy += p1y(j) * b1y;  
			t2yy += p2y(j) * b2y;
			t3yy += p1y(j) * p2y(j);
			// dzy
			t1zy += p1z(j) * b1y;  
			t2zy += p2y(j) * b2z;
			t3zy += p1z(j) * p2y(j);
			// dxy
			t1xy += p1x(j) * b1y;  
			t2xy += p2y(j) * b2x;
			t3xy += p1x(j) * p2y(j);
			for (int k = 0; k < 4; k++ ){
				if ( k != j)
				{
					this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],this->_aLPB->unisolvent_subset_points[k]);
					mpf_class b3(this->_aRBFKernel->basis());
					t4yy += p2y(k) * b3 * p1y(j);
					t4xy += p2y(k) * b3 * p1x(j);
					t4zy += p2y(k) * b3 * p1z(j);
				}
			}
		}
		
		this->_aRBFKernel->set_points(*this->_p1,*this->_p2);
		mpf_class yy = this->_aRBFKernel->dyy() - t1yy - t2yy + t3yy + t4yy;
		mpf_class zy = this->_aRBFKernel->dzy() - t1zy - t2zy + t3zy + t4zy;
		mpf_class xy = this->_aRBFKernel->dxy() - t1xy - t2xy + t3xy + t4xy;

		Tangent *t = static_cast<Tangent*>(this->p1());

		mpf_class value = t->tx()*xy + t->ty()*yy + t->tz()*zy;

		delete t;
		return value.get_d();
	}
	else // fd == DZ
	{
		Matrix <mpf_class, Dynamic, 1> p1x = this->_aLPB->poly_dx(this->p1());
		Matrix <mpf_class, Dynamic, 1> p1y = this->_aLPB->poly_dy(this->p1());
		Matrix <mpf_class, Dynamic, 1> p1z = this->_aLPB->poly_dz(this->p1());
		Matrix <mpf_class, Dynamic, 1> p2z = this->_aLPB->poly_dz(this->p2());

		mpf_class t1zz,t2zz,t3zz,t4zz;
		mpf_class t1xz,t2xz,t3xz,t4xz;
		mpf_class t1yz,t2yz,t3yz,t4yz;

		for (int j = 0; j < 4; j++){
			this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],*this->p2());
			mpf_class b1z(this->_aRBFKernel->dz_p2());
			this->_aRBFKernel->set_points(*this->p1(),this->_aLPB->unisolvent_subset_points[j]);
			mpf_class b2x(this->_aRBFKernel->dx_p1());
			mpf_class b2y(this->_aRBFKernel->dy_p1());
			mpf_class b2z(this->_aRBFKernel->dz_p1());
			// dzz
			t1zz += p1z(j) * b1z;  
			t2zz += p2z(j) * b2z;
			t3zz += p1z(j) * p2z(j);
			// dxz
			t1xz += p1x(j) * b1z;  
			t2xz += p2z(j) * b2x;
			t3xz += p1x(j) * p2z(j);
			// dyz
			t1yz += p1y(j) * b1z;  
			t2yz += p2z(j) * b2y;
			t3yz += p1y(j) * p2z(j);
			for (int k = 0; k < 4; k++ ){
				if ( k != j)
				{
					this->_aRBFKernel->set_points(this->_aLPB->unisolvent_subset_points[j],this->_aLPB->unisolvent_subset_points[k]);
					mpf_class b3(this->_aRBFKernel->basis());
					t4zz += p2z(k) * b3 * p1z(j);
					t4xz += p2z(k) * b3 * p1x(j);
					t4yz += p2z(k) * b3 * p1y(j);
				}
			}
		}
		
		this->_aRBFKernel->set_points(*this->_p1,*this->_p2);
		mpf_class zz = this->_aRBFKernel->dzz() - t1zz - t2zz + t3zz + t4zz;
		mpf_class xz = this->_aRBFKernel->dxz() - t1xz - t2xz + t3xz + t4xz;
		mpf_class yz = this->_aRBFKernel->dyz() - t1yz - t2yz + t3yz + t4yz;

		Tangent *t = static_cast<Tangent*>(this->p1());

		mpf_class value = t->tx()*xz + t->ty()*yz + t->tz()*zz;

		delete t;
		return value.get_d();
	}
}

VectorXd Poly_Zero::basis()
{
	if (!_truncated)
	{
		VectorXd v(1);
		v(0) = 1.0;
		return v;
	}
	else
	{
		VectorXd v;
		return v;
	}
}

VectorXd Poly_Zero::dx()
{
	if (!_truncated)
	{
		VectorXd v(1);
		v(0) = 0.0;
		return v;
	}
	else
	{
		VectorXd v;
		return v;
	}
}

VectorXd Poly_Zero::dy()
{
	if (!_truncated)
	{
		VectorXd v(1);
		v(0) = 0.0;
		return v;
	}
	else
	{
		VectorXd v;
		return v;
	}
}

VectorXd Poly_Zero::dz()
{
	if (!_truncated)
	{
		VectorXd v(1);
		v(0) = 0.0;
		return v;
	}
	else
	{
		VectorXd v;
		return v;
	}
}

VectorXd Poly_First::basis()
{
	if (!_truncated)
	{
		VectorXd v(4);
		v(0) = _p->x();
		v(1) = _p->y();
		v(2) = _p->z();
		v(3) = 1.0;
		return v;
	}
	else
	{
		VectorXd v(3);
		v(0) = _p->x();
		v(1) = _p->y();
		v(2) = _p->z();
		return v;
	}
}

VectorXd  Poly_First::dx()
{
	if (!_truncated)
	{
		VectorXd v(4);
		v(0) = 1.0;
		v(1) = 0.0;
		v(2) = 0.0;
		v(3) = 0.0;
		return v;
	}
	else
	{
		VectorXd v(3);
		v(0) = 1.0;
		v(1) = 0.0;
		v(2) = 0.0;
		return v;
	}
}

VectorXd  Poly_First::dy()
{
	if (!_truncated)
	{
		VectorXd v(4);
		v(0) = 0.0;
		v(1) = 1.0;
		v(2) = 0.0;
		v(3) = 0.0;
		return v;
	}
	else
	{
		VectorXd v(3);
		v(0) = 0.0;
		v(1) = 1.0;
		v(2) = 0.0;
		return v;
	}
}

VectorXd  Poly_First::dz()
{
	if (!_truncated)
	{
		VectorXd v(4);
		v(0) = 0.0;
		v(1) = 0.0;
		v(2) = 1.0;
		v(3) = 0.0;
		return v;
	}
	else
	{
		VectorXd v(3);
		v(0) = 0.0;
		v(1) = 0.0;
		v(2) = 1.0;
		return v;
	}
}

VectorXd  Poly_Second::basis()
{
	if (!_truncated)
	{
		VectorXd v(10);
		v(0) = _p->x() * _p->x();
		v(1) = _p->y() * _p->y();
		v(2) = _p->z() * _p->z();
		v(3) = _p->x() * _p->y();
		v(4) = _p->x() * _p->z();
		v(5) = _p->y() * _p->z();
		v(6) = _p->x();
		v(7) = _p->y();
		v(8) = _p->z();
		v(9) = 1.0;
		return v;
	}
	else
	{
		VectorXd v(9);
		v(0) = _p->x() * _p->x();
		v(1) = _p->y() * _p->y();
		v(2) = _p->z() * _p->z();
		v(3) = _p->x() * _p->y();
		v(4) = _p->x() * _p->z();
		v(5) = _p->y() * _p->z();
		v(6) = _p->x();
		v(7) = _p->y();
		v(8) = _p->z();
		return v;
	}
}

VectorXd Poly_Second::dx()
{
	if (!_truncated)
	{
		VectorXd v(10);
		v(0) = 2.0 * _p->x();
		v(1) = 0.0;
		v(2) = 0.0;
		v(3) = _p->y();
		v(4) = _p->z();
		v(5) = 0.0;
		v(6) = 1.0;
		v(7) = 0.0;
		v(8) = 0.0;
		v(9) = 0.0;
		return v;
	}
	else
	{
		VectorXd v(9);
		v(0) = 2.0 * _p->x();
		v(1) = 0.0;
		v(2) = 0.0;
		v(3) = _p->y();
		v(4) = _p->z();
		v(5) = 0.0;
		v(6) = 1.0;
		v(7) = 0.0;
		v(8) = 0.0;
		return v;
	}
}

VectorXd Poly_Second::dy()
{
	if (!_truncated)
	{
		VectorXd v(10);
		v(0) = 0.0;
		v(1) = 2.0 * _p->y();
		v(2) = 0.0;
		v(3) = _p->x();
		v(4) = 0.0;
		v(5) =_p->z();
		v(6) = 0.0;
		v(7) = 1.0;
		v(8) = 0.0;
		v(9) = 0.0;
		return v;
	}
	else
	{
		VectorXd v(9);
		v(0) = 0.0;
		v(1) = 2.0 * _p->y();
		v(2) = 0.0;
		v(3) = _p->x();
		v(4) = 0.0;
		v(5) =_p->z();
		v(6) = 0.0;
		v(7) = 1.0;
		v(8) = 0.0;
		return v;
	}
}

VectorXd Poly_Second::dz()
{
	if (!_truncated)
	{
		VectorXd v(10);
		v(0) = 0.0;
		v(1) = 0.0;
		v(2) = 2.0 * _p->z();
		v(3) = 0.0;
		v(4) = _p->x();
		v(5) = _p->y();
		v(6) = 0.0;
		v(7) = 0.0;
		v(8) = 1.0;
		v(9) = 0.0;
		return v;
	}
	else
	{
		VectorXd v(9);
		v(0) = 0.0;
		v(1) = 0.0;
		v(2) = 2.0 * _p->z();
		v(3) = 0.0;
		v(4) = _p->x();
		v(5) = _p->y();
		v(6) = 0.0;
		v(7) = 0.0;
		v(8) = 1.0;
		return v;
	}
}
