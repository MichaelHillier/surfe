#include <math_methods.h>
#include <cstdlib>

double Math_methods::RandomDouble( const double &min, const double&max )
{
	double f = (double)rand()/RAND_MAX;
	return min + f*(max - min);
}