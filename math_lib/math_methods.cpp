#include <math_methods.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cstdlib>

void Math_methods::_output_matrix( const std::vector < std::vector < double > > & matrix, const std::string &filename )
{
	std::ofstream out(filename);

	out.precision(15); // sets how many digits for each variable is outputted to the output file


	for (int j=0;j<(int)matrix.size();j++){
		for (int k=0;k<(int)matrix[j].size();k++){
			out<<"  "<<matrix[j][k];
		}
		out<<std::endl;
	}
}

void Math_methods::_output_vector( const std::vector<mpf_class > &vector, const std::string &filename )
{
	std::ofstream out(filename);

	out.precision(15); // sets how many digits for each variable is outputted to the output file


	for (int j=0;j<(int)vector.size();j++){
		out<<vector[j].get_d()<<std::endl;
	}
}

double Math_methods::RandomDouble( const double &min, const double&max )
{
	double f = (double)rand()/RAND_MAX;
	return min + f*(max - min);
}

void Math_methods::_output_vector( const std::vector<double > &vector, const std::string &filename )
{
	std::ofstream out(filename);

	out.precision(15); // sets how many digits for each variable is outputted to the output file


	for (int j=0;j<(int)vector.size();j++){
		out<<vector[j]<<std::endl;
	}
}

void Math_methods::_output_matrix( const std::vector < std::vector < mpf_class > > & matrix, const std::string &filename )
{
	std::ofstream out(filename);

	out.precision(15); // sets how many digits for each variable is outputted to the output file


	for (int j=0;j<(int)matrix.size();j++){
		for (int k=0;k<(int)matrix[j].size();k++){
			out<<"  "<<matrix[j][k].get_d();
		}
		out<<std::endl;
	}
}