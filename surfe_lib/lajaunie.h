// SURFace Estimator(SURFE) - Terms and Conditions of Use

// Unless otherwise noted, computer program source code of the SURFace Estimator(SURFE)
// is covered under Crown Copyright, Government of Canada, and is distributed under the MIT License.

// The Canada wordmark and related graphics associated with this distribution are protected under 
// trademark law and copyright law.No permission is granted to use them outside the parameters of 
// the Government of Canada's corporate identity program. For more information, see
// http://www.tbs-sct.gc.ca/fip-pcim/index-eng.asp

// Copyright title to all 3rd party software distributed with the SURFace Estimator(SURFE) is held 
// by the respective copyright holders as noted in those files.Users are asked to read the 3rd Party
// Licenses referenced with those assets.

// MIT License

// Copyright(c) 2017 Government of Canada

// Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
// and associated documentation files(the "Software"), to deal in the Software without restriction,
// including without limitation the rights to use, copy, modify, merge, publish, distribute, 
// sublicense, and / or sell copies of the Software, and to permit persons to whom the Software is 
// furnished to do so, subject to the following conditions :

// The above copyright notice and this permission notice shall be included in all copies or
// substantial portions of the Software.

//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
// NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef lajaunie_h
#define lajaunie_h

#include <surfe_lib_module.h> // macro for importing / exporting dll

#include "modeling_methods.h"

class SURFE_LIB_EXPORT Lajaunie_Approach : public GRBF_Modelling_Methods {
private:
	int _n_increment_pair;
	bool _get_polynomial_matrix_block(MatrixXd &poly_matrix);
	bool _insert_polynomial_matrix_blocks_in_interpolation_matrix(const MatrixXd &poly_matrix, MatrixXd &interpolation_matrix);
	std::vector < std::vector < Interface > > *_increment_pairs;
	bool _get_increment_pairs();
public:
	// Constructor/Destructor
	Lajaunie_Approach(const model_parameters& m_p, const Basic_input& basic_i);
	~Lajaunie_Approach(){};
	// Methods
	Polynomial_Basis *create_polynomial_basis(const int &poly_order);
	bool get_interpolation_matrix(MatrixXd &interpolation_matrix);
	bool get_equality_values(VectorXd &equality_values);
	bool get_inequality_values(VectorXd &b, VectorXd &r);
	void eval_scalar_interpolant_at_point(Point &p);
	void eval_vector_interpolant_at_point(Point &p);
	bool get_method_parameters();
	bool process_input_data();
	bool setup_system_solver();
	bool get_minimial_and_excluded_input(Basic_input &greedy_input, Basic_input &excluded_input);
	bool measure_residuals(Basic_input &input);
	bool append_greedy_input(Basic_input &input);
	bool convert_modified_kernel_to_rbf_kernel();
	GRBF_Modelling_Methods *clone() { return new Lajaunie_Approach(*this); }
	// Attributes
	Polynomial_Basis *p_basis;
};

#endif