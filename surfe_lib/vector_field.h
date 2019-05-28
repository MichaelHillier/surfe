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

#ifndef vector_field_h
#define vector_field_h

#include <modeling_methods.h>

class Vector_Field : public GRBF_Modelling_Methods {
public:
	// Constructor/Destructor
	Vector_Field(const Parameters &m_params)
	{
		solver = nullptr;
		kernel = nullptr;
		rbf_kernel = nullptr;

		ui_parameters = m_params;

		_iteration = 0;
	}
	~Vector_Field() {};
	// Methods
	bool get_interpolation_matrix(MatrixXd &interpolation_matrix) override;
	bool get_equality_values(VectorXd &equality_values) override;
	void eval_scalar_interpolant_at_point(Point &p) override;
	void eval_vector_interpolant_at_point(Point &p) override;
	void get_method_parameters() override;
	void process_input_data() override {};
	void setup_system_solver() override;
	bool get_minimial_and_excluded_input(Constraints &greedy_input, Constraints &excluded_input) override { return true; }  // TO implement
	bool measure_residuals(Constraints &input) override { return true; }  // TO implement
	bool append_greedy_input(Constraints &input) override { return true; }  // TO implement
	bool convert_modified_kernel_to_rbf_kernel() override { return true; }  // TO IMPLEMENT
	GRBF_Modelling_Methods *clone() override { return new Vector_Field(*this); }
};

#endif
