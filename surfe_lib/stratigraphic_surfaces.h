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

#ifndef stratigraphic_surfaces_h
#define stratigraphic_surfaces_h

#include <modeling_methods.h>

class Stratigraphic_Surfaces : public GRBF_Modelling_Methods {
private:
	// Methods
	bool _get_increment_pairs();
	std::vector<std::vector<Point>> _get_lithostratigraphic_increment_pairs_for_inequality_point(const Inequality &ie_pt);
	double _get_closest_horizon_level_above_given_level(const double &given_level, const std::vector<double> &horizon_levels);
	double _get_closest_horizon_level_below_given_level(const double &given_level, const std::vector<double> &horizon_levels);
	bool _get_polynomial_matrix_block(MatrixXd &poly_matrix);
	bool _insert_polynomial_matrix_blocks_in_interpolation_matrix(const MatrixXd &poly_matrix, MatrixXd &interpolation_matrix);
	// Attributes
	int _n_increment_pairs;
	int _n_sequenced_interface_pairs;
	int _n_sequenced_inequality_pairs;
	int _n_interface_pairs;
	std::vector<std::vector<Interface> > _increment_pairs;

public:
	// Constructor/Destructor
	Stratigraphic_Surfaces(const Parameters &m_params);
	~Stratigraphic_Surfaces() {};
	// Methods
	Polynomial_Basis *create_polynomial_basis(const int &poly_order);
	bool get_interpolation_matrix(MatrixXd &interpolation_matrix) override;
	bool get_equality_values(VectorXd &equality_values) override;
	bool get_inequality_matrix(const MatrixXd &interpolation_matrix, MatrixXd &inequality_matrix);
	bool get_inequality_values(VectorXd &inequality_values);
	bool get_inequality_values(VectorXd &b, VectorXd &r);
	void eval_scalar_interpolant_at_point(Point &p) override;
	void eval_vector_interpolant_at_point(Point &p) override;
	void get_method_parameters() override;
	void process_input_data() override;
	void setup_system_solver() override;
	bool get_minimial_and_excluded_input(Constraints &greedy_input, Constraints &excluded_input) override { return true; } // TO implement
	bool measure_residuals(Constraints &input) override { return true; }  // TO implement
	bool append_greedy_input(Constraints &input) override { return true; }  // TO implement
	bool convert_modified_kernel_to_rbf_kernel() override;
	GRBF_Modelling_Methods *clone() override { return new Stratigraphic_Surfaces(*this); }
	// Attributes
	Polynomial_Basis *p_basis;
};

#endif
