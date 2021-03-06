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

#ifndef modeling_methods_h
#define modeling_methods_h

#include <Eigen/Core>
#include <grbf_exceptions.h>
#include <basis.h>
#include <matrix_solver.h>

#include <omp.h>

using namespace Eigen;

// Abstract base class
class GRBF_Modelling_Methods {
private:
	void _get_distinct_interface_iso_values();
	void _get_interface_points();
	std::vector<double> _get_distinct_inequality_iso_values();
	bool _interface_points_are_coplanar() {
		return true;
	}  // Not implemented yet. should be tested when 2nd order polynomials are
	   // used. Also when unisolvent points are used this should be called.
protected:
	// ATTRIBUTES
	InternalParameters intern_params;  // algorithm parameters
	int _iteration;                 // for greedy progress
	// for interface data
	std::vector<double> interface_iso_values;
	std::vector<std::vector<Interface> > interface_point_lists;
	//std::vector<Interface> interface_test_points;

	bool get_interface_data();
	// fills interface_iso_values, interface_point_lists, interface_test_points data structures
	// validation
	bool check_input_data();

	// METHODS
	bool _update_interface_iso_values();  // this is to prep for output. Is the
	// computed scalar field value using the
	// interpolant @ interface_test_points
	// for iso surface extraction
	bool _output_greedy_debug_objects();
	void _SetIteration(const int &iter) { _iteration = iter; }

public:
	// Destructor
	virtual ~GRBF_Modelling_Methods() {}
	// Methods
	GRBF_Modelling_Methods *get_method(const Parameters &m_parameters);  // factory method to get create the appropriate pointer for given problem
	RBFKernel *create_rbf_kernel(const Parameter_Types::RBF &rbf_type, const bool &anisotropy);
	std::vector<Interface> get_interface_points_ouput() const { return constraints.itrface; }
	Constraints constraints;// algorithm input
	void remove_collocated_constraints(); // cleaning method to ensure valid interpolation matrix
	std::vector<double> get_interface_iso_values() const { return interface_iso_values; }
	void setup_basis_functions();
	bool check_interpolant();
	bool run_greedy_algorithm();
	bool get_equality_matrix(const MatrixXd &interpolation_matrix, MatrixXd &equality_matrix);
	virtual bool get_interpolation_matrix(MatrixXd &interpolation_matrix) = 0;
	virtual bool get_equality_values(VectorXd &equality_values) = 0;
	virtual void eval_scalar_interpolant_at_point(Point &p) = 0;
	virtual void eval_vector_interpolant_at_point(Point &p) = 0;
	virtual void get_method_parameters() = 0;
	virtual void process_input_data() = 0;
	virtual void setup_system_solver() = 0;
	virtual bool get_minimial_and_excluded_input(Constraints &greedy_input, Constraints &excluded_input) = 0;
	virtual bool measure_residuals(Constraints &input) = 0;
	virtual bool append_greedy_input(Constraints &input) = 0;
	virtual bool convert_modified_kernel_to_rbf_kernel() = 0;
	virtual GRBF_Modelling_Methods *clone() = 0;
	// Attributes
	Parameters parameters;  // QT GUI parameters
	System_Solver *solver;
	Kernel *kernel;
	RBFKernel *rbf_kernel;
	std::string error_msg;
	std::vector<Interface> interface_test_points;
};

#endif
