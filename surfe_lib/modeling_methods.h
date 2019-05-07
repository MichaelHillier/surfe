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

#include <surfe_lib_module.h>  // macro for importing / exporting dll

#include <Eigen/Core>
#include <basis.h>
#include <matrix_solver.h>
#include <modelling_input.h>
#include <modelling_parameters.h>

#include <omp.h>

using namespace Eigen;

// Abstract base class
class SURFE_LIB_EXPORT GRBF_Modelling_Methods {
protected:
    // ATTRIBUTES
    model_parameters m_parameters;  // QT GUI parameters
    basic_parameters b_parameters;  // algorithm parameters
    Constraints b_input;            // algorithm input
    int _iteration;                 // for greedy progress
    // METHODS
    bool _update_interface_iso_values();  // this is to prep for output. Is the
    // computed scalar field value using the
    // interpolant @ interface_test_points
    // for iso surface extraction
    void _Progress(char message[], const int &step, const int &total);
    bool _output_greedy_debug_objects();
    void _SetIteration(const int &iter) { _iteration = iter; }

public:
    // Destructor
    virtual ~GRBF_Modelling_Methods() {}
    // Methods
    GRBF_Modelling_Methods *get_method(const model_parameters &m_parameters, const Constraints &input);  // factory method to get create the appropriate pointer for given problem
    RBFKernel *create_rbf_kernel(const Parameter_Types::RBF &rbf_type, const bool &anisotropy);
    std::vector<Evaluation_Point> get_evaluation_points_output() const { return b_input.evaluation_pts;}
    std::vector<Interface> get_interface_points_ouput() const { return b_input.itrface; }
    Constraints get_b_input() const { return b_input; }
    std::vector<double> get_interface_iso_values() const { return b_input.interface_iso_values; }
    bool setup_basis_functions();
    bool check_interpolant();
    bool evaluate_scalar_interpolant();
    bool evaluate_vector_interpolant();
    bool run_algorithm();
    bool run_greedy_algorithm();
    bool get_equality_matrix(const MatrixXd &interpolation_matrix,
                             MatrixXd &equality_matrix);
    virtual bool get_interpolation_matrix(MatrixXd &interpolation_matrix) = 0;
    virtual bool get_equality_values(VectorXd &equality_values) = 0;
    virtual void eval_scalar_interpolant_at_point(Point &p) = 0;
    virtual void eval_vector_interpolant_at_point(Point &p) = 0;
    virtual bool get_method_parameters() = 0;
    virtual bool process_input_data() = 0;
    virtual bool setup_system_solver() = 0;
    virtual bool get_minimial_and_excluded_input(
        Constraints &greedy_input, Constraints &excluded_input) = 0;
    virtual bool measure_residuals(Constraints &input) = 0;
    virtual bool append_greedy_input(Constraints &input) = 0;
    virtual bool convert_modified_kernel_to_rbf_kernel() = 0;
    virtual GRBF_Modelling_Methods *clone() = 0;
    // Attributes
    System_Solver *solver;
    Kernel *kernel;
    RBFKernel *rbf_kernel;
    std::string error_msg;
};

#endif
