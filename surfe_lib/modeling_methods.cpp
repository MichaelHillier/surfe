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

#include <basis.h>
#include <continuous_property.h>
#include <lajaunie.h>
#include <math_methods.h>
#include <matrix_solver.h>
#include <modeling_methods.h>
#include <single_surface.h>
#include <stratigraphic_surfaces.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <time.h>
#include <vector>

double round(double d) { return floor(d + 0.5); }

bool GRBF_Modelling_Methods::_update_interface_iso_values() {
    // this is a messy method.
    // should really only be called if it is a Lajaunie method or Stratigraphic
    // method have put in safe guards to ensure no seg faults function purpose:
    // When using increment approaches we do not know what the scalar field
    // value
    // will be at interface points. Therefore, we wouldn't be able to do a
    // iso-surface extraction. To solve this issue, after the interpolant have
    // been determined we evaluate the interpolant at a interface point (in
    // b_input.interface_test_points) for each interace. Then we will know the
    // right scalar field value for each interface to complete an iso-surface
    // extraction.

    if (constraints.interface_test_points.empty()) return false;

    // evaluate the interpolant at these interface_test_points
    if (solver != nullptr)  // check if we have a valid interpolant first
    {
		for (auto &interface_test_point: constraints.interface_test_points)
            eval_scalar_interpolant_at_point(interface_test_point);
    } 
	else
        return false;

    if (constraints.interface_iso_values.size() != constraints.interface_test_points.size())
        return false;
    // update interface_iso_values to computed scalar field values
    for (int j = 0; j < (int)constraints.interface_iso_values.size(); j++)
        constraints.interface_iso_values[j] = constraints.interface_test_points[j].scalar_field();

    return true;
}

void GRBF_Modelling_Methods::_Progress(char message[], const int &step,
                                       const int &total) {
    // progress width
    const int pwidth = 72;

    // minus label len
    int width = pwidth - strlen(message);
    int pos = (step * width) / total;
    int percent = (step * 100) / total;

    printf("%s[", message);

    // fill progress bar with =
    for (int i = 0; i < pos; i++) printf("%c", '=');

    // fill progress bar with spaces
    printf("% *c", width - pos + 1, ']');
    printf(" %3d%%\r", percent);
}

bool GRBF_Modelling_Methods::setup_basis_functions() {
    rbf_kernel = this->create_rbf_kernel(m_parameters.basis_type, m_parameters.model_global_anisotropy);
    // check RBFKernel pointer
    if (rbf_kernel == nullptr) return false;
    if (b_parameters.modified_basis) {
        if (!constraints.interface_point_lists.empty())
            kernel = new Modified_Kernel(rbf_kernel, constraints.interface_point_lists);
        else
            return false;
    }
	else
        kernel = rbf_kernel;

    return true;
}

bool GRBF_Modelling_Methods::check_interpolant() {

	for (auto &interface_pt: constraints.itrface){
        eval_scalar_interpolant_at_point(interface_pt);
        std::cout << "	Scalar field = " << interface_pt.scalar_field() << std::endl;
    }

	for (auto &planar_pt: constraints.planar){
        eval_vector_interpolant_at_point(planar_pt);
        double vf[3] = {planar_pt.nx_interp(),
						planar_pt.ny_interp(),
						planar_pt.nz_interp()};
    }
	for (auto &tangent_pt: constraints.tangent){
        eval_vector_interpolant_at_point(tangent_pt);
        double vf[3] = {tangent_pt.nx_interp(),
						tangent_pt.ny_interp(),
						tangent_pt.nz_interp()};
    }

    return true;
}

bool GRBF_Modelling_Methods::get_equality_matrix(
    const MatrixXd &interpolation_matrix, MatrixXd &equality_matrix) {
    if (equality_matrix.rows() == 0 ||
        equality_matrix.rows() > interpolation_matrix.rows() ||
        equality_matrix.cols() != interpolation_matrix.cols())
        return false;
    int n_ie = (int)interpolation_matrix.rows() - (int)equality_matrix.rows();
    if (n_ie != b_parameters.n_inequality) return false;

    for (int j = 0; j < equality_matrix.rows(); j++) {
        for (int k = 0; k < equality_matrix.cols(); k++) {
            equality_matrix(j, k) = interpolation_matrix(j + n_ie, k);
        }
    }

    return true;
}

RBFKernel *GRBF_Modelling_Methods::create_rbf_kernel(const Parameter_Types::RBF &rbf_type, const bool &anisotropy) {
    if (anisotropy) {
        if (rbf_type == Parameter_Types::Cubic)
            return new ACubic(constraints.planar);
        else if (rbf_type == Parameter_Types::Gaussian)
            return new AGaussian(m_parameters.shape_parameter, constraints.planar);
        else if (rbf_type == Parameter_Types::IMQ)
            return new AIMQ(m_parameters.shape_parameter, constraints.planar);
        else if (rbf_type == Parameter_Types::MQ)
            return new AMQ(m_parameters.shape_parameter, constraints.planar);
        else if (rbf_type == Parameter_Types::R)
            return new AR(constraints.planar);
        else
            return new ATPS(constraints.planar);
    } else {
        // if (b_input._weights.size() != 0) return new
        // Scaled_Cubic(b_input._weights,b_input._points);
        if (rbf_type == Parameter_Types::Cubic)
            return new Cubic;
        else if (rbf_type == Parameter_Types::Gaussian)
            return new Gaussian(m_parameters.shape_parameter);
        else if (rbf_type == Parameter_Types::IMQ)
            return new IMQ(m_parameters.shape_parameter);
        else if (rbf_type == Parameter_Types::MQ)
            return new MQ(m_parameters.shape_parameter);
        else if (rbf_type == Parameter_Types::R)
            return new R;
        else
            return new TPS;
    }
}

bool GRBF_Modelling_Methods::_output_greedy_debug_objects() {
    if (solver == nullptr)
        return false;
    else {
        if ((int)solver->weights.size() == 0)
            return false;
        else {
//             if (!constraints.evaluation_pts.empty()) TODO
//                 evaluate_scalar_interpolant();
//             else
//                 return false;
        }
    }

    return true;
}

GRBF_Modelling_Methods *GRBF_Modelling_Methods::get_method(
	const model_parameters &m_parameters, const Constraints &input)
{
    if (m_parameters.model_type == Parameter_Types::Single_surface)
        return new Single_Surface(m_parameters, input);
    else if (m_parameters.model_type == Parameter_Types::Lajaunie_approach)
        return new Lajaunie_Approach(m_parameters, input);
    else if (m_parameters.model_type == Parameter_Types::Stratigraphic_horizons)
        return new Stratigraphic_Surfaces(m_parameters, input);
    else
        return new Continuous_Property(m_parameters, input);
}

bool GRBF_Modelling_Methods::run_greedy_algorithm() {
    // check if there are non-zero errors permitted on the data
    if (m_parameters.interface_uncertainty == 0 && m_parameters.angular_uncertainty == 0)
        return false;

    GRBF_Modelling_Methods *greedy_method = get_method(m_parameters, constraints);

    greedy_method->constraints.compute_avg_nn_distances();

    Constraints greedy_input, excluded_input;
    // initialize starting data
    if (!get_minimial_and_excluded_input(greedy_input, excluded_input))
        return false;

    greedy_input.SetInequalityAvgNNDist(greedy_method->constraints.GetInequalityAvgNNDist());
    greedy_input.SetInterfaceAvgNNDist(greedy_method->constraints.GetInterfaceAvgNNDist());
    greedy_input.SetPlanarAvgNNDist(greedy_method->constraints.GetPlanarAvgNNDist());
    greedy_input.SetTangentAvgNNDist(greedy_method->constraints.GetTangentAvgNNDist());
    greedy_method->constraints = greedy_input;

    bool converged = false;
    int iter = 0;
    while (!converged) {
        // run normal algorithm
        if (!greedy_method->process_input_data()) return false;
        if (!greedy_method->get_method_parameters()) return false;
        if (!greedy_method->setup_basis_functions()) return false;
        if (!greedy_method->setup_system_solver()) return false;

        // measure residuals
        if (!greedy_method->measure_residuals(excluded_input)) return false;

        // debug: should output intermediate input constraints and modelled
        // surface
        // using those constraints if (
        // !greedy_method->_output_greedy_debug_objects()) return false;

        // add appropriate data based on residuals
        if (!greedy_method->append_greedy_input(excluded_input))
            converged = true;  // if no input is added convergence is assumed
        iter++;
        greedy_method->_SetIteration(iter);
    }

//     greedy_method->constraints.evaluation_pts = constraints.evaluation_pts;
//     if (!greedy_method->evaluate_scalar_interpolant()) return false;
// 
//     constraints.evaluation_pts = greedy_method->constraints.evaluation_pts;
//     constraints.interface_iso_values.clear();
//     constraints.interface_iso_values =
// 		std::vector<double>(greedy_method->constraints.interface_iso_values);

    return true;
}
