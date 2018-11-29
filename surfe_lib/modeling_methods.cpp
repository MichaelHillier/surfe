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

using namespace Surfe;
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

    if ((int)b_input.interface_test_points->size() == 0) return false;

    // evaluate the interpolant at these interface_test_points
    if (solver != NULL)  // check if we have a valid interpolant first
    {
        for (int j = 0; j < (int)b_input.interface_test_points->size(); j++)
            eval_scalar_interpolant_at_point(
                b_input.interface_test_points->at(j));
    } else
        return false;

    if ((int)b_input.interface_iso_values->size() !=
        (int)b_input.interface_test_points->size())
        return false;
    // update interface_iso_values to computed scalar field values
    for (int j = 0; j < (int)b_input.interface_iso_values->size(); j++)
        b_input.interface_iso_values->at(j) =
            b_input.interface_test_points->at(j).scalar_field();

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
    rbf_kernel = this->create_rbf_kernel(m_parameters.basis_type,
                                         m_parameters.model_global_anisotropy);
    // check RBFKernel pointer
    if (rbf_kernel == NULL) return false;
    if (b_parameters.modified_basis) {
        if ((int)b_input.interface_point_lists->size() != 0)
            kernel =
                new Modified_Kernel(rbf_kernel, *b_input.interface_point_lists);
        else
            return false;
    } else
        kernel = rbf_kernel;

    return true;
}

bool GRBF_Modelling_Methods::check_interpolant() {
    for (int j = 0; j < b_input.itrface->size(); j++) {
        cout << " Interface[" << j << "]: " << endl;
        eval_scalar_interpolant_at_point(b_input.itrface->at(j));
        cout << "	Scalar field = " << b_input.itrface->at(j)
                                                .scalar_field() << endl;
    }

    for (int j = 0; j < b_input.planar->size(); j++) {
        eval_vector_interpolant_at_point(b_input.planar->at(j));
        double vf[3] = {b_input.planar->at(j).nx_interp(),
                        b_input.planar->at(j).ny_interp(),
                        b_input.planar->at(j).nz_interp()};
        cout << " Planar[" << j << "]: " << endl;
        cout << "	Nx = " << b_input.planar->at(j).nx()
             << " Nx interpolated = " << b_input.planar->at(j).nx_interp()
             << endl;
        cout << "	Ny = " << b_input.planar->at(j).ny()
             << " Ny interpolated = " << b_input.planar->at(j).ny_interp()
             << endl;
        cout << "	Nz = " << b_input.planar->at(j).nz()
             << " Nz interpolated = " << b_input.planar->at(j).nz_interp()
             << endl;
    }
    for (int j = 0; j < b_input.tangent->size(); j++) {
        eval_vector_interpolant_at_point(b_input.tangent->at(j));
        double vf[3] = {b_input.tangent->at(j).nx_interp(),
                        b_input.tangent->at(j).ny_interp(),
                        b_input.tangent->at(j).nz_interp()};
        cout << " Tangent[" << j << "]: " << endl;
        cout << "	Tx = " << b_input.tangent->at(j).tx()
             << " Nx interpolated = " << b_input.tangent->at(j).nx_interp()
             << endl;
        cout << "	Ty = " << b_input.tangent->at(j).ty()
             << " Ny interpolated = " << b_input.tangent->at(j).ny_interp()
             << endl;
        cout << "	Tz = " << b_input.tangent->at(j).tz()
             << " Nz interpolated = " << b_input.tangent->at(j).nz_interp()
             << endl;
        cout << " Tx*nx + Ty*ny + Tz*nz = "
             << b_input.tangent->at(j).tx() *
                        b_input.tangent->at(j).nx_interp() +
                    b_input.tangent->at(j).ty() *
                        b_input.tangent->at(j).ny_interp() +
                    b_input.tangent->at(j).tz() *
                        b_input.tangent->at(j).nz_interp() << endl;
    }

    return true;
}

bool GRBF_Modelling_Methods::evaluate_scalar_interpolant() {
    if (solver == NULL)
        return false;
    else {
        if ((int)solver->weights.size() == 0)
            return false;
        else {
            if (b_parameters.modified_basis) {
                if (!convert_modified_kernel_to_rbf_kernel()) {
                    error_msg.append(
                        " QPP solution conversion to Linear Failure.");
                    return false;
                }
            }

            int N = (int)b_input.evaluation_pts->size();
            int add = 0;
            int vv = round((double)N /
                           72.0);  // 72.0 is the width of the progress bar
            double factor = (100.0 * (double)vv) / (double)N;

#pragma omp parallel for schedule(dynamic)
            for (int j = 0; j < N; j++) {
                eval_scalar_interpolant_at_point(b_input.evaluation_pts->at(j));
                eval_vector_interpolant_at_point(b_input.evaluation_pts->at(j));
                if (vv > 0 && j % vv == 0) {
#pragma omp atomic
                    add++;
                    int step = factor * add;
#pragma omp critical
                    _Progress(" Computing Scalar field: ", step, 100);
                }
            }
            cout << endl;
        }
    }
    return true;
}

bool GRBF_Modelling_Methods::run_algorithm() {
    clock_t tstart = clock();

    // set OpenMP parameters
    const int nthreads = omp_get_max_threads();
    omp_set_dynamic(false);
    if (nthreads >= 8)
        omp_set_num_threads(nthreads - 2);
    else
        omp_set_num_threads(nthreads - 1);
    ///////////////////////////////////////

    cout << " Starting SURFE algorithm " << endl;
    cout << " Processing input data...";
    if (!process_input_data()) {
        error_msg = "Error processing input data";
        return false;
    }
    cout << "done!" << endl;
    cout << " Get method parameters...";
    if (!get_method_parameters()) {
        error_msg.append(" Error getting method parameters.");
        return false;
    }
    cout << "done!" << endl;
    cout << " Setup basis functions...";
    if (!setup_basis_functions()) {
        error_msg.append(" Error setting up basis functions.");
        return false;
    }
    cout << "done!" << endl;
    cout << " Solve mathematical problem...";
    if (!setup_system_solver()) {
        cout << "failed" << endl;
        error_msg.append(" Error solving mathematical equations.");
        return false;
    }
    cout << "done!" << endl;
    cout << " Evaluate scalar interpolant at grid nodes...";
    if (!evaluate_scalar_interpolant()) {
        error_msg.append(" Error evaluating interpolant in grid of points.");
        return false;
    }
    cout << "done!" << endl;
    cout << " Total computation time = " << ((double)clock() - tstart) /
                                                CLOCKS_PER_SEC << endl;

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

RBFKernel *GRBF_Modelling_Methods::create_rbf_kernel(
    const Parameter_Types::RBF &rbf_type, const bool &anisotropy) {
    if (anisotropy) {
        if (rbf_type == Parameter_Types::Cubic)
            return new ACubic(*b_input.planar);
        else if (rbf_type == Parameter_Types::Gaussian)
            return new AGaussian(m_parameters.shape_parameter, *b_input.planar);
        else if (rbf_type == Parameter_Types::IMQ)
            return new AIMQ(m_parameters.shape_parameter, *b_input.planar);
        else if (rbf_type == Parameter_Types::MQ)
            return new AMQ(m_parameters.shape_parameter, *b_input.planar);
        else if (rbf_type == Parameter_Types::R)
            return new AR(*b_input.planar);
        else
            return new ATPS(*b_input.planar);
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
    if (solver == NULL)
        return false;
    else {
        if ((int)solver->weights.size() == 0)
            return false;
        else {
            if ((int)b_input.evaluation_pts->size() != 0)
                evaluate_scalar_interpolant();
            else
                return false;
        }
    }

    return true;
}

GRBF_Modelling_Methods *GRBF_Modelling_Methods::get_method(
    const model_parameters &m_parameters, const Basic_input &input) {
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
    if (m_parameters.interface_uncertainty == 0 &&
        m_parameters.angular_uncertainty == 0)
        return false;

    GRBF_Modelling_Methods *greedy_method = get_method(m_parameters, b_input);

    greedy_method->b_input.compute_avg_nn_distances();

    Basic_input greedy_input, excluded_input;
    // initialize starting data
    if (!get_minimial_and_excluded_input(greedy_input, excluded_input))
        return false;
    greedy_input.SetInequalityAvgNNDist(
        greedy_method->b_input.GetInequalityAvgNNDist());
    greedy_input.SetInterfaceAvgNNDist(
        greedy_method->b_input.GetInterfaceAvgNNDist());
    greedy_input.SetPlanarAvgNNDist(
        greedy_method->b_input.GetPlanarAvgNNDist());
    greedy_input.SetTangentAvgNNDist(
        greedy_method->b_input.GetTangentAvgNNDist());
    greedy_method->b_input = greedy_input;

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

    greedy_method->b_input.evaluation_pts = b_input.evaluation_pts;
    if (!greedy_method->evaluate_scalar_interpolant()) return false;

    b_input.evaluation_pts = greedy_method->b_input.evaluation_pts;
    b_input.interface_iso_values->clear();
    b_input.interface_iso_values =
        new std::vector<double>(*greedy_method->b_input.interface_iso_values);

    return true;
}
