#ifndef modeling_methods_h
#define modeling_methods_h

#include <surfe_lib_module.h> // macro for importing / exporting dll

#include <modelling_parameters.h>
#include <modelling_input.h>
#include <matrix_solver.h>
#include <basis.h>
#include <Eigen/Core>

using namespace std;
using namespace Eigen;

// Abstract base class
class SURFE_LIB_EXPORT GRBF_Modelling_Methods {
protected:
	// ATTRIBUTES
	model_parameters m_parameters; // QT GUI parameters
	basic_parameters b_parameters; // algorithm parameters
	Basic_input b_input; // algorithm input
	int _iteration; // for greedy progress
	// METHODS
	bool _update_interface_iso_values(); // this is to prep for output. Is the computed scalar field value using the interpolant @ interface_test_points for iso surface extraction
	void _Progress(char message[], const int &step, const int &total);
	bool _output_greedy_debug_objects();
	void _SetIteration(const int &iter) { _iteration = iter; }
public:
	// Destructor
	virtual ~GRBF_Modelling_Methods(){}
	// Methods
	GRBF_Modelling_Methods *get_method(const model_parameters& m_parameters,const Basic_input& input ); // factory method to get create the appropriate pointer for given problem
	RBFKernel *create_rbf_kernel(const Parameter_Types::RBF &rbf_type, const bool &anisotropy);
	std::vector<Evaluation_Point> *get_evaluation_points_output() const { return b_input.evaluation_pts; }
	std::vector<Interface> *get_interface_points_ouput() const { return b_input.itrface; }
	Basic_input get_b_input() const { return b_input; }
	std::vector<double> *get_interface_iso_values() const { return b_input.interface_iso_values; }
	bool setup_basis_functions();
	bool check_interpolant();
	bool evaluate_scalar_interpolant();
	bool evaluate_vector_interpolant();
	bool run_algorithm();
	bool run_greedy_algorithm();
	bool get_equality_matrix(const MatrixXd &interpolation_matrix, MatrixXd &equality_matrix);
	virtual bool get_interpolation_matrix(MatrixXd &interpolation_matrix) = 0;
	virtual bool get_equality_values(VectorXd &equality_values) = 0;
	virtual void eval_scalar_interpolant_at_point(Point &p) = 0;
	virtual void eval_vector_interpolant_at_point(Point &p) = 0;
	virtual bool get_method_parameters() = 0;
	virtual bool process_input_data() = 0;
	virtual bool setup_system_solver() = 0;
	virtual bool get_minimial_and_excluded_input(Basic_input &greedy_input, Basic_input &excluded_input) = 0;
	virtual bool measure_residuals(Basic_input &input) = 0;
	virtual bool append_greedy_input(Basic_input &input) = 0;
	virtual bool convert_modified_kernel_to_rbf_kernel() = 0;
	virtual GRBF_Modelling_Methods *clone() = 0;
	// Attributes
	System_Solver *solver;
	Kernel *kernel;
	RBFKernel *rbf_kernel;
	std::string error_msg;
};

#endif