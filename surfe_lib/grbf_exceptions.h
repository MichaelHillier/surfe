#ifndef GRBF_EXCEPTIONS_H
#define GRBF_EXCEPTIONS_H

#include <exception>

using namespace std;

class nointerfacedata : public exception {
	const char* what() const throw() override {
		return "No interface data";
	}
};

class nointerfaceincrementpairs : public exception {
	const char* what() const throw() override {
		return "There are no interface increment pairs";
	}
};

class noplanardata : public exception {
	const char* what() const throw() override {
		return "No planar data";
	}
};

class invalidinputdata : public exception {
	const char* what() const throw() override {
		return "Invalid input data as determined by check_input_data()";
	}
};

class failurecomputingglobalanisotropy : public exception {
	const char* what() const throw() override {
		return "Failure computing global anisotropy because there are less than 2 planar constraints";
	}
};


class failurecreatinganisotropickernel : public exception {
	const char* what() const throw() override {
		return "Failure creating an anisotropic kernel";
	}
};

class failuresettingupbasisfunctions : public exception {
	const char* what() const throw() override {
		return "Failure setting up basis functions";
	}
};

class failurecreatingmodifiedkernel: public exception {
	const char* what() const throw() override {
		return "Failure creating modified kernel";
	}
};

class failurecreatinglagrangianpolynomialbasis : public exception {
	const char* what() const throw() override {
		return "Failure creating Lagrangian Polynomial basis";
	}
};

class linearsolverfailure : public exception {
	const char* what() const throw() override {
		return "Eigen's linear solver failed";
	}
};

class pcquadratricsolverfailure : public exception {
	const char* what() const throw() override {
		return "Predictor-Corrector Quadratic Solver failure";
	}
};

class loqoquadratricsolverfailure : public exception {
	const char* what() const throw() override {
		return "LOQO Quadratic Solver failure";
	}
};

class errorcomputinginterpolationmatrix : public exception {
	const char* what() const throw() override {
		return "Error computing interpolation matrix";
	}
};

class errorcomputingequalityvector : public exception {
	const char* what() const throw() override {
		return "Error computing equality vector";
	}
};

class errorcomputinginequalityvector : public exception {
	const char* what() const throw() override {
		return "Error computing inequality vector";
	}
};

class errorupdatinginterfaceisovalues : public exception {
	const char* what() const throw() override {
		return "Error updating interface iso values";
	}
};

class errorcomputinginterpolant : public exception {
	const char* what() const throw() override {
		return "Error computing Interpolant";
	}
};

class SurfeExceptions : public exception {
private:
	std::string errors;
	void append_exceptions(std::string &appended_exception_string, const std::exception& e, int level = 0)
	{
		if (level == 0) {
			appended_exception_string.clear();
			appended_exception_string.append("Exceptions thrown: ");
			appended_exception_string.append(e.what());
		}
		else {
			appended_exception_string.append(", ");
			appended_exception_string.append(e.what());
		}
		try {
			std::rethrow_if_nested(e);
		}
		catch (const std::exception& e) {
			append_exceptions(appended_exception_string, e, level + 1);
		}
		catch (...) {}
	}
public:
	SurfeExceptions(const std::exception &e)
	{
		append_exceptions(errors, e);
	}
	const char* what() const throw() override {
		return errors.c_str();
	}
};


namespace GRBF_Exceptions {
	const nointerfacedata no_iterface_data;
	const nointerfaceincrementpairs no_interface_increment_pairs;
	const noplanardata no_planar_data;
	const invalidinputdata invalid_input_data;
	const failurecomputingglobalanisotropy failure_computing_global_anisotropy;
	const failurecreatinganisotropickernel failure_creating_anisotropic_kernel;
	const failuresettingupbasisfunctions failure_setting_up_basis_functions;
	const failurecreatingmodifiedkernel failure_creating_modified_kernel;
	const failurecreatinglagrangianpolynomialbasis failure_creating_lagrangian_polynomial_basis;
	const linearsolverfailure linear_solver_failure;
	const pcquadratricsolverfailure pc_quadratic_solver_failure;
	const loqoquadratricsolverfailure loqo_quadratic_solver_failure;
	const errorcomputinginterpolationmatrix error_computing_interpolation_matrix;
	const errorcomputingequalityvector error_computing_equality_vector;
	const errorcomputinginequalityvector error_computing_inequality_vector;
	const errorupdatinginterfaceisovalues error_updating_interface_iso_values;
	const errorcomputinginterpolant error_computing_interpolant;
}


#endif // 
