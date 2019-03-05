# README #

This code was development at Natural Resources Canada (Geological Survey of Canada) by Michael Hillier, Eric de Kemp, and Ernst Schetselaar for the purposes of 3D structural geological modelling particularly in sparse data environments.  

### What is this repository for? ###

* This is a library for the SURFE algorithm that implements generalized interpolation using multivariate and scattered structural geologic constraints. It accepts 4 types of constraints: inequalities, interface, planar, and tangent points. It computes an interpolant/approximate for the constraints and evaluates that function at the user supplied points. IF you want to get a surface you will have to generate a list of points that are sampled from a grid/tetrahedral structure and once you have the results of the scalar field at the list of points you will have to put those values back into your grid/tetrahedral structure then perform a marching cube/tetrahedral algorithm.  

### How do I get set up? ###

* Summary of set up
To setup a visual studio project to compile you must get cmake (https://cmake.org/). In cmake you will setup the compiler and dependences.
* Dependencies
GMP (https://gmplib.org/) and Eigen library (http://eigen.tuxfamily.org)
CMake Instructions:
	1. Point to a folder that contains the downloaded source code (e.g. ../surfe_lib)
	2. Point to an output folder for the build results (e.g. ../surfe_lib/_build)
	3. Select configure at bottom left.
	4. Select the version of you want for the vs project (make sure it's 64 bit version) Keep "Use default native compiler" selected.
	5. Then click on finish, some errors will pop, CMake needs to know paths to some libraries. When paths are entered reselect configure option.
	6. Error on EIGEN3_INCLUDE_DIR --> unzip downloaded eigen tar and point to folder that contains eigen cmake file (e.g. …/eigen)
	7. Error on GMP_INCLUDE_DIR --> Point to SURFE downloaded code folder "gmp_lib" (e.g. ../surfe_lib/gmp_lib)
	8. Error on GMP_LIBRARY_C --> Point to SURFE downloaded code library (e.g. ../surfe_lib/gmp_lib/Release/lib/libgmp-3.lib)
	9. Error on GMP_LIBRARY_CPP -> Point to SURFE downloaded code library (e.g. ../surfe_lib/gmp_lib/Release/lib/libgmpxx-3.lib).
 10. When configure is done, click on generate.
 
* What does the algorithm expect as input?
1) Fill the Basic_input (modelling_input.h) data structure with the constraints that you have and their corresponding properties/attributes : 
// input data 
std::vector< Inequality > *inequality;
std::vector< Interface > *itrface;
std::vector< Planar > *planar;
std::vector< Tangent > *tangent;

// evaluation sites in grid
std::vector< Evaluation_Point > *evaluation_pts; // these are the x/y/z locations where the interpolant is going to be evaluated.

2) Fill the model_parameters (modelling_parameters.h) data structure with the appropriate values for your model/data

3) Get the appropriate method
GRBF_Modelling_Methods* method = ... that is appropriate - note I will make a more appropriate method in the base class that will take care of this. 

4) Run the algorithm:  method->run_algorithm() or method->run_greedy_algorithm() (for greedy)

5) Get the results: method->get_evaluation_points_output()

### Who do I talk to? ###

* Michael Hillier - Michael.Hillier@canada.ca
* Other community or team contact
