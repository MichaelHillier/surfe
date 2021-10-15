# README #

This code was development at Natural Resources Canada (Geological Survey of Canada) by Michael Hillier, Eric de Kemp, and Ernst Schetselaar for the purposes of 3D structural geological modelling particularly in sparse data environments.  

Our [paper](https://link.springer.com/article/10.1007/s11004-014-9540-3) associated with this work is:
Hillier MJ, Schetselaar EM, de Kemp EA, Perron G (2014) Three-dimensional modelling of geological surfaces using generalized interpolation with radial basis functions. Math Geosci 46(8):931–953

### What is this repository for? ###

* This is a library for the SURFE algorithm that implements generalized interpolation using multivariate and scattered structural geologic constraints. It accepts 4 types of constraints: inequalities, interface, planar, and tangent points. It computes an interpolant/approximate for the constraints and evaluates that function at the user supplied points.

### How do I get set up? ###

* Summary of set up
	* To setup a visual studio project to compile you must get cmake (https://cmake.org/). In cmake you will setup the compiler and dependences.

* Dependencies
	* [Eigen3](http://eigen.tuxfamily.org)
	* [pybind11](https://github.com/pybind/pybind11)
* Optional Dependencies(for visualization and simple QT interface for parameter selection and data loading)
	* To enable these features check the GEO_BUILDER variable in cmake options
	* [VTK](https://vtk.org/)
	* [Qt5](https://www.qt.io/download)
* CMake Instructions on Windows:
	* Fill in required fields for:
	* EIGEN3_INCLUDE_DIR
	* PYTHON_EXECUTABLE
	* For optional visualization and data loading features:
	* Qt5_DIR
	* VTK_DIR
 
* How can I use this library
	* Use python terminal 
	* Use native C++ code
	See code examples in the test project; main.cpp
	* Link the generated lib/so files in your own application with the corresponding dll's 

* Algorithm Pipeline
	1. Set modelling paramters
	2. Set input constraints: inequality, interface, planar, tangent (Note: Algorithm requires at least 1 planar constraint be given)
	3. Compute interpolant (To obtain a structural model)
	4. Evaluated interpolant at point or set of points.
* Optional features for surface extraction and visualization.
	1. Build a grid upon which the interpolant is evaluated   
	2. Extract Isosurfaces (To obtain surfaces representing iso-contours or structural interfaces)

* How do I specify the modelling parameters and constraints
	1. Manually specifying the individual parameters using the API and the input files containing the data
	2. Optional: Use the Qt gui - called by CreateGRBFInterpolantFromGUIParameters() in Geo_Builder class
![QT GUI](/docs/gui.JPG?raw=true)

* What file types are supported by optional QT interface?
	1. csv files
	2. vtk/vtp files

Developers can also easily expand on the support file types by deriving from the ConstraintFileReader and the corresponding subclasses. See read_input_files for reference.

* Default property names for constraint types
	* Inequality: x, y, z, level
	* Interface: x, y, z, level
	* Planar: x, y, z, dip, strike, azimuth, polarity, normal, nx, ny, nz
	Note: for orientational data not all of the properties have to be specified. 
	Multiple situations: {dip, strike, polarity}, {dip, azimuth, polarity}, {normal (a 3 component vector)}, {nx, ny, nz}
	* Tangent: x, y, z, vector, vx, vy, vz
	Note: for orientational data not all of the properties have to be specified. 
	Two situations:  {vector (a 3 component vector)}, {vx, vy, vz}
	
* Current GUI limitations:
At this time, it is assumed that default properties names are used in the input constraint files. However the ConstraintFileReader class is able to read the property names from the files (not necessarily the default property names) and individual property names can be manually set.

* Optional 3D visualization
API provides users the ability to easily visualize both their data and generated models. VTK is used for the visualization needs.
![3D Visualization](/docs/3dViz.JPG?raw=true)

* Model and Data File export
	* Models and Data can be exported to vtp (VTK File Format).

### Who do I talk to? ###

* Michael Hillier - Michael.Hillier@canada.ca
