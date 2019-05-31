#include <surfe_api.h>
#include <iostream>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

using namespace std;

int main(int argc, char* argv[]) {
	Surfe_API surfe;

// 	try
// 	{
// 		surfe.GetParametersAndConstraints();
// 	}
// 	catch (const std::exception &e)
// 	{
// 		std::cout << "Surfe Exceptions: " << e.what() << " occurred. " << std::endl;
// 		return EXIT_FAILURE;
// 	}

	//surfe.SetInterfaceDataFile("G:/Development/surfe_lib/data/contact_data.csv");
	//surfe.SetPlanarDataFile("G:/Development/surfe_lib/data/planar_data.csv");

// 	surfe.SetGlobalAnisotropy(true);
// 	surfe.SetInterfaceDataFile("G:/Development/surfe_lib/data/anisotropy_test_interface.csv");
// 	surfe.SetPlanarDataFile("G:/Development/surfe_lib/data/anisotropy_test_planar.csv");

	surfe.SetModellingMode(2);
	//surfe.SetGlobalAnisotropy(true);
	surfe.SetPolynomialOrder(2);
	surfe.SetRBFKernel("WendlandC2");
	surfe.SetRBFShapeParameter(100.0);
	surfe.SetInterfaceDataFile("G:/Development/surfe_lib/data/refolded_isoclinal_interface.vtp");
	surfe.SetPlanarDataFile("G:/Development/surfe_lib/data/refolded_isoclinal_normals.vtp");

	//surfe.SetInterfaceDataFile("D:/Development/surfe_lib/data/Overturn_FieldData.vtp");
	//surfe.SetPlanarDataFile("D:/Development/surfe_lib/data/OverturnNormals.vtp");
	try
	{
		surfe.LoadConstraintsFromFiles();
	}
	catch (const std::exception &e)
	{
		std::cout << "Surfe Exceptions: " << e.what() << " occurred. " << std::endl;
		return EXIT_FAILURE;
	}

	try
	{
		surfe.ComputeInterpolant();
	}
	catch (const std::exception&e)
	{
		cout << "Exception: " << e.what() << endl;
		return EXIT_FAILURE;
	}

	//surfe.BuildRegularGrid(-2, 18, -2, 2, -20, 20, 0.125);
	surfe.BuildRegularGrid(10,25);
	//surfe.BuildRegularGrid(25);
	//surfe.BuildRegularGrid(-35000, 35000, 1000);
	//surfe.BuildRegularGrid(-200, 200, 10, 25);

	//surfe.WriteVTKInterfaceConstraints("G:/a_test_itr_pts.vtp");
	//surfe.WriteVTKPlanarConstraints("G:/a_test_planar_pts.vtp");
	surfe.WriteVTKIsoSurfaces("G:/Refold_WC2_100_2nd_poly.vtp");
	//surfe.WriteVTKEvaluationGrid("G:/Refold_WC2_500_2ndpoly.vti");

	surfe.VisualizeVTKData();

	const char *geo_string = surfe.GetEvaluatedVTKGridString();

	vtkSmartPointer<vtkPolyData> iso_surfaces = vtkSmartPointer<vtkPolyData>::New();
	try
	{
		iso_surfaces = surfe.GetIsoSurfaces();
	}
	catch (const std::exception&e)
	{
		cout << "Exception: " << e.what() << endl;
		return EXIT_FAILURE;
	}

	surfe.WriteVTKEvaluationGrid("D:/evaluated_sgrid.vti");
	surfe.WriteVTKIsoSurfaces("D:/iso_surface.vtp");

	return 0;
}