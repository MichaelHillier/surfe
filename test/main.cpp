#include <geo_builder.h>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {

	Geo_Builder model;

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

	//surfe.SetModellingMode(2);
	//surfe.SetGlobalAnisotropy(true);
	//surfe.SetPolynomialOrder(1);
// 	surfe.SetRBFKernel("WendlandC2");
// 	surfe.SetRBFShapeParameter(200.0);
//	surfe.SetRBFKernel("Gaussian");
//	surfe.SetRBFShapeParameter(4.44e-3);
// 	surfe.SetRBFKernel("MaternC4");
// 	surfe.SetRBFShapeParameter(0.02);
	//surfe.SetInterfaceDataFile("G:/Development/surfe_lib/data/refolded_isoclinal_interface.vtp");
	//surfe.SetPlanarDataFile("G:/Development/surfe_lib/data/refolded_isoclinal_normals.vtp");

// 	surfe.SetModellingMode(4);
// 	surfe.SetInequalityDataFile("C:/Development/surfe_lib/data/strat_inequality25.vtp");
// 	surfe.SetInterfaceDataFile("C:/Development/surfe_lib/data/strat_interface25.vtp");
// 	surfe.SetPlanarDataFile("C:/Development/surfe_lib/data/strat_normal25.vtp");

	model.InitializeGRBFInterpolantObject(1);
	model.surfe->SetRestrictedRange(true, 50, 10);
//   	model.surfe->SetRBFKernel("MaternC4");
//   	model.surfe->SetRBFShapeParameter(0.005);
//  	model.surfe->SetPolynomialOrder(1);
	//model.SetInequalityDataFile("C:/Development/surfe_lib/data/conic_ie10.vtp");
	model.SetInterfaceDataFile("../data/Overturn_FieldData.vtp");
	model.SetPlanarDataFile("../data/OverturnNormals.vtp");


	//surfe.SetInterfaceDataFile("C:/Development/surfe_lib/data/disconnection_interface.csv");
	//surfe.SetPlanarDataFile("C:/Development/surfe_lib/data/disconnection_normal.csv");

	//surfe.SetInterfaceDataFile("G:/Development/surfe_lib/data/Overturn_FieldData.vtp");
	//surfe.SetPlanarDataFile("G:/Development/surfe_lib/data/OverturnNormals.vtp");
	try
	{
		model.LoadConstraintsFromFiles();
	}
	catch (const std::exception &e)
	{
		std::cout << "Surfe Exceptions: " << e.what() << " occurred. " << std::endl;
		return EXIT_FAILURE;
	}

	try
	{
		model.surfe->ComputeInterpolant();
	}
	catch (const std::exception&e)
	{
		cout << "Exception: " << e.what() << endl;
		return EXIT_FAILURE;
	}

	model.BuildRegularGrid(10,45);
	//surfe.BuildRegularGrid(-10, 200, -20, 75, -100, 20, 1);
	//surfe.BuildRegularGrid(10,25);
	//surfe.BuildRegularGrid(5,25);
	//surfe.BuildRegularGrid(-35000, 35000, 1000);
	//surfe.BuildRegularGrid(-200, 200, 10, 25);

	//surfe.WriteVTKInterfaceConstraints("C:/Research/SurfeOutput/a_test_itr_pts.vtp");
	//surfe.WriteVTKPlanarConstraints("C:/Research/SurfeOutput/a_test_planar_pts.vtp");
// 	surfe.WriteVTKIsoSurfaces("C:/Research/SurfeOutput/strat_r3.vtp");
// 	surfe.WriteVTKEvaluationGrid("C:/Research/SurfeOutput/strat_r3.vti");
	model.WriteVTKIsoSurfaces("D:/Development/SurfeOutput/overturnR3_50_10.vtp");
	//model.WriteVTKEvaluationGrid("D:/Development/SurfeOutput/overturnAngulargrid25_50.vti");

	model.VisualizeVTKData();

	return 0;
}