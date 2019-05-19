#include <surfe_api.h>
#include <iostream>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

using namespace std;

int main(int argc, char* argv[]) {
	Surfe_API surfe;

	surfe.SetInterfaceDataFile("D:/Development/surfe_lib/data/contact_data.csv");
	surfe.SetPlanarDataFile("D:/Development/surfe_lib/data/planar_data.csv");
	surfe.LoadConstraintsFromFiles();
	// 	try
	// 	{
	// 		surfe.GetUIParametersAndConstraints();
	// 	}
	// 	catch (const std::exception&e)
	// 	{
	// 		cout << "Exception: " << e.what() << endl;
	// 	}

	try
	{
		surfe.ComputeInterpolant();
	}
	catch (const std::exception&e)
	{
		cout << "Exception: " << e.what() << endl;
	}

	surfe.BuildRegularGrid(-35000, 35000, 1000);

	const char *geo_string = surfe.GetEvaluatedVTKGridString();

	vtkSmartPointer<vtkPolyData> iso_surfaces = vtkSmartPointer<vtkPolyData>::New();
	try
	{
		iso_surfaces = surfe.GetIsoSurfaces();
	}
	catch (const std::exception&e)
	{
		cout << "Exception: " << e.what() << endl;
	}

	surfe.WriteVTKEvaluationGrid("D:/evaluated_sgrid.vti");
	surfe.WriteVTKIsoSurfaces("D:/iso_surface.vtp");

	return 0;
}