#include <surfe_api.h>
#include <iostream>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

using namespace std;


int main(int argc, char* argv[]) {

	Surfe_API surfe;

	try
	{
		surfe.GetUIParameters();
	}
	catch (const std::exception&e)
	{
		cout << "Exception: " << e.what() << endl;
	}

	try
	{
		surfe.ComputeInterpolant();
	}
	catch (const std::exception&e)
	{
		cout << "Exception: " << e.what() << endl;
	}

	surfe.ConstructRegularGridOutput(-10, 10, 1);

	vtkSmartPointer<vtkPolyData> iso_surfaces = vtkSmartPointer<vtkPolyData>::New();
	try
	{
		iso_surfaces = surfe.GetIsoSurfacesAsvtkPolyData();
	}
	catch (const std::exception&e)
	{
		cout << "Exception: " << e.what() << endl;
	}


	return 0;
}
