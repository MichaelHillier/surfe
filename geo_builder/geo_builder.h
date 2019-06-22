#ifndef GEO_BUILDER_H
#define GEO_BUILDER_H

#include <surfe_api.h>
#include <modelling_parameters.h>

#include <inputImpl.h>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkImageData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkVertex.h>
#include <vtkMarchingCubes.h>
#include <vtkNew.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetMapper.h>
#include <vtkPointGaussianMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkImageSliceMapper.h>
#include <vtkImageSlice.h>
#include <vtkImageMapToColors.h>
#include <vtkImageProperty.h>
#include <vtkImageActor.h>
#include <vtkLookupTable.h>
#include <vtkArrowSource.h>
#include <vtkAssignAttribute.h>
#include <vtkGlyph3D.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkCellPicker.h>
#include <vtkImagePlaneWidget.h>

class Geo_Builder {
private:
	Surfe_API surfe_;
	InputParameters input_;
	vtkSmartPointer<vtkImageData> grid_;
	InputParameters getGUIParameters();
	void progress(const float &progress_value);
	bool evaluation_completed_;
public:
	Geo_Builder() {
		evaluation_completed_ = false;
	}
	void CreateGRBFInterpolantFromGUIParameters();
	void InitializeGRBFInterpolantObject() { surfe_ = Surfe_API(); }
	void BuildRegularGrid(const double &resolution, const double &xy_percent_padding = 0);
	void BuildRegularGrid(const double &xy_percent_padding = 0);
	void BuildRegularGrid(
		const double &zmin, const double &zmax,
		const double &resolution, const double &xy_percent_padding = 0
	);
	void BuildRegularGrid(
		const double &xmin, const double &xmax,
		const double &ymin, const double &ymax,
		const double &zmin, const double &zmax,
		const double &resolution);
	vtkSmartPointer<vtkImageData> GetEvaluatedGrid();
	vtkSmartPointer<vtkPolyData> GetIsoSurfaces();
	vtkSmartPointer<vtkPolyData> GetInterfaceConstraints();
	vtkSmartPointer<vtkPolyData> GetPlanarConstraints();
	vtkSmartPointer<vtkPolyData> GetTangentConstraints();
	vtkSmartPointer<vtkPolyData> GetInequalityConstraints();
	void WriteVTKInterfaceConstraints(const char *filename);
	void WriteVTKPlanarConstraints(const char *filename);
	void WriteVTKTangentConstraints(const char *filename);
	void WriteVTKInequalityConstraints(const char *filename);
	void WriteVTKEvaluationGrid(const char *filename);
	void WriteVTKIsoSurfaces(const char *filename);
	void VisualizeVTKData();
};

#endif
