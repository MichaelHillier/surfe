#ifndef READ_INPUT_FILES_H
#define READ_INPUT_FILES_H

#include <csv.h>
#include <vector>
#include <modelling_input.h>
#include <grbf_exceptions.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataReader.h>
#include <vtkNew.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

std::string get_file_extension(const char *filename);

std::vector<Interface> build_interface_constraints(const char *interface_file);
std::vector<Planar> build_planar_constraints(const char *planar_file);
std::vector<Tangent> build_tangent_constraints(const char *tangent_file);
std::vector<Inequality> build_inequality_constraints(const char *inequality_file);


class ConstraintFileReader {
protected:
	std::string filename_;
	std::string extension_;
	std::vector<std::string> property_names_;
	void SetFilenameAndExtension(const char *filename);
	virtual void ReadPropertyNames() = 0;
public:
	std::vector<std::string> GetPropertyNames() const { return property_names_; }
};

class CSVConstraintFileReader: public ConstraintFileReader {
private:
	std::string x_name_;
	std::string y_name_;
	std::string z_name_;
protected:
	void SetXName(const char *x) { x_name_ = x; }
	void SetYName(const char *y) { y_name_ = y; }
	void SetZName(const char *z) { z_name_ = z; }
	std::string GetXName() const { return x_name_; }
	std::string GetYName() const { return y_name_; }
	std::string GetZName() const { return z_name_; }
	void ReadPropertyNames() override;
	void SearchForDefaultCoordinatePropertyNames();

};

class VTKConstraintFileReader : public ConstraintFileReader {
protected:
	void ReadPropertyNames() override;
	vtkSmartPointer<vtkPolyData> constraints_vtk;
};

//////////// INTERFACE READERS ////////////
class InterfaceConstraintFileReader : public ConstraintFileReader{
private:
	std::string level_prop_name_;
public:
	void SetLevelPropertyName(const std::string &level_property_name) { level_prop_name_ = level_property_name; }
	std::string GetLevelPropertyName() const { return level_prop_name_; }
	virtual std::vector<Interface> GetConstraints() = 0;
protected:
	void SearchForDefaultPropertyNames();
};
class CSVInterfaceConstraintFileReader : public InterfaceConstraintFileReader, public CSVConstraintFileReader{
public:
	CSVInterfaceConstraintFileReader(const char *interface_filename)
	{
		SetFilenameAndExtension(interface_filename);
		ReadPropertyNames();
		SearchForDefaultCoordinatePropertyNames();
		SearchForDefaultPropertyNames();
	}
	std::vector<Interface> GetConstraints() override;
	void SearchForDefaultPropertyNames() override;
};

class VTKInterfaceConstraintFileReader : public InterfaceConstraintFileReader, public VTKConstraintFileReader {
public:
	VTKInterfaceConstraintFileReader(const char *interface_filename)
	{
		SetFilenameAndExtension(interface_filename);
		ReadPropertyNames();
		SearchForDefaultPropertyNames();
	}
	std::vector<Interface> GetConstraints() override;
	void SearchForDefaultPropertyNames() override;
};
// END OF INTERFACE READERS

//////////// PLANAR READERS ////////////
class PlanarConstraintFileReader{
private:
	std::string normal_prop_name_;
	std::string dip_prop_name_;
	std::string strike_prop_name_;
	std::string azimuth_prop_name_;
	std::string polarity_prop_name_;
	std::string nx_prop_name_;
	std::string ny_prop_name_;
	std::string nz_prop_name_;
protected:
	virtual void SearchForDefaultPropertyNames() = 0;
public:
	void SetNormalVecPropertyName(const std::string &normal_vec_property_name) { normal_prop_name_ = normal_vec_property_name; }
	void SetDipPropertyName(const std::string &dip_property_name) { dip_prop_name_ = dip_property_name; }
	void SetStrikePropertyName(const std::string &strike_property_name) { strike_prop_name_ = strike_property_name; }
	void SetAzimuthPropertyName(const std::string &azimuth_property_name) { azimuth_prop_name_ = azimuth_property_name; }
	void SetPolarityPropertyNam(const std::string &polarity_property_name) { polarity_prop_name_ = polarity_property_name; }
	void SetNxPropertyName(const std::string &nx_property_name) { nx_prop_name_ = nx_property_name; }
	void SetNyPropertyName(const std::string &ny_property_name) { ny_prop_name_ = ny_property_name; }
	void SetNzPropertyName(const std::string &nz_property_name) { nz_prop_name_ = nz_property_name; }
	std::string GetNormalVecPropertyName() const { return normal_prop_name_; }
	std::string GetDipPropertyName() const { return dip_prop_name_; }
	std::string GetStrikePropertyName() const { return strike_prop_name_; }
	std::string GetAzimuthPropertyName() const { return azimuth_prop_name_; }
	std::string GetPolarityPropertyName() const { return polarity_prop_name_; }
	std::string GetNxPropertyName() const { return nx_prop_name_; }
	std::string GetNyPropertyName() const { return ny_prop_name_; }
	std::string GetNzPropertyName() const { return nz_prop_name_; }
	virtual std::vector<Planar> GetConstraints() = 0;
};

class CSVPlanarConstraintFileReader : public PlanarConstraintFileReader, public CSVConstraintFileReader {
public:
	CSVPlanarConstraintFileReader(const char *interface_filename)
	{
		SetFilenameAndExtension(interface_filename);
		ReadPropertyNames();
		SearchForDefaultPropertyNames();
	}
	std::vector<Planar> GetConstraints() override;
	void SearchForDefaultPropertyNames() override;
};

class VTKPlanarConstraintFileReader : public PlanarConstraintFileReader, public VTKConstraintFileReader {
public:
	VTKPlanarConstraintFileReader(const char *interface_filename)
	{
		SetFilenameAndExtension(interface_filename);
		ReadPropertyNames();
		SearchForDefaultPropertyNames();
	}
	std::vector<Planar> GetConstraints() override;
	void SearchForDefaultPropertyNames() override;
};

#endif
