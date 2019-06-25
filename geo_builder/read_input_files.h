#ifndef READ_INPUT_FILES_H
#define READ_INPUT_FILES_H

#include <surfe_api.h>
#include <csv.h>
#include <vector>
#include <string>
#include <algorithm>
#include <modelling_input.h>
#include <grbf_exceptions.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataReader.h>
#include <vtkNew.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

std::string get_file_extension(const char *filename);

std::vector<std::string> ReadPropertyNamesFromCSV(const char *filename);
std::vector<std::string> ReadPropertyNamesFromVTK(const char *filename);
vtkSmartPointer<vtkPolyData> GetvtkPolyDataGeometryFromFile(const char *filename);

class ConstraintFileReader {
private:
 	std::string x_name_;
 	std::string y_name_;
 	std::string z_name_;
protected:
	// members
	Surfe_API *surfe;
	std::string filename_;
	std::string extension_;
	std::vector<std::string> property_names_;
	// Setters/Getters
	void SetFilenameAndExtension(const char *filename);
 	void SetXName(const char *x) { x_name_ = x; }
 	void SetYName(const char *y) { y_name_ = y; }
 	void SetZName(const char *z) { z_name_ = z; }
 	std::string GetXName() const { return x_name_; }
 	std::string GetYName() const { return y_name_; }
 	std::string GetZName() const { return z_name_; }
	std::vector<std::string> GetPropertyNames() const { return property_names_; }
	// Methods
	virtual void SearchForDefaultPropertyNames() = 0;
};

//////////// INTERFACE READERS ////////////
class InterfaceConstraintFileReader : public ConstraintFileReader{
private:
	std::string level_prop_name_;
protected:
	void SearchForDefaultPropertyNames();
	void SetLevelPropertyName(const std::string &level_property_name) { level_prop_name_ = level_property_name; }
	std::string GetLevelPropertyName() const { return level_prop_name_; }
	virtual void GetConstraints() = 0;
};
//////////// PLANAR READERS ////////////
class PlanarConstraintFileReader : public ConstraintFileReader {
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
	bool have_normal;
	bool have_dip_strike_polarity;
	bool have_dip_azimuth_polarity;
	bool have_normal_components;
	void setCaseFromDefaultSearch();
	void SearchForDefaultPropertyNames();
	void SetNormalVecPropertyName(const std::string &normal_vec_property_name) { normal_prop_name_ = normal_vec_property_name; }
	void SetDipPropertyName(const std::string &dip_property_name) { dip_prop_name_ = dip_property_name; }
	void SetStrikePropertyName(const std::string &strike_property_name) { strike_prop_name_ = strike_property_name; }
	void SetAzimuthPropertyName(const std::string &azimuth_property_name) { azimuth_prop_name_ = azimuth_property_name; }
	void SetPolarityPropertyName(const std::string &polarity_property_name) { polarity_prop_name_ = polarity_property_name; }
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
	virtual void GetConstraints() = 0;
};
//////////// TANGENT READERS ////////////
class TangentConstraintFileReader : public ConstraintFileReader {
private:
	std::string vector_prop_name_;
	std::string vx_prop_name_;
	std::string vy_prop_name_;
	std::string vz_prop_name_;
protected:
	bool have_vector;
	bool have_vector_components;
	void setCaseFromDefaultSearch();
	void SearchForDefaultPropertyNames();
	void SetVectorPropertyName(const std::string &vector_property_name) { vector_prop_name_ = vector_property_name; }
	void SetVxPropertyName(const std::string &vx_property_name) { vx_prop_name_ = vx_property_name; }
	void SetVyPropertyName(const std::string &vy_property_name) { vy_prop_name_ = vy_property_name; }
	void SetVzPropertyName(const std::string &vz_property_name) { vz_prop_name_ = vz_property_name; }
	std::string GetVectorPropertyName() const { return vector_prop_name_; }
	std::string GetVxPropertyName() const { return vx_prop_name_; }
	std::string GetVyPropertyName() const { return vy_prop_name_; }
	std::string GetVzPropertyName() const { return vz_prop_name_; }
	virtual void GetConstraints() = 0;
};
//////////// INEQUALITY READERS ////////////
class InequalityConstraintFileReader : public ConstraintFileReader {
private:
	std::string level_prop_name_;
protected:
	void SearchForDefaultPropertyNames();
	void SetLevelPropertyName(const std::string &level_property_name) { level_prop_name_ = level_property_name; }
	std::string GetLevelPropertyName() const { return level_prop_name_; }
	virtual void GetConstraints() = 0;
};

class CSVInterfaceConstraintFileReader : public InterfaceConstraintFileReader{
public:
	CSVInterfaceConstraintFileReader() {}
	static CSVInterfaceConstraintFileReader CreateUsingDefaultPropertyNames(Surfe_API *surfe, const char *interface_filename)
	{
		CSVInterfaceConstraintFileReader reader;
		reader.surfe = surfe;
		reader.SetFilenameAndExtension(interface_filename);
		reader.property_names_ = ReadPropertyNamesFromCSV(interface_filename);
		reader.SearchForDefaultPropertyNames();
		return reader;
	}
	virtual void GetConstraints() override;
};

class VTKInterfaceConstraintFileReader : public InterfaceConstraintFileReader {
public:
	VTKInterfaceConstraintFileReader() {}
	static VTKInterfaceConstraintFileReader CreateUsingDefaultPropertyNames(Surfe_API *surfe, const char *interface_filename)
	{
		VTKInterfaceConstraintFileReader reader;
		reader.surfe = surfe;
		reader.SetFilenameAndExtension(interface_filename);
		reader.property_names_ = ReadPropertyNamesFromVTK(interface_filename);
		reader.SearchForDefaultPropertyNames();
		return reader;
	}
	void GetConstraints() override;
};

class CSVPlanarConstraintFileReader : public PlanarConstraintFileReader{
public:
	CSVPlanarConstraintFileReader(){}
	static CSVPlanarConstraintFileReader CreateUsingDefaultPropertyNames(Surfe_API *surfe, const char *interface_filename)
	{
		CSVPlanarConstraintFileReader reader;
		reader.surfe = surfe;
		reader.SetFilenameAndExtension(interface_filename);
		reader.property_names_ = ReadPropertyNamesFromCSV(interface_filename);
		reader.SearchForDefaultPropertyNames();
		reader.setCaseFromDefaultSearch();
		return reader;
	}
	void GetConstraints() override;
};

class VTKPlanarConstraintFileReader : public PlanarConstraintFileReader {
public:
	VTKPlanarConstraintFileReader() {}
	static VTKPlanarConstraintFileReader CreateUsingDefaultPropertyNames(Surfe_API *surfe, const char *planar_filename)
	{
		VTKPlanarConstraintFileReader reader;
		reader.surfe = surfe;
		reader.SetFilenameAndExtension(planar_filename);
		reader.property_names_ = ReadPropertyNamesFromVTK(planar_filename);
		reader.SearchForDefaultPropertyNames();
		reader.setCaseFromDefaultSearch();
		return reader;
	}
	
	void GetConstraints() override;
};

class CSVTangentConstraintFileReader : public TangentConstraintFileReader {
public:
	CSVTangentConstraintFileReader() {}
	static CSVTangentConstraintFileReader CreateUsingDefaultPropertyNames(Surfe_API *surfe, const char *tangent_filename)
	{
		CSVTangentConstraintFileReader reader;
		reader.surfe = surfe;
		reader.SetFilenameAndExtension(tangent_filename);
		reader.property_names_ = ReadPropertyNamesFromCSV(tangent_filename);
		reader.SearchForDefaultPropertyNames();
		reader.setCaseFromDefaultSearch();
		return reader;
	}
	void GetConstraints() override;
};

class VTKTangentConstraintFileReader : public TangentConstraintFileReader {
public:
	VTKTangentConstraintFileReader() {}
	static VTKTangentConstraintFileReader CreateUsingDefaultPropertyNames(Surfe_API *surfe, const char *tangent_filename)
	{
		VTKTangentConstraintFileReader reader;
		reader.surfe = surfe;
		reader.SetFilenameAndExtension(tangent_filename);
		reader.property_names_ = ReadPropertyNamesFromVTK(tangent_filename);
		reader.SearchForDefaultPropertyNames();
		reader.setCaseFromDefaultSearch();
		return reader;
	}
	void GetConstraints() override;
};

class CSVInequalityConstraintFileReader : public InequalityConstraintFileReader {
public:
	CSVInequalityConstraintFileReader() {}
	static CSVInequalityConstraintFileReader CreateUsingDefaultPropertyNames(Surfe_API *surfe, const char *inequality_filename)
	{
		CSVInequalityConstraintFileReader reader;
		reader.surfe = surfe;
		reader.SetFilenameAndExtension(inequality_filename);
		reader.property_names_ = ReadPropertyNamesFromCSV(inequality_filename);
		reader.SearchForDefaultPropertyNames();
		return reader;
	}
	void GetConstraints() override;
};

class VTKInequalityConstraintFileReader : public InequalityConstraintFileReader {
public:
	VTKInequalityConstraintFileReader() {}
	static VTKInequalityConstraintFileReader CreateUsingDefaultPropertyNames(Surfe_API *surfe, const char *inequality_filename)
	{
		VTKInequalityConstraintFileReader reader;
		reader.surfe = surfe;
		reader.SetFilenameAndExtension(inequality_filename);
		reader.property_names_ = ReadPropertyNamesFromVTK(inequality_filename);
		reader.SearchForDefaultPropertyNames();
		return reader;
	}
	void GetConstraints() override;
};


#endif
