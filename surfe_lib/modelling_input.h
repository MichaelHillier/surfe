// SURFace Estimator(SURFE) - Terms and Conditions of Use

// Unless otherwise noted, computer program source code of the SURFace
// Estimator(SURFE) is covered under Crown Copyright, Government of Canada, and
// is distributed under the MIT License.

// The Canada wordmark and related graphics associated with this distribution
// are protected under trademark law and copyright law.No permission is granted
// to use them outside the parameters of the Government of Canada's corporate
// identity program. For more information, see
// http://www.tbs-sct.gc.ca/fip-pcim/index-eng.asp

// Copyright title to all 3rd party software distributed with the SURFace
// Estimator(SURFE) is held by the respective copyright holders as noted in
// those files.Users are asked to read the 3rd Party Licenses referenced with
// those assets.

// MIT License

// Copyright(c) 2017 Government of Canada

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files(the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and / or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions :

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT
// NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
// PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef _modelling_input_h
#define _modelling_input_h

#include <modelling_parameters.h>
#include <surfe_lib_module.h>
#include <cmath>
#include <vector>
#include <string>
#include <set>

class SURFE_LIB_EXPORT Point {
private:
    double _x;
    double _y;
    double _z;
    double _c;
    double _scalar_field;
    double _field_normal[3];

public:
    Point(const double &x_coord, const double &y_coord, const double &z_coord,
          const double &c_coord = 0)
        : _x(x_coord), _y(y_coord), _z(z_coord), _c(c_coord) {
        _scalar_field = 0;
        for (int j = 0; j < 3; j++) _field_normal[j] = 0;
    }
	// Setters
    void set_x(const double &x_coord) { _x = x_coord; }
    void set_y(const double &y_coord) { _y = y_coord; }
    void set_z(const double &z_coord) { _z = z_coord; }
    void set_c(const double &c_coord) { _c = c_coord; }
    double scalar_field() const { return _scalar_field; }
    void set_scalar_field(const double &scalar_field_value) { _scalar_field = scalar_field_value; }
	void set_vector_field(const double &nx, const double &ny, const double &nz) { 
		_field_normal[0] = nx;
		_field_normal[1] = ny;
		_field_normal[2] = nz;
	}
	// Getters
	double x() const { return _x; }
	double y() const { return _y; }
	double z() const { return _z; }
	double c() const { return _c; }
    double nx_interp() const { return _field_normal[0]; }
    double ny_interp() const { return _field_normal[1]; }
    double nz_interp() const { return _field_normal[2]; }
};

class SURFE_LIB_EXPORT Evaluation_Point : public Point {
public:
    Evaluation_Point(const double &x_coord, const double &y_coord,
                     const double &z_coord, const double &c_coord = 0)
        : Point(x_coord, y_coord, z_coord, c_coord) {}
};

class SURFE_LIB_EXPORT Interface : public Point {
private:
    double _level;
    double _residual;
    double _level_bound[2];

public:
    Interface(const double &x_coord, const double &y_coord,
              const double &z_coord, const double &lvl,
              const double &c_coord = 0)
        : Point(x_coord, y_coord, z_coord, c_coord), _level(lvl) {
        _residual = 0.0;
        _level_bound[0] = 0.0;
        _level_bound[1] = 0.0;
    }
    double level() const { return _level; }
    double residual() const { return _residual; }
    double level_lower_bound() const { return _level_bound[0]; }
    double level_upper_bound() const { return _level_bound[1]; }
    void setResidual(const double &res) { _residual = res; }
    void setLevel(const double &v) { _level = v; }
    void setLevelBounds(const double &level_uncertainty) {
        _level_bound[0] = -1.0 * level_uncertainty;
        _level_bound[1] = level_uncertainty;
    }
};

class SURFE_LIB_EXPORT Inequality : public Point {
private:
    double _inequality_level;
    bool _residual;

public:
    Inequality(const double &x_coord, const double &y_coord,
               const double &z_coord, const double &lvl,
               const double &c_coord = 0)
        : Point(x_coord, y_coord, z_coord, c_coord), _inequality_level(lvl) {
        _residual = true;
    }
    double level() const { return _inequality_level; }
    bool residual() const { return _residual; }
    void setResidual(const bool &res) { _residual = res; }
};

class SURFE_LIB_EXPORT Planar : public Point {
private:
    double _dip;
    double _strike;
    int _polarity;
    double _normal[3];
    double _residual;
    double _normal_bound[3][2];
    bool _compute_strike_dip_polarity_from_normal();
    bool _compute_normal_from_strike_dip_polarity();

public:
    Planar(const double &x_coord, const double &y_coord, const double &z_coord,
           const double &nx, const double &ny, const double &nz,
           const double &c_coord = 0)
        : Point(x_coord, y_coord, z_coord, c_coord) {
        _normal[0] = nx;
        _normal[1] = ny;
        _normal[2] = nz;
        _residual = 0.0;
        // compute strike, dip, and polarity
        _compute_strike_dip_polarity_from_normal();
    }
    Planar(const double &x_coord, const double &y_coord, const double &z_coord,
           const double &dip, const double &strike, const int &polarity,
           const double &c_coord = 0)
        : Point(x_coord, y_coord, z_coord, c_coord),
          _dip(dip),
          _strike(strike),
          _polarity(polarity) {
        // compute normal
        _compute_normal_from_strike_dip_polarity();
        // residual default
        _residual = 0.0;
    }
    bool getDipVector(double (&vector)[3]);
    bool getStrikeVector(double (&vector)[3]);
    double dip() const { return _dip; }
    double strike() const { return _strike; }
    int polarity() const { return _polarity; }
    double nx() const { return _normal[0]; }
    double ny() const { return _normal[1]; }
    double nz() const { return _normal[2]; }
    double nx_lower_bound() const { return _normal_bound[0][0]; }
    double nx_upper_bound() const { return _normal_bound[0][1]; }
    double ny_lower_bound() const { return _normal_bound[1][0]; }
    double ny_upper_bound() const { return _normal_bound[1][1]; }
    double nz_lower_bound() const { return _normal_bound[2][0]; }
    double nz_upper_bound() const { return _normal_bound[2][1]; }
    void setNormalBounds(const double &delta_strike, const double &delta_dip);
    double residual() const { return _residual; }
    void setResidual(const double &res) { _residual = res; }
    void setNormal(const double &nx, const double &ny, const double &nz) {
        _normal[0] = nx;
        _normal[1] = ny;
        _normal[2] = nz;
    }
};

class SURFE_LIB_EXPORT Tangent : public Point {
private:
    double _tangent[3];
    double _residual;
    double _angle_bound[2];
    double _inner_product_constraint;

public:
    Tangent(const double &x_coord, const double &y_coord, const double &z_coord,
            const double &tx, const double &ty, const double &tz,
            const double &c_coord = 0)
        : Point(x_coord, y_coord, z_coord, c_coord) {
        _tangent[0] = tx;
        _tangent[1] = ty;
        _tangent[2] = tz;
        _residual = 0.0;
        _inner_product_constraint = 0.0;  // default value 0.0 means that the
                                          // angle b/t the gradient of the
        // scalar field and tangent vector is 90 degrees.
    }
    double tx() const { return _tangent[0]; }
    double ty() const { return _tangent[1]; }
    double tz() const { return _tangent[2]; }
    double residual() const { return _residual; }
    double angle_lower_bound() const { return _angle_bound[0]; }
    double angle_upper_bound() const { return _angle_bound[1]; }
    double inner_product_constraint() const {
        return _inner_product_constraint;
    }
    void setResidual(const double &res) { _residual = res; }
    void setAngleBounds(const double &angle) {
        // t . del s = Cos(ϴ)*||t||*||del s||
        // ||t|| = 1
        // 0 <= ||del s|| <= + inf (but in reality ~ 2)
        double a = cos((90.0 - angle) * D2R) * 2.0;
        if (a < 0) {
            _angle_bound[0] = a;
            _angle_bound[1] = 0;
        } else {
            _angle_bound[0] = 0;
            _angle_bound[1] = a;
        }
    }
    void setInnerProductConstraint(const double &ip_constraint) {
        _inner_product_constraint = ip_constraint;
    }
};

class SURFE_LIB_EXPORT Constraints {
private:
    // Attributes
    double _avg_nn_dist_ie;
    double _avg_nn_dist_itr;
    double _avg_nn_dist_p;
    double _avg_nn_dist_t;
    // interface
    void _get_distinct_interface_iso_values();
    void _get_interface_points();
    bool _interface_points_are_coplanar() {
        return true;
    }  // Not implemented yet. should be tested when 2nd order polynomials are
       // used. Also when unisolvent points are used this should be called.
    // inequality
    std::set<double> _get_distinct_inequality_iso_values();

public:
	Constraints() {
		_avg_nn_dist_ie = -99999.0;   // no data value
		_avg_nn_dist_itr = -99999.0;  // no data value
		_avg_nn_dist_p = -99999.0;    // no data value
		_avg_nn_dist_t = -99999.0;    // no data value
	}
	//////////////////
	// Data members //
	//////////////////
	// 4 constraint types
	std::vector<Inequality> inequality;
	std::vector<Interface> itrface;
	std::vector<Planar> planar;
	std::vector<Tangent> tangent;

	// for interface data
	std::vector<double> interface_iso_values;
	std::vector<std::vector<Interface> > interface_point_lists;
	std::vector<Interface> interface_test_points;

	// evaluation points - where interpolant is evaluated
	std::vector<Evaluation_Point> evaluation_pts;

	/////////////
	// Methods //
	/////////////
	// Appending constraints
	void add_interface_constraint(const Point& pt) { itrface.emplace_back(pt); }
	void add_planar_constraint(const Planar& planar_pt) { planar.emplace_back(planar_pt); }
	void add_tangent_constraint(const Tangent& tangent_pt) { tangent.emplace_back(tangent_pt); }
	void add_inequality_constraint(const Inequality& inequality_pt) { inequality.emplace_back(inequality_pt); }

	// Clear current (if exists) evaluation_pts vector 
	void reset_evaluation_points() { evaluation_pts.clear(); };

    bool get_interface_data();  // fills interface_iso_values,
                                // interface_point_lists,
                                // interface_test_points data structures

    // validation
    bool check_input_data();

    // spatial analysis
    double compute_inequality_avg_nn_distance();
    double compute_interface_avg_nn_distance();
    double compute_planar_avg_nn_distance();
    double compute_tangent_avg_nn_distance();
    void compute_avg_nn_distances();
    double GetInequalityAvgNNDist() const { return _avg_nn_dist_ie; }
    double GetInterfaceAvgNNDist() const { return _avg_nn_dist_itr; }
    double GetPlanarAvgNNDist() const { return _avg_nn_dist_p; }
    double GetTangentAvgNNDist() const { return _avg_nn_dist_t; }
	void SetInequalityAvgNNDist(const double &dist) { _avg_nn_dist_ie = dist; }
	void SetInterfaceAvgNNDist(const double &dist) { _avg_nn_dist_itr = dist; }
	void SetPlanarAvgNNDist(const double &dist) { _avg_nn_dist_p = dist; }
	void SetTangentAvgNNDist(const double &dist) { _avg_nn_dist_t = dist; }
};


double distance_btw_pts(const Point &p1, const Point &p2);
int nearest_neighbour_index(const Point &p, const std::vector<Point> &pts);
std::vector<int> get_n_nearest_neighbours_to_point(const int &n, const Point &p, const std::vector<Point> &pts);
int furtherest_neighbour_index(const Point &p, const std::vector<Point> &pts);
int furtherest_neighbour_index(const std::vector<Point> &pts1, const std::vector<Point> &pts2);
double avg_nn_distance(const std::vector<Point> &pts);
bool Find_STL_Vector_Indices_FurtherestTwoPoints(const std::vector<Point> &pts, int (&TwoIndexes)[2]);
int Find_STL_Vector_Index_ofPointClosestToOtherPointWithinDistance(const Point &p, const std::vector<Point> &pts, const double &dist);
void calculate_bounds(const std::vector<Point> &pts, double (&bounds)[6]);
std::vector<int> get_extremal_point_data_indices_from_points(const std::vector<Point> &pts);
bool is_index_in_list(const int &index, const std::vector<int> &list);
bool find_fill_distance(const Constraints &input, double &fill_distance);
// The below functions will intelligently* get the indices within the STL vector
// of points that have large residuals Intelligently* : Doesn't blindly capture
// all points with large residuals
//                 - Considers the magnitude of the residuals
//                 - Considers the distance to other large residual points
//                 - Considers the variability with close large residual points
std::vector<int> Get_Inequality_STL_Vector_Indices_With_Large_Residuals(const std::vector<Inequality> &inequality, const double &avg_nn_distance);
std::vector<int> Get_Interface_STL_Vector_Indices_With_Large_Residuals(const std::vector<Interface> &itrface, const double &itrface_uncertainty,const double &avg_nn_distance);
std::vector<int> Get_Planar_STL_Vector_Indices_With_Large_Residuals(const std::vector<Planar> &planar, const double &angular_uncertainty,const double &avg_nn_distance);
std::vector<int> Get_Tangent_STL_Vector_Indices_With_Large_Residuals(const std::vector<Tangent> &tangent, const double &angular_uncertainty,const double &avg_nn_distance);
bool get_maximal_axial_variability_order(const double (&bounds)[6], Parameter_Types::AXIS (&axis_order)[3]);

#endif
