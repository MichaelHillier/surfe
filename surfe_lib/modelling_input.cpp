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

#include <modelling_parameters.h>
#include <continuous_property.h>
#include <math_methods.h>
#include <modeling_methods.h>
#include <modelling_input.h>
#include <algorithm>
#include <functional>

bool Constraints::get_interface_data()
{
	if (itrface.size() < 2) return false;
	
	// flush existing interface data containers - useful for greedy methods
	interface_iso_values.clear();
	interface_point_lists.clear();
	interface_test_points.clear();

	_get_distinct_interface_iso_values();
	_get_interface_points();

	return true;
}

bool Constraints::check_input_data()
{
	// check interface data...

	// check planar data ...

	// check tangent data ...

	// check inequality data ...

	// if using inequality data, check the level property data to ensure it is
	// consistent with interface data level property ...
	if (!inequality.empty()) 
	{
		std::vector<double> inequality_iso_values = _get_distinct_inequality_iso_values();
		if (inequality_iso_values.empty()) 
			return false;
		for (const auto &ineql_iso_value: inequality_iso_values){
			// if one of the inequality iso values is the same as the
			// itrface iso values data is not properly conditioned
			for (const auto &iter_iso_value: interface_iso_values){
				if (ineql_iso_value == iter_iso_value) 
					return false;
			}
		}
	}

	return true;
}

double Constraints::compute_inequality_avg_nn_distance()
{
	std::vector<Point> pts(inequality.begin(), inequality.end());
	return avg_nn_distance(pts);
}

double Constraints::compute_interface_avg_nn_distance()
{
	std::vector<Point> pts(itrface.begin(), itrface.end());
	return avg_nn_distance(pts);
}

double Constraints::compute_planar_avg_nn_distance()
{
	std::vector<Point> pts(planar.begin(), planar.end());
	return avg_nn_distance(pts);
}

double Constraints::compute_tangent_avg_nn_distance()
{
	std::vector<Point> pts(tangent.begin(), tangent.end());
	return avg_nn_distance(pts);
}

void Constraints::compute_avg_nn_distances() 
{
#pragma omp parallel sections
	{
#pragma omp section
		{ _avg_nn_dist_ie = compute_inequality_avg_nn_distance(); }
#pragma omp section
		{ _avg_nn_dist_itr = compute_interface_avg_nn_distance(); }
#pragma omp section
		{ _avg_nn_dist_p = compute_planar_avg_nn_distance(); }
#pragma omp section
		{ _avg_nn_dist_t = compute_tangent_avg_nn_distance(); }
	}
}

void Constraints::_get_distinct_interface_iso_values()
{
	std::set<double> distinct_iso_values;
	for (const auto &interface_constraint : itrface) 
		distinct_iso_values.insert(interface_constraint.level());
	// sort the vector (largest to smallest) - done for convenience and for functional reasons
	std::vector<double> distinct_iso_values_vec(distinct_iso_values.begin(), distinct_iso_values.end());
	std::sort(distinct_iso_values_vec.begin(), distinct_iso_values_vec.end(), std::greater<double>());
	for (const auto& value : distinct_iso_values_vec)
		interface_iso_values.push_back(value);
}

void Constraints::_get_interface_points()
{
	// interface[0][0,1,2,3,....] points 0,1,2,3,.... belong to the 0th
	// interface
	// ...
	// interface[m = interface_iso_values.size()][76,45,43,4,.....] points
	// 76,45,43,4,..... belong to the mth interface
	interface_point_lists.resize(interface_iso_values.size());
	for (int j = 0; j < (int)interface_iso_values.size(); j++) {
		for (const auto &interface_pt: itrface){
			if (interface_pt.level() == interface_iso_values.at(j))
			{
				// add to 2D vector
				interface_point_lists[j].push_back(interface_pt);
			}
		}
	}

	for (int j = 0; j < (int)interface_point_lists.size(); j++) {
		// set the test_interface_points
		interface_test_points.push_back(interface_point_lists[j][0]);
		if ((int)interface_point_lists.at(j).size() == 1)
		{
			// need to have at least 2 points per interface
			// remove this interface from the list
			interface_point_lists.erase(interface_point_lists.begin() + j);
			j--;
		}
	}
}

std::vector<double> Constraints::_get_distinct_inequality_iso_values()
{
	std::set<double> distinct_iso_values;
	for (const auto &inequality_constraint : inequality) distinct_iso_values.insert(inequality_constraint.level());
	// sort the vector (largest to smallest) - done for convenience and for functional reasons
	std::vector<double> distinct_iso_values_vec(distinct_iso_values.begin(), distinct_iso_values.end());
	std::sort(distinct_iso_values_vec.begin(), distinct_iso_values_vec.end(), std::greater<double>());
	return distinct_iso_values_vec;
}

bool Planar::_compute_strike_dip_polarity_from_normal()
{
	if (_normal[2] < 0) 
		_polarity = 1;
	else 
		_polarity = 0;

	// get dip first
	_dip = acos(_normal[2]) * R2D;  // could do better. e.g. for overturn cases
									// puts _dip > 90. but there
	// formula's for getting normals works so sticking with it for now
	// get dip_direction
	double dip_direction = atan2(_normal[1], _normal[0]) * R2D;

	// if negative azimuth get positive angle
	if (dip_direction < 0) 
		dip_direction += 360;

	// get strike
	_strike = 360 - dip_direction;

	return true;  // check this computation
}

bool Planar::_compute_normal_from_strike_dip_polarity()
{
	// Get down dip vector - v
	double vx = cos(-1.0 * (_strike * D2R)) * cos(-1.0 * (_dip * D2R));
	double vy = sin(-1.0 * (_strike * D2R)) * cos(-1.0 * (_dip * D2R));
	double vz = sin(-1.0 * (_dip * D2R));

	// Get strike vector - vp
	double vpx = -1.0 * vy;
	double vpy = vx;
	double vpz = 0;

	// Compute cross product of v x vp
	double Nx = (vy * vpz) - (vz * vpy);
	double Ny = (vz * vpx) - (vx * vpz);
	double Nz = (vx * vpy) - (vy * vpx);

	// Compute Length
	double Length = sqrt(Nx * Nx + Ny * Ny + Nz * Nz);

	// Normalize
	Nx = Nx / Length;
	Ny = Ny / Length;
	Nz = Nz / Length;

	// Check polarity of computed normal is consistent with the property
	// "polarity" polarity == 1 . Is overturned: Must point downward polarity
	// ==
	// 0 . Is upright: Must point upward polarity != 0/1 . Is unknown: Can
	// point
	// upward or downward, it doesn't matter.
	//                    Retains its initialized polarity set above.
	if ((_polarity == 1 && Nz > 0) || (_polarity == 0 && Nz < 0)) {
		// Flip vector
		Nx = -Nx;
		Ny = -Ny;
		Nz = -Nz;
	}

	// assign normal
	_normal[0] = Nx;
	_normal[1] = Ny;
	_normal[2] = Nz;

	return true;
}

bool Planar::getDipVector(double(&vector)[3])
{
	// Get down dip vector - v
	double vx = cos(-1.0 * (_strike * D2R)) * cos(-1.0 * (_dip * D2R));
	double vy = sin(-1.0 * (_strike * D2R)) * cos(-1.0 * (_dip * D2R));
	double vz = sin(-1.0 * (_dip * D2R));

	// normalize
	double length = sqrt(vx * vx + vy * vy + vz * vz);
	vx /= length;
	vy /= length;
	vz /= length;

	vector[0] = vx;
	vector[1] = vy;
	vector[2] = vz;

	return true;
}

bool Planar::getStrikeVector(double(&vector)[3])
{
	// Get strike vector
	double vx = -sin(-1.0 * (_strike * D2R));
	double vy = cos(-1.0 * (_strike * D2R));
	double vz = 0;

	vector[0] = vx;
	vector[1] = vy;
	vector[2] = vz;

	return true;
}

void Planar::setNormalBounds(const double &delta_strike,const double &delta_dip)
{
	// theta = strike
	// phi = dip
	// delta_strike = error on strike measurement
	// delta_dip    = error on dip measurement

	// given the errors on the strike/dip angles a matrix will be computed to
	// hold
	// the lower and upper bounds on the computed normal components: nx_lower <
	// nx
	// < nx_upper ny_lower < ny < ny_upper nz_lower < nz < nz_upper

	double theta = _strike * D2R;
	double phi = _dip * D2R;
	double dtheta = delta_strike * D2R;
	double dphi = delta_dip * D2R;

	double eqn1x = cos(dtheta + theta) * sin(dphi + phi);
	double eqn1y = -1.0 * sin(dtheta + theta) * sin(dphi + phi);
	double eqn1z = cos(dphi + phi);

	double eqn2x = -1.0 * cos(dtheta - theta) * sin(dphi - phi);
	double eqn2y = -1.0 * sin(dtheta - theta) * sin(dphi - phi);
	double eqn2z = cos(dphi - phi);

	double eqn3x = -1.0 * cos(dtheta + theta) * sin(dphi - phi);
	double eqn3y = sin(dtheta + theta) * sin(dphi - phi);
	double eqn3z = cos(dphi - phi);

	double eqn4x = cos(dtheta - theta) * sin(dphi + phi);
	double eqn4y = sin(dtheta - theta) * sin(dphi + phi);
	double eqn4z = cos(dphi + phi);

	double nx_lower = eqn1x;
	if (eqn2x < nx_lower) nx_lower = eqn2x;
	if (eqn3x < nx_lower) nx_lower = eqn3x;
	if (eqn4x < nx_lower) nx_lower = eqn4x;
	double ny_lower = eqn1y;
	if (eqn2y < ny_lower) ny_lower = eqn2y;
	if (eqn3y < ny_lower) ny_lower = eqn3y;
	if (eqn4y < ny_lower) ny_lower = eqn4y;
	double nz_lower = eqn1z;
	if (eqn2z < nz_lower) nz_lower = eqn2z;
	if (eqn3z < nz_lower) nz_lower = eqn3z;
	if (eqn4z < nz_lower) nz_lower = eqn4z;

	double nx_upper = eqn1x;
	if (eqn2x > nx_upper) nx_upper = eqn2x;
	if (eqn3x > nx_upper) nx_upper = eqn3x;
	if (eqn4x > nx_upper) nx_upper = eqn4x;
	double ny_upper = eqn1y;
	if (eqn2y > ny_upper) ny_upper = eqn2y;
	if (eqn3y > ny_upper) ny_upper = eqn3y;
	if (eqn4y > ny_upper) ny_upper = eqn4y;
	double nz_upper = eqn1z;
	if (eqn2z > nz_upper) nz_upper = eqn2z;
	if (eqn3z > nz_upper) nz_upper = eqn3z;
	if (eqn4z > nz_upper) nz_upper = eqn4z;

	double l_lower = sqrt(nx_lower * nx_lower + ny_lower * ny_lower + nz_lower * nz_lower);  // length of lower bound vector
	double l_upper = sqrt(nx_upper * nx_upper + ny_upper * ny_upper + nz_upper * nz_upper);  // length of upper bound vector

	_normal_bound[0][0] = nx_lower;
	_normal_bound[0][1] = nx_upper;
	_normal_bound[1][0] = ny_lower;
	_normal_bound[1][1] = ny_upper;
	_normal_bound[2][0] = nz_lower;
	_normal_bound[2][1] = nz_upper;
}

std::vector<Point> convert_constraints_to_points(const Constraints& constraints)
{
	std::vector<Point> pts;
	pts.insert(std::end(pts), std::begin(constraints.inequality), std::end(constraints.inequality));
	pts.insert(std::end(pts), std::begin(constraints.itrface), std::end(constraints.itrface));
	pts.insert(std::end(pts), std::begin(constraints.planar), std::end(constraints.planar));
	pts.insert(std::end(pts), std::begin(constraints.tangent), std::end(constraints.tangent));
	return pts;
}

double distance_btw_pts(const Point &p1, const Point &p2) 
{
	double dx = p1.x() - p2.x();
	double dy = p1.y() - p2.y();
	double dz = p1.z() - p2.z();
	double dc = p1.c() - p2.c();

	return sqrt(dx * dx + dy * dy + dz * dz + dc * dc);
}

int nearest_neighbour_index(const Point &p,	const std::vector<Point> &pts)
{
	double min_distance = DBL_MAX;
	int index =	-1;  // if this is actually returned by this function a seg fault will result

	// calculate all non-zero distances
	for (int j = 0; j < (int)pts.size(); j++) {
		double distance = distance_btw_pts(p, pts[j]);
		if (distance != 0)  // in case p is in the set of pts
		{
			if (distance < min_distance) 
			{
				min_distance = distance;
				index = j;
			}
		}
	}

	return index;
}
std::vector<int> get_n_nearest_neighbours_to_point(const int &n, const Point &p, const std::vector<Point> &pts) {
	std::vector<int> nn_indices;

	int n_pts = (int)pts.size();

	// calculate all non-zero distances + create index array
	std::vector<double> dist_arr;
	std::vector<int> index_arr;
	for (int j = 0; j < (int)pts.size(); j++) {
		double distance = distance_btw_pts(p, pts[j]);
		if (distance != 0) 
		{
			dist_arr.push_back(distance);
			index_arr.push_back(j);
		}
	}
	// sort distances and corresponding indices
	Math_methods::sort_vector_w_index(dist_arr,	index_arr);  // smallest to largest

	if (n <= n_pts)
		for (int j = 0; j < n; j++) 
			nn_indices.push_back(index_arr[j]);
	else
		for (int j = 0; j < n_pts; j++)
			nn_indices.push_back(index_arr[j]);  // will return less than n nn index. best option vs getting seg fault

	return nn_indices;
}

int furtherest_neighbour_index(const Point &p, const std::vector<Point> &pts) 
{
	int index = 0;
	double largest_distance = distance_btw_pts(p, pts[0]);
	for (int j = 1; j < (int)pts.size(); j++) {
		double distance_j = distance_btw_pts(p, pts[j]);
		if (distance_j > largest_distance)
		{
			largest_distance = distance_j;
			index = j;
		}
	}
	return index;
}

int furtherest_neighbour_index(const std::vector<Point> &pts1, const std::vector<Point> &pts2) {
	// find a point in PTS1 dataset that is furthest away from PTS2 dataset

	int index = 0;
	double largest_distance = 0.0;
	for (int j = 0; j < (int)pts1.size(); j++) {
		for (int k = 0; k < (int)pts2.size(); k++) {
			double distance_jk = distance_btw_pts(pts1[j], pts2[k]);
			if (distance_jk > largest_distance) {
				largest_distance = distance_jk;
				index = j;
			}
		}
	}
	return index;
}

double avg_nn_distance(const std::vector<Point> &pts)
{
	double average_nn_distance = 0.0;
	int n = (int)pts.size();
	for (int j = 0; j < n; j++) {
		double min_dist = DBL_MAX;
		for (int k = 0; k < n; k++) {
			if (k != j) {
				double dist = distance_btw_pts(pts.at(j), pts.at(k));
				if (dist < min_dist) min_dist = dist;
			}
		}
		if (n == 1) 
			min_dist = 0;  // trap this edge case
		average_nn_distance += min_dist;
	}
	if (n != 0) 
		average_nn_distance /= n;
	return average_nn_distance;
}

bool Find_STL_Vector_Indices_FurtherestTwoPoints(const std::vector<Point> &pts, int(&TwoIndexes)[2]) 
{
	if (pts.size() < 2) return false;

	double largest_distance = -DBL_MAX;
	for (int j = 0; j < (int)pts.size(); j++) {
		for (int k = 0; k < (int)pts.size(); k++) {
			double distance_jk = distance_btw_pts(pts[j], pts[k]);
			if (distance_jk > largest_distance) {
				largest_distance = distance_jk;
				TwoIndexes[0] = j;
				TwoIndexes[1] = k;
			}
		}
	}

	return true;
}

int Find_STL_Vector_Index_ofPointClosestToOtherPointWithinDistance(const Point &p, const std::vector<Point> &pts, const double &dist)
{
	double smallest_residual = DBL_MAX;
	int index = -1;
	for (int j = 0; j < (int)pts.size(); j++) {
		double distance = distance_btw_pts(p, pts[j]);
		double residual = abs(distance - dist);
		if (residual < smallest_residual) {
			smallest_residual = residual;
			index = j;
		}
	}

	return index;
}

void calculate_bounds(const std::vector<Point> &pts, double(&bounds)[6])
{
	// bounds[6] = { xmin, xmax, ymin, ymax, zmin, zmax }
	bounds[0] = pts[0].x();
	bounds[1] = pts[0].x();
	bounds[2] = pts[0].y();
	bounds[3] = pts[0].y();
	bounds[4] = pts[0].z();
	bounds[5] = pts[0].z();

	for (int j = 1; j < (int)pts.size(); j++) {
		if (pts[j].x() < bounds[0]) bounds[0] = pts[j].x();  // x-min
		if (pts[j].x() > bounds[1]) bounds[1] = pts[j].x();  // x-max
		if (pts[j].y() < bounds[2]) bounds[2] = pts[j].y();  // y-min
		if (pts[j].y() > bounds[3]) bounds[3] = pts[j].y();  // y-max
		if (pts[j].z() < bounds[4]) bounds[4] = pts[j].z();  // z-min
		if (pts[j].z() > bounds[5]) bounds[5] = pts[j].z();  // z-max
	}
}

std::vector<int> get_extremal_point_data_indices_from_points(const std::vector<Point> &pts) 
{
	// return a vector of indices from pts[] that maximally samples the data
	// space
	int n = (int)pts.size();

	// sort x-coord, y-coord, and z-coord data and corresponding indices
	std::vector<double> x_coord;
	std::vector<double> y_coord;
	std::vector<double> z_coord;
	std::vector<int> x_idx;
	std::vector<int> y_idx;
	std::vector<int> z_idx;
	for (int j = 0; j < n; j++) {
		x_coord.push_back(pts[j].x());
		y_coord.push_back(pts[j].y());
		z_coord.push_back(pts[j].z());
		x_idx.push_back(j);
		y_idx.push_back(j);
		z_idx.push_back(j);
	}
	Math_methods::sort_vector_w_index(x_coord, x_idx);
	Math_methods::sort_vector_w_index(y_coord, y_idx);
	Math_methods::sort_vector_w_index(z_coord, z_idx);

	// compute the order of axial sampling variability, largest to smallest
	std::vector<double> ranges;
	std::vector<int> axis;
	ranges.push_back(x_coord[n - 1] - x_coord[0]);  // x-range
	ranges.push_back(y_coord[n - 1] - y_coord[0]);  // y-range
	ranges.push_back(z_coord[n - 1] - z_coord[0]);  // z-range
	for (int j = 0; j < 3; j++) axis.push_back(j);  // x-axis = 0, y-axis = 1, z-axis = 2
	Math_methods::sort_vector_w_index(ranges, axis);  // sort smallest to largest
	int axis_order[3] = { axis[2], axis[1], axis[0] };  // re-arrange largest to smallest

	std::vector<int> data_indices;  // result container
	for (int j = 0; j < 3; j++) {
		if (axis_order[j] == 0)  // x-axis
		{
			if (!is_index_in_list(x_idx[0], data_indices))
				data_indices.push_back(x_idx[0]);
			else {
				for (int k = 1; k < n; k++) {
					if (!is_index_in_list(x_idx[k], data_indices)) {
						data_indices.push_back(x_idx[k]);
						break;
					}
				}
			}

			if (!is_index_in_list(x_idx[n - 1], data_indices))
				data_indices.push_back(x_idx[n - 1]);
			else {
				for (int k = 1; k < n; k++) {
					if (!is_index_in_list(x_idx[n - k - 1], data_indices)) {
						data_indices.push_back(x_idx[n - k - 1]);
						break;
					}
				}
			}
		}
		if (axis_order[j] == 1)  // y-axis
		{
			if (!is_index_in_list(y_idx[0], data_indices))
				data_indices.push_back(y_idx[0]);
			else {
				for (int k = 1; k < n; k++) {
					if (!is_index_in_list(y_idx[k], data_indices)) {
						data_indices.push_back(y_idx[k]);
						break;
					}
				}
			}

			if (!is_index_in_list(y_idx[n - 1], data_indices))
				data_indices.push_back(y_idx[n - 1]);
			else {
				for (int k = 1; k < n; k++) {
					if (!is_index_in_list(y_idx[n - k - 1], data_indices)) {
						data_indices.push_back(y_idx[n - k - 1]);
						break;
					}
				}
			}
		}
		if (axis_order[j] == 2)  // z-axis
		{
			if (!is_index_in_list(z_idx[0], data_indices))
				data_indices.push_back(z_idx[0]);
			else {
				for (int k = 1; k < n; k++) {
					if (!is_index_in_list(z_idx[k], data_indices)) {
						data_indices.push_back(z_idx[k]);
						break;
					}
				}
			}

			if (!is_index_in_list(z_idx[n - 1], data_indices))
				data_indices.push_back(z_idx[n - 1]);
			else {
				for (int k = 1; k < n; k++) {
					if (!is_index_in_list(z_idx[n - k - 1], data_indices)) {
						data_indices.push_back(z_idx[n - k - 1]);
						break;
					}
				}
			}
		}
	}

	return data_indices;
}

bool get_maximal_axial_variability_order(const double(&bounds)[6], Parameter_Types::AXIS(&axis_order)[3])
{
	std::vector<double> ranges;
	std::vector<int> axis;
	for (int j = 0; j < 3; j++) {
		ranges.push_back(abs(bounds[2 * j] - bounds[2 * j + 1]));
		axis.push_back(j);
	}
	Math_methods::sort_vector_w_index(ranges, axis);

	int idx = 0;
	for (int j = 2; j >= 0; j--) {
		if (axis[j] == 0) axis_order[idx] = Parameter_Types::AXIS::Xaxis;
		if (axis[j] == 1) axis_order[idx] = Parameter_Types::AXIS::Yaxis;
		if (axis[j] == 2) axis_order[idx] = Parameter_Types::AXIS::Zaxis;
		idx++;
	}

	return true;
}

bool is_index_in_list(const int &index, const std::vector<int> &list) 
{
	for (const auto &list_index : list){
		if (list_index == index) 
			return true;
	}
	return false;
}

double get_largest_distance_between_points(const std::vector<Point> &pts)
{
	double largest_distance = 0.0;
	for (const auto& pt_a : pts) {
		for (const auto& pt_b : pts) {
			double d = distance_btw_pts(pt_a, pt_b);
			if (d > largest_distance)
				largest_distance = d;
		}
	}
	return largest_distance;
}

std::vector<int> Get_Inequality_STL_Vector_Indices_With_Large_Residuals(
	const std::vector<Inequality> &inequality, const double &avg_nn_distance)
{
	// Function will intelligently* get the indices within the STL vector of
	// Inequality points that have large residuals Intelligently* : Doesn't
	// blindly capture all points with large residuals
	//                 - Considers the distance to other large residual points
	//               These are a bit special
	std::vector<int> inequality_indices_to_include;  // what we are going to
	// function on function exit

	double Large_distance = DBL_MAX;

	std::vector<int> inequality_residuals_indices;
	for (int j = 0; j < (int)inequality.size(); j++) {
		if (!inequality.at(j).residual())
			inequality_residuals_indices.push_back(j);
	}
	if (!inequality_residuals_indices.empty()) {
		// Will always accept the first residual over the threshold
		inequality_indices_to_include.push_back(inequality_residuals_indices[0]);
		inequality_residuals_indices.pop_back();

		for (const auto &index : inequality_residuals_indices){
			// find closest point currently in inequality_indices_to_include[]
			// to
			// index
			double nn_dist = Large_distance;
			int nn_index = -1;  // should always be overwritten. Will get a seg fault if this isn't the case
			for (const auto &ineql_index_to_include: inequality_indices_to_include){
				double dist = distance_btw_pts(inequality[index],inequality[ineql_index_to_include]);
				if (dist < nn_dist) 
				{
					nn_dist = dist;
					nn_index = ineql_index_to_include;
				}
			}
			// Distance to other Large Residual Points Condition
			if (nn_dist > avg_nn_distance)
				inequality_indices_to_include.push_back(index);
		}
	}

	std::sort(inequality_indices_to_include.begin(),inequality_indices_to_include.end());

	return inequality_indices_to_include;
}

std::vector<int> Get_Interface_STL_Vector_Indices_With_Large_Residuals(
	const std::vector<Interface> &itrface, const double &itrface_uncertainty, const double &avg_nn_distance) 
{
	// Function will intelligently* get the indices within the STL vector of
	// Interface points that have large residuals Intelligently* : Doesn't
	// blindly
	// capture all points with large residuals
	//                 - Considers the magnitude of the residuals
	//                 - Considers the distance to other large residual points
	//                 - Considers the variability with close large residual
	//                 points

	std::vector<int> itrface_indices_to_include;  // what we are going to
												  // function
												  // on function exit

	double Large_distance = DBL_MAX;

	double largest_residual = 0;
	int ref_index = -1;
	std::vector<double> large_itrface_residuals;
	std::vector<int> large_itrface_residuals_indices;
	for (int j = 0; j < (int)itrface.size(); j++) {
		double error = itrface.at(j).residual();
		if (error > itrface_uncertainty) {
			if (error > largest_residual) {
				largest_residual = error;
				ref_index = j;
			}
			large_itrface_residuals.push_back(error);
			large_itrface_residuals_indices.push_back(j);
		}
	}
	// 	if (ref_index != -1 )
	// 	{
	// 		itrface_indices_to_include.push_back(ref_index);
	// 	}
	if (!large_itrface_residuals.empty())
	{
		// Residual Magnitude Condition
		// sort all residuals over the threshold(angular_uncertainty) in
		// smallest to
		// largest - along with linked indices
		Math_methods::sort_vector_w_index(large_itrface_residuals,large_itrface_residuals_indices);
		// Will always accept largest residual over the threshold
		itrface_indices_to_include.push_back(large_itrface_residuals_indices[(int)large_itrface_residuals_indices.size() - 1]);
		large_itrface_residuals.pop_back();  // probably don't need to do this since we never use it again
		large_itrface_residuals_indices.pop_back();

		for (int j = 0; j < (int)large_itrface_residuals_indices.size(); j++) {
			// current index being analyzed (next largest residual in list):
			int index = large_itrface_residuals_indices[(int)large_itrface_residuals_indices.size() - j - 1];
			// find closest point currently in itrface_indices_to_include[] to
			// index
			double nn_dist = Large_distance;
			int nn_index = -1;  // should always be overwritten. Will get a seg fault if this isn't the case
			for (const auto &itr_index_to_include : itrface_indices_to_include){
				double dist = distance_btw_pts(itrface[index],itrface[itr_index_to_include]);
				if (dist < nn_dist) {
					nn_dist = dist;
					nn_index = itr_index_to_include;
				}
			}
			// Distance to other Large Residual Points Condition
			if (nn_dist > avg_nn_distance)
				itrface_indices_to_include.push_back(index);
		}
	}

	std::sort(itrface_indices_to_include.begin(),itrface_indices_to_include.end());

	return itrface_indices_to_include;
}

std::vector<int> Get_Planar_STL_Vector_Indices_With_Large_Residuals(
	const std::vector<Planar> &planar, const double &angular_uncertainty, const double &avg_nn_distance) {
	// Function will intelligently* get the indices within the STL vector of
	// Planar points that have large residuals Intelligently* : Doesn't blindly
	// capture all points with large residuals
	//                 - Considers the magnitude of the residuals
	//                 - Considers the distance to other large residual points
	//                 - Considers the variability with close large residual
	//                 points -- This has been taken out for now.. needs more
	//                 testing

	std::vector<int> planar_indices_to_include;  // what we are going to function on function exit

	double Large_distance = DBL_MAX;

	double largest_residual = 0;
	int ref_index = -1;
	std::vector<double> large_planar_residuals;
	std::vector<int> large_planar_residuals_indices;
	for (int j = 0; j < (int)planar.size(); j++) {
		double grad_err = planar.at(j).residual() * R2D;
		if (grad_err > angular_uncertainty) {
			if (grad_err > largest_residual) {
				largest_residual = grad_err;
				ref_index = j;
			}
			large_planar_residuals.push_back(grad_err);
			large_planar_residuals_indices.push_back(j);
		}
	}
	// 	if (ref_index != -1 )
	// 	{
	// 		planar_indices_to_include.push_back(ref_index);
	// 	}
	if (!large_planar_residuals.empty()) 
	{
		// Residual Magnitude Condition
		// sort all residuals over the threshold(angular_uncertainty) in
		// smallest to
		// largest - along with linked indices
		Math_methods::sort_vector_w_index(large_planar_residuals,large_planar_residuals_indices);
		// Will always accept largest residual over the threshold
		planar_indices_to_include.push_back(large_planar_residuals_indices[(int)large_planar_residuals_indices.size() - 1]);
		large_planar_residuals.pop_back();  // probably don't need to do this since we never use it again
		large_planar_residuals_indices.pop_back();

		for (int j = 0; j < (int)large_planar_residuals_indices.size(); j++) {
			// current index being analyzed (next largest residual in list):
			int index = large_planar_residuals_indices[(int)large_planar_residuals_indices.size() - j - 1];
			// find closest point currently in planar_indices_to_include[] to
			// index
			double nn_dist = Large_distance;
			int nn_index = -1;  // should always be overwritten. Will get a seg fault if this isn't the case
			for (const auto &planar_index_to_include : planar_indices_to_include){
				double dist = distance_btw_pts(planar[index], planar[planar_index_to_include]);
				if (dist < nn_dist) 
				{
					nn_dist = dist;
					nn_index = planar_index_to_include;
				}
			}
			// Distance to other Large Residual Points Condition
			if (nn_dist > avg_nn_distance)
				planar_indices_to_include.push_back(index);
		}
	}

	std::sort(planar_indices_to_include.begin(),planar_indices_to_include.end());

	return planar_indices_to_include;
}

std::vector<int> Get_Tangent_STL_Vector_Indices_With_Large_Residuals(
	const std::vector<Tangent> &tangent, const double &angular_uncertainty, const double &avg_nn_distance)
{
	// Function will intelligently* get the indices within the STL vector of
	// tangent points that have large residuals Intelligently* : Doesn't blindly
	// capture all points with large residuals
	//                 - Considers the magnitude of the residuals
	//                 - Considers the distance to other large residual points
	//                 - Considers the variability with close large residual
	//                 points

	std::vector<int> tangent_indices_to_include;  // what we are going to function on function exit

	double Large_distance = DBL_MAX;

	std::vector<double> large_tangent_residuals;
	std::vector<int> large_tangent_residuals_indices;
	for (int j = 0; j < (int)tangent.size(); j++) {
		double grad_err = tangent.at(j).residual() * R2D;
		if (grad_err > angular_uncertainty)
		{
			large_tangent_residuals.push_back(grad_err);
			large_tangent_residuals_indices.push_back(j);
		}
	}
	if (!large_tangent_residuals.empty()) {
		// Residual Magnitude Condition
		// sort all residuals over the threshold(angular_uncertainty) in
		// smallest to
		// largest - along with linked indices
		Math_methods::sort_vector_w_index(large_tangent_residuals, large_tangent_residuals_indices);
		// Will always accept largest residual over the threshold
		tangent_indices_to_include.push_back(large_tangent_residuals_indices[(int)large_tangent_residuals_indices.size() - 1]);
		large_tangent_residuals.pop_back();  // probably don't need to do this since we never use it again
		large_tangent_residuals_indices.pop_back();

		for (int j = 0; j < (int)large_tangent_residuals_indices.size(); j++) {
			// current index being analyzed (next largest residual in list):
			int index = large_tangent_residuals_indices[(int)large_tangent_residuals_indices.size() - j - 1];
			// find closest point currently in tangent_indices_to_include[] to
			// index
			double nn_dist = Large_distance;
			int nn_index = -1;  // should always be overwritten. Will get a seg fault if this isn't the case
			for (const auto &tangent_index_to_include : tangent_indices_to_include) {
				for (int k = 0; k < (int)tangent_indices_to_include.size(); k++) {
					double dist = distance_btw_pts(tangent[index], tangent[tangent_index_to_include]);
					if (dist < nn_dist)
					{
						nn_dist = dist;
						nn_index = tangent_index_to_include;
					}
				}
				// Distance to other Large Residual Points Condition
				if (nn_dist > avg_nn_distance)
					tangent_indices_to_include.push_back(index);
			}
		}

		std::sort(tangent_indices_to_include.begin(), tangent_indices_to_include.end());

		return tangent_indices_to_include;
	}
}