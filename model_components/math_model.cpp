
//This file contains mathematical functions

#include "math_structures.cpp"
#include <cmath>


Points ls_intersections(Line line, Sphere sphere) {
    //This function calculates the intersection points between a line and a sphere
    //The function returns a vector of points

    Points intersection_points;

    double dx = line.dx;
    double dy = line.dy;
    double dz = line.dz;

    double cx = line.point.x;
    double cy = line.point.y;
    double cz = line.point.z;

    double h = sphere.centre.x;
    double k = sphere.centre.y;
    double l = sphere.centre.z;

    double r = sphere.radius;

    //calculate the coefficients of the quadratic equation
    double A = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
    double B = 2 * dx * (cx - h) + 2 * dy * (cy - k) + 2 * dz * (cz - l);
    double C = pow(cx - h, 2) + pow(cy - k, 2) + pow(cz - l, 2) - pow(r, 2);

    //calculate the discriminant
    double discriminant = pow(B, 2) - 4 * A * C;

    //if the discriminant is negative, there are no real solutions
    if (discriminant < 0) {
        cout << "No real solutions" << endl;
        return intersection_points;
    }

    //if the discriminant is zero, there is one real solution
    if (discriminant == 0) {
        double t = -B / (2 * A);

        Point intersection_point;
        intersection_point.x = dx * t + cx;
        intersection_point.y = dy * t + cy;
        intersection_point.z = dz * t + cz;

        intersection_points.push_back(intersection_point);

        return intersection_points;
    };

    //if the discriminant is positive, there are two real solutions
    if (discriminant > 0) {
        double t1 = (-B + sqrt(discriminant)) / (2 * A);
        double t2 = (-B - sqrt(discriminant)) / (2 * A);

        Point intersection_point1;
        intersection_point1.x = dx * t1 + cx;
        intersection_point1.y = dy * t1 + cy;
        intersection_point1.z = dz * t1 + cz;

        Point intersection_point2;
        intersection_point2.x = dx * t2 + cx;
        intersection_point2.y = dy * t2 + cy;
        intersection_point2.z = dz * t2 + cz;

        intersection_points.push_back(intersection_point1);
        intersection_points.push_back(intersection_point2);

        

    };

    return intersection_points;

};


Magnitudes find_comp_contributions(Velocity i, Velocity j) {

    //This function calculates the contributions of each component of the velocity vector of i to the velocity vector of j
    //The function returns a vector of magnitudes

    Magnitudes contributions;

    //1. unpack the absolute velocities of i and j

    double i_vx = sqrt(pow(i.vx, 2));
    double i_vy = sqrt(pow(i.vy, 2));
    double i_vz = sqrt(pow(i.vz, 2));

    double j_vx = sqrt(pow(j.vx, 2));
    double j_vy = sqrt(pow(j.vy, 2));
    double j_vz = sqrt(pow(j.vz, 2));

    //2. Find the total magnitude of the system by taking the sum of squares

    double total_magnitude = i_vx + i_vy + i_vz + j_vx + j_vy + j_vz;

    //3. Calculate the contributions of each component of i and i to the total magnitude
    Magnitude i_contributions;
    Magnitude j_contributions;

    //push them to the contributions vector of Magnitude

    contributions.magnitudes.push_back(i_contributions);
    contributions.magnitudes.push_back(j_contributions);


    contributions.magnitudes[0].x_mag = (i_vx / total_magnitude);
    contributions.magnitudes[0].y_mag = (i_vy / total_magnitude);
    contributions.magnitudes[0].z_mag = (i_vz / total_magnitude);

    contributions.magnitudes[1].x_mag =(j_vx / total_magnitude);
    contributions.magnitudes[1].y_mag =(j_vy / total_magnitude);
    contributions.magnitudes[1].z_mag = (j_vz / total_magnitude);

    return contributions;
};


PPoints polar_converter(Points cartesian_points){
    //This function converts a vector of cartesian points to polar points and scales them.

    PPoints polar_points;

    //1. Convert the cartesian points to polar points



    for (int i = 0; i < cartesian_points.points.size(); i++) {
        Point cartesian_point = cartesian_points.points[i];


        
        double r = sqrt(pow(cartesian_point.x, 2) + pow(cartesian_point.y, 2) + pow(cartesian_point.z, 2));
        double theta;
        double phi;

        if (r == 0) {
            theta = 0;
            phi = 0;
            
        } else {
            theta = atan2(cartesian_point.y,cartesian_point.x);
            phi = acos(cartesian_point.z / r);
        };

        PPoint polar_point;
        polar_point.r = r;
        polar_point.theta = theta;
        polar_point.phi = phi;

        polar_points.push_back(polar_point);

    };

    //4. Return the polar points

    return polar_points;
    
};

Points cartesian_converter(PPoints polar_points) {



    //This function converts a vector of polar points to cartesian points and scales them back to the original scale.


    Points cartesian_points;

    //1. Convert the polar points to cartesian points

    for (int i = 0; i < polar_points.ppoints.size(); i++) {
        PPoint polar_point = polar_points.ppoints[i];

        Point cartesian_point;

        cartesian_point.x = polar_point.r * cos(polar_point.theta) * sin(polar_point.phi);
        cartesian_point.y = polar_point.r * sin(polar_point.theta) * sin(polar_point.phi);
        cartesian_point.z = polar_point.r * cos(polar_point.phi);

        cartesian_points.push_back(cartesian_point);

    };
    
    return cartesian_points;

};

_Points rotate_cartesians(Points cartesian_points) {
    //This function rotates a vector of cartesian points to a new coordinate system.

    _Points _rotated_points;
    Points rotated_points = cartesian_points;

    // Extract the first two points to calculate the midpoint
    Points ij_points;
    ij_points.points.push_back(cartesian_points.points[0]);
    ij_points.points.push_back(cartesian_points.points[1]);

    // Calculate the midpoint
    Point midpoint_ij = ij_points.calculate_midpoint();

    // Shift the first two points to the origin

    rotated_points.points[0] = rotated_points.points[0].subtract(midpoint_ij);
    rotated_points.points[1] = rotated_points.points[1].subtract(midpoint_ij);

    // Calculate rotation angles based on shifted points
    auto delta = rotated_points.points[1].subtract(rotated_points.points[0]);
    double theta_z = atan2(delta.y, delta.x);
    double theta_y = atan2(delta.z, sqrt(pow(delta.x, 2) + pow(delta.y, 2)));



    // Rotate points around the origin
    rotated_points.rotate_points(-theta_z, theta_y); 

    

    // Store the results
    _rotated_points.points = rotated_points.points;
    _rotated_points.rotation_angles = {theta_z, theta_y};
    _rotated_points.rotation_intercept = midpoint_ij;

    return _rotated_points;
}

Points unrotate_cartesians(const std::vector<Point>& rotated_points, const Rotator& angles, const Point& intercept) {
    Points unrotated_points;
    unrotated_points.points = rotated_points;

    // Rotate the points back by the given angles
    unrotated_points.undo_rotate_points(-angles.theta, -angles.phi);

    // Add the intercept back to the position points (first two elements)
    unrotated_points.points[0] = unrotated_points.points[0].add(intercept);
    unrotated_points.points[1] = unrotated_points.points[1].add(intercept);

    return unrotated_points;
}

void compare_points(const Points& original_points, const Points& restored_points) {

    double aggregate_difference = 0;
    for (size_t i = 0; i < original_points.points.size(); ++i) {
        double dx = std::abs(original_points.points[i].x - restored_points.points[i].x);
        double dy = std::abs(original_points.points[i].y - restored_points.points[i].y);
        double dz = std::abs(original_points.points[i].z - restored_points.points[i].z);

        double total_difference = dx + dy + dz;
        aggregate_difference += total_difference;

    

        
    }

    cout << "The aggregate difference between the original and restored points is: " << aggregate_difference << endl;
}


double custom_max(double a, double b, double c) {
    return max({isnan(a) ? -numeric_limits<double>::infinity() : a, 
                     isnan(b) ? -numeric_limits<double>::infinity() : b, 
                     isnan(c) ? -numeric_limits<double>::infinity() : c});
}


double viewer_distance(Particle particle, Point viewer_position) {
    //This function will calculate the distance between a point and the viewer

    double x = particle.x - viewer_position.x;
    double y = particle.y - viewer_position.y;
    double z = particle.z - viewer_position.z;

    double distance = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

    return distance;
    
    
}

Point calculate_viewer_position(double azimuth, double elevation, double distance) {
    double azimuth_rad = azimuth * 3.1415 / 180.0;
    double elevation_rad = elevation * 3.1415 / 180.0;

    double x = distance * cos(elevation_rad) * cos(azimuth_rad);
    double y = distance * cos(elevation_rad) * sin(azimuth_rad);
    double z = distance * sin(elevation_rad);

    return {x, y, z};
}

double calc_sys_momentum(const Particles& particles) {
    double total_momentum_x = 0.0;
    double total_momentum_y = 0.0;
    double total_momentum_z = 0.0;

    // Calculate the vector sum of momentum for each particle
    for (const auto& p : particles.particles) {
        total_momentum_x += p.vx;
        total_momentum_y += p.vy;
        total_momentum_z += p.vz;
    }

    // Calculate the magnitude of the total momentum vector
    double total_momentum = sqrt(total_momentum_x * total_momentum_x +
                                      total_momentum_y * total_momentum_y +
                                      total_momentum_z * total_momentum_z);

    return total_momentum;
}

double calc_sys_kinetic_energy(const Particles& particles) {
    double total_kinetic_energy = 0.0;

    // Calculate the kinetic energy for each particle and sum them
    for (const auto& p : particles.particles) {
        double velocity_squared = p.vx * p.vx + p.vy * p.vy + p.vz * p.vz;
        total_kinetic_energy += 0.5 * velocity_squared;  // Assuming mass = 1 for each particle
    }

    return total_kinetic_energy;
}