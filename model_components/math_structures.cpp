//This file contains mathematical structures.


struct Point {
    // This struct will store the location parameters of a point
    double x;
    double y;
    double z;

    // Method to add another point to this point
    Point add(const Point& other) {
        return Point{x + other.x, y + other.y, z + other.z};
    }

    // Method to subtract another point from this point
    Point subtract(const Point& other) {
        return Point{x - other.x, y - other.y, z - other.z};
    }
};


struct PPoint {
    //This struct will store a 3D point using polar coordinates
    double r;
    double theta; //angle in the xy plane
    double phi; //angle in the xz plane
    
};

struct Rotator {
    //This struct will store the rotation parameters of a point
    double theta; //to rotate in the xy plane
    double phi;  //to rotate in the xz plane

};



struct Velocity {
    //This struct will store the velocity parameters of a point
    double vx;
    double vy;
    double vz;

    double calc_magnitude() {
        //calculate the magnitude of the velocity vector
        return sqrt(pow(vx, 2) + pow(vy, 2) + pow(vz, 2));
    }

    bool check_zero() {
        //check if all velocities are equal to zero
        if (vx == 0 && vy == 0 && vz == 0) {
            return true;
        }
        return false;
    }
};

struct Points {
    //This struct will store a vector of points
    vector<Point> points;


    void push_back(const Point point) {
        points.push_back(point);
    }

    //to calculate the midpoint between the n points stored in the vector:

    Point calculate_midpoint() const {

        return {(points[0].x + points[1].x) / 2, (points[0].y + points[1].y) / 2, (points[0].z + points[1].z) / 2};
    }

    //to calculate the distance between two points:

    double calculate_distance() {
        //check size of vector, if not 2, return error
        if (points.size() != 2) {
            cout << "Error: vector does not contain 2 points" << endl;
            return -1;
        }

        double x_diff = points[0].x - points[1].x;
        double y_diff = points[0].y - points[1].y;
        double z_diff = points[0].z - points[1].z;

        return sqrt(pow(x_diff, 2) + pow(y_diff, 2) + pow(z_diff, 2));
    }

    Velocity calculate_distance_components() {
        //check size of vector, if not 2, return error
        

        double x_diff = points[0].x - points[1].x;
        double y_diff = points[0].y - points[1].y;
        double z_diff = points[0].z - points[1].z;

        Velocity distance_components;

        distance_components.vx = x_diff;
        distance_components.vy = y_diff;
        distance_components.vz = z_diff;

        return distance_components;
    }

    Point find_closest(Point point) {
        //calculate the closest point in the vector to a given point
        
        Point closest_point;

        //1 calculate the distances between the given point and all points in the vector

        vector<double> distances;

        for (int i = 0; i < points.size(); i++) {
            double x_diff = point.x - points[i].x;
            double y_diff = point.y - points[i].y;
            double z_diff = point.z - points[i].z;

            distances.push_back(sqrt(pow(x_diff, 2) + pow(y_diff, 2) + pow(z_diff, 2)));
        }

        //2 find the index of the minimum distance

        int min_index = 0;
        double min_distance = distances[0];

        for (int i = 1; i < distances.size(); i++) {
            if (distances[i] < min_distance) {
                min_distance = distances[i];
                min_index = i;
            }
        }

        //3 return the closest point

        closest_point = points[min_index];

        return closest_point;
    }

    void rotate_points(double theta, double phi) {
        for (Point& point : points) {
            // Store original coordinates
            double x = point.x;
            double y = point.y;
            double z = point.z;

            // Rotate in the azimuthal plane (around z-axis)
            double newX = x * cos(theta) - y * sin(theta);
            double newY = x * sin(theta) + y * cos(theta);

            // Update the point's x and y with azimuthal rotation
            point.x = newX;
            point.y = newY;

            // Store the updated coordinates for the second rotation
            x = newX;
            y = newY;

            // Rotate in the altitude plane (around y-axis)
            newX = x * cos(phi) + z * sin(phi);
            double newZ = -x * sin(phi) + z * cos(phi);

            // Update the point's x and z with altitude rotation
            point.x = newX;
            point.z = newZ;
        }
    }

    void undo_rotate_points(double theta_z, double theta_y) {
        for (Point& point : points) {
            // Rotate back around y-axis
            double x_new = point.x * cos(theta_y) + point.z * sin(theta_y);
            double z_new = -point.x * sin(theta_y) + point.z * cos(theta_y);
            point.x = x_new;
            point.z = z_new;

            // Rotate back around z-axis
            double y_new = point.y * cos(theta_z) - point.x * sin(theta_z);
            x_new = point.y * sin(theta_z) + point.x * cos(theta_z);
            point.x = x_new;
            point.y = y_new;
        }
    }




    void subtract_intercept(Point intercept) {
        //subtract the intercept from all points in the vector

        for (int i = 0; i < points.size(); i++) {
            points[i].x -= intercept.x;
            points[i].y -= intercept.y;
            points[i].z -= intercept.z;
        }
    }

    void add_intercept(const Point& intercept) {
        for (auto& point : points) {
            point.x += intercept.x;
            point.y += intercept.y;
            point.z += intercept.z;
        }
    }
};

struct PPoints {
    //This struct will store a vector of polar points
    vector<PPoint> ppoints;


    void push_back(const PPoint point) {
        ppoints.push_back(point);
    }
};


struct _Points {
    //This struct will store a vector of polar points, rotation_angle and rotation_intercept
    vector<Point> points;
    Rotator rotation_angles;
    Point rotation_intercept;

};




struct Scaler {
    //This struct will store the rotation parameters of a point and the scaling parameters of a point
    double theta;
    double phi;

    double r;

};



struct PPoints_plus_scaler {
    //This struct will store both a PPoints struct and a rotation struct

    PPoints ppoints;
    Scaler polar_scaler;

};

struct Points_plus_scaler {
    //This struct will store both a Points struct and a rotation struct

    Points points;
    Scaler scaler;

};



struct Velocities {
    //This struct will store a vector of velocities
    vector<Velocity> velocities;


    void push_back(const Velocity velocity) {
        velocities.push_back(velocity);
    }

    

    bool check_zero() {
        //check if all velocities are equal to zero
        for (int i = 0; i < velocities.size(); i++) {
            if (velocities[i].vx != 0 || velocities[i].vy != 0 || velocities[i].vz != 0) {
                return false;
            }
        }
        return true;
    }

    void solve_zero() {
        //if all velocities are zero, assign a velocity of 0.1 to 1 of the Velocity components

        if (check_zero()) {
           
            velocities[0].vx = 0.1;
            
        }
    };

};

struct Gradient {
    //This struct will store the gradient parameters of a point. 3d space has only 3

    double gradient_xy;
    double gradient_xz;
    
    
};

struct Intercept {
    //This struct will store the intercept parameters of a point
    double intercept_x;
    double intercept_y;
    double intercept_z;

    void calc_intercept(Point point, Gradient gradient) {
        //calculate the intercept of a point given a gradient
        intercept_x = point.x - gradient.gradient_xy * point.y;
        intercept_y = point.y - gradient.gradient_xy * point.x;
        intercept_z = point.z - gradient.gradient_xz * point.x;
    }

};

struct Line {
    Point point;
    double dx;
    double dy;
    double dz;

    // Define the line in 3D space given a point and a velocity
    void define_line(Point point1, Velocity velocity) {
        point = point1;
        dx = velocity.vx;
        dy = velocity.vy;
        dz = velocity.vz;
    }

    // Define the line in 3D space given two points
    void define_line2(Point point1, Point point2) {
        point = point1;
        dx = point2.x - point1.x;
        dy = point2.y - point1.y;
        dz = point2.z - point1.z;
    }

    // Define the perpendicular line in 3D space given two points
    void define_perp(Point point1, Point point2) {
        point = point1;
        dx = point2.x - point1.x;
        dy = point2.y - point1.y;
        dz = point2.z - point1.z;

        // Calculate the perpendicular direction
        double temp = dx;
        dx = -dy;
        dy = temp;
        dz = -dz;
    }
};

struct Sphere {
    //This struct will store the parameters of a  sphere in 3D space

    Point centre;
    double radius;
};

struct Magnitude {
    //This struct will store the 3 coordinate magnitudes in a vector
    double x_mag;
    double y_mag;
    double z_mag;

};


struct Magnitudes {
    //This struct will store the magnitude of a vector
    vector<Magnitude> magnitudes;

};

struct Particles_plus_skipped {
    //This struct will store both a Particles struct and a double for the timestep scaling factor

    Particles particles;
    double timestep_adjustment;

};



Point normalize(Point vector) {
    // Helper function to normalize a 3D vector
    double length = sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
    return {vector.x / length, vector.y / length, vector.z / length};
}


double dot_product(Point vector1, Point vector2) {
    // Helper function to calculate the dot product of two 3D vectors
    return vector1.x * vector2.x + vector1.y * vector2.y + vector1.z * vector2.z;
}


Point scalar_multiply(double scalar, Point vector) {
    // Helper function to multiply a vector by a scalar
    return {scalar * vector.x, scalar * vector.y, scalar * vector.z};
}


Point vector_subtract(Point vector1, Point vector2) {
    // Helper function to subtract one vector from another
    return {vector1.x - vector2.x, vector1.y - vector2.y, vector1.z - vector2.z};
}

Point vector_add(Point vector1, Point vector2) {
    
    // Helper function to add two vectors together
    return {vector1.x + vector2.x, vector1.y + vector2.y, vector1.z + vector2.z};
}

