#include <cmath>
#include <vector>
#include <iostream>

struct Point {
    double x, y, z;
};

struct Points {
    std::vector<Point> points;
};

Point normalize(Point vector) {
    double length = std::sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
    return {vector.x / length, vector.y / length, vector.z / length};
}

double dot_product(Point vector1, Point vector2) {
    return vector1.x * vector2.x + vector1.y * vector2.y + vector1.z * vector2.z;
}

Point scalar_multiply(double scalar, Point vector) {
    return {scalar * vector.x, scalar * vector.y, scalar * vector.z};
}

Point vector_subtract(Point vector1, Point vector2) {
    return {vector1.x - vector2.x, vector1.y - vector2.y, vector1.z - vector2.z};
}

Point vector_add(Point vector1, Point vector2) {
    return {vector1.x + vector2.x, vector1.y + vector2.y, vector1.z + vector2.z};
}

Points resolve_obliques(Point i, Point j, Point i_v, Point j_v, double restitution_param) {
    Points collided_points;

    Point normal_vector = {j.x - i.x, j.y - i.y, j.z - i.z};
    Point normal_unit_vector = normalize(normal_vector);

    double v1_normal_initial = dot_product(i_v, normal_unit_vector);
    double v2_normal_initial = dot_product(j_v, normal_unit_vector);

    Point v1_tangent_initial = vector_subtract(i_v, scalar_multiply(v1_normal_initial, normal_unit_vector));
    Point v2_tangent_initial = vector_subtract(j_v, scalar_multiply(v2_normal_initial, normal_unit_vector));

    double v1_normal_final = (v1_normal_initial * (1 - restitution_param) + v2_normal_initial * (1 + restitution_param)) / 2;
    double v2_normal_final = (v2_normal_initial * (1 - restitution_param) + v1_normal_initial * (1 + restitution_param)) / 2;

    Point v1_normal_vector_final = scalar_multiply(v1_normal_final, normal_unit_vector);
    Point v2_normal_vector_final = scalar_multiply(v2_normal_final, normal_unit_vector);

    Point velocity_i_final = vector_add(v1_normal_vector_final, v1_tangent_initial);
    Point velocity_j_final = vector_add(v2_normal_vector_final, v2_tangent_initial);

    collided_points.points = {i, j, velocity_i_final, velocity_j_final};

    return collided_points;
}

int main() {
    Point i = {0, 0, 0};
    Point j = {1, 0, 0};
    Point i_v = {1, 0, 0};
    Point j_v = {-1, 0, 0};
    double restitution_param = 0.8;

    Points result = resolve_obliques(i, j, i_v, j_v, restitution_param);

    double initial_magnitude = std::sqrt(i_v.x * i_v.x + i_v.y * i_v.y + i_v.z * i_v.z) +
                               std::sqrt(j_v.x * j_v.x + j_v.y * j_v.y + j_v.z * j_v.z);

    Point velocity_i_final = result.points[2];
    Point velocity_j_final = result.points[3];
    double final_magnitude = std::sqrt(velocity_i_final.x * velocity_i_final.x + velocity_i_final.y * velocity_i_final.y + velocity_i_final.z * velocity_i_final.z) +
                             std::sqrt(velocity_j_final.x * velocity_j_final.x + velocity_j_final.y * velocity_j_final.y + velocity_j_final.z * velocity_j_final.z);

    Point initial_momentum = vector_add(i_v, j_v);
    Point final_momentum = vector_add(velocity_i_final, velocity_j_final);

    std::cout << "Initial total magnitude: " << initial_magnitude << std::endl;
    std::cout << "Final total magnitude: " << final_magnitude << std::endl;

    std::cout << "Initial momentum: (" << initial_momentum.x << ", " << initial_momentum.y << ", " << initial_momentum.z << ")" << std::endl;
    std::cout << "Final momentum: (" << final_momentum.x << ", " << final_momentum.y << ", " << final_momentum.z << ")" << std::endl;

    return 0;
}
