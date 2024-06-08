//This file contains all the functions that are used for the PlotUniverse Tool


//libraries
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <math.h>
#include <time.h>
#include <windows.h>

#include <random> //include a normal distribution random number generator
#include <fstream>
#include <sstream>
#include <map>
#include <variant>

//other files in project
#include "math_model.cpp"
//#include "math_structures.cpp"






//namespaces
//using namespace std;

Particles sort_to_viewer(Particles particles, Point viewer) {
    // Calculate distances and store them with their indices
    vector<pair<float, int>> distance_indices;
    for (int i = 0; i < particles.particles.size(); i++) {
        float distance = viewer_distance(particles.particles[i], viewer);
        distance_indices.push_back({distance, i});
    }

    // Sort the indices based on the distances
    sort(distance_indices.begin(), distance_indices.end());

    // Create sorted particles based on sorted indices
    Particles sorted_particles;
    for (auto &di : distance_indices) {
        sorted_particles.particles.push_back(particles.particles[di.second]);
    }

    return sorted_particles;
}

void plot_particles(Particles particles, FILE* gnuplotPipe, double radius = 5) {
    //This function will plot the particles using GNUplot

    //1. Loop through the particles and plot each particle

    
    //to use radius instead of point size, use the following command

    //fprintf(gnuplotPipe, "set view equal xyz\n");
    
    fprintf(gnuplotPipe, "unset key\n");
    

    fprintf(gnuplotPipe, "splot '-' with points pt 7 ps %f lc rgb variable\n", radius);

    //2. Send the particles to the gnuplot pipe

    for (int i = 0; i < particles.particles.size(); i++) {
        
        fprintf(gnuplotPipe, "%f %f %f %d\n", particles.particles[i].x, particles.particles[i].y, particles.particles[i].z, particles.particles[i].rgb);

    }
    fprintf(gnuplotPipe, "e\n");



    //3. refresh the plot and ensure the window is the correct size compared to the particle radius


    fprintf(gnuplotPipe, "refresh\n");

    

}

Particles kinetic_interact_old(Particle particle_i, Particle particle_j){

    //This function will change the velocities of the 2 particles based on the kinetic interaction

    //1. Loop through the particles and change the velocities of the particles for each dimension

    float m1 = 1;
    float m2 = 1;

    float u1 = particle_i.vx;
    float u2 = particle_j.vx;

    float v1 = (m1 - m2) / (m1 + m2) * u1 + 2 * m2 / (m1 + m2) * u2;

    float v2 = 2 * m1 / (m1 + m2) * u1 + (m2 - m1) / (m1 + m2) * u2;

    particle_i.vx = v1;
    particle_j.vx = v2;

    u1 = particle_i.vy;
    u2 = particle_j.vy;

    v1 = (m1 - m2) / (m1 + m2) * u1 + 2 * m2 / (m1 + m2) * u2;

    v2 = 2 * m1 / (m1 + m2) * u1 + (m2 - m1) / (m1 + m2) * u2;

    particle_i.vy = v1;
    particle_j.vy = v2;

    u1 = particle_i.vz;
    u2 = particle_j.vz;

    v1 = (m1 - m2) / (m1 + m2) * u1 + 2 * m2 / (m1 + m2) * u2;

    v2 = 2 * m1 / (m1 + m2) * u1 + (m2 - m1) / (m1 + m2) * u2;

    particle_i.vz = v1;
    particle_j.vz = v2;

    Particles particles;

    particles.push_back(particle_i);
    particles.push_back(particle_j);

    return particles;

}


Points resolve_obliques(Point i, Point j, Point i_v, Point j_v, double restitution_param) {
    /**
 * Resolves the collision between two particles and calculates the new velocities.
 *
 * Steps:
 * 1. Calculate the normal vector (line of centers) between the two particles.
 * 2. Normalize the normal vector.
 * 3. Decompose the velocities of the particles into normal and tangent components.
 * 4. Apply conservation of momentum and Newton's law of restitution to calculate the final normal components.
 * 5. Combine the normal and tangent components to get the final velocities of the particles.
 * https://www.savemyexams.com/a-level/further-maths_further-mechanics-1/edexcel/17/revision-notes/elastic-collisions-in-2d/elastic-collisions-in-2d/oblique-collisions-of-two-spheres/

 */
    //calculate the magntiude of both vectors

    Points collided_points;

    // Print before collision
    //cout << "Before collision:" << endl;
    //cout << "Location i: (" << i.x << ", " << i.y << ", " << i.z << "), Velocity i: (" << i_v.x << ", " << i_v.y << ", " << i_v.z << ")" << endl;
    //cout << "Location j: (" << j.x << ", " << j.y << ", " << j.z << "), Velocity j: (" << j_v.x << ", " << j_v.y << ", " << j_v.z << ")" << endl;

    // Step 1: Calculate the normal vector (line of centres)
    Point normal_vector = {j.x - i.x, j.y - i.y, j.z - i.z};

    // Step 2: Normalize the normal vector
    Point normal_unit_vector = normalize(normal_vector);

    // Step 3: Decompose velocities into normal and tangent components
    double v1_normal_initial = dot_product(i_v, normal_unit_vector);
    double v2_normal_initial = dot_product(j_v, normal_unit_vector);

    Point v1_tangent_initial = vector_subtract(i_v, scalar_multiply(v1_normal_initial, normal_unit_vector));
    Point v2_tangent_initial = vector_subtract(j_v, scalar_multiply(v2_normal_initial, normal_unit_vector));

    // Step 4.1: Apply conservation of momentum in the normal direction
    double total_mass = 2;  // Since both spheres have equal mass
    double v1_normal_final = (v1_normal_initial * (1 - restitution_param) + v2_normal_initial * (1 + restitution_param)) / total_mass;
    double v2_normal_final = (v2_normal_initial * (1 - restitution_param) + v1_normal_initial * (1 + restitution_param)) / total_mass;

    // Step 4.2: Calculate new normal components based on Newton's law of restitution
    v1_normal_final = (v1_normal_initial * (1 - restitution_param) + v2_normal_initial * (1 + restitution_param)) / 2;
    v2_normal_final = (v2_normal_initial * (1 - restitution_param) + v1_normal_initial * (1 + restitution_param)) / 2;

    // Step 5: Combine components to get the final velocities
    Point v1_normal_vector_final = scalar_multiply(v1_normal_final, normal_unit_vector);
    Point v2_normal_vector_final = scalar_multiply(v2_normal_final, normal_unit_vector);

    Point velocity_i_final = vector_add(v1_normal_vector_final, v1_tangent_initial);
    Point velocity_j_final = vector_add(v2_normal_vector_final, v2_tangent_initial);

    collided_points.points = {i, j, velocity_i_final, velocity_j_final};

    return collided_points;
}


Points resolve_collision(Point i, Point j, Point i_v, Point j_v, double restitution_param) {
    //This function will resolve the collision between 2 particles in the new coordinate system

    //1. inititalise the resolved points

    Points resolved_points;

    //2. resolve the kinetic interaction

    resolved_points = resolve_obliques(i, j, i_v, j_v, restitution_param);

    //further actions can be added here

    //3.return the resolved points

    return resolved_points;

}







Particles kinetic_interact_new(Particle particle_i, Particle particle_j, double restitution_param = 1) {
    //This function will change the velocities of the 2 particles based on the kinetic interaction

    //1. Store the particle location and velocities in a structure
    Points cartesian_points;


    Point i_cartesian = { particle_i.x, particle_i.y, particle_i.z };
    Point j_cartesian = { particle_j.x, particle_j.y, particle_j.z };
    Point i_v_cartesian = { particle_i.vx, particle_i.vy, particle_i.vz };
    Point j_v_cartesian = { particle_j.vx, particle_j.vy, particle_j.vz };

    cartesian_points.points = { i_cartesian, j_cartesian, i_v_cartesian, j_v_cartesian };


    //2. rotate the coordinate system and set the origin to the midpoint between i and j

    _Points _rotated_cartesian_points;
    
    _rotated_cartesian_points = rotate_cartesians(cartesian_points);

    Rotator angles = _rotated_cartesian_points.rotation_angles;
    //float azimuth_rotation = angles.theta;
    //float altitude_rotation = angles.phi;

    Point intercept = _rotated_cartesian_points.rotation_intercept;

    vector<Point> rotated_cartesian_points_raw = _rotated_cartesian_points.points; // Convert Points to std::vector<Point>

    Points rotated_cartesian_points;

    rotated_cartesian_points.points = { rotated_cartesian_points_raw[0], rotated_cartesian_points_raw[1], rotated_cartesian_points_raw[2], rotated_cartesian_points_raw[3] };


    //3.resolve the collission

    Point i = { rotated_cartesian_points.points[0].x, rotated_cartesian_points.points[0].y, rotated_cartesian_points.points[0].z };
    Point j = { rotated_cartesian_points.points[1].x, rotated_cartesian_points.points[1].y, rotated_cartesian_points.points[1].z };

    Point i_v = { rotated_cartesian_points.points[2].x, rotated_cartesian_points.points[2].y, rotated_cartesian_points.points[2].z };
    Point j_v = { rotated_cartesian_points.points[3].x, rotated_cartesian_points.points[3].y, rotated_cartesian_points.points[3].z };


    Points new_cartesian_points = resolve_collision(i, j, i_v, j_v, restitution_param);

    //4. revert the coordinate system rotation

    Points unrotated_cartesian_points;

    unrotated_cartesian_points = unrotate_cartesians(new_cartesian_points.points, angles, intercept);

    //check the difference between the unrotated_cartesian_points and the original cartesian_points
    //compare_points(cartesian_points, unrotated_cartesian_points);


    //5. Update the particle locations and velocities

    particle_i.x = unrotated_cartesian_points.points[0].x;
    particle_i.y = unrotated_cartesian_points.points[0].y;
    particle_i.z = unrotated_cartesian_points.points[0].z;

    particle_i.vx = unrotated_cartesian_points.points[2].x;
    particle_i.vy = unrotated_cartesian_points.points[2].y;
    particle_i.vz = unrotated_cartesian_points.points[2].z;

    particle_j.x = unrotated_cartesian_points.points[1].x;
    particle_j.y = unrotated_cartesian_points.points[1].y;
    particle_j.z = unrotated_cartesian_points.points[1].z;

    particle_j.vx = unrotated_cartesian_points.points[3].x;
    particle_j.vy = unrotated_cartesian_points.points[3].y;
    particle_j.vz = unrotated_cartesian_points.points[3].z;


    Particles particles_post;

    particles_post.push_back(particle_i);
    particles_post.push_back(particle_j);

    //6. return the updated particles

    return particles_post;

}

Particles backtrack_collision(Particle particle_i, Particle particle_j, double radius) {
    //This function will calculate the exact point of collision between 2 particles

    //store the original particles

    Particles original_particles;

    original_particles.push_back(particle_i);
    original_particles.push_back(particle_j);

    //1.  Store the particles in a structure

    Point i = { particle_i.x, particle_i.y, particle_i.z };
    Point j = { particle_j.x, particle_j.y, particle_j.z };
    Velocity i_v = { particle_i.vx, particle_i.vy, particle_i.vz };
    Velocity j_v = { particle_j.vx, particle_j.vy, particle_j.vz };


    //store i and j and i_v and j_v in the Points and Velocities structures

    Points ij_points;


    ij_points.points = { i, j };

    //make a copy of the points

    Points ij_points_old = ij_points;

    Velocities ij_velocities;
    
    ij_velocities.velocities = { i_v, j_v };




    //2. Calculate the distance between the particles

    float distance_pre = ij_points.calculate_distance();

    //cout << "Distance between particle i and j before collission: " << distance_pre << endl;


    //1. calculate the midpoint between the 2 particles

    Point midpoint_ij = ij_points.calculate_midpoint();


    //2 Find the slope parameters of velocity vectors of particle i, from its origin //PER DEFINITION THE VELOCITY VECTOR IS THE SLOPE OF THE LINE


    //check if all velocity of particle i is zero, if so add a small value to the velocity

    if (ij_velocities.check_zero()) {
        ij_velocities.solve_zero();
    }

    //find the velocity vector magnitudes

    double magnitude_i = ij_velocities.velocities[0].calc_magnitude();
    double magnitude_j = ij_velocities.velocities[1].calc_magnitude();


    //find the relative magnitudes

    double relative_magnitude_i = 2*(magnitude_i / (magnitude_i + magnitude_j)); //times 2 because moving from radius to two radii
    double relative_magnitude_j = 2*(magnitude_j / (magnitude_i + magnitude_j));

    Magnitudes relative_contributions;

    relative_contributions = find_comp_contributions(ij_velocities.velocities[0], ij_velocities.velocities[1]);

    //separate ito relative_contributions_i and relative_contributions_j

    Magnitude relative_contributions_i = relative_contributions.magnitudes[0];
    Magnitude relative_contributions_j = relative_contributions.magnitudes[1];

    //to multiply each element in relative_contributions_i *2

    relative_contributions_i.x_mag *= 2;
    relative_contributions_i.y_mag *= 2;
    relative_contributions_i.z_mag *= 2;

    relative_contributions_j.x_mag *= 2;
    relative_contributions_j.y_mag *= 2;
    relative_contributions_j.z_mag *= 2;

     


    //check if relative_magnitude_i is not equal to any of {0,1,2}
    
    //vector <float> check_values = { 0, 1, 2 };

    //check if none of the check_values are equal to relative_magnitude_i



    //alternatively
    //bool found = false;
    //for (float f : check_values) {
    //    if (f == relative_magnitude_i) {
    //        found = true;
    //        break;
    //    }
    //}

    //if (!found) {
    //    cout << "Error: relative_magnitude_i is not equal to any of {0,1,2}" << endl;
    //}

    //3. find the intersection point between:

    //a. the line from the midpoint in the direction of the velocity vector of particle i
    //b. the SPHERE defined as the origin of particle j with radius = radius

    Velocity non_zero_v;

    //check if all velocities of particle i are zero, if so define non_zero_v as velocity of particle j

    if (ij_velocities.velocities[0].check_zero()) {
        non_zero_v = ij_velocities.velocities[1];
    }
    else {
        non_zero_v = ij_velocities.velocities[0];
    }

    Line midpoint_i_line; 
    midpoint_i_line.define_line(midpoint_ij, non_zero_v);

    Sphere sphere_j;
    sphere_j.centre = j;
    sphere_j.radius = radius;

    Points intersections_iside = ls_intersections(midpoint_i_line, sphere_j);

    Line midpoint_j_line;
    
    midpoint_j_line.define_line(midpoint_ij, ij_velocities.velocities[1]);

    Sphere sphere_i;
    sphere_i.centre = i;
    sphere_i.radius = radius;

    Points intersections_jside = ls_intersections(midpoint_j_line, sphere_i);

    //4. Find the intersection point that is closest to the midpoint

    Point intersection_i = intersections_iside.find_closest(midpoint_ij);
    Point intersection_j = intersections_jside.find_closest(midpoint_ij);

    //5. Calculate the distance between the midpoint and the intersection points

    Points midpoint_i;
    midpoint_i.points = { midpoint_ij, intersection_i};
    Points midpoint_j;
    midpoint_j.points = { midpoint_ij, intersection_j};


    Velocity distance_midpoint_i = midpoint_i.calculate_distance_components();
    Velocity distance_midpoint_j = midpoint_j.calculate_distance_components();
    

   //6a. change location of point i by the distance_midpoint_i and point j by the distance_midpoint_j

    //if intersection_i is not nan, then change the location of point i by the distance_midpoint_i

    if (!isnan(intersection_i.x)) {
        i = { i.x - distance_midpoint_i.vx*relative_contributions_i.x_mag, i.y - distance_midpoint_i.vy*relative_contributions_i.y_mag, i.z - distance_midpoint_i.vz*relative_contributions_i.z_mag };
    }

    if (!isnan(intersection_j.x)) {
        j = { j.x - distance_midpoint_j.vx*relative_contributions_j.x_mag, j.y - distance_midpoint_j.vy*relative_contributions_j.y_mag, j.z - distance_midpoint_j.vz*relative_contributions_j.z_mag };
    }



    //6c. multiply the velocity of point i by the relative_contributions of the velocities of point i and j



    //float transfered = distance_midpoint_j.vy*relative_contributions_j.y_mag;

    //cout << "transfered: " << transfered << endl;

    //i_v = { i_v.vx + distance_midpoint_i.vx*relative_contributions_i.x_mag, i_v.vy + distance_midpoint_i.vy*relative_contributions_i.y_mag, i_v.vz + distance_midpoint_i.vz*relative_contributions_i.z_mag }; Should this be incorporated now or as a correction to the velocity after the collision?
    //j_v = { j_v.vx + distance_midpoint_j.vx*relative_contributions_j.x_mag, j_v.vy + distance_midpoint_j.vy*relative_contributions_j.y_mag, j_v.vz + distance_midpoint_j.vz*relative_contributions_j.z_mag };

    Velocity i_v_skipped = {0,0,0};
    Velocity j_v_skipped = {0,0,0};
    
    if (!isnan(intersection_i.x)) {
        i_v_skipped = {distance_midpoint_i.vx*relative_contributions_i.x_mag, distance_midpoint_i.vy*relative_contributions_i.y_mag, distance_midpoint_i.vz*relative_contributions_i.z_mag };
    }
    
    if (!isnan(intersection_j.x)) {
        j_v_skipped = {distance_midpoint_j.vx*relative_contributions_j.x_mag, distance_midpoint_j.vy*relative_contributions_j.y_mag, distance_midpoint_j.vz*relative_contributions_j.z_mag };
    }

    //calculate the time step adjustment for the next iteration

    double magnitude_i_v_skipped = i_v_skipped.calc_magnitude();
    double timestep_adjustment_i = (magnitude_i_v_skipped/2 + magnitude_i)/magnitude_i;

    double magnitude_j_v_skipped = j_v_skipped.calc_magnitude();
    double timestep_adjustment_j = (magnitude_j_v_skipped/2 + magnitude_j)/magnitude_j;

    //get max of timestep_adjustment_i and timestep_adjustment_j

    double min_timestep = 1;

    double timestep_adjustment = custom_max(timestep_adjustment_i, timestep_adjustment_j,min_timestep);
    

    //7. Validation: find the distance between particle i and j, it should be equal to 2*radius

    ij_points.points[0] = i;
    ij_points.points[1] = j;

    ij_velocities.velocities[0] = i_v;
    ij_velocities.velocities[1] = j_v;

    double distance_post = ij_points.calculate_distance();

    //check if distance_post is not equal to 2*radius with a tolerance of 0.001

    //print the distance_post

    //cout << "Distance between particle i and j after collission: " << distance_post << endl;

    //if (abs(distance_post - 2*radius) > 0.001) {
    //    cout << "Error: distance between particle i and j after collission is not equal to 2*radius" << endl;
    //}

    //8. Update the particles

    particle_i.x = ij_points.points[0].x;
    particle_i.y = ij_points.points[0].y;
    particle_i.z = ij_points.points[0].z;

    particle_i.time_scaler = timestep_adjustment;


    particle_i.vx = ij_velocities.velocities[0].vx;
    particle_i.vy = ij_velocities.velocities[0].vy;
    particle_i.vz = ij_velocities.velocities[0].vz;

    particle_j.x = ij_points.points[1].x;
    particle_j.y = ij_points.points[1].y;
    particle_j.z = ij_points.points[1].z;

    particle_j.time_scaler = timestep_adjustment;

    particle_j.vx = ij_velocities.velocities[1].vx;
    particle_j.vy = ij_velocities.velocities[1].vy;
    particle_j.vz = ij_velocities.velocities[1].vz;

    Particles pre_collission_particles;
 
    pre_collission_particles.push_back(particle_i);
    pre_collission_particles.push_back(particle_j);


    //check if momentum is conserved

    double sys_momentum_pre = calc_sys_momentum(original_particles);
    double sys_momentum_post = calc_sys_momentum(pre_collission_particles);

    if (sys_momentum_pre != sys_momentum_post) {
        cout << "System momentum not conserved after backtrack" << endl;
    }

    return pre_collission_particles;





};

Particles kinetic_backtrack_interact(Particle particle_i, Particle particle_j, double radius, double restitution_param) {
    //This function will change the velocities of the 2 particles based on the kinetic interaction

    //1. Define the exact point of collision

    Particles pre_backtrack_particles;

    pre_backtrack_particles.particles = { particle_i, particle_j };


    Particles pre_collission_particles = backtrack_collision(particle_i, particle_j, radius);


        
    particle_i = pre_collission_particles.particles[0];
    particle_j = pre_collission_particles.particles[1];

    //2. interact kinetically at the point of collision



    Particles particles_post = kinetic_interact_new(particle_i, particle_j, restitution_param);




    particle_i = particles_post.particles[0];
    particle_j = particles_post.particles[1];

    //3. Return the updated particles

    Particles updated_particles;

    updated_particles.push_back(particle_i);
    updated_particles.push_back(particle_j);


    return updated_particles;

    

}

Particles interact_pair(Particle particle_i, Particle particle_j, double radius, double restitution_param, string interaction = "none") {
    //This function will perform an interaction between particles


    Particles particles;

    //1. Check the interaction type

    if (interaction == "none") {

        particles.push_back(particle_i);
        particles.push_back(particle_j);


        return particles;
    }
    else if (interaction == "kinetic") {
        //2. Perform the interaction

        particles = kinetic_interact_old(particle_i, particle_j);





    } else if (interaction == "kinetic_backtrack")
    {
        particles = kinetic_backtrack_interact(particle_i, particle_j, radius, restitution_param);

    }
    

    else {
        cout << "Invalid interaction type" << endl;}
        
    
    return particles;
}

Particles detect_collission(Particles particles, double radius, string interaction, double restitution_param) {
    //This function will check for collisions between particles, perform the interaction and store the collided particles in a vector with Ids

    //Particles particles_ij;

    vector<int> collided_particles;

    Particles particles_ij;


    //1. Loop through the particles and check for collisions

    for (int i = 0; i < particles.particles.size(); i++) {
        for (int j = i+1; j < particles.particles.size(); j++) {
            
            Particle particle_i = particles.particles[i];
            Particle particle_j = particles.particles[j];

            double distance = sqrt(pow(particle_i.x - particle_j.x, 2) + pow(particle_i.y - particle_j.y, 2) + pow(particle_i.z - particle_j.z, 2));
            if (distance < (radius + radius)) {
                collided_particles.push_back(particle_i.ID);
                collided_particles.push_back(particle_j.ID);

                //calculate the pre interaction momentum
                Particles pre_interaction_ij;

                pre_interaction_ij.push_back(particle_i);
                pre_interaction_ij.push_back(particle_j);

                //double pre_interaction_momentum = calc_sys_momentum(pre_interaction_ij);

                

                //2. Perform the interaction between the particles

                particles_ij = interact_pair(particle_i, particle_j, radius, restitution_param, interaction);

                //check momentum after collission

                //double post_interaction_momentum = calc_sys_momentum(particles_ij);

                //3. Check if the system momentum is conserved

                //if (pre_interaction_momentum != post_interaction_momentum) {
                //    cout << "System momentum not conserved" << endl;
                //}

                //4. update particle i and j

                particles.particles[i] = particles_ij.particles[0];
                particles.particles[j] = particles_ij.particles[1];

            }
        }
    }


    return particles;
}

Particles update_particles(Particles particles, double radius, double restitution_param,  string interaction = "none") {
    //This function will update the particles

    //1. Check for collisions between particles, store the collided particles in a vector with Ids

    Particles updated_particles;
    Particles collided_particles;
    Particle particle;

    
    collided_particles = detect_collission(particles, radius, interaction, restitution_param);

    //double sys_momentum_pre = calc_sys_momentum(collided_particles);

    //2.  Add the velocities to the coordinates of the particles

    for (int i = 0; i < collided_particles.particles.size(); i++) {
        Particle particle = collided_particles.particles[i];


        particle.x += particle.vx*particle.time_scaler;
        particle.y += particle.vy*particle.time_scaler;
        particle.z += particle.vz*particle.time_scaler;

        updated_particles.push_back(particle);
    }


    //double sys_momentum_post = calc_sys_momentum(updated_particles);

    //3. Check if the system momentum is conserved

    //if (sys_momentum_pre != sys_momentum_post) {
    //    cout << "System momentum before and after step 2 not conserved" << endl;
    //}




    return updated_particles;
}

Particles handle_overlap(Particles particles, double radius = 5) {
    //This function will handle overlapping particles

    //1. Loop through the particles and check for overlapping particles

    int removed_particles = 0;
    for (int i = 0; i < particles.particles.size(); i++) {
        for (int j = i+1; j < particles.particles.size(); j++) {
            
            Particle particle_i = particles.particles[i];
            Particle particle_j = particles.particles[j];

            double distance = sqrt(pow(particle_i.x - particle_j.x, 2) + pow(particle_i.y - particle_j.y, 2) + pow(particle_i.z - particle_j.z, 2));
            if (distance < (radius + radius)) { //check if distance < sum of radii
                //2. Delete particle j

                //remove particle at index j

                particles.particles.erase(particles.particles.begin() + j);

                removed_particles++;

                
            }
        }

    
    }
    cout << "Removed " << removed_particles << " overlapping particles" << endl;

    return particles;
}




void simulate_particles(Particles particles, Scenario selected_scenario_params) {
    //This function will simulate the particles using GNUplot. It will initialise the plot, populate the particles, update the plot, refresh the plot, and close the plot.


    // Step 0: extract the parameters from the scenario

    string title;
    string interaction;

    int steps;
    int speed;
    double radius;
    double volume;
    double restitution_param = 1;


    title = selected_scenario_params.name;

    for (const auto& map : selected_scenario_params.params) {

        // //get the title of the scenario

        // if (map.first == "scenario") {
        //     title = get<string>(map.second[0]);
        // }
        
        //get the interaction type

        if (map.first == "interaction") {
            interaction = get<string>(map.second[0]);
        }

        //get the number of steps

        if (map.first == "steps") {
            if (holds_alternative<float>(map.second[0])){

                steps = get<float>(map.second[0]);
            
            }
        }

        //get the speed of the simulation

        if (map.first == "speed") {
            
            if (holds_alternative<float>(map.second[0])){

                speed = get<float>(map.second[0]);
            }
        }

        //get the radius of the particles

        if (map.first == "radius") {
            radius = get<float>(map.second[0]);
        }

        //get the restituion parameter

        if (map.first == "restitution_param") {
            restitution_param = get<float>(map.second[0]);
        }

    }

    int step = 0;

    // Step 1: Handle Overlapping Particles
    
    particles = handle_overlap(particles, radius);
    int n = particles.particles.size();

    // Step 2: Calculate the initial momentum of the system

    double sys_momentum = 0;
    double sys_kinetic_energy = 0;

    sys_momentum = calc_sys_momentum(particles);
    sys_kinetic_energy = calc_sys_kinetic_energy(particles);




    // Step 3: Open GNUplot Pipe

    cout << "Connecting to GNUplot..." << endl;


    FILE* gnuplotPipe = popen("gnuplot -persistent", "w"); 

    // Step 4: Set Plot Parameters



    //set plot parameters

    


    fprintf(gnuplotPipe, "set title '%s' font 'Arial Bold,16'\n", title.c_str());
    fprintf(gnuplotPipe, "set label 1 'Number of Particles: %d' at screen 0.02,0.95 font 'Arial Bold,12'\n", n);
    fprintf(gnuplotPipe, "set label 2 'Timestep: %d' at screen 0.02,0.90 font 'Arial Bold,12'\n", step);
    fprintf(gnuplotPipe, "set label 3 'Restitution Parameter: %.2f' at screen 0.02,0.85 font 'Arial Bold,12'\n", restitution_param);
    fprintf(gnuplotPipe, "set label 4 'System momentum: %.2f' at screen 0.02,0.80 font 'Arial Bold,12'\n", sys_momentum);
    fprintf(gnuplotPipe, "set label 5 'System kinetic energy: %.2f' at screen 0.02,0.75 font 'Arial Bold,12'\n", sys_kinetic_energy);
    fprintf(gnuplotPipe, "set label 6 'FPS: %d' at screen 0.02,0.70 font 'Arial Bold,12'\n", 0);
    
    


    //set axis limitsfprintf(gnuplotPipe, "set term wxt position x, y\n");

    fprintf(gnuplotPipe, "set xrange [-100:100]\n");
    fprintf(gnuplotPipe, "set yrange [-100:100]\n");
    fprintf(gnuplotPipe, "set zrange [-100:100]\n");

    


    //set gridlines for all 3 axes
    
    fprintf(gnuplotPipe, "set grid xtics ytics ztics\n");

    //set gridlines for xy plane

    fprintf(gnuplotPipe, "set xyplane at 0\n");

    //set axis titles

    fprintf(gnuplotPipe, "set xlabel 'X' font 'Arial Bold,12'\n");
    fprintf(gnuplotPipe, "set ylabel 'Y' font 'Arial Bold,12'\n");
    fprintf(gnuplotPipe, "set zlabel 'Z' font 'Arial Bold,12'\n");

    //define the point the viewer sees the plot from


    //double azimuth = -30; 
    //double elevation = 60;
    //double distance = 100;

    //Point viewer = calculate_viewer_position(azimuth, elevation, distance);

    //set view equal to the azimuth and elevation

    //fprintf(gnuplotPipe, "set view %f, %f, 1, 1\n", azimuth, elevation);

    fprintf(gnuplotPipe, "set view equal xyz\n"); //otherwise axes will be scaled differently





    
    

    //5. Plot Particles
    //particles = sort_to_viewer(particles, viewer); //ensure rendering as per distance to viewer

    plot_particles(particles, gnuplotPipe, radius);

    fflush(gnuplotPipe);

    //Ask user to start simulation
    cout << "Press any key to start simulation..." << endl;
    cin.ignore();
    cin.get();

    
    //Update Particles and refresh plot
    
    double t_start_time;
    double t_end_time;
    double fps;
    
    for (int i = 0; i < steps; i++) {
        step++; 

        //timestamp

        t_start_time = clock();
        


        fprintf(gnuplotPipe, "set label 2 'Timestep: %d' at screen 0.02,0.90 font 'Arial Bold,12'\n", step);

        particles = update_particles(particles, radius, restitution_param, interaction);

        //particles = sort_to_viewer(particles, viewer); //ensure rendering as per distance to viewer
        



        plot_particles(particles, gnuplotPipe, radius);

        sys_momentum = calc_sys_momentum(particles);

        fprintf(gnuplotPipe, "set label 4 'System momentum: %.2f' at screen 0.02,0.80 font 'Arial Bold,12'\n", sys_momentum);

        sys_kinetic_energy = calc_sys_kinetic_energy(particles);

        fprintf(gnuplotPipe, "set label 5 'System kinetic energy: %.2f' at screen 0.02,0.75 font 'Arial Bold,12'\n", sys_kinetic_energy);

        t_end_time = clock();

        fps = 1 / ((t_end_time - t_start_time) / CLOCKS_PER_SEC);

        fprintf(gnuplotPipe, "set label 6 'FPS: %.2f' at screen 0.02,0.70 font 'Arial Bold,12'\n", fps);








        fflush(gnuplotPipe);


    }

    //6. Terminate simulation

    cout << "Press any key to terminate the simulation..." << endl;
    cin.ignore();
    cin.get();

    cout << "Closing GNUplot..." << endl;

    pclose(gnuplotPipe);

}
