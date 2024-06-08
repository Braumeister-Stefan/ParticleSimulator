//This file contains the functions needed to initialise, save and retrieve complex objects for the simulation

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
#include <iostream>
#include <sstream>
#include <map>
#include <variant>
#include <filesystem>



//namespaces
using namespace std;

//structures
#include "structures.cpp"



Object sample_point(vector<vector<float>> limits) {
    //This function will sample a point within the limits of the chosen object

    Object object_point;

    //1. Define a uniform distribution that samples 1 x,y,z coordinate within the cube limits, store in cube_point

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<float> dis_x(limits[0][0], limits[0][1]);
    uniform_real_distribution<float> dis_y(limits[1][0], limits[1][1]);
    uniform_real_distribution<float> dis_z(limits[2][0], limits[2][1]);

    float sample_x = dis_x(gen);
    float sample_y = dis_y(gen);
    float sample_z = dis_z(gen);

    //2. Store the x,y,z coordinates in a map called loc_parameters

    vector<variant<float, string>> loc_parameters = {sample_x, sample_y, sample_z};

    object_point.params["loc_parameters"] = loc_parameters;

    return object_point;
}

vector<float> calculate_distances(Objects cube_points,Object cube_point,int particle_radius) {
    //This function will calculate the distance between the new point and all the points in cube_points

    vector<float> distances;

    //1. Retrieve the x, y, z coordinates of the new point

    float x = get<float>(cube_point.params["loc_parameters"][0]);
    float y = get<float>(cube_point.params["loc_parameters"][1]);
    float z = get<float>(cube_point.params["loc_parameters"][2]);

    //2. Calculate the distance between the new point and all the points in cube_points

    for (int i = 0; i < cube_points.objects.size(); i++) {
        float distance = sqrt(pow(get<float>(cube_points.objects[i].params["loc_parameters"][0]) - x, 2) + pow(get<float>(cube_points.objects[i].params["loc_parameters"][1]) - y, 2) + pow(get<float>(cube_points.objects[i].params["loc_parameters"][2]) - z, 2));
        distances.push_back(distance);

        }

    return distances;
}

Objects sample_cube(int n, int diagonal, vector<variant<float, string>> loc_parameters, int particle_radius) {
    //This function will define a cube with n points and a diagonal of length diagonal

    Objects cube_points;
    Object cube_point;

    //1. Starting from the origin, and the diagonal, define the x,y,z limits of the cube.

    //define length of the cube 

    float length = diagonal / sqrt(3); //follows pythagoras theorem in 3D


    //2a. Retrieve the x, y, z coordinates of the object as the first 3 values of the vector loc_parameters

    float x = get<float>(loc_parameters[0]);
    float y = get<float>(loc_parameters[1]);
    float z = get<float>(loc_parameters[2]);

    //2b define limits of the cube and store in a vector of dimension 3*2

    vector<vector<float>> cube_limits = {{x - length/2, x + length/2}, {y - length/2, y + length/2}, {z - length/2, z + length/2}};


    for (int i = 0; i < n*5; i++) {

        if (cube_points.objects.size() == n) { //once we have the chosen amount of points, break the loop. 
            break;
        }

        //3 sample point within the cube limits
        
        cube_point = sample_point(cube_limits);

        //3. Check the distance between the new point and all the points in cube_points. If the distance is greater than the particle radius, store the point in cube_points, else, repeat step 2

        vector<float> distances;


        

        if (cube_points.objects.size() == 0) {
            cube_points.push_back(cube_point);
        } else {
            //calculate the distance between cube_point and all the points in cube_points
            distances = calculate_distances(cube_points, cube_point,particle_radius);


            //check if any of the distances is less than the particle radius. If not, store the point in cube_points

            if (all_of(distances.begin(), distances.end(), [particle_radius](float i) {return i > (particle_radius + particle_radius);})) {
                cube_points.push_back(cube_point);
                
                //if cube_points.objects.size() is a multiple of 10, print the number of points stored

                if (cube_points.objects.size() % 100 == 0) {
                    cout << "Points stored: " << cube_points.objects.size() << endl;
                }

            } else {
                //cout << "Point not stored" << endl;
                
            }
            
        }

        distances.clear();
        

    }


    //5. Return cube_points

    return cube_points;
    
}

bool sphere_func(float x, float y, float z, float radius, vector<variant<float, string>> loc_parameters) {
    //This function will check if a point is within a sphere

    float center_x = get<float>(loc_parameters[0]);
    float center_y = get<float>(loc_parameters[1]);
    float center_z = get<float>(loc_parameters[2]);

    float LH_side = pow(x - center_x, 2) + pow(y - center_y, 2) + pow(z - center_z, 2);
    float RH_side = pow(radius, 2);

    bool in_sphere = LH_side  <= RH_side;


    return in_sphere;
}



Objects sample_sphere(int n, int diagonal, vector<variant<float, string>> loc_parameters, int particle_radius) {
    //This function will define a sphere with n points and a diagonal of length diagonal

    Objects sphere_points;
    Object sphere_point;

    //1. Starting from the origin, and the diagonal, define the x,y,z limits of the sphere

    float radius = diagonal / 2;


    //2a. Retrieve the x, y, z coordinates of the object as the first 3 values of the vector loc_parameters

    float x = get<float>(loc_parameters[0]);
    float y = get<float>(loc_parameters[1]);
    float z = get<float>(loc_parameters[2]);

    //2b define the limits of the sphere as the cube that contains the sphere

    float cube_length = diagonal; // / sqrt(3); //follows pythagoras theorem in 3D.this would make the cube the largest cube inside the sphere

    vector<vector<float>> cube_limits = {{x - cube_length/2, x + cube_length/2}, {y - cube_length/2, y + cube_length/2}, {z - cube_length/2, z + cube_length/2}};


    for (int i = 0; i < n*5; i++) {

        if (sphere_points.objects.size() == n) { //once we have the chosen amount of points, break the loop. 
            break;
        }

        //3a sample point within the cube limits
        
        sphere_point = sample_point(cube_limits);

        //3b check if the point is within the sphere. If not, go to next iteration

        bool in_sphere = sphere_func(get<float>(sphere_point.params["loc_parameters"][0]), get<float>(sphere_point.params["loc_parameters"][1]), get<float>(sphere_point.params["loc_parameters"][2]), radius, loc_parameters);

        if (!in_sphere) {

            //cout << "Point not inside sphere" << endl;
            continue;
        }


        //4. Check the distance between the new point and all the points in sphere_points. If the distance is greater than the particle radius, store the point in sphere_points, else, repeat step 2

        vector<float> distances;


    
        if (sphere_points.objects.size() == 0) {
            sphere_points.push_back(sphere_point);
        } else {
            //calculate the distance between cube_point and all the points in cube_points
            distances = calculate_distances(sphere_points, sphere_point,particle_radius);


            //check if any of the distances is less than the sum of particle radii.

            if (all_of(distances.begin(), distances.end(), [particle_radius](float i) {return i > (particle_radius +particle_radius);})) {
                sphere_points.push_back(sphere_point);
                
                //if sphere_points.objects.size() is a multiple of 10, print the number of points stored

                if (sphere_points.objects.size() % 100 == 0) {
                    cout << "Points stored: " << sphere_points.objects.size() << endl;
                }

            } else {
                //cout << "Point not stored" << endl;
                
            }
            
        }

        distances.clear();
        

    }


    //5. Return cube_points

    return sphere_points;

}

void store_flattened_obj(Objects objects, string file_name) {
    //This function will save the flattened object to a text file. 

    Object object_i;
    ofstream file;

    string name_i;
    string line_i;
    string pair_i;



    file.open(file_name);

    //1. Loop through each object in objects

    for (int i = 0; i < objects.objects.size(); i++) {

        object_i = objects.objects[i];

        //2. Select the name of the object and store in name_i parameter

        name_i = "object(" + object_i.name + ")";
        line_i = name_i + ";";


        //3. Loop through each key value pair in object_i.params, store as string pair_i with convention n(1);rgb_parameters(255,0,0)

        for (const auto& pair : object_i.params) {
            pair_i = pair.first + "(";

            for (int j = 0; j < pair.second.size(); j++) {
                pair_i += to_string(get<float>(pair.second[j]));

                if (j != pair.second.size() - 1) {
                    pair_i += ",";
                }
            }

            pair_i += ")";

            line_i += pair_i;

            //3.1. If not the last pair, add a semicolon to the line_i

            if (pair.first != object_i.params.rbegin()->first) {
                line_i += ";";
            }
    
        }

        //4. Write line_i to the text file
        file << line_i << endl;
        line_i.clear();
        pair_i.clear();
        name_i.clear();

    }


    file.close();

    cout << "Object saved to " << file_name << endl;

}

Objects retrieve_flattened_obj(string file_name) {
    //This function will load the flattened object from a text file and should reverse the process of store_flattened_obj

    Object object_i;
    Objects objects;
    ofstream file;

    string name_i;
    string line_i;
    string pair_i;

    ifstream file_read(file_name);

    //1. Loop through each line in the text file. a line has the format object(cube_alpha_0);linear_parameters(-0.990000,-0.500000,0.000000);loc_parameters(23.705654,16.346561,-2.187499);rgb_parameters(0.000000,30.000000,33.000000)

    while (getline(file_read, line_i)) {

        //2. Split the line by the delimiter ";"

        vector<string> line_split;
        stringstream ss(line_i);
        string token;

        while (getline(ss, token, ';')) {
            line_split.push_back(token);
        }

        //3. Loop through each element in line_split

        for (int i = 0; i < line_split.size(); i++) {

            //4. Split the element by the delimiter "(", remove from the second element ")"

            vector<string> element_split;
            stringstream ss_element(line_split[i]);
            string token_element;

            while (getline(ss_element, token_element, '(')) {
                element_split.push_back(token_element);

            }

            element_split[1].pop_back(); //remove the last character from element_split[1]





            //5. If element_split[0] is equal to "object", store the name of the object in name_i

            if (element_split[0] == "object") {
                object_i.name = element_split[1];
            } else {

                //6. Store the first element in element_split as the key of the object

                string key = element_split[0];

                //7. Split the element by the delimiter ","

                vector<string> element_split_2;
                stringstream ss_element_2(element_split[1]);
                string token_element_2;

                while (getline(ss_element_2, token_element_2, ',')) {
                    element_split_2.push_back(token_element_2);
                }

                //8. Store the remaining elements in element_split_2 as the values of the object

                vector<variant<float, string>> values;

                for (int j = 0; j < element_split_2.size(); j++) {
                    values.push_back(stof(element_split_2[j]));
                }

                //9. Store the key value pair in object_i.params

                object_i.params[key] = values;

                     
            } 
            

        }

        //11. Store the object in objects

        objects.push_back(object_i);

        //12. Clear object_i and relevant variables

        object_i.params.clear();
        object_i.name.clear();
        line_i.clear();
        pair_i.clear();
        name_i.clear();


    }

    file_read.close();

    cout << "Object retrieved from " << file_name << endl;

    return objects;
            
}


Objects unfold_object(Object object, Scenario selected_scenario_params) {
    //This function will unfold an object into n points, e.g. a cube, a sphere etc

    Objects flattened_object;
    int particle_radius; //this is not optimal. Will eventually be part of object params, instead of scenario_params. Subsequently it can be removed
    bool refresh_objects = false;

    for (const auto& map : selected_scenario_params.params) {
        if (map.first == "radius") {
            particle_radius = get<float>(map.second[0]);
            
        }

        if (map.first == "refresh_objects") {
            refresh_objects = get<string>(map.second[0]) == "true";


        }
    }
    
    //1. If refresh_objects is false, check if the object has already been defined. If so, return the object

    string file_name = "inputs\\defined_objects\\" + object.name + ".txt";

    //check if a file exists by check how many lines are in the file. If there are more than 0, the file exists
    bool file_exists = false;

    ifstream file_read(file_name);
    int line_count = 0;
    string line;
    while (getline(file_read, line)) {
        line_count++;
    }

    if (line_count > 0) {
        file_exists = true;
    }

    if (!refresh_objects && file_exists) {

        ifstream file(file_name);
        if (file.good()) {
            cout << "Object " << object.name << " already defined, importing..." << endl;
            return retrieve_flattened_obj(file_name);
        }

    } else {
        cout << "Object " << object.name << " not defined yet/user request to redefine. Defining..." << endl;

    }
    

    //2. Retrieve the value of the key "shape" and store in shape, n in int n. Also store the key equal to diagonal in "diagonal"

    string shape;
    int n;
    int diagonal;

    vector<variant<float, string>> loc_parameters;
    vector<variant<float, string>> linear_parameters;
    vector<variant<float, string>> rgb_parameters;

    //loop through keys in selected_scenario_params and store the value of the key radius in particle_radius


    for (const auto& pair : object.params) {

        if (pair.first == "shape") {
            shape = get<string>(pair.second[0]);
            }
        if (pair.first == "n") {
            n = get<float>(pair.second[0]);
        }
        if (pair.first == "diagonal") {
            diagonal = get<float>(pair.second[0]);
        }
        if (pair.first == "loc_parameters") {
            // Retrieve the x, y, z coordinates of the particle as the first 3 values of the vector
            loc_parameters = pair.second;
        }
        if (pair.first == "linear_parameters") {
            // Retrieve the vx, vy, vz velocities of the particle as the next 3 values of the vector
            linear_parameters = pair.second;
        }
        if (pair.first == "rgb_parameters") {
            // Retrieve the r, g, b values of the particle as the next 3 values of the vector
            rgb_parameters = pair.second;
        }
        }
       

    //3. Check the shape parameter. If equal to "cube", call the function define_cube, if equal to "sphere", call the function define_sphere

    if (shape == "cube") {

        flattened_object = sample_cube(n, diagonal, loc_parameters, particle_radius);
    } else if (shape == "sphere") {
        flattened_object = sample_sphere(n, diagonal, loc_parameters, particle_radius);
    } else {
        cout << "Shape parameter: " << shape << " not implemented" << endl;

        return flattened_object;

    }

    //4. add the name (object name + i), linear and rgb parameters to each point in the flattened object

    for (int i = 0; i < flattened_object.objects.size(); i++) {
        flattened_object.objects[i].name = (object.name + "_" + to_string(i)).c_str();
        flattened_object.objects[i].params["linear_parameters"] = linear_parameters;
        flattened_object.objects[i].params["rgb_parameters"] = rgb_parameters;
    }

    


    //6. Export the flattened object to a text file


    store_flattened_obj(flattened_object, file_name);

    return flattened_object;

}



