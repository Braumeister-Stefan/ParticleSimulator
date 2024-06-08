//this file contains all functions related to loading the simulation inputs from the input file

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


//namespaces
using namespace std;



//other files in project
#include "object_loader_model.cpp"




Scenarios load_scenarios(const string& filename) {
    //This function will load a scenario for each line in the file.

    ifstream file(filename);

    string line;
    
    vector<string> line_parts;
    string part;
    string name_part;
    string value_part;

    vector<variant<float, string>> parameter_values;
    Parameter single_parameter;
    Scenario scenario;
    Scenarios scenarios;


    //1. Loop through each line


    while (getline(file, line)) {
        stringstream ss(line);

        //check if the line is empty, if so, skip it
       if (line.empty()) {
            continue;
        }
        //2. Loop through each part of the line, each part is separated by a semi colon. do the following:
        

        while (getline(ss, part, ';')) {
            stringstream ss_part(part);

            //2a. Extract the name of the parameter

            getline(ss_part, name_part, '(');

            //2b. Extract the values of the parameter

            getline(ss_part, value_part, ')');

            //2c. Store the comma separated values in value_part in a vector parameter_values

            stringstream ss_values(value_part);
            float temp_float;
            string temp_string;

            while (getline(ss_values, value_part, ',')) {
                stringstream ss_value(value_part);
                if (ss_value >> temp_float) {
                    parameter_values.push_back(temp_float);
                } else {
                    parameter_values.push_back(value_part);
                }
            }

            //2c. Store the name of the parameter and its values in a map called single_parameter

            single_parameter.name = name_part;
            single_parameter.values = parameter_values;
            parameter_values.clear();


            //2d. Store the single_parameter in scenario. if the name of the parameter is "scenario", store the value in the name key of the scenario map

            if (single_parameter.name == "scenario") {
                scenario.name = get<string>(single_parameter.values[0]);
            } else {
                scenario.params[single_parameter.name] = single_parameter.values;
            }

            single_parameter.values.clear();
            single_parameter.name.clear();

            
            

        }

        //store the name of the scenario in the name key of the scenario map

        

        //2e. Store the scenario in scenarios

        
        scenarios.push_back(scenario);

        //3. Clear the parameter_values vector

        parameter_values.clear();
        

    }

    //3. Notify the user that the type has been loaded

    cout << "Loaded scenarios from file: " << filename << endl;

    return scenarios;
  
}

Objects load_objects(const string& filename) {
    //This function will load an object for each line in the file.

    ifstream file(filename);

    string line;
    
    vector<string> line_parts;
    string part;
    string name_part;
    string value_part;

    vector<variant<float, string>> parameter_values;
    Parameter single_parameter;
    Object object;
    Objects objects;


    //1. Loop through each line


    while (getline(file, line)) {
        stringstream ss(line);

        //check if the line is empty, if so, skip it
       if (line.empty()) {
            continue;
        }
        //2. Loop through each part of the line, each part is separated by a semi colon. do the following:
        
        
        //cout << "parameter values size: " << parameter_values.size() << endl;

        while (getline(ss, part, ';')) {
            
            
            name_part.clear();
            value_part.clear();

            stringstream ss_part(part);

            //2a. Extract the name of the parameter

            getline(ss_part, name_part, '(');

            //2b. Extract the values of the parameter

            getline(ss_part, value_part, ')');

            //2c. Store the comma separated values in value_part in a vector parameter_values

            stringstream ss_values(value_part);
            float temp_float;
            string temp_string;

            while (getline(ss_values, value_part, ',')) {
                


                stringstream ss_value(value_part);
                if (ss_value >> temp_float) {
                    parameter_values.push_back(temp_float);
                } else {
                    parameter_values.push_back(value_part);
                }
            }


            //2c. Store the name of the parameter and its values in a map called single_parameter
            
            single_parameter.values.clear();
            single_parameter.name.clear();

            single_parameter.name = name_part;
            single_parameter.values = parameter_values;
            parameter_values.clear();

            

            


            //2d. Store the single_parameter in a vector called multiple_parameters

            if (single_parameter.name == "object") {
                object.name = get<string>(single_parameter.values[0]);
            } else {
                object.params[single_parameter.name] = single_parameter.values;
            }
   

        }

        //2e. Store the vector multiple_parameters in a vector called structured_line

        objects.push_back(object);

        //3. Clear the parameter_values vector

        
        object.params.clear();
        object.name.clear();

    }

    //3. Notify the user that the type has been loaded

    cout << "Loaded objects from file: " << filename << endl;

    return objects;
  
}

string select_scenario(Scenarios scenarios) {
    // This function will allow the user to select a scenario from a list of scenarios

    //1. Extract the scenario names and store in a vector. The scenario name is the value of the first key of each map in the vector of maps

    vector<string> scenario_names;
    string key_str;
    string val_str;


    //print the size of the scenarios vector
    //cout << "The size of the scenarios vector is: " << scenarios.size() << endl;


    for (int i = 0; i < scenarios.scenarios.size(); i++) {
        //select element i of scenarios and store in scenario_i

        Scenario scenario_i = scenarios.scenarios[i];

        //store the name of the scenario in val_str

        val_str = scenario_i.name;
    
        scenario_names.push_back(val_str);

    }

    //2. Print the scenario names to the user, give each scenario a number

    cout << "The following scenarios are available: " << endl;

    for (int i = 0; i < scenario_names.size(); i++) {
        cout << i+1 << ": " << scenario_names[i] << endl;
    }

    //3. Choose the scenario for the simulation

    int selection;

    //if there is only one scenario, select it automatically

    if (scenario_names.size() == 1) {
        selection = 1;
        cout << "Automatically selected scenario: " << scenario_names[0] << endl;
    } else {
        cout << "Please enter the number of the scenario you would like to select: ";
        cin >> selection;

        cout << "Selected scenario: " << selection << ": " << scenario_names[selection-1] << endl;

    }

    //4. Return the selected scenario name

    return scenario_names[selection-1];
    
}

Scenario retrieve_scenario(Scenarios scenarios, string selection) {
    //This function will retrieve the selected scenario from the list of scenarios, and broadcast the parameters to the user

    Scenario selected_scenario_params;

    //1. Select the scenario from the list of scenarios

    for (int i = 0; i < scenarios.scenarios.size(); i++) {
        Scenario scenario_i = scenarios.scenarios[i];
        string scenario_i_name = scenario_i.name;


        if (scenario_i_name == selection) {
            selected_scenario_params = scenario_i;
        }
    }

    //2. Print the selected scenario parameters to the user
    cout << endl;
    cout << "The selected scenario has the following parameters: " << endl;

    cout << endl;

    for (const auto& map : selected_scenario_params.params) {
        cout << map.first << ": ";
        for (const auto& variant : map.second) {
            std::visit([](auto&& arg) {
                cout << arg << " ";
            }, variant);
        }
        cout << endl;
    }

    return selected_scenario_params;
}

Objects retrieve_objects(Objects all_objects, Scenario selected_scenario_params) {
    //This function will retrieve the objects associated with the selected scenario, and broadcast the parameters to the user

    Objects selected_objects;

    
    //1. extract the values of the key "objects" in selected_scenario_params, and store them in a vector called requested_objects

    vector<string> requested_objects;

    for (const auto& map : selected_scenario_params.params) {
        if (map.first == "objects") {
            for (const auto& variant : map.second) {
                requested_objects.push_back(get<string>(variant));
            }
        }
    }

    //2. Loop through all_objects and if the value of the first key is in requested_objects, store the entire object in selected_objects

    for (int i = 0; i < all_objects.objects.size(); i++) {
        Object object_i = all_objects.objects[i];
        string object_i_name = object_i.name;

        if (find(requested_objects.begin(), requested_objects.end(), object_i_name) != requested_objects.end()) {
            selected_objects.push_back(object_i);
        }
    }

    //3. Print the selected objects to the user

    cout << endl;
    cout << "The selected objects are: " << endl;
    cout << endl;

    for (int i = 0; i < selected_objects.objects.size(); i++) {
        Object object_i = selected_objects.objects[i];
        cout << "Object " << i+1 << ": " << object_i.name << endl;
        for (const auto& map : object_i.params) {
            cout << map.first << ": ";
            for (const auto& variant : map.second) {
                std::visit([](auto&& arg) {
                    cout << arg << " ";
                }, variant);
            }
            cout << endl;
        }
        cout << endl;
    }
    


    return selected_objects;

}


Particles generate_particles(Objects selected_objects,Scenario selected_scenario_params) {
    //This function will convert the objects to a multidimensional vector of particles. Each particle should have an ID, <3 coordinates>, <3 velocities>, <3 rgb>. Use no variants, but use map

    Particle particle;
    Particles particles;

    //1. Some objects might have complexity, e.g. a cube, a sphere etc instead of being a single point. If this is the case, the object would have key "n" with a value greater than 1 and a shape key with the name of the shape. 
    //Check below is n>1, if so, call the appropriate function to split the object into n points

    Object object;
    Objects flattened_object, flattened_objects;



    for (int i = 0; i < selected_objects.objects.size(); i++) {
        //for each object, check if n < 2, if so, store the object in flattened_objects

        object = selected_objects.objects[i];

        int complex_object_count = 0;

        //object.params is a map structured as follows :map<string, vector<variant<float, string>>>

        //check if object.params has a key "n" and if the value of the key is greater than 1

        if (object.params.count("n") > 0) {
            if (get<float>(object.params["n"][0]) > 1) {

                //this loop will deal with complex object
                complex_object_count++;

                cout << "Complex object found" << endl;
                //print all the keys in the object

                //below, I want to feed the initial object to split_object, not just the current map

                flattened_object = unfold_object(object, selected_scenario_params);

                //loop through each item in flattened_object and store in flattened_objects

                for (int i = 0; i < flattened_object.objects.size(); i++) {
                    flattened_objects.push_back(flattened_object.objects[i]);
                }
                

            } else {
                

                //this loop will deal with simple object

                flattened_objects.push_back(object);
                

                }
            object.params.clear();

            //cout << complex_object_count << " complex objects have been flattened" << endl;
            

            }
        }

    
    cout << "The number of flattened objects is: " << flattened_objects.objects.size() << endl;

    //2. Loop through the flattened objects and convert each object to a particle

    
    for (int i = 0; i < flattened_objects.objects.size(); i++) {
        
        object = flattened_objects.objects[i];
        Particle particle;


        //keep track of parameter names that could not be processed and store them in a vector
        vector<string> unprocessed_keys;

        //2a generate an ID for the particle, based on i and the number of previous particles generated for i.

        particle.ID = i+1;
               
        //Store the coordinates of the particle

        //retrieve the x, y, z coordinates of the particle as the first 3 values of the vector of the pair with key loc_parameters

        for (const auto& pair : object.params) {

            if (pair.first == "loc_parameters") {
                vector<variant<float, string>> loc_parameters = pair.second;

                particle.x = get<float>(loc_parameters[0]);
                particle.y = get<float>(loc_parameters[1]);
                particle.z = get<float>(loc_parameters[2]);

            } else if (pair.first == "linear_parameters") {
                vector<variant<float, string>> linear_parameters = pair.second;

                particle.vx = get<float>(linear_parameters[0]);
                particle.vy = get<float>(linear_parameters[1]);
                particle.vz = get<float>(linear_parameters[2]);

            } else if (pair.first == "rgb_parameters") {
                vector<variant<float, string>> rgb_parameters = pair.second;

                particle.r = get<float>(rgb_parameters[0]);
                particle.g = get<float>(rgb_parameters[1]);
                particle.b = get<float>(rgb_parameters[2]);

                int rgb = ((int)(get<float>(rgb_parameters[0]) * 255) << 16) | ((int)(get<float>(rgb_parameters[1]) * 255) << 8) | (int)(get<float>(rgb_parameters[2]) * 255);

                particle.rgb = rgb;

            } else {
                unprocessed_keys.push_back(pair.first);  

            }
            
        }
        //2e. Store the particle in the vector of particles
        particles.push_back(particle);
        
    }

    cout << "The number of particles generated is: " << particles.particles.size() << endl;

    return particles;
}