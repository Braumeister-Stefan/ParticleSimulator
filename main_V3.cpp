//This project allows the user to plot n points on a 2D plane given the x,y coordinates of each point. 

//other files in project
#include "model_components/initialisation_model.cpp"
#include "model_components/simulation_model.cpp"




//namespaces
using namespace std;


int main() {
    

    //1. Load run parameters from text file6
    Scenarios scenarios = load_scenarios("inputs\\scenario_params.txt");
    Objects all_objects = load_objects("inputs\\object_params.txt");
    

    //2. Allow user to select a scenario and list scenario parameters and object parameters

    string selection = select_scenario(scenarios);

    //3. Import the user defined scenario and associated objects

    Scenario selected_scenario = retrieve_scenario(scenarios, selection);
    Objects selected_objects = retrieve_objects(all_objects, selected_scenario);

    //4. Convert the objects to a multidimensional vector of particles. Each particle should have an ID, <3 coordinates>, <3 velocities>, <3 rgb>
    Particles particles = generate_particles(selected_objects, selected_scenario);

    //5. Start the simulation
    
    simulate_particles(particles, selected_scenario);

    return 0;

}



