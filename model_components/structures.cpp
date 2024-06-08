//This file contains all structures needed for the PlotUniverse tool

#include <vector>
#include <map>
#include <variant>
#include <string>

using namespace std;

struct Particle {
    //This struct will store the properties of a particle
    int ID;
    float x;
    float y;
    float z;
    float vx;
    float vy;
    float vz;
    float r;
    float g;
    float b;
    int rgb;
    double time_scaler = 1;

};

struct Particles {
    //This struct will store a vector of particles
    vector<Particle> particles;

    void push_back(const Particle particle) {
        particles.push_back(particle);
    }
};

struct Scenario {
    //This struct will store the parameters of a scenario
    string name;
    map<string, vector<variant<float, string>>> params;


};

struct Scenarios {
    //This struct will store a vector of scenarios
    vector<Scenario> scenarios;

    void push_back(const Scenario scenario) {
        scenarios.push_back(scenario);
    }
};

struct Object {
    //This struct will store the parameters of an object
    string name;
    map<string, vector<variant<float, string>>> params;
};

struct Objects {
    //This struct will store a vector of objects
    vector<Object> objects;

    void push_back(const Object object) {
        objects.push_back(object);
    }
};

struct Parameter {
    //This struct will store the parameters of a parameter
    string name;
    vector<variant<float, string>> values;
};





