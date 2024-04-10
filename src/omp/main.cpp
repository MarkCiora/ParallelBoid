#include <iostream>
#include <random>
#include <string>
#include <chrono>

#include "vec3.h"
#include "boid.h"
#include "ini.h"

int main(int argv, char **argc){
    // read from parameter file
    auto currentTime = std::chrono::system_clock::now();
    auto duration = currentTime.time_since_epoch();
    double currentTimeInSeconds = std::chrono::duration_cast<std::chrono::duration<int>>(duration).count();
    srand (currentTimeInSeconds);
    mINI::INIFile file("parameters.ini");
    mINI::INIStructure ini;
    file.read(ini);

    float time = 0;

    if (ini.has("boids")) {
        if (ini["boids"].has("nboids")) boid::nboids = stoi(ini.get("boids").get("nboids"));
        else boid::nboids = 10;

        if (ini["boids"].has("time")) time = stof(ini.get("boids").get("time"));
        else time = 10.0;
    }

    boid::new_boids_random();
    boid::run(time);

    return 0;
}