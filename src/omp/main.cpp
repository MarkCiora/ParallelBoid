#include <iostream>
#include <random>
#include <string>
#include <chrono>

#include "vec3.h"
#include "boid.h"
#include "ini.h"

#include <omp.h>

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

    if (ini.has("weights")) {
        if (ini["weights"].has("w_collision")) boid::w_collision = stof(ini.get("weights").get("w_collision"));
        else boid::w_collision = 0.3;

        if (ini["weights"].has("w_alignment")) boid::w_alignment = stof(ini.get("weights").get("w_alignment"));
        else boid::w_collision = 0.4;

        if (ini["weights"].has("w_centering")) boid::w_centering = stof(ini.get("weights").get("w_centering"));
        else boid::w_centering = 0.3;
    }

    boid::new_boids_random();

    double start = 0.0;
    double end = 0.0;
    double total = 0.0;

    start = omp_get_wtime();
    boid::run(time);
    end = omp_get_wtime();
    total = end - start;
    
    std::cout << "Execution Time: " << total << std::endl;

    return 0;
}