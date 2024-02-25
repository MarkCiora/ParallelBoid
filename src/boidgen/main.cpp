#include <iostream>
#include <random>
#include <string>

#include "vec3.h"
#include "boid.h"

int main(int argv, char **argc){

    if (argv == 1){
        boid::new_boids_random();
    } else if (argv == 2){
        boid::nboids = std::stoi(argc[1]);
        boid::new_boids_random();
    }

    boid::run(1.0);

    return 0;
}