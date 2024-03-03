#include <iostream>
#include <random>
#include <string>

#include "vec3.h"
#include "boid.h"

int main(int argv, char **argc){

    float time = 10.0;
    if (argv == 1){
        boid::new_boids_random();
    } else if (argv == 2){
        boid::nboids = std::stoi(argc[1]);
        boid::new_boids_random();
    }else if (argv == 3){
        boid::nboids = std::stoi(argc[1]);
        time = std::stof(argc[2]);
        boid::new_boids_random();
    }


    boid::run(time);

    return 0;
}