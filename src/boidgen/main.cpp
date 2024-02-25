#include <iostream>
#include <random>

#include "vec3.h"
#include "boid.h"

int main(int argv, char **argc){

    boid::new_boids_random(10);
    boid::print_boids();

    return 0;
}