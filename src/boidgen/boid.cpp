#include "boid.h"

#include <random>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>

#include "vec3.h"

vec3* boid::pos = nullptr;
vec3* boid::vel = nullptr;
vec3* boid::acc = nullptr;
vec3* boid::sim_boids = nullptr;
vec3 boid::dim_low = vec3(-10., -10., -10.);
vec3 boid::dim_high = vec3(10., 10., 10.);
vec3 boid::vel_low = vec3(-1., -1., -1.);
vec3 boid::vel_high = vec3(1., 1., 1.);
int boid::nboids = 0;
float boid::dt = 1. / (float)(60);
float boid::time = 0.0;

void boid::new_boids_random(int n){
    kill();
    nboids = n;
    pos = new vec3[n];
    vel = new vec3[n];
    acc = new vec3[n];
    vec3 dim_diff = dim_high - dim_low;
    for (int i = 0; i < 3*n; i++){
        ((float*)pos)[i] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        ((float*)pos)[i] *= ((float*)&dim_diff)[i%3];
        ((float*)pos)[i] += ((float*)&dim_low)[i%3];
        ((float*)vel)[i] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        ((float*)vel)[i] *= ((float*)&vel_low)[i%3];
        ((float*)vel)[i] += ((float*)&vel_high)[i%3];
        ((float*)acc)[i] = 0;
    }
}

void boid::kill(){
    delete [] pos;
    delete [] vel;
    delete [] acc;
    delete [] sim_boids;
}

void boid::step_sim(){
    calc_acc();
    physics_update();
    time += dt;
}

void boid::run(float time){
    int steps = static_cast<int>(time / dt) + 1;
    sim_boids = new vec3[steps * nboids];
    for (int i = 0; i < steps; i++){

    }
}

void boid::print_boids(){
    std::cout << std::fixed << std::setprecision(3);
    for (int i = 0; i < nboids; i++){
        std::cout << i << ": " << pos[i] << " + " << dt << "*" << vel[i] << std::endl;
    }
}

void boid::write_sim_boids(std::string outfile){

}

void boid::calc_acc(){
    //eileen code boid rules here
}

void boid::physics_update(){
    for(int i = 0; i < nboids; i++){
        vel[i] += acc[i] * dt;
        pos[i] += vel[i] * dt;
    }
}