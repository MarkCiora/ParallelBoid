#pragma once

#include <string>

#include "vec3.h"

class boid{
public:
    static vec3* pos;
    static vec3* vel;
    static vec3* acc;
    static vec3* sim_boids;
    static vec3 dim_low, dim_high;
    static vec3 vel_low, vel_high;
    static vec3 acc_low, acc_high;
    static int nboids;
    static float dt;
    static float time;

    static void new_boids_random();
    static void kill();
    static void step_sim();
    static void run(float time = 10.0);
    
    static void print_boids();
    static void write_sim_boids(std::string outfile);

private:
    static void calc_acc();
    static void physics_update();
};