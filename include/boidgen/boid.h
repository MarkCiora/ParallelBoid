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
    static vec3 center;
    static int nboids;
    static int steps;
    static float dt;
    static float time;

    static void new_boids_random();
    static void kill();
    static void step_sim();
    static void run(float time = 10.0);
    
    static void print_boids();

private:
    static void set_center_all();

    static void calc_acc_all();
    static void physics_update();
    static void write_sim_boids();

    static float get_dist(vec3, vec3);
};