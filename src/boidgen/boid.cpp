#include "boid.h"

#include <algorithm>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <fstream>

#include "vec3.h"

vec3* boid::pos = nullptr;
vec3* boid::vel = nullptr;
vec3* boid::acc = nullptr;
vec3* boid::sim_boids = nullptr;
vec3 boid::dim_low = vec3(-5., -5., -5.);
vec3 boid::dim_high = vec3(5., 5., 5.);
vec3 boid::vel_low = vec3(-1., -1., -1.);
vec3 boid::vel_high = vec3(1., 1., 1.);
vec3 boid::center = vec3(0., 0., 0.);
int boid::nboids = 2;
int boid::steps = 0;
float boid::dt = 1. / (float)(60);
float boid::time = 0.0;

float gtfo_distance = 0.5;
float w_collision = 0.4;
float w_alignment = 0.4;
float w_centering = 0.3;

void boid::new_boids_random(){
    kill();
    pos = new vec3[nboids];
    vel = new vec3[nboids];
    acc = new vec3[nboids];
    vec3 dim_diff = dim_high - dim_low;
    for (int i = 0; i < 3*nboids; i++){
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
    calc_acc_all();
    physics_update();
    time += dt;
}

void boid::run(float time){
    steps = static_cast<int>(time / dt) + 1;
    int sim_boids_index = 0;
    sim_boids = new vec3[steps * nboids];
    std::cout << "Step 0: " << std::endl;
    for (int j = 0; j < nboids; j++){
        sim_boids[sim_boids_index] = pos[j];
        sim_boids_index++;
    }
    //print_boids();
    for (int i = 1; i < steps; i++){
        std::cout << "Step " << i + 1 << ": " << std::endl;
        step_sim();
        for (int j = 0; j < nboids; j++){
            sim_boids[sim_boids_index] = pos[j];
            sim_boids_index++;
        }
        print_boids();
    }
    write_sim_boids();
}

void boid::print_boids(){
    std::cout << std::fixed << std::setprecision(3);
    for (int i = 0; i < nboids; i++){
        std::cout << i << ": " << pos[i] << " + " << dt << "*" << vel[i] << std::endl;
    }
    std::cout << "center: " << center << std::endl;
}

void boid::write_sim_boids(){
    //write
    // int nboids
    // int steps
    // float dt
    // float time
    // float array 3*nboids*steps
    std::ofstream file("boid_data.bin", std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        exit(-1);
    }
    file.write(reinterpret_cast<char*>(&nboids), sizeof(nboids));
    file.write(reinterpret_cast<char*>(&steps), sizeof(steps));
    file.write(reinterpret_cast<char*>(&dt), sizeof(dt));
    file.write(reinterpret_cast<char*>(&time), sizeof(time));
    file.write(reinterpret_cast<char*>(sim_boids), sizeof(vec3) * nboids * steps);
    file.close();
}

// calculate acceleration using all boids with each other
void boid::calc_acc_all(){
    vec3* collision = nullptr;
    vec3* alignment = nullptr;
    vec3* centering = nullptr;
    int counter = 0;

    // collision
    vec3 avg_diff = vec3(0,0,0);
    collision = new vec3[nboids];

    for (int i = 0; i < nboids; i++){
        counter = 0;
        avg_diff.clear();
        for (int j = 0; j < nboids; j++){
            if (j == i) continue;
            vec3 diff = (pos[i] - pos[j]);
            float distsq = diff.normsqrd();
            if (distsq < gtfo_distance * gtfo_distance){
                avg_diff += (diff / distsq - diff / diff.norm()) * gtfo_distance;
                counter++;
            }
        }
        if (counter > 0) {
            avg_diff /= counter;
            collision[i] = avg_diff;
        }
    }

    // alignment
    vec3 avg_vel = vec3(0,0,0);
    alignment = new vec3[nboids];

    for (int i = 0; i < nboids; i++){
        avg_vel += vel[i];
    }
    avg_vel /= nboids;

    for (int i = 0; i < nboids; i++){
        alignment[i] = avg_vel - vel[i];
    }

    // centering
    vec3 avg_pos = vec3(0,0,0);
    centering = new vec3[nboids];

    for (int i = 0; i < nboids; i++){
        avg_pos += pos[i];
    }
    avg_pos /= nboids;

    for (int i = 0; i < nboids; i++){
        centering[i] = avg_pos - pos[i];
    }

    // calculate acceleration for all boids
    for (int i = 0; i < nboids; i++){
        acc[i] = ((collision[i] * w_collision) + (alignment[i] * w_alignment) + (centering[i] * w_centering)).normalized();
    }

    delete [] alignment;
    delete [] collision;
    delete [] centering;
}

void boid::set_center_all(){
    center.clear();
    for (int i = 0; i < nboids; i++){
        center += pos[i];
    }
    center /= nboids;
}

void boid::physics_update(){
    for(int i = 0; i < nboids; i++){
        vel[i] += acc[i] * dt;
        // if (vel[i].normsqrd() >= 1){
            vel[i] /= vel[i].norm();
        // }
        pos[i] += vel[i] * dt;
    }
}