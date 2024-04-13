#include "boid.h"

#include <algorithm>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <fstream>
#include <cmath>

#include "vec3.h"

vec3* boid::pos = nullptr;
vec3* boid::vel = nullptr;
vec3* boid::acc = nullptr;
vec3* boid::sim_boids = nullptr;
vec3 boid::dim_low = vec3(-10., -10., -10.);
vec3 boid::dim_high = vec3(10., 10., 10.);
vec3 boid::vel_low = vec3(-1., -1., -1.);
vec3 boid::vel_high = vec3(1., 1., 1.);
vec3 boid::center = vec3(0., 0., 0.);
int boid::nboids = 2;
int boid::steps = 0;
float boid::dt = 1. / (float)(60);
float boid::time = 0.0;
float boid::centering_distance = 1.5;
float boid::alignment_distance = 1.;

float gtfo_distance = 1;
float boid::w_collision = 0.4;
float boid::w_alignment = 0.4;
float boid::w_centering = 0.3;

void boid::new_boids_random(){
    kill();
    pos = new vec3[nboids];
    vel = new vec3[nboids];
    acc = new vec3[nboids];
    vec3 dim_diff = dim_high - dim_low;
    vec3 vel_diff = vel_high - vel_low;
    for (int i = 0; i < 3*nboids; i++){
        ((float*)pos)[i] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        ((float*)pos)[i] *= ((float*)&dim_diff)[i%3];
        ((float*)pos)[i] += ((float*)&dim_low)[i%3];
        ((float*)vel)[i] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        ((float*)vel)[i] *= ((float*)&vel_diff)[i%3];
        ((float*)vel)[i] += ((float*)&vel_low)[i%3];
        ((float*)acc)[i] = 0;
    }
    for (int i = 0; i < nboids; i++) vel[i].normalize();
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

__global__
void sim_kernel(int nboids, int steps, float time, float dt,
                float centering_distance, float alignment_distance,
                float w_col, float w_ali, float w_cen,
                float* dpos, float *dvel, float *dacc,
                float* ddim_low, float* ddim_high,
                float* dvel_low, float* dvelhigh,
                float* dcenter){

}

void boid::run(float time){
    steps = static_cast<int>(time / dt) + 1;
    int sim_boids_index = 0;
    sim_boids = new vec3[steps * nboids];

    for (int j = 0; j < nboids; j++){
        sim_boids[sim_boids_index] = pos[j];
        sim_boids_index++;
    }

    float *dpos = nullptr;
    float *dvel = nullptr;
    float *dacc = nullptr;
    float *ddim_low = nullptr;
    float *ddim_high = nullptr;
    float *dvel_low = nullptr;
    float *dvel_high = nullptr;
    float *dcenter = nullptr;

    sim_kernel<<<1,1>>>( nboids, steps, time, dt,
                centering_distance, alignment_distance,
                w_collision, w_alignment, w_centering,
                dpos, dvel, dacc,
                ddim_low, ddim_high,
                dvel_low, dvel_high,
                dcenter);

    // for (int i = 1; i < steps; i++){
    //     step_sim();
    //     for (int j = 0; j < nboids; j++){
    //         sim_boids[sim_boids_index] = pos[j];
    //         sim_boids_index++;
    //     }
    // }

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
    vec3* wall_avoidance_high = nullptr;
    vec3* wall_avoidance_low = nullptr;
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
                avg_diff += (diff * gtfo_distance / distsq - diff / diff.norm()) * gtfo_distance;
                counter++;
            }
        }
        if (counter > 0) {
            avg_diff /= counter;
            collision[i] = avg_diff;
        }
    }

    // alignment
    vec3* avg_vel = new vec3[nboids];
    alignment = new vec3[nboids];

    for (int i = 0; i < nboids; i++){
        int alignment_counter = 1;
        avg_vel[i] = vec3(0,0,0);
        for (int j = 0; j < nboids; j++){
            if ((pos[i] - pos[j]).norm() <= alignment_distance){
                alignment_counter++;
                avg_vel[i] = avg_vel[i] + vel[j];
            }
        }
        avg_vel[i] = avg_vel[i] / ((float)alignment_counter);
    }

    for (int i = 0; i < nboids; i++){
        alignment[i] = avg_vel[i] - vel[i];
    }

    // centering
    vec3* avg_pos = new vec3[nboids];
    centering = new vec3[nboids];

    for (int i = 0; i < nboids; i++){
        int centering_counter = 1;
        avg_pos[i] = vec3(0,0,0);
        for (int j = 0; j < nboids; j++){
            if ((pos[i] - pos[j]).norm() <= centering_distance){
                centering_counter++;
                avg_pos[i] = avg_pos[i] + pos[j];
            }
        }
        avg_pos[i] = avg_pos[i] / ((float)centering_counter);
    }

    for (int i = 0; i < nboids; i++){
        centering[i] = avg_pos[i] - pos[i];
    }

    // walls
    wall_avoidance_high = new vec3[nboids];
    wall_avoidance_low = new vec3[nboids];

    for (int i = 0; i < nboids; i++){
        float diff_top = std::fabs(pos[i].x - dim_high.x);
        float diff_bottom = std::fabs(pos[i].x - dim_low.x);

        float diff_right = std::fabs(pos[i].y - dim_high.y);
        float diff_left = std::fabs(pos[i].y - dim_low.y);

        float diff_front = std::fabs(pos[i].z - dim_high.z);
        float diff_back = std::fabs(pos[i].z - dim_low.z);

        if (diff_top < gtfo_distance){
            wall_avoidance_high[i].x = -(gtfo_distance / diff_top - 1) * gtfo_distance;
        }
        
        if (diff_right < gtfo_distance){
            wall_avoidance_high[i].y = -(gtfo_distance / diff_right - 1) * gtfo_distance;
        }

        if (diff_front < gtfo_distance){
            wall_avoidance_high[i].z = -(gtfo_distance / diff_front - 1) * gtfo_distance;
        }

        if (diff_bottom < gtfo_distance){
            wall_avoidance_low[i].x = (gtfo_distance / diff_bottom - 1) * gtfo_distance;
        }

        if (diff_left < gtfo_distance){
            wall_avoidance_low[i].y = (gtfo_distance / diff_left - 1) * gtfo_distance;
        }

        if (diff_back < gtfo_distance){
            wall_avoidance_low[i].z = (gtfo_distance / diff_back - 1) * gtfo_distance;
        }
    }


    // calculate acceleration for all boids
    for (int i = 0; i < nboids; i++){
        acc[i] = (collision[i] * w_collision) + (alignment[i] * w_alignment) + (centering[i] * w_centering) + wall_avoidance_high[i] + wall_avoidance_low[i];
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
        vec3 dir = vel[i].normalized();
        vec3 acc_dir = acc[i].normalized();
        float dot_product = dir.x * acc_dir.x + dir.y * acc_dir.y + dir.z * acc_dir.z;
        vec3 true_acc_dir = (acc_dir - dir * dot_product).normalized();
        float true_acc_mag = acc[i].norm();
        // if (true_acc_mag >= 3) true_acc_mag = 3;
        vel[i] += true_acc_dir * true_acc_mag * dt;
        vel[i].normalize();
        pos[i] += vel[i] * dt;
        // vel[i] += acc[i] * dt;
        // if (vel[i].normsqrd() >= 1){
        //     vel[i] /= vel[i].norm();
        // }
        // pos[i] += vel[i] * dt;
    }
}