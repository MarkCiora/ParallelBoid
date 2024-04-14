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
void sim_kernel(int nboids, int steps, float time, float dt, float gtfo_distance,
                float centering_distance, float alignment_distance,
                float w_col, float w_ali, float w_cen,
                float* dpos, float *dvel, float *dacc,
                float* ddim_low, float* ddim_high,
                float* dvel_low, float* dvelhigh,
                float* data_array){
//stuff
    int index = blockIdx.x*blockDim.x + threadIdx.x;
    //FIRST DATA_ARRAY ENTRY
    if (index < nboids){
        data_array[3*index] = dpos[3*index];
        data_array[3*index+1] = dpos[3*index+1];
        data_array[3*index+2] = dpos[3*index+2];
        // printf("%i, %i, %f, %f, %f \n",
        //             index,
        //             0,
        //             data_array[0*nboids*3 + 3*index],
        //             data_array[0*nboids*3 + 3*index+1],
        //             data_array[0*nboids*3 + 3*index+2]);
    }
    __syncthreads();

    float x, y, z, vx, vy, vz, acx, acy, acz;
    for(int time_step = 1; time_step < steps; time_step++){
        if (index < nboids){
                        
        //calc acc
            float cox=0, coy=0, coz=0;
            float alx=0, aly=0, alz=0;
            float cex=0, cey=0, cez=0;
            // vec3 collision = vec3(0,0,0);
            // vec3 alignment = vec3(0,0,0);
            // vec3 centering = vec3(0,0,0);
            float wax=0, way=0, waz=0;
            // vec3 wa_high = vec3(0,0,0);
            // vec3 wa_low = vec3(0,0,0);
            int collision_counter = 0;
            int alignment_counter = 1;
            int centering_counter = 1;
            x = dpos[3*index];
            y = dpos[3*index+1];
            z = dpos[3*index+2];
            vx = dvel[3*index];
            vy = dvel[3*index+1];
            vz = dvel[3*index+2];
            float avg_diffx=0, avg_diffy=0, avg_diffz=0;
            float avg_velx=vx, avg_vely=vy, avg_velz=vz;
            float avg_posx=x, avg_posy=y, avg_posz=z;
            // vec3 avg_diff = vec3(0,0,0);
            // vec3 avg_vel = vec3(vx,vy,vz);
            // vec3 avg_pos = vec3(x,y,z);

            //all in 1
            for (int j = 1; j < nboids; j++){
                int i = (j + index) % nboids;
                float jx = dpos[3*i];
                float jy = dpos[3*i+1];
                float jz = dpos[3*i+2];
                float dx = jx - x;
                float dy = jy - y;
                float dz = jz - z;
                float vjx = dvel[3*i];
                float vjy = dvel[3*i+1];
                float vjz = dvel[3*i+2];

                //collision
                float distsq = dx*dx + dy*dy + dz*dz;
                if (distsq < gtfo_distance * gtfo_distance){
                    float diff_factor = (gtfo_distance / distsq - rsqrtf(distsq)) * gtfo_distance;
                    avg_diffx -= diff_factor * dx;
                    avg_diffy -= diff_factor * dy;
                    avg_diffz -= diff_factor * dz;
                    collision_counter++;
                }
                if (distsq < alignment_distance * alignment_distance){
                    avg_velx += vjx;
                    avg_vely += vjy;
                    avg_velz += vjz;
                    alignment_counter++;
                }
                if (distsq < centering_distance * centering_distance){
                    avg_posx += jx;
                    avg_posy += jy;
                    avg_posz += jz;
                    centering_counter++;
                }
            }
            if (collision_counter > 0) {
                cox = avg_diffx / ((float)collision_counter);
                coy = avg_diffy / ((float)collision_counter);
                coz = avg_diffz / ((float)collision_counter);
            }
            alx = avg_velx / ((float)alignment_counter) - vx;
            aly = avg_vely / ((float)alignment_counter) - vy;
            alz = avg_velz / ((float)alignment_counter) - vz;
            cex = avg_posx / ((float)centering_counter) - x;
            cey = avg_posy / ((float)centering_counter) - y;
            cez = avg_posz / ((float)centering_counter) - z;

            float diff_top = fabs(10.0 - x);
            float diff_bottom = fabs(10.0 + x);
            float diff_right = fabs(10.0 - y);
            float diff_left = fabs(10.0 + y);
            float diff_front = fabs(10.0 - z);
            float diff_back = fabs(10.0 + z);
            if (diff_top < gtfo_distance) wax = (gtfo_distance / diff_top - 1) * gtfo_distance;
            if (diff_right < gtfo_distance) way = (gtfo_distance / diff_right - 1) * gtfo_distance;
            if (diff_front < gtfo_distance) waz = (gtfo_distance / diff_front - 1) * gtfo_distance;
            if (diff_bottom < gtfo_distance) wax = -(gtfo_distance / diff_bottom - 1) * gtfo_distance;
            if (diff_left < gtfo_distance) way = -(gtfo_distance / diff_left - 1) * gtfo_distance;
            if (diff_back < gtfo_distance) waz = -(gtfo_distance / diff_back - 1) * gtfo_distance;

            acx = (cox)*w_col + alx*w_ali + cex*w_cen - wax;
            acy = (coy)*w_col + aly*w_ali + cey*w_cen - way;
            acz = (coz)*w_col + alz*w_ali + cez*w_cen - waz;
            // printf("pos %i, %i, %f, %f, %f \n",
            //         index,
            //         time_step,
            //         x,
            //         y,
            //         z);
            // printf("vel %i, %i, %f, %f, %f \n",
            //         index,
            //         time_step,
            //         vx,
            //         vy,
            //         vz);
            // printf("acc %i, %i, %f, %f, %f \n",
            //         index,
            //         time_step,
            //         acx,
            //         acy,
            //         acz);
            // printf("co %i, %i, %f, %f, %f \n",
            //         index,
            //         time_step,
            //         cox,
            //         coy,
            //         coz);
            // printf("wa %i, %i, %f, %f, %f \n",
            //         index,
            //         time_step,
            //         wax,
            //         way,
            //         waz);
            // printf("al %i, %i, %f, %f, %f \n",
            //         index,
            //         time_step,
            //         alx,
            //         aly,
            //         alz);
            // printf("ce %i, %i, %f, %f, %f \n",
            //         index,
            //         time_step,
            //         cex,
            //         cey,
            //         cez);
        }
        __syncthreads();

        //physics update
        if (index < nboids){
            float speed = sqrtf(vx*vx + vy*vy + vz*vz);
            float dirx = vx / speed;
            float diry = vy / speed;
            float dirz = vz / speed;
            float acc_mag = sqrtf(acx*acx + acy*acy + acz*acz);
            // printf("%f, %f, %f, %f\n", acx, acy, acz, acc_mag);
            if (acc_mag > 0.00001){
                float adirx = acx / acc_mag;
                float adiry = acy / acc_mag;
                float adirz = acz / acc_mag;
                float dir_dots = dirx*adirx + diry*adiry + dirz*adirz;
                float acx_true = adirx - dirx * dir_dots;
                float acy_true = adiry - diry * dir_dots;
                float acz_true = adirz - dirz * dir_dots;
                float true_mag = sqrtf(acx_true*acx_true + acy_true*acy_true + acz_true*acz_true);
                acx_true = acx_true * acc_mag / true_mag;
                acy_true = acy_true * acc_mag / true_mag;
                acz_true = acz_true * acc_mag / true_mag;
                vx += acx_true * dt;
                vy += acy_true * dt;
                vz += acz_true * dt;
                speed = sqrtf(vx*vx + vy*vy + vz*vz);
                vx /= speed;
                vy /= speed;
                vz /= speed;
            }
            x += vx * dt;
            y += vy * dt;
            z += vz * dt;
            dpos[3*index] = x;
            dpos[3*index+1] = y;
            dpos[3*index+2] = z;
            dvel[3*index] = vx;
            dvel[3*index+1] = vy;
            dvel[3*index+2] = vz;

            //DATA_ARRAY ENTRY
            data_array[time_step*nboids*3 + 3*index] = x;
            data_array[time_step*nboids*3 + 3*index+1] = y;
            data_array[time_step*nboids*3 + 3*index+2] = z;
            // printf("%i, %i, %f, %f, %f \n",
            //         index,
            //         time_step,
            //         data_array[time_step*nboids*3 + 3*index],
            //         data_array[time_step*nboids*3 + 3*index+1],
            //         data_array[time_step*nboids*3 + 3*index+2]);
            // printf("\n");
        }
        __syncthreads();
    }

    // if (index < nboids)
    // for (int time_step = 0; time_step < steps; time_step++){
    //         printf("%i, %i, %f, %f, %f \n",
    //                 index,
    //                 time_step,
    //                 data_array[time_step*nboids*3 + 3*index],
    //                 data_array[time_step*nboids*3 + 3*index+1],
    //                 data_array[time_step*nboids*3 + 3*index+2]);
    // }
    
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
    float *data_array = nullptr;

    cudaMalloc(&dpos, nboids*3*sizeof(float));
    cudaMalloc(&dvel, nboids*3*sizeof(float));
    cudaMalloc(&dacc, nboids*3*sizeof(float));
    cudaMalloc(&ddim_low, 3*sizeof(float));
    cudaMalloc(&ddim_high, 3*sizeof(float));
    cudaMalloc(&dvel_low, 3*sizeof(float));
    cudaMalloc(&dvel_high, 3*sizeof(float));
    cudaMalloc(&data_array, 3*nboids*steps*sizeof(float));

    cudaMemcpy(dpos, (float*)pos, nboids*3*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dvel, (float*)vel, nboids*3*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(ddim_low, (float*)(&dim_low), 3*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(ddim_high, (float*)(&dim_high), 3*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dvel_low, (float*)(&vel_low), 3*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dvel_high, (float*)(&vel_high), 3*sizeof(float), cudaMemcpyHostToDevice);

    sim_kernel<<<(nboids + 255)/256, 256>>>(
                nboids, steps, time, dt, gtfo_distance,
                centering_distance, alignment_distance,
                w_collision, w_alignment, w_centering,
                dpos, dvel, dacc,
                ddim_low, ddim_high,
                dvel_low, dvel_high,
                data_array);

    cudaDeviceSynchronize();

    cudaMemcpy((float*)(sim_boids), data_array, 3*nboids*steps*sizeof(float), cudaMemcpyDeviceToHost);

    cudaDeviceSynchronize();
    cudaFree(dpos);
    cudaFree(dvel);
    cudaFree(dacc);
    cudaFree(ddim_high);
    cudaFree(ddim_low);
    cudaFree(dvel_high);
    cudaFree(dvel_low);
    cudaFree(data_array);
    
    // printf("bozo\n");
    // for (int time_step = 0; time_step < steps; time_step++){
    //         printf("%i, %i, %f, %f, %f \n",
    //                 0,
    //                 time_step,
    //                 ((float*)(sim_boids))[time_step*nboids*3 + 3*0],
    //                 ((float*)(sim_boids))[time_step*nboids*3 + 3*0+1],
    //                 ((float*)(sim_boids))[time_step*nboids*3 + 3*0+2]);
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