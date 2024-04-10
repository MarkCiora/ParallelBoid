#include "obstacle.h"

class wall: public obstacle {
    vec3 top_left;
    vec3 top_right;
    vec3 bottom_left;
    vec3 bottom_right;

    void generate(vec3 dim_low, vec3 dim_high) {
        int x = 0;
    }
};