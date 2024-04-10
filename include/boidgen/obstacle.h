#pragma once

#include "vec3.h"

class obstacle {
public:
    static vec3* pos;
    virtual void generate();
    virtual void kill();
};