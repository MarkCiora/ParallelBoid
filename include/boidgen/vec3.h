#pragma once

#include <iostream>

class vec3{
public:
    float x,y,z;
    vec3();
    vec3(float x_, float y_, float z_);
    void clear();
    float norm();
    float normsqrd();
    static vec3 cross(const vec3& v1, const vec3& v2);
    vec3 normalized();
    vec3& normalize();
    vec3 operator+(const vec3& other) const;
    vec3& operator+=(const vec3& other);
    vec3 operator-(const vec3& other) const;
    vec3& operator-=(const vec3& other);
    vec3 operator*(const float& other) const;
    vec3& operator*=(const float& other);
    vec3 operator/(const float& other) const;
    vec3& operator/=(const float& other);
    friend std::ostream& operator<<(std::ostream& os, const vec3& obj);
};