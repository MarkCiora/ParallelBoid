#include "vec3.h"

#include <iostream>
#include <iomanip>
#include <cmath>

vec3::vec3(){
    x = 0; y = 0; z = 0;
}

vec3::vec3(float x_, float y_, float z_){
    x = x_;
    y = y_;
    z = z_;
}

void vec3::clear(){
    x = 0; y = 0; z = 0;
}

float vec3::norm(){
    return std::sqrt(x*x + y*y + z*z);
}

float vec3::normsqrd(){
    return x*x + y*y + z*z;
}

vec3 vec3::normalized(){
    return *this / norm();
}

vec3& vec3::normalize(){
    *this = *this / norm();
    return *this;
}

vec3 vec3::cross(const vec3& v1, const vec3& v2){
    float x = v1.y * v2.z - v1.z * v2.y;
    float y = v1.z * v2.x - v1.x * v2.z;
    float z = v1.x * v2.y - v1.y * v2.x;
    return vec3(x,y,z);
}

vec3 vec3::operator+(const vec3& other) const{
    return vec3(x + other.x, y + other.y, z + other.z);
}

vec3& vec3::operator+=(const vec3& other){
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}

vec3 vec3::operator-(const vec3& other) const{
    return vec3(x - other.x, y - other.y, z - other.z);
}

vec3& vec3::operator-=(const vec3& other){
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
}

vec3 vec3::operator*(const float& other) const{
    return vec3(x * other, y * other, z * other);
}

vec3& vec3::operator*=(const float& other){
    x *= other;
    y *= other;
    z *= other;
    return *this;
}

vec3 vec3::operator/(const float& other) const{
    return vec3(x / other, y / other, z / other);
}

vec3& vec3::operator/=(const float& other){
    x /= other;
    y /= other;
    z /= other;
    return *this;
}

std::ostream& operator<<(std::ostream& os, const vec3& obj) {
    os << "(" << 
        std::setw(6) << obj.x << ", " << 
        std::setw(6) << obj.y << ", " << 
        std::setw(6) << obj.z << ")";
    return os;
}