#ifndef PLANE_H
#define PLANE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>

// Plane representing an infinite plane
// function: ax+by+cz=d
// choose your representation , add more fields and fill in the functions

class Plane : public Object3D {
public:
    Plane() {
        this->n = Vector3f::ZERO;
        this->d = 0.0;
    }

    Plane(const Vector3f &normal, float d, Material *m) : Object3D(m) {
        this->n = normal;
        this->d = d;
    }

    ~Plane() override = default;

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        Vector3f Rd = r.getDirection();
        Vector3f R0 = r.getOrigin();
        float a = Vector3f::dot(n, Rd);
        float b = Vector3f::dot(n, R0);
        if (a == 0) return false;
        float t = (d - b) / a;
        if (t > tmin && t < h.getT()){
            h.set(t, material, n);
            return true;
        }
        return false;
    }

protected:
    Vector3f n;
    float d;

};

#endif //PLANE_H
		

