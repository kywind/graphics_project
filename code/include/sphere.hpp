#ifndef SPHERE_H
#define SPHERE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>


class Sphere : public Object3D {
public:
    Sphere() {
        // unit ball at the center
        this->center = Vector3f::ZERO;
        this->radius = 1.0;
    }

    Sphere(const Vector3f &center, float radius, Material *material) : Object3D(material) {
        // 
        this->center = center;
        this->radius = radius;
    }

    ~Sphere() override = default;

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        //
        Vector3f Rd = r.getDirection().normalized();
        float len = r.getDirection().length();
        Vector3f R0 = r.getOrigin();
        Vector3f l = center - R0;
        float l2 = l.squaredLength();
        float r2 = radius * radius;
        float tp = Vector3f::dot(l, Rd);
        float d2 = l2 - tp * tp;
        if (d2 > r2) return false;
        float tpr = std::sqrt(r2 - d2);
        float t;
        if ((tp - tpr) / len > tmin) t = (tp - tpr) / len;
        else if ((tp + tpr) / len > tmin) t = (tp + tpr) / len;
        else return false;
        if (t < h.getT()) {
            h.set(t, material, (r.pointAtParameter(t) - center).normalized());
            return true;
        }
        return false;
    }

protected:
    Vector3f center;
    float radius;

};


#endif
