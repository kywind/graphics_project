#ifndef SFLIGHT_H
#define SFLIGHT_H

#include <Vector3f.h>
#include <vector>
#include "object3d.hpp"
#include "hit.hpp"
#include "ray.hpp"
#include "group.hpp"
#include "mesh.hpp"

class SurfaceLight {
    
public:
    SurfaceLight() = delete;

    SurfaceLight(const Vector3f &p, const Vector3f &c, const Vector3f n, const Vector3f &x, const float &a, const float &b) {
        position = p;
        color = c;
        norm = n.normalized();
        x_axis = x.normalized();
        x_len = a;
        y_len = b;
        y_axis = Vector3f::cross(norm, x_axis);
        offset = Vector3f::dot(position, norm);
    }
    
    ~SurfaceLight() = default;

    bool intersect(const Ray &r, float &min_t, Vector3f &c, float tmin) {
        Vector3f rayDir = r.getDirection().normalized();
        Vector3f rayOri = r.getOrigin();
        float ray_cosine = Vector3f::dot(norm, rayDir);
        float ray_offset = Vector3f::dot(norm, rayOri);
        if (ray_cosine == 0) return false;
        float t = (offset - ray_offset) / ray_cosine;

        Vector3f hitCoor = rayOri + t * rayDir - position;
        float hitX = Vector3f::dot(hitCoor, x_axis);
        float hitY = Vector3f::dot(hitCoor, y_axis);

        if (t > tmin && t < min_t && hitX > -x_len / 2 && hitX < x_len / 2 && hitY > -y_len / 2 && hitY < y_len / 2) {
            min_t = t;
            if (ray_cosine < 0) c = color;  // uniform light at every angle
            return true;
        }
        return false;
    }

    Vector3f getPosition() {
        return position;
    }

    Vector3f getSample(float rand1, float rand2) {  // uniform in range (0, 1)
        rand1 = rand1 - 0.5;
        rand2 = rand2 - 0.5;
        return position + x_axis * rand1 * x_len + y_axis * rand2 * y_len;
    }

    float getArea(Vector3f dir) {
        return x_len * y_len * abs(Vector3f::dot(norm, dir.normalized()));
    }

protected:

    Vector3f position;
    float offset;
    Vector3f color;
    Vector3f norm;
    Vector3f x_axis;
    Vector3f y_axis;
    float x_len;
    float y_len;

};

#endif // SFLIGHT_H

