#ifndef CAMERA_H
#define CAMERA_H

#include "ray.hpp"
#include <vecmath.h>
#include <float.h>
#include <cmath>


class Camera {
public:
    Camera(const Vector3f &center, const Vector3f &direction, const Vector3f &up, int imgW, int imgH) {
        this->center = center;
        this->direction = direction.normalized();
        this->horizontal = Vector3f::cross(this->direction, up).normalized();
        this->up = Vector3f::cross(this->horizontal, this->direction);
        this->width = imgW;
        this->height = imgH;
    }

    // Generate rays for each screen-space coordinate
    virtual Ray generateRay(const Vector2f &point) = 0;

    virtual ~Camera() = default;

    virtual void generateOrigin(const float theta, const float r) = 0;

    int getWidth() const { return width; }
    int getHeight() const { return height; }
    bool use_aperture;

protected:
    // Extrinsic parameters
    Vector3f center;
    Vector3f direction;
    Vector3f up;
    Vector3f horizontal;
    // Intrinsic parameters
    int width;
    int height;
};

// Perspective camera
// You can add new functions or variables whenever needed.
class PerspectiveCamera : public Camera {

public:
    PerspectiveCamera(const Vector3f &center, const Vector3f &direction,
            const Vector3f &up, int imgW, int imgH, float angle) : Camera(center, direction, up, imgW, imgH) {
        // angle is in radian.
        this->fx = height / 2.0 / std::tan(angle / 2.0);
        this->fy = fx;
        use_aperture = false;
    }

    Ray generateRay(const Vector2f &point) override {
        // 
        Vector3f v((point[0] - width / 2.0) / fx, (height / 2.0 - point[1]) / fy, 1);
        return Ray(center, Matrix3f(horizontal, -up, direction, true) * v.normalized());
    }

    void generateOrigin(const float theta, const float r) override {};

protected:
    float fx;
    float fy;

};

class ApertureCamera : public Camera {

public:
    ApertureCamera(const Vector3f &center, const Vector3f &direction,
            const Vector3f &up, int imgW, int imgH, float angle, float f_a, float u_a) : Camera(center, direction, up, imgW, imgH) {
        // angle is in radian.
        this->fx = height / 2.0 / std::tan(angle / 2.0);
        this->fy = fx;
        f_aperture = f_a;
        u_aperture = u_a;
        use_aperture = true;
    }

    Ray generateRay(const Vector2f &point) override {
        Vector3f v((point[0] - width / 2.0) / fx, (height / 2.0 - point[1]) / fy, 1);
        Vector3f rayDir = Matrix3f(horizontal, -up, direction, true) * v;
        Vector3f u_pos = center + u_aperture * rayDir;
        return Ray(origin, (u_pos - origin).normalized());
    }

    void generateOrigin(const float theta, float r) override {
        origin[0] = 1 / f_aperture * r * sin(theta * 2 * M_PI);
        origin[1] = 1 / f_aperture * r * cos(theta * 2 * M_PI);
        origin[2] = 0;
        origin = Matrix3f(horizontal, -up, direction, true) * origin + center;
    }

protected:
    float fx;
    float fy;
    float f_aperture;
    float u_aperture;
    Vector3f origin;
};

#endif //CAMERA_H
