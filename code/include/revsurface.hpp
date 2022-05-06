#ifndef REVSURFACE_HPP
#define REVSURFACE_HPP

#include "object3d.hpp"
#include "curve.hpp"
#include <tuple>

class RevSurface : public Object3D {

    Curve *pCurve;

public:
    RevSurface(Curve *pCurve, Material* material) : pCurve(pCurve), Object3D(material) {
        // Check flat.
        for (const auto &cp : pCurve->getControls()) {
            if (cp[2] != 0.0) {
                printf("Profile of revSurface taust be flat on xy plane.\n");
                exit(0);
            }
        }
    }

    ~RevSurface() override {
        delete pCurve;
    }

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        float t, theta, tau, t_in, t_out;
        if (!bbox_intersect(r, t_in, t_out) || t_in > h.getT()) return false;
        Vector3f rayDir = r.getDirection();
        Vector3f rayOri = r.getOrigin();
        bool flag = false;
        int T = 5;
        for (int repeat = 0; repeat <= T; repeat ++) {
            t = (t_out - t_in) * repeat * 1. / T + t_in;
            Vector3f pt(rayOri + rayDir * t);
            theta = atan2(-pt[2], pt[0]) + M_PI;
            tau = (pCurve->y_max - pt[1]) / (pCurve->y_max - pCurve->y_min);
            Vector3f normal;

            if (solve(r, t, theta, tau, normal))
                if (!isnan(tau) && !isnan(theta) && !isnan(t))
                if (!isinf(tau) && !isinf(theta) && !isinf(t))
                if (t > tmin && tau >= pCurve->t_min && tau <= pCurve->t_max && t < h.getT()) {
                    h.set(t, material, normal.normalized());
                    flag = true;
                }
        }
        return flag;
    }

    bool solve(const Ray &r, float &t, float &theta, float &tau, Vector3f &normal) {
        Vector3f dV, dT, point, ray_point;
        Vector3f rayDir = r.getDirection();
        Vector3f rayOri = r.getOrigin();
        Matrix3f rot;
        CurvePoint cp;
        float d, loss, dot_prod;
        
        for (int i = 0; i < 10; i++) {
            theta = fmod(theta, 2 * M_PI);
            if (theta < 0) theta += 2 * M_PI;
            if (tau >= pCurve->t_max) tau = pCurve->t_max;
            if (tau <= pCurve->t_min) tau = pCurve->t_min;

            rot = Matrix3f::rotateY(theta);
            cp = pCurve->get_point(tau);
            point = rot * cp.V;
            dV = rot * cp.T;
            dT = Vector3f(cp.V[0] * sin(theta), 0, cp.V[0] * cos(theta));
            normal = Vector3f::cross(dT, dV);

            ray_point = rayOri + rayDir * t;
            Vector3f f = ray_point - point;
            loss = f.squaredLength();
            if (loss < 1e-7) return true;
            d = Vector3f::dot(rayDir, normal);
            t -= Vector3f::dot(Vector3f::cross(f, dT), dV) / d;
            tau -= Vector3f::dot(Vector3f::cross(f, dT), rayDir) / d;
            theta -= Vector3f::dot(Vector3f::cross(f, dV), rayDir) / d;

            //t -= 0.5 * loss / Vector3f::dot(rayDir, f);
            //tau += 0.5 * loss / Vector3f::dot(cp.T, rot.transposed() * f);
            //dot_prod = -sin(theta) * ray_point[0] * cp.V[0] + cos(theta) * ray_point[0] * cp.V[2] + \
            //           -cos(theta) * ray_point[2] * cp.V[0] - sin(theta) * ray_point[2] * cp.V[2];
            //theta += 0.5 * loss / dot_prod;
        }
        return false;
    }


    bool bbox_intersect(const Ray &r, float & t_in, float & t_out) {

        Vector3f min_point(-pCurve->x_max, pCurve->y_min, -pCurve->x_max);
        Vector3f max_point(pCurve->x_max, pCurve->y_max, pCurve->x_max);

        Vector3f rayDir = r.getDirection();
        Vector3f rayOri = r.getOrigin();
        float t_enter_x, t_enter_y, t_enter_z, t_exit_x, t_exit_y, t_exit_z;

        t_enter_x = (min_point[0] - rayOri[0]) / rayDir[0];
        t_enter_y = (min_point[1] - rayOri[1]) / rayDir[1];
        t_enter_z = (min_point[2] - rayOri[2]) / rayDir[2];

        t_exit_x = (max_point[0] - rayOri[0]) / rayDir[0];
        t_exit_y = (max_point[1] - rayOri[1]) / rayDir[1];
        t_exit_z = (max_point[2] - rayOri[2]) / rayDir[2];

        //if (!isfinite(t_enter_x) || !isfinite(t_exit_x) || isnan(t_enter_x) || isnan(t_exit_x)) std::cout << "error" << std::endl;
        //if (!isfinite(t_enter_y) || !isfinite(t_exit_y) || isnan(t_enter_y) || isnan(t_exit_y)) std::cout << "error" << std::endl;
        //if (!isfinite(t_enter_z) || !isfinite(t_exit_z) || isnan(t_enter_z) || isnan(t_exit_z)) std::cout << "error" << std::endl;
        //if (t_enter_x == t_exit_x) std::cout << "error" << std::endl;
        //if (t_enter_y == t_exit_y) std::cout << "error" << std::endl;
        //if (t_enter_z == t_exit_z) std::cout << "error" << std::endl;

        if (t_enter_x > t_exit_x) {
            t_exit_x = t_enter_x + t_exit_x;
            t_enter_x = t_exit_x - t_enter_x;
            t_exit_x = t_exit_x - t_enter_x;
        }
        if (t_enter_y > t_exit_y) {
            t_exit_y = t_enter_y + t_exit_y;
            t_enter_y = t_exit_y - t_enter_y;
            t_exit_y = t_exit_y - t_enter_y;
        }
        if (t_enter_z > t_exit_z) {
            t_exit_z = t_enter_z + t_exit_z;
            t_enter_z = t_exit_z - t_enter_z;
            t_exit_z = t_exit_z - t_enter_z;
        }

        t_enter_x = max(max(t_enter_x, t_enter_y), t_enter_z);
        t_exit_x = min(min(t_exit_x, t_exit_y), t_exit_z);

        if (t_enter_x > t_exit_x) return false;
        
        t_in = t_enter_x;
        t_out = t_exit_x;
        return true;
    }

};

#endif //REVSURFACE_HPP
