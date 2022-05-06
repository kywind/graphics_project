#ifndef CURVE_HPP
#define CURVE_HPP

#include "object3d.hpp"
#include <vecmath.h>
#include <vector>
#include <utility>

#include <algorithm>

// Bernstein class to compute spline basis function.

// The CurvePoint object stores information about a point on a curve
// after it has been tesselated: the vertex (V) and the tangent (T)
// It is the responsiblility of functions that create these objects to fill in all the data.
struct CurvePoint {
    Vector3f V; // Vertex
    Vector3f T; // Tangent  (unit)
};

class Curve : public Object3D {
protected:
    std::vector<Vector3f> controls;
public:
    explicit Curve(std::vector<Vector3f> points) : controls(std::move(points)) {
        y_min = controls[0][1];
        y_max = controls[0][1];
        x_max = controls[0][0];
        for (int i=0; i<(int)controls.size(); i++) {
            y_min = min(controls[i][1], y_min);
            y_max = max(controls[i][1], y_max);
            x_max = max(abs(controls[i][0]), x_max);
        }
    }

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        return false;
    }

    std::vector<Vector3f> &getControls() {
        return controls;
    }

    virtual CurvePoint get_point(float t) = 0;

    virtual void discretize(int resolution, std::vector<CurvePoint>& data) = 0;

    float t_min, t_max, y_min, y_max, x_max;
};


class BezierCurve : public Curve {
public:
    explicit BezierCurve(const std::vector<Vector3f> &points) : Curve(points) {
        if (points.size() < 4 || points.size() % 3 != 1) {
            printf("Number of control points of BezierCurve must be 3n+1!\n");
            exit(0);
        }
        n = controls.size() - 1;
        mem = new float[2*(n+2)];
        t_min = 0 + 1e-3;
        t_max = 1 - 1e-3;
    }

    CurvePoint get_point(float t) override {
        CurvePoint point;
        calcBaseFunc(n, t);
        for (int i = 0; i <= n; i ++) {
            point.V += mem[n+2+i] * controls[i];
            point.T += n * (mem[i] - mem[i+1]) * controls[i];
        }
        return point;
    }

    void discretize(int resolution, std::vector<CurvePoint>& data) override {
        data.clear();
        float tc = 0;
        for (int j = 0; j < resolution; j ++) {
            CurvePoint point;
            calcBaseFunc(n, tc);
            for (int i = 0; i <= n; i ++) {
                point.V += mem[n+2+i] * controls[i];
                point.T += n * (mem[i] - mem[i+1]) * controls[i];
            }
            point.T.normalize();
            data.push_back(point);
            tc += 1.0 / resolution;
        }
    }

protected:

    int n;
    float* mem;

    float pow(float n, int a) {
        if (a <= 0) {
            return 1.0;
        }
        return n * pow(n, a-1);
    } 
    
    void calcBaseFunc(int n, float t) {  // only calculate the basefunc of n and n-1
        mem[0] = 0.;
        mem[1] = pow(1-t, n-1);
        for (int i = 1; i <= n-1; i ++) {
            mem[i+1] = mem[i] * (n-i) * t / i / (1-t);
        }
        mem[n+1] = 0.;
        mem[n+2] = pow(1-t, n);
        for (int i = 1; i <= n; i ++) {
            mem[i+n+2] = mem[i+n+1] * (n-i+1) * t / i / (1-t);
        }
    }
};


class BsplineCurve : public Curve {
public:
    BsplineCurve(const std::vector<Vector3f> &points) : Curve(points) {
        if (points.size() < 4) {
            printf("Number of control points of BspineCurve must be more than 4!\n");
            exit(0);
        }
        k = 3;
        n = controls.size() - 1;
        knots = new float[n+k+2];
        mem = new float[(n+k+1)*(k+1)];
        for (int i = 0; i < n+k+2; i ++) knots[i] = 1. * i / (n+k+1);
        t_min = knots[k] + 1e-3;
        t_max = knots[n+1] - 1e-3;
    }

    CurvePoint get_point(float t) override {
        int i = k;
        for (int j = 0; j < n+k+1; j ++) { // for row 0
            if (t >= knots[j] && t < knots[j+1]) i = j;
        }
        //if (i < k || i > n) std::cout << "error: " << i << " ( " << k << ", " << n << " ) " << std::endl;
        calcBaseFunc(i, t);
        CurvePoint point;
        for (int kk = k; kk >= 0; kk --) { // kk: index of control point
            point.V += mem[i-kk+k*(n+k+1)] * controls[i-kk];
            point.T += k * (mem[i-kk+(k-1)*(n+k+1)] / (knots[i-kk+k] - knots[i-kk]) - \
                    mem[i-kk+1+(k-1)*(n+k+1)] / (knots[i-kk+k+1] - knots[i-kk+1])) * controls[i-kk];
        }
        return point;
    }

    void discretize(int resolution, std::vector<CurvePoint>& data) override {
        data.clear();
        for (int i = k; i <= n; i ++) { // iterate through control points; i: the control point 
            float tc = knots[i];
            calcBaseFunc(i, tc);
            for (int j = 0; j < resolution; j ++) { // iterate through sample points (#=resolution per control point)
                CurvePoint point;
                for (int kk = k; kk >= 0; kk --) { // kk: index of control point
                    point.V += mem[i-kk+k*(n+k+1)] * controls[i-kk];
                    point.T += k * (mem[i-kk+(k-1)*(n+k+1)] / (knots[i-kk+k] - knots[i-kk]) - \
                            mem[i-kk+1+(k-1)*(n+k+1)] / (knots[i-kk+k+1] - knots[i-kk+1])) * controls[i-kk];
                }
                point.T.normalize();
                data.push_back(point);
                tc += (knots[i+1] - knots[i]) / resolution;
            }
        }
    }

protected:
    float *mem;
    float *knots;
    int k, n;
    void calcBaseFunc(int i, float tc) {
        for (int j = 0; j < n+k+1; j ++) { // for row 0
            if (tc >= knots[j] && tc < knots[j+1]) mem[j] = 1.;
            else mem[j] = 0.;
        }
        for (int p = 1; p <= k; p ++) { // for every nonzero row
            for (int j = i-k; j <= i+k-p; j ++) { // for every entry on row p
                mem[j+p*(n+k+1)] = \
                        (tc - knots[j]) / (knots[j+p] - knots[j]) * mem[j+(p-1)*(n+k+1)] + \
                        (knots[j+p+1] - tc) / (knots[j+p+1] - knots[j+1]) * mem[j+1+(p-1)*(n+k+1)];
            }
        }
    }
};

#endif // CURVE_HPP
