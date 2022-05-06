#ifndef MESH_H
#define MESH_H

#include <vector>
#include "object3d.hpp"
#include "meshtriangle.hpp"
#include "material.hpp"
#include "texture.hpp"
#include "triangle.hpp"
#include "Vector2f.h"
#include "Vector3f.h"


class Mesh : public Object3D {

public:

    Mesh() = default;
    Mesh(const char *filename, const bool interpolate, Material *m);

    Mesh(Material *m): Object3D(m) {}

	Texture* texture;
    bool use_texture;

    std::vector<Vector3f> *v;
    std::vector<MeshTriangle> *t;
    std::vector<Vector3f> *n;
    std::vector<Vector2f> *tex;

    std::vector<int> v_id;
    std::vector<int> t_id;

    Vector3f minPos;
    Vector3f maxPos;
    int divide_dir; // 0: x, 1: y, 2: z
    bool is_leaf;
    Mesh* left;
    Mesh* right;
    bool interpolate;

    bool intersect(const Ray &r, Hit &h, float tmin) override;

    void set_texture(Texture* text) {
        texture = text;
        use_texture = true;
        for (int triId = 0; triId < (int) t->size(); ++triId) {
            (*t)[triId].set_texture(text);
        }
    }

private:

    // Normal can be used for light estimation
    void computeNormal();

    void computeKdTree(int stop_iters);

    static bool v_sort_x (Vector3f a, Vector3f b) { 
        return a[0] < b[0]; 
    }
    static bool v_sort_y (Vector3f a, Vector3f b) { 
        return a[1] < b[1]; 
    }
    static bool v_sort_z (Vector3f a, Vector3f b) { 
        return a[2] < b[2]; 
    }
    static bool testIntersection(Vector3f a, Vector3f b, Vector3f c, Vector3f minP, Vector3f maxP) {
        float l0 = (maxP[0]-minP[0])/2, l1 = (maxP[1]-minP[1])/2, l2 = (maxP[2]-minP[2])/2;
        if (l0 <= 0 || l1 <= 0 || l2 <= 0) return false;
        Vector3f center((minP[0]+maxP[0])/2, (minP[1]+maxP[1])/2, (minP[2]+maxP[2])/2);
        Vector3f p0, p1, p2;
        Vector3f u0(1, 0, 0), u1(0, 1, 0), u2(0, 0, 1);
        p0 = a - center; 
        p1 = b - center;
        p2 = c - center;
        Vector3f e0 = p1 - p0;
        if (!testAxisEdge(p0, p1, p2, e0, l0, l1, l2)) return false;
        Vector3f e1 = p2 - p1;
        if (!testAxisEdge(p1, p2, p0, e1, l0, l1, l2)) return false;
        Vector3f e2 = p0 - p2;
        if (!testAxisEdge(p2, p0, p1, e2, l0, l1, l2)) return false;
        if (min(min(p0[0], p1[0]), p2[0]) > l0 || max(max(p0[0], p1[0]), p2[0]) < -l0) return false;
        if (min(min(p0[1], p1[1]), p2[1]) > l1 || max(max(p0[1], p1[1]), p2[1]) < -l1) return false;
        if (min(min(p0[2], p1[2]), p2[2]) > l2 || max(max(p0[2], p1[2]), p2[2]) < -l2) return false;
        Vector3f normal = Vector3f::cross(p1-p0, p2-p0).normalized();
        float d = -Vector3f::dot(normal, p0);
        Vector3f pmin, pmax;
        if (normal[0] > 0) { pmin[0] = -l0; pmax[0] = l0; } else { pmin[0] = l0; pmax[0] = -l0; }
        if (normal[1] > 0) { pmin[1] = -l1; pmax[1] = l1; } else { pmin[1] = l1; pmax[1] = -l1; }
        if (normal[2] > 0) { pmin[2] = -l2; pmax[2] = l2; } else { pmin[2] = l2; pmax[2] = -l2; }
        if (Vector3f::dot(normal, pmin) + d > 0 || Vector3f::dot(normal, pmax) + d < 0) return false;
        return true;
    }
    static bool testAxisEdge(Vector3f v0, Vector3f v1, Vector3f v2, Vector3f e, float l0, float l1, float l2) {
        float x = abs(e[0]), y = abs(e[1]), z = abs(e[2]);
        float s0 = -e[2] * v0[1] + e[1] * v0[2];
        float s2 = -e[2] * v2[1] + e[1] * v2[2];
        float r = z * l1 + y * l2;
        if (max(s0, s2) < -r || min(s0, s2) > r) return false;
        s0 = e[2] * v0[0] - e[0] * v0[2];
        s2 = e[2] * v2[0] - e[0] * v2[2];
        r = z * l0 + x * l2;
        if (max(s0, s2) < -r || min(s0, s2) > r) return false;
        s0 = -e[1] * v0[0] + e[0] * v0[1];
        s2 = -e[1] * v2[0] + e[0] * v2[1];
        r = y * l0 + x * l1;
        if (max(s0, s2) < -r || min(s0, s2) > r) return false;
        return true;
    }

};

#endif


