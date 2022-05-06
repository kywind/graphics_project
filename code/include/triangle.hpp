#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>
#include <iostream>
using namespace std;

// TODO: implement this class and add more fields as necessary,
class Triangle: public Object3D {

public:
	Triangle() = delete;

    // a b c are three vertex positions of the triangle
	Triangle( const Vector3f& a, const Vector3f& b, const Vector3f& c, Material* m) : Object3D(m) {
		normal = Vector3f::cross(b-a, c-a).normalized();
		vertices[0] = a;
		vertices[1] = b;
		vertices[2] = c;
	}

	Triangle( const Vector3f& a, const Vector3f& b, const Vector3f& c) {
		normal = Vector3f::cross(b-a, c-a).normalized();
		vertices[0] = a;
		vertices[1] = b;
		vertices[2] = c;
	}

	bool intersect( const Ray& ray,  Hit& hit , float tmin) override {
		Vector3f Rd = ray.getDirection();
        Vector3f R0 = ray.getOrigin();
		Vector3f E1 = vertices[0] - vertices[1];
		Vector3f E2 = vertices[0] - vertices[2];
		Vector3f S = vertices[0] - R0;
		float det0 = Matrix3f(Rd, E1, E2, true).determinant();
		if (det0 == 0) return false;
		float det1 = Matrix3f(S, E1, E2, true).determinant();
		float det2 = Matrix3f(Rd, S, E2, true).determinant();
		float det3 = Matrix3f(Rd, E1, S, true).determinant();
		Vector3f res = 1.0 / det0 * Vector3f(det1, det2, det3);
		if (res[0] > tmin && res[0] < hit.getT() && res[1] >= 0 && res[2] >= 0 && res[1] + res[2] <= 1){
			hit.set(res[0], material, normal);
            return true;
		}
        return false;
	}
	Vector3f normal;
	Vector3f vertices[3];
protected:

};

#endif //TRIANGLE_H
