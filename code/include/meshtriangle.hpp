#ifndef MESHTRIANGLE_H
#define MESHTRIANGLE_H

#include "object3d.hpp"
#include "texture.hpp"
#include <vecmath.h>
#include <cmath>
#include <iostream>
#include <cstdlib>
using namespace std;


class MeshTriangle: public Object3D {

public:
	MeshTriangle() = default;

	MeshTriangle(Material* m) : Object3D(m) {};

	int vertex_id[3];
	int norm_id[3];
	int uv_id[3];

	bool interpolate;
	bool use_texture;
	bool use_uv;
	Texture* texture;

	Vector3f vertex[3], norm[3];
	Vector2f uv[3];

	void initiate(std::vector<Vector3f> *vertices, \
				  std::vector<Vector2f> *uvs, \
				  std::vector<Vector3f> *norms, \
				  bool interpolate) {
		vertex[0] = (*vertices)[vertex_id[0]];
		vertex[1] = (*vertices)[vertex_id[1]];
		vertex[2] = (*vertices)[vertex_id[2]];
		norm[0] = (*norms)[norm_id[0]].normalized();
		norm[1] = (*norms)[norm_id[1]].normalized();
		norm[2] = (*norms)[norm_id[2]].normalized();
		if (uv_id[0] == -1 || uv_id[1] == -1 || uv_id[2] == -1 || uvs->size() <= 0) {
			uv[0] = 0;
			uv[1] = 0;
			uv[2] = 0;
			use_uv = false;
		} else {
			uv[0] = (*uvs)[uv_id[0]];
			uv[1] = (*uvs)[uv_id[1]];
			uv[2] = (*uvs)[uv_id[2]];
			use_uv = true;
		}
        this->interpolate = interpolate;
		this->use_texture = false;
	}

	void set_texture(Texture* text) {
		texture = text;
		use_texture = use_uv;
	}

	bool intersect( const Ray& ray,  Hit& hit , float tmin) override {
		Vector3f Rd = ray.getDirection();
        Vector3f R0 = ray.getOrigin();
		Vector3f E1 = vertex[0] - vertex[1];
		Vector3f E2 = vertex[0] - vertex[2];
		Vector3f S = vertex[0] - R0;
		float det0 = Matrix3f(Rd, E1, E2, true).determinant();
		if (det0 == 0) return false;
		float det1 = Matrix3f(S, E1, E2, true).determinant();
		float det2 = Matrix3f(Rd, S, E2, true).determinant();
		float det3 = Matrix3f(Rd, E1, S, true).determinant();
		Vector3f res = 1.0 / det0 * Vector3f(det1, det2, det3);

		if (res[0] > tmin && res[0] < hit.getT() && res[1] >= 0 && res[2] >= 0 && res[1] + res[2] <= 1){
    		Vector3f norm_hit = (norm[0] + norm[1] + norm[2]) / 3;
			Vector3f pos_hit = (1 - res[1] - res[2]) * vertex[0] + res[1] * vertex[1] + res[2] * vertex[2];
			Vector2f uv_hit = (1 - res[1] - res[2]) * uv[0] + res[1] * uv[1] + res[2] * uv[2];
			Material* m_hit = material;
			
			if (use_texture) m_hit = texture->get_material(uv_hit);
			if (interpolate) norm_hit = (1 - res[1] - res[2]) * norm[0] + res[1] * norm[1] + res[2] * norm[2]; 

			// compute bump
			if (use_texture && texture->use_bump_map) {
				Vector3f diff_i(uv[1][0] - uv_hit[0], uv[2][0] - uv_hit[0], 0);
				Vector3f diff_j(uv[1][1] - uv_hit[1], uv[2][1] - uv_hit[1], 0);
				Matrix3f A = Matrix3f(vertex[1] - pos_hit, vertex[2] - pos_hit, norm_hit, false).inverse();
				Vector3f tangent = (A * diff_i).normalized();
				Vector3f bitangent = (A * diff_j).normalized();
				Vector3f uv_norm = texture->get_bump(uv_hit);
				norm_hit = (uv_norm[0] * tangent + uv_norm[1] * bitangent +  uv_norm[2] * norm_hit).normalized();
			}
			hit.set(res[0], m_hit, norm_hit);
            return true;
		}
        return false;
	}

protected:

};

#endif //MESHTRIANGLE_H
