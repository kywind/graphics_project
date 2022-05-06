#ifndef MATERIAL_H
#define MATERIAL_H

#include <cassert>
#include <vecmath.h>

#include "ray.hpp"
#include "hit.hpp"
#include <iostream>

// Shade function that computes Phong introduced in class.
class Material {

public:

    explicit Material(const Vector3f &d_color = Vector3f::ZERO, \
                      const Vector3f &s_color = Vector3f::ZERO, \
                      float s = 0, \
                      Vector3f a_color = Vector3f::ZERO, \
                      float w_d = 1, \
                      float w_r = 0, \
                      float w_t = 0, \
                      float n = 1, \
                      float a_o = 1, \
                      bool emi = false, \
                      Vector3f emi_color = Vector3f::ZERO): \
                      diffuseColor(d_color), \
                      specularColor(s_color), \
                      shininess(s), \
                      ambientColor(a_color), \
                      w_diffusion(w_d), \
                      w_reflection(w_r), \
                      w_transmission(w_t), \
                      ref_index(n), \
                      ambient_occlusion(a_o), \
                      emission(emi), \
                      emissionColor(emi_color) {}

    virtual ~Material() = default;

    Vector3f Shade(const Ray &ray, const Hit &hit, const Vector3f &dirToLight, const Vector3f &lightColor) {
        Vector3f shaded = Vector3f::ZERO;
        Vector3f normal = hit.getNormal();
        float diffuse = Vector3f::dot(dirToLight, normal);
        Vector3f R = 2 * diffuse * normal - dirToLight;
        float specular = Vector3f::dot(R, -ray.getDirection());
        if (diffuse > 0) shaded += lightColor * diffuseColor * diffuse;
        if (specular > 0) shaded += lightColor * specularColor * pow(specular, shininess);
        return shaded;
    }

    Vector3f diffuseColor;
    Vector3f specularColor;
    Vector3f ambientColor;
    float shininess;
    float w_diffusion;
    float w_reflection;
    float w_transmission;
    float ref_index;
    float ambient_occlusion;
    bool emission;
    Vector3f emissionColor;

protected:
    
};


#endif // MATERIAL_H
