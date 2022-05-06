#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <random>
#include <fstream>

#include "scene_parser.hpp"
#include "image.hpp"
#include "camera.hpp"
#include "group.hpp"
#include "light.hpp"
#include "surfacelight.hpp"

#include <string>

using namespace std;


std::default_random_engine random_engine(time(NULL)); 
std::uniform_real_distribution<float> rand_nums(0.0, 1.0);


Vector3f random_sample(Vector3f norm, float angle) { // sample from the hemisphere
    Vector3f dx, dy;
    float hit_proj = sqrt(norm[0] * norm[0] + norm[1] * norm[1]);
    if (hit_proj > 0.1) {
        dx = Vector3f(norm[1] / hit_proj, -norm[0] / hit_proj, 0);
        dx.normalize();
    }
    else {
        hit_proj = sqrt(norm[0] * norm[0] + norm[2] * norm[2]);
        dx = Vector3f(norm[2] / hit_proj, 0, -norm[0] / hit_proj);
        dx.normalize();
    }
    dy = Vector3f::cross(norm, dx);

    float xi_1 = 1 - rand_nums(random_engine) * (1 - cos(angle));
    float xi_2 = rand_nums(random_engine);

    Vector3f dir = dx * cos(2 * M_PI * xi_2) * sqrt(1 - xi_1 * xi_1) + dy * sin(2 * M_PI * xi_2) * sqrt(1 - xi_1 * xi_1) + norm * xi_1;
    return dir;      
}


void path_trace(SceneParser &scene, Group* baseGroup, Ray ray, int depth, float weight, float tmin, \
                Vector3f &color, bool no_light) {
    color = Vector3f::ZERO;
    if (weight < 0.01) return;
    Hit hit;
    bool isIntersect = baseGroup->intersect(ray, hit, tmin);
    float t_light = 99999;
    bool intersectLight = false;
    Vector3f lightColor;
    for (int li=0; li<scene.getNumSurfaceLights(); li++) {
        SurfaceLight* light = scene.getSurfaceLight(li);
        // record the nearest intersect with light
        intersectLight = intersectLight || light -> intersect(ray, t_light, lightColor, tmin);
    }
    // compare the intersection with the nearest object and the nearest light
    // if intersect nothing, return backgound
    // if intersect light first, return the BRDF result of light color
    if (!intersectLight && !isIntersect) {
        color = scene.getBackgroundColor();
        return;
    }
    if (intersectLight && t_light < hit.getT()) {
        if (!no_light) color += lightColor;
        return;
    }
    if (depth <= 0) return;

    Vector3f hitPoint = ray.pointAtParameter(hit.getT());
    if ((*hit.getMaterial()).emission) color += (*hit.getMaterial()).emissionColor;
        
    float w_r = (*hit.getMaterial()).w_reflection;
    float w_t = (*hit.getMaterial()).w_transmission;
    float w_d = (*hit.getMaterial()).w_diffusion;

    // goal: sample a new ray w.r.t. w_r, w_t and w_d (to reflect, transmiss or to diffuse)
    
    Vector3f hitNorm = hit.getNormal().normalized();  // norm of the hit point
    Vector3f rayNorm = ray.getDirection().normalized();  // norm of the ray
    float cos_theta_i = -Vector3f::dot(hitNorm, rayNorm);
    Vector3f tDir;
    Vector3f rDir = 2 * hitNorm * cos_theta_i + rayNorm;  // get reflect direction

    if (w_t > 1e-5) { // simulate a transmission first (to get the real w_r and w_t and the reflect + transmiss ray) 
        float n = (*hit.getMaterial()).ref_index;
        float cos_theta_t, w_decay;
        if (cos_theta_i < 0) { // inner surface
            n = 1 / n;
            hitNorm = -hitNorm;
            cos_theta_i = -cos_theta_i;
            tDir = (rayNorm + cos_theta_i * hitNorm) / n;
            cos_theta_t = sqrt(1 - tDir.length() * tDir.length());
            w_decay = pow(1 - cos_theta_t, 5);
        } else {
            tDir = (rayNorm + cos_theta_i * hitNorm) / n;
            cos_theta_t = sqrt(1 - tDir.length() * tDir.length());
            w_decay = pow(1 - cos_theta_i, 5);
        }
        if (tDir.length() >= 1) { // fully reflective
            w_r = w_r + w_t;  // tDir is no longer important
            w_t = 0;
        } else {
            w_t -= w_t * w_decay;
            w_r += w_t * w_decay;
            tDir = tDir - hitNorm * cos_theta_t;  
        }
    }

    Ray tRay(hitPoint, tDir);  // get the true transmission ray
    Ray rRay(hitPoint, rDir); // get the true reflection ray
    //bool reflect = false, transmiss = false, diffuse = false;  // select only one path to reduce complexity (not used)
    int sample_repeat = 1;
    for (int i=0; i<sample_repeat; i++) {
        // mode selection (not used)
        // float choice = rand_nums(random_engine);
        //  if (choice < w_d) diffuse = 1;
        //  else {
        //      if (w_t > 1e-5) transmiss = 1;
        //      if (w_r > 1e-5) reflect = 1;
        //  }
        if (w_d > 1e-5) {  // if (diffuse) for "mode selection" mode
            // sample uniformly from the light
            int light_repeat = 1;
            for (int lri=0; lri<light_repeat; lri++) {
                for (int li=0; li<scene.getNumSurfaceLights(); li++) {
                    SurfaceLight* light = scene.getSurfaceLight(li);
                    // get light direction
                    Vector3f lightDir = light -> getPosition() - hitPoint;
                    float lightDist = lightDir.length();
                    //float lightAngle = atan(light -> getRadius() / lightDist);
                    lightDir.normalize();
                    Vector3f dDir;
                    hitNorm = hit.getNormal().normalized();

                    // sample a point on the light source
                    float rand_1, rand_2;
                    rand_1 = rand_nums(random_engine);
                    rand_2 = rand_nums(random_engine);
                    Vector3f lightPosition = light -> getSample(rand_1, rand_2);
                    float area = light -> getArea(lightDir);
                    float area_weight = area / (2 * M_PI * lightDist * lightDist);
                    dDir = lightPosition - hitPoint;
                    dDir.normalize();

                    Vector3f shaded = Vector3f::ZERO; // the shaded color
                    if (Vector3f::dot(dDir, hitNorm) > 0) {
                        Ray dRay(hitPoint, dDir);
                        Vector3f newColor;
                        path_trace(scene, baseGroup, dRay, 0, weight * w_d, tmin, newColor, false);

                        // get diffuse colors 
                        Vector3f aColor = (*hit.getMaterial()).ambientColor;
                        Vector3f dColor = (*hit.getMaterial()).diffuseColor;
                        Vector3f sColor = (*hit.getMaterial()).specularColor;
                        float shininess = (*hit.getMaterial()).shininess;
                        float a_o = (*hit.getMaterial()).ambient_occlusion;

                        // start shading (BRDF)
                        float cos_theta_t = Vector3f::dot(dDir, hitNorm);
                        Vector3f R = 2 * cos_theta_t * hitNorm - dDir;
                        float cos_theta_specular = Vector3f::dot(R, -rayNorm);
                            
                        shaded += aColor * a_o * scene.getAmbientColor();
                        if (cos_theta_t > 0) shaded += newColor * dColor * cos_theta_t;
                        if (cos_theta_specular > 0) shaded += newColor * sColor * pow(cos_theta_specular, shininess);
                    }
                    color += shaded * 1.0 * w_d * area_weight / light_repeat / sample_repeat;
                }
            }
            // sample for a ray with light turned off
            Vector3f dDir, randomDir, lightDir;
            float lightAngle;
            hitNorm = hit.getNormal().normalized();
            Vector3f shaded = Vector3f::ZERO; // the shaded color
            dDir = random_sample(hitNorm, M_PI / 2);
            Ray dRay(hitPoint, dDir);
            Vector3f newColor;
            path_trace(scene, baseGroup, dRay, depth-1, weight * w_d, tmin, newColor, true);

            // get diffuse colors 
            Vector3f aColor = (*hit.getMaterial()).ambientColor;
            Vector3f dColor = (*hit.getMaterial()).diffuseColor;
            Vector3f sColor = (*hit.getMaterial()).specularColor;
            float shininess = (*hit.getMaterial()).shininess;
            float a_o = (*hit.getMaterial()).ambient_occlusion;

            // start shading (BRDF)
            float cos_theta_t = Vector3f::dot(dDir, hitNorm);
            Vector3f R = 2 * cos_theta_t * hitNorm - dDir;
            float cos_theta_specular = Vector3f::dot(R, -rayNorm);

            shaded += aColor * a_o * scene.getAmbientColor();
            if (cos_theta_t > 0) shaded += newColor * dColor * cos_theta_t;
            if (cos_theta_specular > 0) shaded += newColor * sColor * pow(cos_theta_specular, shininess);

            color += shaded * 1.0 * w_d / sample_repeat;
        } 
        if (w_t > 1e-5) { // if (transmiss) for "mode selection" mode
            Vector3f newColor;
            path_trace(scene, baseGroup, tRay, depth-1, weight * w_t, tmin, newColor, no_light);
            color += newColor * 1.0 / sample_repeat * w_t;
        } 
        if (w_r > 1e-5) { // if (reflect) for "mode selection" mode
            Vector3f newColor;
            path_trace(scene, baseGroup, rRay, depth-1, weight * w_r, tmin, newColor, no_light);
            color += newColor * 1.0 / sample_repeat * w_r;
        }
    }
}


void path_tracing_main(string inputFile, string outputFile) {
    // Main RayCasting Logic
    // First, parse the scene using SceneParser.

    SceneParser scene(inputFile.c_str());
    Group* baseGroup = scene.getGroup();
    Camera* camera = scene.getCamera();
    Image img(camera->getWidth(), camera->getHeight());

    // Then loop over each pixel in the image, shooting a ray
    // through that pixel and finding its intersection with
    // the scene.  Write the color at the intersection to that
    // pixel in your output image.
    
    int depth = 5;  // path sampling depth
    float weight = 1.0;  // path sampling weight
    float tmin = 1e-3; // path sampling min intersection t
    float repeat_x[5] = {0,  0.25, -0.25, 0.25, -0.25};  // anti-aliasing
    float repeat_y[5] = {0, -0.25,  0.25, 0.25, -0.25};  // anti-aliasing
    int repeat = 5; // repetition for anti-aliasing
    int repeat_aperture = 1; // repetition for aperture sampling
    int round = 5; // repetition for path tracing main cycle (to reduce noise)

    for (int r=1; r<=round; r++) {
        for (int x=0; x<camera->getWidth(); x++) {
            for (int y=0; y<camera->getHeight(); y++) {
                Vector3f finalColor;
                for (int t=0; t<repeat; t++) {
                    if (camera->use_aperture) {
                        for (int t_a=0; t_a<repeat_aperture; t_a++) {
                            float rand_1 = rand_nums(random_engine);
                            float rand_2 = rand_nums(random_engine);
                            camera->generateOrigin(rand_1, rand_2);
                            Ray camRay = camera->generateRay(Vector2f(x * 1.0 + repeat_x[t], y * 1.0 + repeat_y[t]));
                            Vector3f color = Vector3f::ZERO;
                            path_trace(scene, baseGroup, camRay, depth, weight, tmin, color, false);
                            finalColor += 1.0 / repeat / repeat_aperture * color;
                            //std::cout << "round " << r << ", point " << x << " " << y << ", repeat (anti-aliase) " << t << ", repeat (aperture) " << t_a << std::endl;
                        }
                    } else {
                        Ray camRay = camera->generateRay(Vector2f(x * 1.0 + repeat_x[t], y * 1.0 + repeat_y[t]));
                        Vector3f color = Vector3f::ZERO;
                        path_trace(scene, baseGroup, camRay, depth, weight, tmin, color, false);
                        finalColor += 1.0 / repeat * color;
                    }
                }
                img.AddPixel(x, y, finalColor * 1.0 / round);
            }
            //std::cout << "round " << r << ", point " << x << std::endl;
        }
        std::cout << "round " << r << " completed. " << std::endl;
        img.SaveImage(outputFile.c_str());
    }
}


void ray_trace(SceneParser &scene, Group* baseGroup, Ray ray, int depth, float weight, float tmin, Vector3f &color) {
    color = Vector3f::ZERO;
    if (weight < 0.01) return;
    Hit hit;
    bool isIntersect = baseGroup->intersect(ray, hit, tmin);
    if (!isIntersect) {
        color = scene.getBackgroundColor();
        return;
    }
    Vector3f hitPoint = ray.pointAtParameter(hit.getT());
    if ((*hit.getMaterial()).w_diffusion > 1e-5) {
        for (int li=0; li<scene.getNumLights(); li++) {
            Light* light = scene.getLight(li);
            Vector3f lightDir, lightColor;
            light->getIllumination(hitPoint, lightDir, lightColor); // lightDir: point to light source
            bool isShadow = light->isShadow(baseGroup, lightDir, hitPoint, tmin);
            if (!isShadow) color += (*hit.getMaterial()).w_diffusion * hit.getMaterial()->Shade(ray, hit, lightDir, lightColor);
        }
        color += (*hit.getMaterial()).w_diffusion * (*hit.getMaterial()).ambientColor * (*hit.getMaterial()).ambient_occlusion * scene.getAmbientColor();
    }
    if (depth > 1) {
        float w_r = (*hit.getMaterial()).w_reflection;
        float w_t = (*hit.getMaterial()).w_transmission;
        if (w_t > 1e-5) { // transmission
            Vector3f hitNorm = hit.getNormal().normalized();
            Vector3f rayNorm = ray.getDirection().normalized();
            float cos_theta_i = -Vector3f::dot(hitNorm, rayNorm);
            float n = (*hit.getMaterial()).ref_index;
            float cos_theta_t, w_decay;
            Vector3f newDir;
            if (cos_theta_i < 0) { // inner surface
                n = 1 / n;
                hitNorm = -hitNorm;
                cos_theta_i = -cos_theta_i;
                newDir = (rayNorm + cos_theta_i * hitNorm) / n;
                cos_theta_t = sqrt(1 - newDir.length() * newDir.length());
                w_decay = pow(1 - cos_theta_t, 5);
            } else {
                newDir = (rayNorm + cos_theta_i * hitNorm) / n;
                cos_theta_t = sqrt(1 - newDir.length() * newDir.length());
                w_decay = pow(1 - cos_theta_i, 5);
            }
            if (newDir.length() >= 1) { // fully reflective
                w_r = w_r + w_t;
            } else {
                w_t -= w_t * w_decay;
                w_r += w_t * w_decay;
                newDir = newDir - hitNorm * cos_theta_t;
                Ray refRay(hitPoint, newDir);
                Vector3f newColor;
                ray_trace(scene, baseGroup, refRay, depth-1, weight * w_t, tmin, newColor);
                color += w_t * newColor;
            }
        }
        if (w_r > 1e-5) { // reflection
            Vector3f hitNorm = hit.getNormal().normalized();
            Vector3f rayNorm = ray.getDirection().normalized();
            float product = -Vector3f::dot(hitNorm, rayNorm);
            if (product < 0) { // inner surface
                hitNorm = -hitNorm;
                product = -product;
            }
            Vector3f newDir = 2 * hitNorm * product + rayNorm;
            Ray refRay(hitPoint, newDir);
            Vector3f newColor;
            ray_trace(scene, baseGroup, refRay, depth-1, weight * w_r, tmin, newColor);
            color += w_r * newColor;
        }
    }
}


void ray_tracing_main(string &inputFile, string &outputFile) {
    // Main RayCasting Logic
    // First, parse the scene using SceneParser.
    SceneParser scene(inputFile.c_str());
    Group* baseGroup = scene.getGroup();
    Camera* camera = scene.getCamera();
    Image img(camera->getWidth(), camera->getHeight());

    // Then loop over each pixel in the image, shooting a ray
    // through that pixel and finding its intersection with
    // the scene.  Write the color at the intersection to that
    // pixel in your output image.
    
    int depth = 10;  // depth
    float weight = 1.0;
    float tmin = 1e-3;
    float repeat_x[5] = {0,  0.25, -0.25, 0.25, -0.25};  // anti-aliasing
    float repeat_y[5] = {0, -0.25,  0.25, 0.25, -0.25};  // anti-aliasing
    int repeat = 5; // repetition for anti-aliasing
    int repeat_aperture = 10; // repetition for aperture sampling

    for (int x=0; x<camera->getWidth(); x++) {
        for (int y=0; y<camera->getHeight(); y++) {
            Vector3f finalColor;
            for (int t=0; t<repeat; t++) {
                if (camera->use_aperture) {
                    for (int t_a=0; t_a<repeat_aperture; t_a++) {
                        float rand_1 = rand_nums(random_engine);
                        float rand_2 = rand_nums(random_engine);
                        camera->generateOrigin(rand_1, rand_2);
                        Ray camRay = camera->generateRay(Vector2f(x * 1.0 + repeat_x[t], y * 1.0 + repeat_y[t]));
                        Vector3f color = Vector3f::ZERO;
                        ray_trace(scene, baseGroup, camRay, depth, weight, tmin, color);
                        finalColor += 1.0 / repeat / repeat_aperture * color;
                    }
                } else {
                    Ray camRay = camera->generateRay(Vector2f(x * 1.0 + repeat_x[t], y * 1.0 + repeat_y[t]));
                    Vector3f color = Vector3f::ZERO;
                    ray_trace(scene, baseGroup, camRay, depth, weight, tmin, color);
                    finalColor += 1.0 / repeat * color;
                }
            }
            img.SetPixel(x, y, finalColor);
        }
        //std::cout << "point " << x << std::endl;
    }
    img.SaveImage(outputFile.c_str());
}


int main(int argc, char *argv[]) {
    for (int argNum = 1; argNum < argc; ++argNum) {
        std::cout << "Argument " << argNum << " is: " << argv[argNum] << std::endl;
    }
    if (argc != 4) {
        cout << "Usage: ./bin/FinalProject <input scene file> <output bmp file> <trace type>" << endl;
        return 1;
    }
    string inputFile = argv[1];
    string outputFile = argv[2];  // only bmp is allowed.
    string traceType = argv[3];
    string ray("ray"), path("path");
    if (traceType == ray) ray_tracing_main(inputFile, outputFile);
    else if (traceType == path) path_tracing_main(inputFile, outputFile);
    return 0;
}

