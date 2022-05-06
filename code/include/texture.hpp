#ifndef TEXTURE_H
#define TEXTURE_H

//#include <opencv2/highgui/highgui.hpp>
//#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include "material.hpp"

class Texture {

public: 
    Texture() = default;

    Texture(Material* m) {
        material = m;
        use_bump_map = false;
        use_diffuse_map = false;
        use_specular_map = false;
        use_ambient_map = false;
        use_emission_map = false;
    }

    Material* material;  // the base material
    cv::Mat bump_map;
    cv::Mat diffuse_map;
    cv::Mat specular_map;
    cv::Mat ambient_map;
    cv::Mat emission_map;
    bool use_bump_map;
    bool use_diffuse_map;
    bool use_specular_map;
    bool use_ambient_map;
    bool use_emission_map;

    void read_bump_map(const char *filename) {
        bump_map = cv::imread(filename);
    }

    void read_diffuse_map(const char *filename) {
        diffuse_map = cv::imread(filename);
    }

    void read_specular_map(const char *filename) {
        specular_map = cv::imread(filename);
    }

    void read_ambient_map(const char *filename) {
        ambient_map = cv::imread(filename);
    }

    void read_emission_map(const char *filename) {
        emission_map = cv::imread(filename);
    }

    Vector3f get_bump(Vector2f uv) {
        float x, y;
        x = (1 - uv[1]) * bump_map.rows;
        y = uv[0] * bump_map.cols;
        int x0, y0; float x1, y1;
        x0 = (int) x; y0 = (int) y; x1 = x - x0; y1 = y - y0;
        cv::Scalar c[4];
        c[0] = bump_map.at<cv::Vec3b>(x0, y0);
        c[1] = bump_map.at<cv::Vec3b>(x0+1, y0);
        c[2] = bump_map.at<cv::Vec3b>(x0, y0+1);
        c[3] = bump_map.at<cv::Vec3b>(x0+1, y0+1);
        Vector3f color[4];
        color[0][0] = c[0](2) * 1.0 / 255; color[0][1] = c[0](1) * 1.0 / 255; color[0][2] = c[0](0)* 1.0 / 255; 
        color[1][0] = c[1](2) * 1.0 / 255; color[1][1] = c[1](1) * 1.0 / 255; color[1][2] = c[1](0)* 1.0 / 255; 
        color[2][0] = c[2](2) * 1.0 / 255; color[2][1] = c[2](1) * 1.0 / 255; color[2][2] = c[2](0)* 1.0 / 255; 
        color[3][0] = c[3](2) * 1.0 / 255; color[3][1] = c[3](1) * 1.0 / 255; color[3][2] = c[3](0)* 1.0 / 255; 
        color[0] = color[0] * (1 - x1) + color[1] * x1;
        color[2] = color[2] * (1 - x1) + color[3] * x1;
        color[0] = color[0] * (1 - y1) + color[2] * y1;
        Vector3f mean(0.5, 0.5, 0.5);
        color[0] = (color[0] - mean) * 2;
        return color[0];
    }

    Vector3f read_pixel(cv::Mat* map, Vector2f uv) {
        float x, y;
        x = (1 - uv[1]) * map->rows;
        y = uv[0] * map->cols;
        int x0, y0; float x1, y1;
        x0 = (int) x; y0 = (int) y; x1 = x - x0; y1 = y - y0;
        cv::Scalar c[4];
        c[0] = map->at<cv::Vec3b>(x0, y0);
        c[1] = map->at<cv::Vec3b>(x0+1, y0);
        c[2] = map->at<cv::Vec3b>(x0, y0+1);
        c[3] = map->at<cv::Vec3b>(x0+1, y0+1);
        Vector3f color[4];
        color[0][0] = c[0](2) * 1.0 / 255; color[0][1] = c[0](1) * 1.0 / 255; color[0][2] = c[0](0)* 1.0 / 255; 
        color[1][0] = c[1](2) * 1.0 / 255; color[1][1] = c[1](1) * 1.0 / 255; color[1][2] = c[1](0)* 1.0 / 255; 
        color[2][0] = c[2](2) * 1.0 / 255; color[2][1] = c[2](1) * 1.0 / 255; color[2][2] = c[2](0)* 1.0 / 255; 
        color[3][0] = c[3](2) * 1.0 / 255; color[3][1] = c[3](1) * 1.0 / 255; color[3][2] = c[3](0)* 1.0 / 255; 
        color[0] = color[0] * (1 - x1) + color[1] * x1;
        color[2] = color[2] * (1 - x1) + color[3] * x1;
        color[0] = color[0] * (1 - y1) + color[2] * y1;
        return color[0];
    }

    Material* get_material(Vector2f uv) {
        if (use_diffuse_map) {
            Vector3f color = read_pixel(&diffuse_map, uv);
            material->diffuseColor = color;
            material->ambientColor = color;
        }
        if (use_specular_map) {  // specular map (roughness map) in Phong model indicates shininess
            Vector3f color = read_pixel(&specular_map, uv);
            float color_mean = (color[0] + color[1] + color[2]) / 3;   
            //material->shininess = 20 / (color_mean + 0.1);
            material->w_reflection = 0.1 - 0.1 * color_mean;  // sometimes reflection can also be regarded as roughness
            material->w_diffusion = 0.9 + 0.1 * color_mean;
        }
        if (use_ambient_map) {
            Vector3f color = read_pixel(&ambient_map, uv);
            material->ambient_occlusion = (color[0] + color[1] + color[2]) / 3;
        }
        if (use_emission_map) {
            Vector3f color = read_pixel(&emission_map, uv);
            if (color[0] > 0.02 || color[1] > 0.02 || color[2] > 0.02) {
                material->emission = true;
                material->emissionColor = color * 2;
            } else {
                material->emission = false;
            }
        }
        return material;
    }


private: 

};

#endif