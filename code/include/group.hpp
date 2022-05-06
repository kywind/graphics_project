#ifndef GROUP_H
#define GROUP_H


#include "object3d.hpp"
#include "ray.hpp"
#include "hit.hpp"
#include <iostream>
#include <vector>


// Group - add data structure to store a list of Object*
class Group : public Object3D {

public:

    Group() {
        std::vector <Object3D*> objects;
    }

    explicit Group (int num_objects) {
        for (int i=0; i<num_objects; i++){
            objects.push_back(nullptr);
        }
    }

    ~Group() override {
        int num_objects = getGroupSize();
        for (int i=0; i<num_objects; i++){
            delete objects[i];
        }
    }

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        int num_objects = getGroupSize();
        bool flag = false;
        for (int i=0; i<num_objects; i++){
            if (objects[i]->intersect(r, h, tmin)) flag = true;
        }
        return flag;
    }

    void addObject(int index, Object3D *obj) {
        objects[index] = obj;
    }

    int getGroupSize() {
        return (int) objects.size();
    }

private:
    std::vector <Object3D*> objects;
};

#endif
	
