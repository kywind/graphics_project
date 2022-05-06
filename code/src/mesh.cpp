#include "mesh.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <utility>
#include <sstream>


bool Mesh::intersect(const Ray &ray, Hit &hit, float tmin) {
    //if ((maxPos[0] - minPos[0]) * (maxPos[1] - minPos[1]) * (maxPos[2] - minPos[2]) <= 0) return false;
    Vector3f rayDir = ray.getDirection().normalized();
    Vector3f rayOri = ray.getOrigin();
    float t_enter_x, t_enter_y, t_enter_z, t_exit_x, t_exit_y, t_exit_z;

    t_enter_x = (minPos[0] - rayOri[0]) / rayDir[0];
    t_enter_y = (minPos[1] - rayOri[1]) / rayDir[1];
    t_enter_z = (minPos[2] - rayOri[2]) / rayDir[2];

    t_exit_x = (maxPos[0] - rayOri[0]) / rayDir[0];
    t_exit_y = (maxPos[1] - rayOri[1]) / rayDir[1];
    t_exit_z = (maxPos[2] - rayOri[2]) / rayDir[2];

    //if (!isfinite(t_enter_x) || !isfinite(t_exit_x) || isnan(t_enter_x) || isnan(t_exit_x)) std::cout << "error" << std::endl;
    //if (!isfinite(t_enter_y) || !isfinite(t_exit_y) || isnan(t_enter_y) || isnan(t_exit_y)) std::cout << "error" << std::endl;
    //if (!isfinite(t_enter_z) || !isfinite(t_exit_z) || isnan(t_enter_z) || isnan(t_exit_z)) std::cout << "error" << std::endl;
    //if (t_enter_x == t_exit_x) std::cout << "error" << std::endl;
    //if (t_enter_y == t_exit_y) std::cout << "error" << std::endl;
    //if (t_enter_z == t_exit_z) std::cout << "error" << std::endl;

    bool flip = false;
    if (t_enter_x > t_exit_x) {
        if (divide_dir == 0) flip = true;
        t_exit_x = t_enter_x + t_exit_x;
        t_enter_x = t_exit_x - t_enter_x;
        t_exit_x = t_exit_x - t_enter_x;
    }
    if (t_enter_y > t_exit_y) {
        if (divide_dir == 1) flip = true;
        t_exit_y = t_enter_y + t_exit_y;
        t_enter_y = t_exit_y - t_enter_y;
        t_exit_y = t_exit_y - t_enter_y;
    }
    if (t_enter_z > t_exit_z) {
        if (divide_dir == 2) flip = true;
        t_exit_z = t_enter_z + t_exit_z;
        t_enter_z = t_exit_z - t_enter_z;
        t_exit_z = t_exit_z - t_enter_z;
    }

    t_enter_x = max(max(t_enter_x, t_enter_y), t_enter_z);
    t_exit_x = min(min(t_exit_x, t_exit_y), t_exit_z);

    
    if (t_enter_x > t_exit_x || t_exit_x < 0) return false;

    // have intersection with bounding box
    if (is_leaf) {
        bool result = false;
        for (int ti = 0; ti < (int) t_id.size(); ti++) {
            MeshTriangle& tri = (*t)[t_id[ti]];
            if (tri.intersect(ray, hit, tmin)) result = true;
        }
        //if (hit.getT() < t_enter_x && hit.getT() > t_exit_x) result = false;
        return result;
    }
    
    //if (flip) {
    //    return right -> intersect(ray, hit, tmin) || left->intersect(ray, hit, tmin);
    //}
    //else {
    //    return left -> intersect(ray, hit, tmin) || right->intersect(ray, hit, tmin);
    //}   
    
    bool ans = left -> intersect(ray, hit, tmin);
    ans = right -> intersect(ray, hit, tmin) || ans;
    return ans;
}


Mesh::Mesh(const char *filename, const bool interpolate, Material *material) : Object3D(material) {

    this->interpolate = interpolate;
    n = new std::vector<Vector3f>;
    tex = new std::vector<Vector2f>;
    v = new std::vector<Vector3f>;
    t = new std::vector<MeshTriangle>;

    std::ifstream f;
    f.open(filename);
    if (!f.is_open()) {
        std::cout << "Cannot open " << filename << "\n";
        return;
    }
    std::string line;
    std::string vTok("v");
    std::string fTok("f");
    std::string texTok("vt");
    std::string nTok("vn");
    char bslash = '/', space = ' ';
    std::string tok;
    bool mode = false;

    while (true) {
        std::getline(f, line);
        if (f.eof()) break;
        if (line.size() < 3) continue;
        if (line.at(0) == '#') continue;
        if (line.at(0) == 'o') continue;
        if (line.at(0) == 's') continue;
        std::stringstream ss(line);
        ss >> tok;

        if (tok == vTok) {
            Vector3f vec;
            ss >> vec[0] >> vec[1] >> vec[2];
            v->push_back(vec);
        } else if (tok == fTok) {
            if (line.find(bslash) != std::string::npos) {
                mode = true;
                bool is_triangle_mesh = true;
                string s[4];
                std::stringstream facess(line);
                facess >> tok >> s[0] >> s[1] >> s[2];
                if (facess >> s[3]) is_triangle_mesh = false;
                int a[4], b[4], c[4];

                for (int index = 0; index < 4; index ++) {
                    int k = 0;
                    a[index] = 0;
                    b[index] = 0;
                    c[index] = 0;
                    for (; s[index][k] != '/'; k++) {
			        	a[index] = a[index] * 10 + (s[index][k] - 48);
			        }
                    k++;
			        for (; s[index][k] != '/'; k++) {
			        	b[index] = b[index] * 10 + (s[index][k] - 48);
			        }
                    k++;
			        for (; k < s[index].length(); k++) {
			        	c[index] = c[index] * 10 + (s[index][k] - 48);
			        }
                }

                MeshTriangle trig1(material);
                trig1.vertex_id[0] = a[0]-1; trig1.vertex_id[1] = a[1]-1; trig1.vertex_id[2] = a[2]-1;
                trig1.uv_id[0] = b[0]-1;     trig1.uv_id[1] = b[1]-1;     trig1.uv_id[2] = b[2]-1;
                trig1.norm_id[0] = c[0]-1;   trig1.norm_id[1] = c[1]-1;   trig1.norm_id[2] = c[2]-1;
                t->push_back(trig1);

                if (!is_triangle_mesh) {
                    MeshTriangle trig2(material);
                    trig2.vertex_id[0] = a[0]-1; trig2.vertex_id[1] = a[2]-1; trig2.vertex_id[2] = a[3]-1;
                    trig2.uv_id[0] = b[0]-1;     trig2.uv_id[1] = b[2]-1;     trig2.uv_id[2] = b[3]-1;
                    trig2.norm_id[0] = c[0]-1;   trig2.norm_id[1] = c[2]-1;   trig2.norm_id[2] = c[3]-1;
                    t->push_back(trig2);
                }
            } else {
                assert (!mode);
                bool is_triangle_mesh = true;
                string s[4];
                std::stringstream facess(line);
                facess >> tok >> s[0] >> s[1] >> s[2];
                if (facess >> s[3]) is_triangle_mesh = false;

                int a[4];
                for (int index = 0; index < 4; index ++) {
                    a[index] = 0;
                    for (int k=0; k < s[index].length(); k++) {
			        	a[index] = a[index] * 10 + (s[index][k] - 48);
			        }
                }
                MeshTriangle trig1(material);
                trig1.vertex_id[0] = a[0]-1; trig1.vertex_id[1] = a[1]-1; trig1.vertex_id[2] = a[2]-1;
                t->push_back(trig1);

                if (!is_triangle_mesh) {
                    MeshTriangle trig2(material);
                    trig2.vertex_id[0] = a[0]-1; trig2.vertex_id[1] = a[2]-1; trig2.vertex_id[2] = a[3]-1;
                    t->push_back(trig2);
                }
            }
        } else if (tok == texTok) {
            assert (mode);
            Vector2f texcoord;
            ss >> texcoord[0];
            ss >> texcoord[1];
            tex->push_back(texcoord);
        }
        else if (tok == nTok) {
            assert (mode);
            Vector3f norm;
            ss >> norm[0] >> norm[1] >> norm[2];
            n->push_back(norm);
        }
    }
    f.close();
    
    // this is the root node, need to compute minPos and maxPos
    for (int vi=0; vi < (int) v->size(); vi++) {
        v_id.push_back(vi);
        minPos[0] = min((*v)[vi][0], minPos[0]);
        minPos[1] = min((*v)[vi][1], minPos[1]);
        minPos[2] = min((*v)[vi][2], minPos[2]);
        maxPos[0] = max((*v)[vi][0], maxPos[0]);
        maxPos[1] = max((*v)[vi][1], maxPos[1]);
        maxPos[2] = max((*v)[vi][2], maxPos[2]);
    }
    minPos -= Vector3f(1e-3);
    maxPos += Vector3f(1e-3);

    divide_dir = 0;
    if (mode == 0) computeNormal(); // save norm data into t and n
    
    for (int triId = 0; triId < (int) t->size(); ++triId) {
        t_id.push_back(triId);
        (*t)[triId].initiate(v, tex, n, interpolate);
    }
    computeKdTree(0);
}


void Mesh::computeNormal() {
    if (!interpolate) {
        for (int triId = 0; triId < (int) t->size(); ++triId) {
            MeshTriangle& tri = (*t)[triId];
            Vector3f a = (*v)[tri.vertex_id[1]] - (*v)[tri.vertex_id[0]];
            Vector3f b = (*v)[tri.vertex_id[2]] - (*v)[tri.vertex_id[0]];
            b = Vector3f::cross(a, b);
            n->push_back(b);
            (*t)[triId].norm_id[0] = triId;
            (*t)[triId].norm_id[1] = triId;
            (*t)[triId].norm_id[2] = triId;
        }
    } else {
        for (int vi = 0; vi < (int) v->size(); ++vi) {
            n->push_back(Vector3f::ZERO);
        }
        for (int triId = 0; triId < (int) t->size(); ++triId) {
            MeshTriangle& tri = (*t)[triId];
            Vector3f a = (*v)[tri.vertex_id[1]] - (*v)[tri.vertex_id[0]];
            Vector3f b = (*v)[tri.vertex_id[2]] - (*v)[tri.vertex_id[0]];
            b = Vector3f::cross(a, b);
            (*n)[tri.vertex_id[0]] += b;
            (*n)[tri.vertex_id[1]] += b;
            (*n)[tri.vertex_id[2]] += b;
            (*t)[triId].norm_id[0] = tri.vertex_id[0];
            (*t)[triId].norm_id[1] = tri.vertex_id[1];
            (*t)[triId].norm_id[2] = tri.vertex_id[2];
            
        }
        for (int vi = 0; vi < (int) v->size(); ++vi) {
            (*n)[vi].normalize();
        }
    }
    
}


void Mesh::computeKdTree(int stop_iters) {

    std::vector<Vector3f> v_temp;
    for (int vi=0; vi < (int) v_id.size(); vi++) v_temp.push_back((*v)[v_id[vi]]);
    bool stop = ((int) v_temp.size() < 3) || stop_iters >= 3;

    if (stop) is_leaf = true;

    else {
        is_leaf = false;
        //  get left_max and right_min using the median of triangles (v and t)
        Vector3f left_max = maxPos;
        Vector3f right_min = minPos;

        if (divide_dir == 0) sort(v_temp.begin(), v_temp.end(), v_sort_x);
        else if (divide_dir == 1) sort(v_temp.begin(), v_temp.end(), v_sort_y);
        else sort(v_temp.begin(), v_temp.end(), v_sort_z);

        int query = (int) v_temp.size() / 2;
        float new_mid = v_temp[query][divide_dir];
        left_max[divide_dir] = new_mid+1e-3;
        right_min[divide_dir] = new_mid-1e-3;

        std::vector<int> v_id_l, v_id_r;
        std::vector<int> t_id_l, t_id_r;

        // test intersection of every triangle with the box of left and right
        int left_count = 0, right_count = 0;
        for (int ti = 0; ti < (int) t_id.size(); ti++) {
            MeshTriangle& tri = (*t)[t_id[ti]];
            bool inter_left = testIntersection((*v)[tri.vertex_id[0]], (*v)[tri.vertex_id[1]], (*v)[tri.vertex_id[2]], minPos, left_max);
            bool inter_right = testIntersection((*v)[tri.vertex_id[0]], (*v)[tri.vertex_id[1]], (*v)[tri.vertex_id[2]], right_min, maxPos);
            if (inter_left) {
                left_count ++;
                t_id_l.push_back(t_id[ti]);
            }
            if (inter_right) {
                right_count ++;
                t_id_r.push_back(t_id[ti]);
            }
        }

        for (int vi=0; vi < (int) v_id.size(); vi++) {
            if ((*v)[v_id[vi]][divide_dir] < new_mid) v_id_l.push_back(v_id[vi]);
            else v_id_r.push_back(v_id[vi]);
        }

        if (max(t_id_l.size(), t_id_r.size()) == t_id.size()) stop_iters ++;
        else stop_iters = 0;

        left = new Mesh(material);
        right = new Mesh(material);

        left->t = this->t;
        left->v = this->v;
        right->t = this->t;
        right->v = this->v;

        left->t_id = t_id_l;
        left->v_id = v_id_l;
        left->divide_dir = (divide_dir + 1) % 3;
        left->minPos = this->minPos;
        left->maxPos = left_max;

        right->t_id = t_id_r;
        right->v_id = v_id_r;
        right->divide_dir = (divide_dir + 1) % 3;
        right->minPos = right_min;
        right->maxPos = this->maxPos;

        left->computeKdTree(stop_iters);
        right->computeKdTree(stop_iters);

    }

}
