#include <fstream>

#include "dupin.hh"
#include <cmath>

Dupin::Dupin(std::string filename) : Object(filename) {
  reload();
}

Dupin::~Dupin() {
}

void Dupin::drawWithNames(const Visualization &vis) const {
  if (!vis.show_wireframe)
    return;
  for (auto v : mesh.vertices()) {
    glPushName(v.idx());
    glRasterPos3dv(mesh.point(v).data());
    glPopName();
  }
}

Vector Dupin::postSelection(int selected) {
  return mesh.point(BaseMesh::VertexHandle(selected));
}

void Dupin::movement(int selected, const Vector &pos) {
  mesh.set_point(BaseMesh::VertexHandle(selected), pos);
}

void Dupin::updateBaseMesh() {
    mesh.clear();
    std::vector<BaseMesh::VertexHandle> handles, tri;
    for (size_t i = 0; i < resolution.first; ++i) {
        float u = range.first + (range.second - range.first) * i / resolution.first;
        for (size_t j = 0; j < resolution.second; ++j) {
            float v = range.first + (range.second - range.first) * j / resolution.second;
            Vector p(calculateX(u, v), calculateY(u, v), calculateZ(u, v));
            handles.push_back(mesh.add_vertex(p));
        }
        // Adding the starter point to the end
        handles.push_back(handles[i * (resolution.second + 1)]);
    }
    //Adding the starter circle to the end
    for(size_t i = 0; i < resolution.second + 1; ++i) {
        handles.push_back(handles[i]);
    }
    for (size_t i = 0; i < resolution.first; ++i) {
        for (size_t j = 0; j < resolution.second; ++j) {
            tri.clear();
            tri.push_back(handles[i * (resolution.second + 1) + j]);
            tri.push_back(handles[(i + 1) * (resolution.second + 1) + j + 1]);
            tri.push_back(handles[i * (resolution.second + 1) + j + 1]);
            mesh.add_face(tri);
            tri.clear();
            tri.push_back(handles[i * (resolution.second + 1) + j]);
            tri.push_back(handles[(i + 1) * (resolution.second + 1) + j]);
            tri.push_back(handles[(i + 1) * (resolution.second + 1) + j + 1]);
            mesh.add_face(tri);
        }
    }
    Object::updateBaseMesh(false, false);
}

bool Dupin::reload() {
    float read_a, read_b, read_c, read_d = 0.0f;
    try {
        std::ifstream f(filename);
        f.exceptions(std::ios::failbit | std::ios::badbit);
        f >> read_a >> read_b >> read_c >> read_d;
        a = read_a; b = read_b; c = read_c; d = read_d;

    } catch (std::ios_base::failure &e) {
        return false;
    }
    updateBaseMesh();
    return true;
}

void Dupin::setA(float value) {
    a = value;
    updateBaseMesh();
}

void Dupin::setB(float value) {
    b = value;
    updateBaseMesh();
}

void Dupin::setC(float value) {
    c = value;
    updateBaseMesh();
}

void Dupin::setD(float value) {
    d = value;
    updateBaseMesh();
}

float Dupin::calculateX(float u, float v) const {
    float numerator = d * (c - a * std::cos(u) * std::cos(v)) +  b * b * std::cos(u);
    float denominator = a - c * std::cos(u) * std::cos(v);
    return numerator / denominator;
}

float Dupin::calculateY(float u, float v) const {
    float numerator = b * std::sin(u) * (a - d * std::cos(v));
    float denominator = a - c * std::cos(u) * std::cos(v);
    return numerator / denominator;
}

float Dupin::calculateZ(float u, float v) const {
    float numerator = b * std::sin(v) * (c * std::cos(u) - d);
    float denominator = a - c * std::cos(u) * std::cos(v);
    return numerator / denominator;
}
