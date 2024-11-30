#include <fstream>

#include "dupin.hh"
#include <cmath>

Dupin::Dupin(std::string filename) : Object(filename) {
  reload();
}

Dupin::~Dupin() {
}

void Dupin::draw(const Visualization &vis) const {
    Object::draw(vis);
    patch->draw(vis);
    if(vis.show_control_points){
        glDisable(GL_LIGHTING);
        glColor3d(0.3, 0.3, 1.0);
        glPointSize(18.0);
        glBegin(GL_POINTS);
        for(const auto &p : controlPoints){
          glVertex3dv(p.data());
        }
        glEnd();
        glPointSize(1.0);
        glEnable(GL_LIGHTING);
    }
}

void Dupin::drawWithNames(const Visualization &vis) const {
  if (!vis.show_control_points)
    return;
  int name = 0;
  for (auto v : controlPoints) {
    glPushName(name++);
    glRasterPos3dv(v.data());
    glPopName();
  }
}

Vector Dupin::postSelection(int selected) {
  return controlPoints[selected];
}

void Dupin::movement(int selected, const Vector &pos) {
  controlPoints[selected] = Vector(pos[0], controlPoints[selected][1], controlPoints[selected][2]);
  updateParameters();
}

void Dupin::updateParameters() {
    a = 0.25 * (controlPoints[0][0] + controlPoints[1][0] - controlPoints[2][0] - controlPoints[3][0]);
    d = 0.25 * (controlPoints[0][0] - controlPoints[1][0] + controlPoints[2][0] - controlPoints[3][0]);
    c = -0.25 * (controlPoints[0][0] - controlPoints[1][0] - controlPoints[2][0] + controlPoints[3][0]);
    b = std::sqrt(a*a - c*c);
    updateBaseMesh();
}

void Dupin::updateBaseMesh() {
    mesh.clear();
    std::vector<BaseMesh::VertexHandle> handles, tri;
    Vector translation = Vector(0.0f, 0.0f, 0.0f);
    if(controlPoints.size() == 4) {
        translation = 0.25 * (controlPoints[0] + controlPoints[1] + controlPoints[2] + controlPoints[3]);
    }
    for (size_t i = 0; i < resolution.first; ++i) {
        float u = range.first + (range.second - range.first) * i / resolution.first;
        for (size_t j = 0; j < resolution.second; ++j) {
            float v = range.first + (range.second - range.first) * j / resolution.second;
            Vector p(calculateX(u, v), calculateY(u, v), calculateZ(u, v));
            handles.push_back(mesh.add_vertex(p+translation));
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

    std::vector<float> x;
    x.push_back(a-c+d);
    x.push_back(a+c-d);
    x.push_back(-a+c+d);
    x.push_back(-a-c-d);
    for(auto i : x){
        controlPoints.push_back(Vector(i, 0.0f, 0.0f));
    }
    patch = new BezierDupinPatch(a, b, c, d, std::make_pair(0.1, 0.5), std::make_pair(2.3, 3.4));
    return true;
}

void Dupin::setA(float value) {
    a = value;
    c = std::sqrt(a*a - b*b);
    updateBaseMesh();
}

void Dupin::setB(float value) {
    b = value;
    c = std::sqrt(a*a - b*b);
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
