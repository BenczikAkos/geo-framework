#include <fstream>
#include <cmath>
#include "tubular.hh"

Tubular::Tubular(std::string filename) : Object(filename) {
    reload();
}

Tubular::~Tubular() {}

void Tubular::draw(const Visualization &vis) const {
    if (vis.show_control_points) {
        c0.draw();
        c1.draw();
    }
    Object::draw(vis);
}

void Tubular::drawWithNames(const Visualization &vis) const {
    if (!vis.show_control_points)
        return;
    for (size_t i = 0; i < c0.cp.size(); ++i) {
        const auto &p = c0.cp[i];
        glPushName(i);
        glRasterPos3dv(p.data());
        glPopName();
    }
    for (size_t i = 0; i < c1.cp.size(); ++i) {
        const auto &p = c1.cp[i];
        glPushName(i + c0.cp.size());
        glRasterPos3dv(p.data());
        glPopName();
    }
}

Vector Tubular::postSelection(int selected) {
    if(selected < c0.cp.size())
        return c0.cp[selected];
    selected -= c0.cp.size();
    return c1.cp[selected];
}

void Tubular::movement(int selected, const Vector &pos) {
    if(selected < c0.cp.size()){
        c0.cp[selected] = pos;
    }
    else {
        selected -= c0.cp.size();
        c1.cp[selected] = pos;
    }
    updateBaseMesh();
}

void Tubular::updateBaseMesh() {
    mesh.clear();
    mesh.request_vertex_normals();
    std::vector<BaseMesh::VertexHandle> handles, tri;

    for(int i = 0; i < v_resolution; ++i) {
        double v = i / (double)v_resolution;
        for(int j = 0; j < u_resolution; ++j) {
            double u = j / (double)u_resolution;
            VectorVector c0_der, c1_der;
            auto vertex = evaluate(u, v, c0_der, c1_der);
            vertices.push_back(vertex);
            auto v_handle = mesh.add_vertex(vertex);
            mesh.set_normal(v_handle, ru(u, v)%rv(u, v));
            handles.push_back(v_handle);
        }
        handles.push_back(handles[i * u_resolution + i]); // Close the loop
    }
    for(int i = 0; i < v_resolution-1; ++i){
        for(int j = 0; j < u_resolution; ++j){
            tri.clear();
            tri.push_back(handles[i * (u_resolution + 1) + j]);
            tri.push_back(handles[i * (u_resolution + 1) + j + 1]);
            tri.push_back(handles[(i + 1) * (u_resolution + 1) + j]);
            mesh.add_face(tri);
            tri.clear();
            tri.push_back(handles[(i + 1) * (u_resolution + 1) + j]);
            tri.push_back(handles[i * (u_resolution + 1) + j + 1]);
            tri.push_back(handles[(i + 1) * (u_resolution + 1) + j + 1]);
            mesh.add_face(tri);
        }
    }
    Object::updateBaseMesh(false, false); //TODO: exact normal and curvature calculator
}

bool Tubular::reload() {
    try{
        std::ifstream file(filename);
        file.exceptions(std::ios::failbit | std::ios::badbit);
        readBSpline(file, c0);
        readBSpline(file, c1);
    }
    catch (std::ifstream::failure &) {
        return false;
    }

    updateBaseMesh();
    return true;
}

void Tubular::updateMu(double newMu) {
    mu = newMu;
    updateBaseMesh();
}


void Tubular::readBSpline(std::ifstream &input, BSpline &result) {
    int degree, num_of_control_points;
    input >> degree >> num_of_control_points;
    result.p = degree;
    result.n = num_of_control_points - 1;
    for(int i = 0; i < num_of_control_points+1; ++i) {
        double knotValue = 0.0;
        input >> knotValue;
        result.knots.push_back(knotValue);
    }
    Vector cp;
    for(int i = 0; i < num_of_control_points; ++i) {
        input>> cp[0] >> cp[1] >> cp[2];
        result.cp.push_back(cp);
    }
}

Point Tubular::evaluate(double u, double v, VectorVector &c0_der, VectorVector &c1_der) const {
    c0_der.resize(3);
    c1_der.resize(3);
    u = std::lerp(c0.knots.front(), c0.knots.back(), u);
    Point C0_u = c0.derivatives(u, 2, c0_der);
    Point C1_u = c1.derivatives(u, 2, c1_der);
    Point B0 = c0_der[1]%c0_der[2];
    Point B1 = c1_der[1]%c1_der[2];
    Point vertex = C0_u * F0(v) + mu * B0 * G0(v) + C1_u * F1(v) + mu * B1 * G1(v);
    return vertex;
}

// Vector Tubular::normal(BaseMesh::VertexHandle vh) const {
//     return Vector(0.0); // TODO
// }

// double Tubular::meanCurvature(BaseMesh::VertexHandle vh) const {
//     return 0.0; // TODO
// }

#pragma region HelperFuncs
double Tubular::F0(double v) const {
    return 2.0 * v * v * v - 3.0 * v * v + 1.0;
}

double Tubular::dF0(double v) const {
    return 6.0 * v * v - 6.0 * v;
}

double Tubular::F1(double v) const {
    return -2.0 * v * v * v + 3.0 * v * v;
}

double Tubular::dF1(double v) const {
    return -6.0 * v * v + 6.0 * v;
}

double Tubular::G0(double v) const {
    return v * v * v - 2.0 * v * v + v;
}

double Tubular::dG0(double v) const {
    return 3.0 * v * v - 4.0 * v + 1.0;
}

double Tubular::G1(double v) const {
    return v * v * v - v * v;
}

double Tubular::dG1(double v) const {
    return 3.0 * v * v - 2.0 * v;
}
#pragma endregion

Vector Tubular::ru(double u, double v) const {
    VectorVector c0_der, c1_der;
    c0_der.resize(4); c0.derivatives(u, 3, c0_der);
    c1_der.resize(4); c1.derivatives(u, 3, c1_der);
    auto dB0 = c0_der[1]%c0_der[3];
    auto dB1 = c1_der[1]%c1_der[3];
    return c0_der[1] * F0(v) + mu * dB0 * G0(v) + c1_der[1] * F1(v) + mu * dB1 * G1(v);
}

Vector Tubular::rv(double u, double v) const {
    VectorVector c0_der, c1_der;
    c0_der.resize(3); 
    Point C0_u = c0.derivatives(u, 2, c0_der);
    c1_der.resize(3);
    Point C1_u =  c1.derivatives(u, 2, c1_der);
    auto B0 = c0_der[1]%c0_der[2];
    auto B1 = c1_der[1]%c1_der[2];
    return C0_u * dF0(v) + mu * B0 * dG0(v) + C1_u * dF1(v) + mu * B1 * dG1(v);
}
