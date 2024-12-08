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
            u = std::lerp(c0.knots.front(), c0.knots.back(), u);
            VectorVector c0_der, c1_der;
            auto vertex = evaluate(u, v, c0_der, c1_der);
            vertices.push_back(vertex);
            auto v_handle = mesh.add_vertex(vertex);
            mesh.set_normal(v_handle, ru(u, v)%rv(u, v));
            mesh.data(v_handle).mean = calculateMeanCurvature(u,v); //TODO: fix
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
    Object::updateBaseMesh(true, true); //TODO: exact normal and curvature calculator
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
    Point C0_u = c0.derivatives(u, 2, c0_der);
    Point C1_u = c1.derivatives(u, 2, c1_der);
    Point B0 = c0_der[1]%c0_der[2];
    Point B1 = c1_der[1]%c1_der[2];
    Point vertex = C0_u * F0(v) + mu * B0 * G0(v) + C1_u * F1(v) + mu * B1 * G1(v);
    return vertex;
}

Vector Tubular::normal(BaseMesh::VertexHandle vh) const {
    return mesh.normal(vh).normalized();
}

double Tubular::meanCurvature(BaseMesh::VertexHandle vh) const {
    return mesh.data(vh).mean;
}

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
    //u = std::lerp(c0.knots.front(), c0.knots.back(), u);
    c0_der.resize(4); c0.derivatives(u, 3, c0_der);
    c1_der.resize(4); c1.derivatives(u, 3, c1_der);
    Vector dB0 = c0_der[1]%c0_der[3]; // + c0_der[2]%c0_der[2] = 0
    Vector dB1 = c1_der[1]%c1_der[3]; // + c1_der[2]%c1_der[2] = 0
    auto result = c0_der[1] * F0(v) + mu * dB0 * G0(v) + c1_der[1] * F1(v) + mu * dB1 * G1(v);
    return c0_der[1] * F0(v) + mu * dB0 * G0(v) + c1_der[1] * F1(v) + mu * dB1 * G1(v);
}

Vector Tubular::rv(double u, double v) const {
    VectorVector c0_der, c1_der;
    c0_der.resize(3); 
    Point C0_u = c0.derivatives(u, 2, c0_der);
    c1_der.resize(3);
    Point C1_u =  c1.derivatives(u, 2, c1_der);
    Vector B0 = c0_der[1]%c0_der[2];
    Vector B1 = c1_der[1]%c1_der[2];
    return C0_u * dF0(v) + mu * B0 * dG0(v) + C1_u * dF1(v) + mu * B1 * dG1(v);
}

Vector Tubular::ruu(double u, double v) const {
    VectorVector c0_der, c1_der;
    c0_der.resize(5); c0.derivatives(u, 4, c0_der);
    c1_der.resize(5); c1.derivatives(u, 4, c1_der);
    auto d2B0 = c0_der[1]%c0_der[4] + c0_der[2]%c0_der[3];
    auto d2B1 = c1_der[1]%c1_der[4] + c1_der[2]%c1_der[3];
    return c0_der[2] * F0(v) + mu * d2B0 * G0(v) + c1_der[2] * F1(v) + mu * d2B1 * G1(v);
}

Vector Tubular::rvv(double u, double v) const {
    VectorVector c0_der, c1_der;
    c0_der.resize(3); 
    Point C0_u = c0.derivatives(u, 2, c0_der);
    c1_der.resize(3);
    Point C1_u =  c1.derivatives(u, 2, c1_der);
    auto B0 = c0_der[1]%c0_der[2];
    auto B1 = c1_der[1]%c1_der[2];
    return C0_u * (12.0 * v - 6.0)
        + mu * B0 * (6.0 * v - 4.0) 
        + C1_u * (-12.0 * v - 6.0) 
        + mu * B1 * (6.0 * v - 2.0);
}

Vector Tubular::ruv(double u, double v) const {
    VectorVector c0_der, c1_der;
    c0_der.resize(4); c0.derivatives(u, 3, c0_der);
    c1_der.resize(4); c1.derivatives(u, 3, c1_der);
    auto dB0 = c0_der[1]%c0_der[3];
    auto dB1 = c1_der[1]%c1_der[3];
    return c0_der[1] * dF0(v) + mu * dB0 * dG0(v) + c1_der[1] * dF1(v) + mu * dB1 * dG1(v);
}

double Tubular::calculateMeanCurvature(double u, double v) const {
    double E = ru(u, v) | ru(u, v);
    double F = ru(u, v) | rv(u, v);
    double G = rv(u, v) | rv(u, v);
    double L = ruu(u, v) | (ru(u,v) % rv(u, v));
    double M = ruv(u, v) | (ru(u,v) % rv(u, v));
    double N = rvv(u, v) | (ru(u,v) % rv(u, v));
    return (N * E - 2.0 * M * F + L * G) / (2.0 * (E * G - F * F));
}

