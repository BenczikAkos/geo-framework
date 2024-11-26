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
    return Vector(0.0); // TODO
}

void Tubular::movement(int selected, const Vector &pos) {
    // TODO
}

void Tubular::updateBaseMesh() {
    mesh.clear();
    std::vector<BaseMesh::VertexHandle> handles, tri;

    for(int i = 0; i < v_resolution; ++i) {
        double v = i / (double)v_resolution;
        for(int j = 0; j < u_resolution; ++j) {
            double u = j / (double)u_resolution;
            VectorVector c0_der, c1_der;
            auto vertex = evaluate(u, v, c0_der, c1_der);
            vertices.push_back(vertex);
            handles.push_back(mesh.add_vertex(vertex));
        }
    }
    for(int i = 0; i < u_resolution-1; ++i) {
        for(int j = 0; j < v_resolution-1; ++j) {
            tri.clear();
            tri.push_back(handles[i * u_resolution + j]);
            tri.push_back(handles[i * u_resolution + j + 1]);
            tri.push_back(handles[(i + 1) * u_resolution + j]);
            mesh.add_face(tri);
            tri.clear();
            tri.push_back(handles[(i + 1) * u_resolution + j]);
            tri.push_back(handles[i * u_resolution + j + 1]);
            tri.push_back(handles[(i + 1) * u_resolution + j + 1]);
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

double Tubular::F1(double v) const {
    return -2.0 * v * v * v + 3.0 * v * v;
}

double Tubular::G0(double v) const {
    return v * v * v - 2.0 * v * v + v;
}

double Tubular::G1(double v) const {
    return v * v * v - v * v;
}

Vector Tubular::ru(double u, double v) const {
    double h = 0.001;
    return Vector(0.0); // TODO
}

Vector Tubular::rv(double u, double v) const {
    double h = 0.001;
    return Vector(0.0); // TODO
}
#pragma endregion
