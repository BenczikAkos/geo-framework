#include <fstream>
#include <cmath>
#include "tubular.hh"

Tubular::Tubular(std::string filename) : Object(filename) {
    reload();
}

Tubular::~Tubular() {}

void Tubular::drawWithNames(const Visualization &vis) const {
    if (!vis.show_control_points)
        return;
    //c0.drawWithNames(vis);
    //c1.drawWithNames(vis);
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
            auto vertex = evaluate(u, v);
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

Vector Tubular::evaluate(double u, double v) const {
    Vector vertex(0.0);
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
    return evaluate(u-h, v) - evaluate(u+h, v);
}

Vector Tubular::rv(double u, double v) const {
    double h = 0.001;
    return evaluate(u, v-h) - evaluate(u, v+h);
}
#pragma endregion
