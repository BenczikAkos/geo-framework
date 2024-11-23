#include <fstream>
#include <cmath>
#include "tubular.hh"

Tubular::Tubular(std::string filename) : Object(filename) {
    reload();
}

Tubular::~Tubular() {}

void Tubular::draw(const Visualization &vis) const {
    Object::draw(vis);
}

void Tubular::drawWithNames(const Visualization &vis) const {
    if (!vis.show_control_points)
        return;
    c0.value().drawWithNames(vis);
    c1.value().drawWithNames(vis);
}

Vector Tubular::postSelection(int selected) {
    return Vector(0.0); // TODO
}

void Tubular::movement(int selected, const Vector &pos) {
    // TODO
}

void Tubular::updateBaseMesh() {
    mesh.clear();

    // Object::updateBaseMesh(true, true); TODO: exact normal and curvature calculator
}

bool Tubular::reload() {
    try{
        std::ifstream file(filename);
        file.exceptions(std::ios::failbit | std::ios::badbit);
        readBSpline(file, c0);
        readBSpline(file, c1);
        std::cout<<"C0"<<std::endl;
        c0.value().printControlPoints();
        c0.value().printKnots();
        std::cout<<"C1"<<std::endl;
        c1.value().printControlPoints();
        c1.value().printKnots();

    }
    catch (std::ifstream::failure &){
        return false;
    }

    updateBaseMesh();
    return true;
}

void Tubular::readBSpline(std::ifstream &input, std::optional<BSpline> &result) {
    int degree, num_of_control_points;
    input >> degree >> num_of_control_points;
    std::vector<double> knots;
    for(int i = 0; i < num_of_control_points+1; ++i) {
        double knotValue = 0.0;
        input >> knotValue;
        knots.push_back(knotValue);
    }
    std::vector<Vector> control_points(num_of_control_points);
    for(int i = 0; i < num_of_control_points; ++i){
        input>> control_points[i][0] >> control_points[i][1] >> control_points[i][2];
    }
    result.emplace(degree, num_of_control_points, knots, control_points);
}

#pragma region HelperFuncs
Vector Tubular::normal(BaseMesh::VertexHandle vh) const {
    return Vector(0.0); // TODO
}

double Tubular::meanCurvature(BaseMesh::VertexHandle vh) const {
    return 0.0; // TODO
}

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
#pragma endregion
