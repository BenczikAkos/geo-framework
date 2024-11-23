#include <fstream>
#include <cmath>
#include "tubular.hh"

Tubular::Tubular(std::string filename) : Object(filename) {
    reload();
}

Tubular::~Tubular() {}

void Tubular::draw(const Visualization &vis) const {
    Object::draw(vis);
    glDisable(GL_LIGHTING);
    glColor3d(1.0, 0.0, 1.0);
    glPointSize(3.0);
    glBegin(GL_POINTS);
    for (const auto &p : vertices)
      glVertex3dv(p.data());
    glEnd();

    glPointSize(1.0);
    glEnable(GL_LIGHTING);
}

void Tubular::drawWithNames(const Visualization &vis) const {
    // if (!vis.show_control_points)
    //     return;
    // c0.value().drawWithNames(vis);
    // c1.value().drawWithNames(vis);
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
            Vector B0 = c0.value().dEvaluateBSplineNormalizedInput(u)%c0.value().ddEvaluateBSplineNormalizedInput(u);
            Vector B1 = c1.value().dEvaluateBSplineNormalizedInput(u)%c1.value().ddEvaluateBSplineNormalizedInput(u);
            Vector vertex = c0.value().evaluateBSplineNormalizedInput(u) * F0(v)
                            + mu * B0 * G0(v)
                            + c1.value().evaluateBSplineNormalizedInput(u) * F1(v)
                            + mu * B1 * G1(v);
            vertices.push_back(vertex);
            handles.push_back(mesh.add_vertex(vertex));
        }
    }
    // for(const auto &v: vertices){
    //     std::cout<<v<<std::endl;
    // }
    for(int i = 0; i < u_resolution - 1; ++i) {
        for(int j = 0; j < v_resolution - 1; ++j) {
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
    for(int i = 0; i < num_of_control_points; ++i) {
        input>> control_points[i][0] >> control_points[i][1] >> control_points[i][2];
    }
    result.emplace(degree, num_of_control_points, knots, control_points);
}

Vector Tubular::normal(BaseMesh::VertexHandle vh) const {
    return Vector(0.0); // TODO
}

double Tubular::meanCurvature(BaseMesh::VertexHandle vh) const {
    return 0.0; // TODO
}

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
#pragma endregion
