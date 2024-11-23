#include <fstream>
#include <cmath>
#include "bspline.hh"

BSpline::BSpline(std::string filename) : Object(filename) {
  reload();
}

BSpline::BSpline(size_t degree, size_t num_control_points, std::vector<double> &knotsInput, std::vector<Vector> &controlPointsInput) {
  this->degree = degree;
  this->num_control_points = num_control_points;
  std::cout<<this->degree<<"   "<<this->num_control_points<<std::endl;
  //making circular knot values
  knots.clear();
  knots.push_back(knotsInput[0]);
  for(size_t i = 0; i < degree; ++i){
    double circularKnot = knots[0] - (knotsInput[knotsInput.size() - 1 - i] - knotsInput[knotsInput.size() - 2 - i]);
    knots.insert(knots.begin(), circularKnot);
  }
  for(size_t i = 1; i < knotsInput.size(); ++i) {
    knots.push_back(knotsInput[i]);
  }
  for(size_t i = 0; i < degree; ++i){
    double circularKnot = knots[knots.size()-1] + (knotsInput[i+1] - knotsInput[i]);
    knots.push_back(circularKnot);
  }
  //circular control points
  control_points.clear();
  for(int i = degree-1; i >= 0; --i) {
    control_points.push_back(controlPointsInput[controlPointsInput.size() - 1 - i]);
  }
  for(size_t i = 0; i < controlPointsInput.size(); ++i){
    control_points.push_back(controlPointsInput[i]);
  }
  for(size_t i = 0; i <= degree; ++i){
    control_points.push_back(controlPointsInput[i]);
  }
}

BSpline::~BSpline() {}

void BSpline::draw(const Visualization &vis) const {
  Object::draw(vis);
  glDisable(GL_LIGHTING);
  glLineWidth(3.0);
  glColor3d(1.0, 1.0, 1.0);

  // Draw the curve
  glBegin(GL_LINE_STRIP);
  for (const auto &p : curve_points)
    glVertex3dv(p.data());
  glEnd();


  if (vis.show_control_points) {
    glDisable(GL_LIGHTING);
    glLineWidth(3.0);
    glColor3d(0.3, 0.3, 1.0);

    // Draw control polygon
    glBegin(GL_LINE_STRIP);
    for(int i = degree; i < degree+num_control_points; ++i){
      glVertex3dv(control_points[i].data());
    }
    glVertex3dv(control_points[degree].data()); // Close the polygon
    glEnd();

    glLineWidth(1.0);
    glPointSize(8.0);
    glColor3d(1.0, 0.0, 1.0);

    // Draw control points
    glBegin(GL_POINTS);
    for (const auto &p : control_points)
      glVertex3dv(p.data());
    glEnd();

    glPointSize(1.0);
    glEnable(GL_LIGHTING);
  }
}

void BSpline::drawWithNames(const Visualization &vis) const {
  if (!vis.show_control_points)
    return;

  for (size_t i = 0; i < control_points.size(); ++i) {
    const auto &p = control_points[i];
    glPushName(i);
    glRasterPos3dv(p.data());
    glPopName();
  }
}

Vector BSpline::postSelection(int selected) {
  return control_points[selected];
}

void BSpline::movement(int selected, const Vector &pos) {
  control_points[selected] = pos;
  int prev_iteration_CP = selected - num_control_points;
  if(0 <= prev_iteration_CP && prev_iteration_CP < control_points.size()){
    control_points[prev_iteration_CP] = pos;
  }
  int next_iteration_CP = selected + num_control_points;  
  if(0 <= next_iteration_CP && next_iteration_CP < control_points.size()){
    control_points[next_iteration_CP] = pos;
  }
  updateBaseMesh();
}

double BSpline::basisFunction(size_t i, size_t k, double t) const {
  if (k == 0) {
    return (t >= knots[i] && t < knots[i + 1]) ? 1.0 : 0.0;
  }
  double denom1 = knots[i + k] - knots[i];
  double denom2 = knots[i + k + 1] - knots[i + 1];
  double term1 = (denom1 > 0.0) ? (t - knots[i]) / denom1 * basisFunction(i, k - 1, t) : 0.0;
  double term2 = (denom2 > 0.0) ? (knots[i + k + 1] - t) / denom2 * basisFunction(i + 1, k - 1, t) : 0.0;
  return term1 + term2;
}

Vector BSpline::evaluateBSpline(double t) const {
  Vector result(0.0, 0.0, 0.0);
  for (size_t i = 0; i < control_points.size(); ++i) {
    result += control_points[i] * basisFunction(i, degree, t);
  }
  return result;
}

void BSpline::updateBaseMesh() {
  mesh.clear();
  curve_points.clear();
  size_t resolution = 100;
  for(size_t i = 0; i < resolution; ++i) {
    double t = (knots[knots.size()-1-degree] - knots[degree]) * i / resolution + knots[degree];
    auto bSplinePoint = evaluateBSpline(t);
    mesh.add_vertex(bSplinePoint);
    curve_points.push_back(bSplinePoint);
  }
  curve_points.push_back(curve_points[0]);
  Object::updateBaseMesh();
}


bool BSpline::reload() {
  try {
    std::ifstream file(filename);
    file.exceptions(std::ios::failbit | std::ios::badbit);

    file >> degree >> num_control_points;

    knots.clear();
    control_points.clear();

    // Read the given knot intervals and construct periodic knots
    std::vector<double> knot_intervals(num_control_points+1);
    for (double &interval : knot_intervals)
      file >> interval;
    knots.push_back(knot_intervals[0]);
    for(size_t i = 0; i < degree; ++i){
      double circularKnot = knots[0] - (knot_intervals[knot_intervals.size() - 1 - i] - knot_intervals[knot_intervals.size() - 2 - i]);
      knots.insert(knots.begin(), circularKnot);
    }
    for(size_t i = 1; i < knot_intervals.size(); ++i) {
      knots.push_back(knot_intervals[i]);
    }
    for(size_t i = 0; i < degree; ++i){
      double circularKnot = knots[knots.size()-1] + (knot_intervals[i+1] - knot_intervals[i]);
      knots.push_back(circularKnot);
    }
    // Read control points
    std::vector<Vector> file_control_points(num_control_points);
    control_points.clear();
    for (size_t i = 0; i < num_control_points; ++i) {
      file >> file_control_points[i][0] >> file_control_points[i][1] >> file_control_points[i][2];
    }
    for(int i = degree-1; i >= 0; --i) {
        control_points.push_back(file_control_points[file_control_points.size() - 1 - i]);
    }
    for(size_t i = 0; i < file_control_points.size(); ++i){
      control_points.push_back(file_control_points[i]);
    }
    for(size_t i = 0; i <= degree; ++i){
      control_points.push_back(file_control_points[i]);
    }

  } catch (std::ifstream::failure &) {
    return false;
  }

  updateBaseMesh();
  return true;
}

void BSpline::printKnots() const {
  std::cout<<"Printing knots "<<std::endl;
  for(const double& knot: knots){
    std::cout << knot<<std::endl;
  }
}

void BSpline::printControlPoints() const {
  std::cout<<"Printing CPs"<<std::endl;
  for(int i = 0; i < control_points.size(); ++i){
    std::cout<< i << " " << control_points[i] <<std::endl;
  }
  std::cout<<std::endl;
}