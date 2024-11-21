#include <fstream>
#include <cmath>
#include "bspline.hh"

BSpline::BSpline(std::string filename) : Object(filename) {
  reload();
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
    for (const auto &p : control_points)
      glVertex3dv(p.data());
    glVertex3dv(control_points[0].data()); // Close the polygon
    glEnd();

    glLineWidth(1.0);
    glPointSize(4.0);
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
    std::cout<<"t: " <<t << " basisF: " << basisFunction(i, degree, t) <<std::endl;
    result += control_points[i] * basisFunction(i, degree, t);
  }
  return result;
}

void BSpline::updateBaseMesh() {
  mesh.clear();
  curve_points.clear();
  size_t resolution = 100;
  for(size_t i = 0; i < resolution; ++i) {
    double t = (knots[knots.size()-1] - knots[0]) * i / resolution + knots[0];
    auto bSplinePoint = evaluateBSpline(t);
    mesh.add_vertex(bSplinePoint);
    curve_points.push_back(bSplinePoint);
    std::cout<<"value: "<<t<<" "<<bSplinePoint<<std::endl;
  }
  Object::updateBaseMesh();
}


bool BSpline::reload() {
  try {
    std::ifstream file(filename);
    file.exceptions(std::ios::failbit | std::ios::badbit);

    size_t num_control_points;
    file >> degree >> num_control_points;

    knots.clear();
    control_points.clear();

    // Read the given knot intervals and construct periodic knots
    std::vector<double> knot_intervals(num_control_points-1);
    for (double &interval : knot_intervals)
      file >> interval;
    knots.push_back(knot_intervals[0]);
    knots.push_back(knot_intervals[0]);
    for(size_t i = 0; i < num_control_points-1; ++i) {
      knots.push_back(knot_intervals[i]);
    }
    knots.push_back(knot_intervals[knot_intervals.size() - 1]);
    knots.push_back(knot_intervals[knot_intervals.size() - 1]);
    for(size_t i = 0; i < knots.size(); ++i) {
      std::cout<<"knots: "<<knots[i]<<std::endl;
    }
    // for(size_t i = 0; i < degree*2; ++i) {
    //     knots.push_back(knots[knots.size() - 1] + (knot_intervals[i + 1] - knot_intervals[i]));
    // }
    // Read control points
    control_points.resize(num_control_points);
    for (size_t i = 0; i < num_control_points; ++i) {
      file >> control_points[i][0] >> control_points[i][1] >> control_points[i][2];
    }
    // for(size_t i = 0; i < degree*2; ++i) {
    //     control_points.push_back(control_points[i]);
    // }

  } catch (std::ifstream::failure &) {
    return false;
  }

  updateBaseMesh();
  return true;
}


