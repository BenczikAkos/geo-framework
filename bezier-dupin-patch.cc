#include "bezier-dupin-patch.hh"

BezierDupinPatch::BezierDupinPatch(float a, float b, float c, float d,
                                  std::pair<double, double> uRange,
                                  std::pair<double, double> vRange,
                                  std::vector<Vector> dupin_control_points) : Object("") {
    control_points.resize(9);
    weights.resize(9);
    std::vector<double> g(3);
    std::vector<double> h(3);
    std::vector<double> G(3);
    std::vector<double> H(3);
    g[0] = std::tan(0.5 * uRange.first);
    g[2] = std::tan(0.5 * uRange.second);
    g[1] = 0.5 * (g[0] + g[2]);
    h[0] = std::tan(0.5 * vRange.first);
    h[2] = std::tan(0.5 * vRange.second);
    h[1] = 0.5 * (h[0] + h[2]);
    G[0] = g[0] * g[0];
    G[1] = g[0] * g[2];
    G[2] = g[2] * g[2];
    H[0] = h[0] * h[0];
    H[1] = h[0] * h[2];
    H[2] = h[2] * h[2];
    Vector translation = Vector(0.0f, 0.0f, 0.0f);
    if(dupin_control_points.size() == 4) {
        translation = 0.25 * (dupin_control_points[0] + dupin_control_points[1] + dupin_control_points[2] + dupin_control_points[3]);
    }
    for(size_t i = 0; i < 3; ++i){
        for(size_t j = 0; j < 3; ++j){
            size_t index = i * 3 + j;
            weights[index] = a*(1.0+G[i])*(1.0+H[j]) - c*(1.0-G[i])*(1.0-H[j]);
            double CPdenominator = 1.0 / (a - c * ((1.0-G[i])/(1.0+G[i])) * ((1.0-H[j])/(1.0+H[j])));
            double CPx = b*b*(1.0-G[i])/(1.0+G[i]) - d*a*(1.0-G[i])/(1.0+G[i])*(1.0-H[j])/(1.0+H[j]) + d*c;
            double CPy = 2.0*b* (g[i]/(1.0+G[i])) * (a-d*(1.0-H[j])/(1.0+H[j]));
            double CPz = 2.0*b*(c*((1.0-G[i])/(1.0+G[i])) - d) * (h[j]/(1+H[j]));
            control_points[index] = CPdenominator * Vector(CPx, CPy, CPz) + translation;
        }
    }
    updateBaseMesh();
}


BezierDupinPatch::~BezierDupinPatch() {
}

void BezierDupinPatch::draw(const Visualization &vis) const {
  Object::draw(vis);
  if (vis.show_control_points) {
    glDisable(GL_LIGHTING);
    glLineWidth(3.0);
    glColor3d(0.3, 0.3, 1.0);
    for (size_t k = 0; k < 2; ++k)
      for (size_t i = 0; i <= 2; ++i) {
        glBegin(GL_LINE_STRIP);
        for (size_t j = 0; j <= 2; ++j) {
          size_t const index = k ? j * 3 + i : i * 3 + j;
          const auto &p = control_points[index];
          glVertex3dv(p.data());
        }
        glEnd();
      }
    glLineWidth(1.0);
    glPointSize(8.0);
    glColor3d(1.0, 0.0, 1.0);
    glBegin(GL_POINTS);
    for (const auto &p : control_points)
      glVertex3dv(p.data());
    glEnd();
    glPointSize(1.0);
    glEnable(GL_LIGHTING);
  }
}

void BezierDupinPatch::drawWithNames(const Visualization &vis) const {
  if (!vis.show_control_points)
    return;
  for (size_t i = 0; i < control_points.size(); ++i) {
    const auto &p = control_points[i];
    glPushName(i);
    glRasterPos3dv(p.data());
    glPopName();
  }
}

Vector BezierDupinPatch::postSelection(int selected) {
  return control_points[selected];
}

void BezierDupinPatch::movement(int selected, const Vector &pos) {
  control_points[selected] = pos;
}

static void bernstein(size_t n, double u, std::vector<double> &coeff) {
  coeff.clear(); coeff.reserve(n + 1);
  coeff.push_back(1.0);
  double u1 = 1.0 - u;
  for (size_t j = 1; j <= n; ++j) {
    double saved = 0.0;
    for (size_t k = 0; k < j; ++k) {
      double tmp = coeff[k];
      coeff[k] = saved + tmp * u1;
      saved = tmp * u;
    }
    coeff.push_back(saved);
  }
}

void BezierDupinPatch::updateBaseMesh() {
  size_t resolution = 50;

  mesh.clear();
  std::vector<BaseMesh::VertexHandle> handles, tri;

  std::vector<double> coeff_u, coeff_v;
  for (size_t i = 0; i < resolution; ++i) {
    double u = (double)i / (double)(resolution - 1);
    bernstein(2, u, coeff_u);
    for (size_t j = 0; j < resolution; ++j) {
      double v = (double)j / (double)(resolution - 1);
      bernstein(2, v, coeff_v);
      Vector p(0.0, 0.0, 0.0);
      double denominator = 0.0;
      for (size_t k = 0, index = 0; k <= 2; ++k) {
        for (size_t l = 0; l <= 2; ++l, ++index) {
            denominator += weights[index] * coeff_u[k] * coeff_v[l];
            p += control_points[index] * weights[index] * coeff_u[k] * coeff_v[l];
        }
      }
      handles.push_back(mesh.add_vertex(p/denominator));
    }
  }
  for (size_t i = 0; i < resolution - 1; ++i)
    for (size_t j = 0; j < resolution - 1; ++j) {
      tri.clear();
      tri.push_back(handles[i * resolution + j]);
      tri.push_back(handles[(i + 1) * resolution + j]);
      tri.push_back(handles[i * resolution + j + 1]);
      mesh.add_face(tri);
      tri.clear();
      tri.push_back(handles[(i + 1) * resolution + j]);
      tri.push_back(handles[(i + 1) * resolution + j + 1]);
      tri.push_back(handles[i * resolution + j + 1]);
      mesh.add_face(tri);
    }
  Object::updateBaseMesh(false, false);
}

bool BezierDupinPatch::reload() {
  return true;
}
