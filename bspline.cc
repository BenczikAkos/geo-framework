#include "bspline.hh"

#include <algorithm>
#include <cassert>
#include <limits>

double BSpline::getKnot(int i) const {
  if (i >= (int)n + 2)
    return knots.back() + (knots[i-n-1] - knots.front());
  if (i >= 0)
    return knots[i];
  return knots.front() - (knots.back() - knots[n+i+1]);
}

Point BSpline::getCP(int i) const {
  return cp[((i % (n + 1)) + (n + 1)) % (n + 1)];
}

size_t BSpline::findSpan(double u) const
{
  return (std::upper_bound(knots.begin(), knots.end(), u) - knots.begin()) - 1;
}

void BSpline::basisFunctionDerivatives(size_t i, double u, size_t d, DoubleMatrix &der) const
{
  der.clear(); der.resize(d + 1);
  DoubleVector left(p + 1), right(p + 1), a[2];
  a[0].resize(p + 1); a[1].resize(p + 1);
  DoubleMatrix ndu(p + 1);
  ndu[0].resize(p + 1); ndu[0][0] = 1.0;
  for (size_t j = 1; j <= p; ++j) {
    ndu[j].resize(p + 1);
    left[j] = u - getKnot(i+1-j);
    right[j] = getKnot(i+j) - u;
    double saved = 0.0;
    for (size_t r = 0; r < j; ++r) {
      // lower triangle
      ndu[j][r] = right[r+1] + left[j-r];
      double tmp = ndu[r][j-1] / ndu[j][r];
      // upper triangle
      ndu[r][j] = saved + tmp * right[r+1];
      saved = tmp * left[j-r];
    }
    ndu[j][j] = saved;
  }
  for (size_t j = 0; j <= p; ++j)
    der[0].push_back(ndu[j][p]);
  for (size_t r = 0; r <= p; ++r) {
    size_t s1 = 0, s2 = 1;
    a[0][0] = 1.0;
    for (size_t k = 1; k <= d; ++k) {
      double dd = 0.0;
      int rk = r - k;
      int pk = p - k;
      if (r >= k) {
        a[s2][0] = a[s1][0] / ndu[pk+1][rk];
        dd = a[s2][0] * ndu[rk][pk];
      }
      size_t j1 = rk >= -1 ? 1 : -rk;
      size_t j2 = (int)r - 1 <= pk ? k - 1 : p - r;
      for (size_t j = j1; j <= j2; ++j) {
        a[s2][j] = (a[s1][j] - a[s1][j-1]) / ndu[pk+1][rk+j];
        dd += a[s2][j] * ndu[rk+j][pk];
      }
      if (r <= (size_t)pk) {
        a[s2][k] = -a[s1][k-1] / ndu[pk+1][r];
        dd += a[s2][k] * ndu[r][pk];
      }
      der[k].push_back(dd);
      std::swap(s1, s2);
    }
  }
  size_t r = p;
  for (size_t k = 1; k <= d; ++k) {
    for (size_t j = 0; j <= p; ++j)
      der[k][j] *= r;
    r *= p - k;
  }
}

Point BSpline::derivatives(double u, size_t d, VectorVector &der) const
{
  if (u >= knots.back()) // >= instead of == because of numerical problems
    u = knots.front();
  size_t du = std::min(d, p);
  der.clear();
  size_t span = findSpan(u);
  DoubleMatrix nder; basisFunctionDerivatives(span, u, du, nder);
  for (size_t k = 0; k <= du; ++k) {
    der.emplace_back(0.0, 0.0, 0.0);
    for (size_t j = 0; j <= p; ++j)
      der[k] += getCP(span-p+j) * nder[k][j];
  }
  for (size_t k = p + 1; k <= d; ++k)
    der.emplace_back(0.0, 0.0, 0.0);
  return der[0];
}

void BSpline::draw() const {
  glDisable(GL_LIGHTING);
  glLineWidth(3.0);
  glColor3d(0.3, 0.3, 1.0);
  // Draw control polygon
  glBegin(GL_LINE_STRIP);
  for(int i = 0; i < cp.size(); ++i){
    glVertex3dv(cp[i].data());
  }
  glVertex3dv(cp.front().data()); // Close the polygon
  glEnd();
  glLineWidth(1.0);
  glPointSize(8.0);
  glColor3d(1.0, 0.0, 1.0);
  // Draw control points
  glBegin(GL_POINTS);
  for (const auto &p : cp)
    glVertex3dv(p.data());
  glEnd();
  glPointSize(1.0);
  glEnable(GL_LIGHTING);
}
