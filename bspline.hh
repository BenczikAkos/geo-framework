#include "object.hh"

#include <cstddef>

using Point        = Vector;
using PointVector  = std::vector<Point>;
using VectorVector = std::vector<Vector>;
using DoubleVector = std::vector<double>;
using DoubleMatrix = std::vector<DoubleVector>;

struct BSpline
{
  size_t p;                     // degree
  size_t n;                     // n + 1 = cp.size()
  DoubleVector knots;           // n + 2 values
  PointVector cp;               // n + 1 points

  double getKnot(int i) const;
  Point getCP(int i) const;

  size_t findSpan(double u) const;
  void basisFunctionDerivatives(size_t i, double u, size_t d, DoubleMatrix &der) const;
  Point derivatives(double u, size_t d, VectorVector &der) const;
};
