#pragma once

#include "object.hh"

class BSpline : public Object {
public:
  BSpline(std::string filename);
  BSpline(size_t degree, size_t num_control_points, std::vector<double> &knotsInput, std::vector<Vector> &controlPointsInput);
  virtual ~BSpline();
  virtual void draw(const Visualization &vis) const override;
  virtual void drawWithNames(const Visualization &vis) const override;
  virtual Vector postSelection(int selected) override;
  virtual void movement(int selected, const Vector &pos) override;
  virtual void updateBaseMesh() override;
  virtual bool reload() override;
  void printKnots() const;
  void printControlPoints() const;
  Vector evaluateBSplineNormalizedInput(double t) const;
  Vector dEvaluateBSplineNormalizedInput(double t) const;
  Vector ddEvaluateBSplineNormalizedInput(double t) const;

private:
  size_t degree;
  size_t num_control_points;
  std::vector<double> knots;
  std::vector<Vector> control_points;
  std::vector<Vector> curve_points;

  double basisFunction(size_t i, size_t k, double t) const;
  Vector evaluateBSpline(double t) const;
  Vector dEvaluateBSpline(double t) const;
  Vector ddEvaluateBSpline(double t) const;
};
