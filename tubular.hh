#pragma once

#include "bspline.hh"

class Tubular : public Object {
public:
  Tubular(std::string filename);
  virtual ~Tubular();
  virtual void draw(const Visualization &vis) const override;
  virtual void drawWithNames(const Visualization &vis) const override;
  virtual Vector postSelection(int selected) override;
  virtual void movement(int selected, const Vector &pos) override;
  virtual bool reload() override;

  void updateMu(double newMu);
protected:
  virtual void updateBaseMesh() override;
  virtual Vector normal(BaseMesh::VertexHandle vh) const override;
  virtual double meanCurvature(BaseMesh::VertexHandle vh) const override;
private:
    void readBSpline(std::ifstream &input, BSpline &result);
    std::vector<Point> vertices;
    Point evaluate(double u, double v, VectorVector &c0_der, VectorVector &c1_der) const;
    double F0(double v) const;
    double dF0(double v) const;
    double F1(double v) const;
    double dF1(double v) const;
    double G0(double v) const;
    double dG0(double v) const;
    double G1(double v) const;
    double dG1(double v) const;
    Vector ru(double u, double v) const;
    Vector rv(double u, double v) const;
    Vector ruu(double u, double v) const;
    Vector rvv(double u, double v) const;
    Vector ruv(double u, double v) const;
    double calculateMeanCurvature(double u, double v) const;
    BSpline c0;
    BSpline c1;
    double mu = 1.0;
    int u_resolution = 100;
    int v_resolution = 100;
};
