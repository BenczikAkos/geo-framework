#pragma once

#include <optional>
#include "object.hh"
#include "bspline.hh"

class Tubular : public Object {
public:
  Tubular(std::string filename);
  virtual ~Tubular();
  virtual void drawWithNames(const Visualization &vis) const override;
  virtual Vector postSelection(int selected) override;
  virtual void movement(int selected, const Vector &pos) override;
  virtual void updateBaseMesh() override;
  virtual bool reload() override;

  void updateMu(double newMu);
protected:
  //virtual Vector normal(BaseMesh::VertexHandle vh) const override;
  //virtual double meanCurvature(BaseMesh::VertexHandle vh) const override;
private:
    void readBSpline(std::ifstream &input, BSpline &result);
    std::vector<Vector> vertices;
    Vector evaluate(double u, double v) const;
    double F0(double v) const;
    double F1(double v) const;
    double G0(double v) const;
    double G1(double v) const;
    Vector ru(double u, double v) const;
    Vector rv(double u, double v) const;
    BSpline c0;
    BSpline c1;
    double mu = 1.0;
    int u_resolution = 100;
    int v_resolution = 100;
};
