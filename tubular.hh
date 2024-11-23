#pragma once

#include <optional>
#include "object.hh"
#include "bspline.hh"

class Tubular : public Object {
public:
  Tubular(std::string filename);
  virtual ~Tubular();
  virtual void draw(const Visualization &vis) const override;
  virtual void drawWithNames(const Visualization &vis) const override;
  virtual Vector postSelection(int selected) override;
  virtual void movement(int selected, const Vector &pos) override;
  virtual void updateBaseMesh() override;
  virtual bool reload() override;
protected:
  virtual Vector normal(BaseMesh::VertexHandle vh) const override;
  virtual double meanCurvature(BaseMesh::VertexHandle vh) const override;
private:
    void readBSpline(std::ifstream &input, std::optional<BSpline> &result);
    double F0(double v) const;
    double F1(double v) const;
    double G0(double v) const;
    double G1(double v) const;
    std::optional<BSpline> c0;
    std::optional<BSpline> c1;
    double mu = 1;
    int u_resolution = 100;
    int v_resolution = 100;
};
