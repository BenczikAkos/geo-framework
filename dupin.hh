#pragma once

#include "object.hh"
#include "window.hh"
#include "bezier-dupin-patch.hh"
#include <memory>

class Dupin : public Object {
public:
  Dupin(std::string filename, QWidget *parent);
  virtual ~Dupin();
  virtual void draw(const Visualization &vis) const override;
  virtual void drawWithNames(const Visualization &vis) const override;
  virtual Vector postSelection(int selected) override;
  virtual void movement(int selected, const Vector &pos) override;
  virtual void updateBaseMesh() override;
  virtual bool reload() override;
  void modifyA(float value);
  void modifyB(float value);
  void modifyD(float value);  
private:
  Window *window;
  float a, b, c, d = 0.0f;
  std::pair<float, float> range = {0.0f, std::numbers::pi * 2.0f};
  std::pair<size_t, size_t> resolution = {200, 50}; //first is bigCircle, second is channel resolution
  std::vector<Vector> controlPoints;
  std::unique_ptr<BezierDupinPatch> patch;
  void updateParameters();
  void updateControlPoints();
  float calculateX(float u, float v) const;
  float calculateY(float u, float v) const;
  float calculateZ(float u, float v) const;
};
