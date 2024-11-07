#pragma once

#include "object.hh"
#include <numbers>

class Dupin : public Object {
public:
  Dupin(std::string filename);
  virtual ~Dupin();
  virtual void drawWithNames(const Visualization &vis) const override;
  virtual Vector postSelection(int selected) override;
  virtual void movement(int selected, const Vector &pos) override;
  virtual void updateBaseMesh() override;
  virtual bool reload() override;
private:
  float a, b, c, d = 0.0f;
  std::pair<float, float> range = {0.0f, std::numbers::pi * 2.0f};
  std::pair<size_t, size_t> resolution = {200, 50}; //first is bigCircle, second is channel resolution
  float calculateX(float u, float v) const;
  float calculateY(float u, float v) const;
  float calculateZ(float u, float v) const;
};
