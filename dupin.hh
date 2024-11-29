#pragma once

#include "object.hh"
#include <numbers>

class Dupin : public Object {
public:
  Dupin(std::string filename);
  virtual ~Dupin();
  virtual void draw(const Visualization &vis) const override;
  virtual void drawWithNames(const Visualization &vis) const override;
  virtual Vector postSelection(int selected) override;
  virtual void movement(int selected, const Vector &pos) override;
  virtual void updateBaseMesh() override;
  virtual bool reload() override;
  void setA(float value);
  void setB(float value);
  void setD(float value);  
private:
  float a, b, c, d = 0.0f;
  std::pair<float, float> range = {-10.0f, 10.0f};
  std::pair<size_t, size_t> resolution = {200, 50}; //first is bigCircle, second is channel resolution
  std::vector<Vector> controlPoints;
  void updateParameters();
  float calculateX(float u, float v) const;
  float calculateY(float u, float v) const;
  float calculateZ(float u, float v) const;
  //helper functions:
  float t0(float t) const;
  float t1(float t) const;
  float t2(float t) const;
};
