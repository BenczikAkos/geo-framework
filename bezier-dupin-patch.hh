#pragma once

#include "object.hh"

class BezierDupinPatch : public Object {
public:
    BezierDupinPatch(float a, float b, float c, float d, std::pair<double, double> uRange, std::pair<double, double> vRange);
    virtual ~BezierDupinPatch();
    virtual void draw(const Visualization &vis) const override;
    virtual void drawWithNames(const Visualization &vis) const override;
    virtual Vector postSelection(int selected) override;
    virtual void movement(int selected, const Vector &pos) override;
    virtual void updateBaseMesh() override;
    virtual bool reload() override;
private:
    std::vector<Vector> control_points;
    std::vector<double> weights;
};
