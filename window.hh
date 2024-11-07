#pragma once

#include <QtWidgets/QMainWindow>
#include <QtWidgets/QSlider>

#include "viewer.hh"

class QApplication;
class QProgressBar;

class Window : public QMainWindow {
  Q_OBJECT

public:
  explicit Window(QApplication *parent);

private slots:
  void open(bool clear_others);
  void setCutoff();
  void setRange();
  void setSlicing();
  void setupUI();

private:
  QApplication *parent;
  Viewer *viewer;
  QSlider *sliderA;
  QSlider *sliderB;
  QSlider *sliderC;
  QSlider *sliderD;
  QString last_directory;
};
