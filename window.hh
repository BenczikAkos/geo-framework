#pragma once

#include <QtWidgets/QMainWindow>
#include <QtWidgets/QLabel>

#include "viewer.hh"

class QApplication;
class QProgressBar;

class Window : public QMainWindow {
  Q_OBJECT

public:
  explicit Window(QApplication *parent);
  void setParameterSpinBoxes(double a, double b, double d);

private slots:
  void open(bool clear_others);
  void setCutoff();
  void setRange();
  void setSlicing();
  void setupParameterEditor();

private:
  QApplication *parent;
  Viewer *viewer;
  QLabel *labelA;
  QLabel *labelB;
  QLabel *labelC;
  QLabel *labelD;
  QString last_directory;
};
