#pragma once

#include <QtWidgets/QMainWindow>
#include <QtWidgets/QSpinBox>

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
  void setupParameterEditor();

private:
  QApplication *parent;
  Viewer *viewer;
  QDoubleSpinBox *spinBoxA;
  QDoubleSpinBox *spinBoxB;
  QDoubleSpinBox *spinBoxD;
  QString last_directory;
};
