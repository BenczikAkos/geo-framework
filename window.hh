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
  void setupUI();
  void setDupinA(int value);
  void setDupinB(int value);
  void setDupinC(int value);
  void setDupinD(int value);

private:
  QApplication *parent;
  Viewer *viewer;
  QSpinBox *spinBoxA;
  QSpinBox *spinBoxB;
  QSpinBox *spinBoxC;
  QSpinBox *spinBoxD;
  QString last_directory;
};
