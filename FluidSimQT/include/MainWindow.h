#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "NGLScene.h"


namespace Ui {
    class MainWindow;
}


class MainWindow : public QMainWindow
{
    Q_OBJECT

public slots:
    void updateSimulationParameters();

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    double currentPressure;
    double currentViscosity;
    double currentParticleSpacing;



private:
    BoundingBox m_boundingBox;
    std::unique_ptr<Physics> m_physics;
    std::unique_ptr<Emitter> m_emitter;

    Ui::MainWindow *m_ui;
		/// @brief our openGL widget
    NGLScene *m_gl;
};

#endif // MAINWINDOW_H
