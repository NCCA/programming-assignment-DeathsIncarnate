// //
// // Created by s5616052 on 14/05/25.
// //
//
// #ifndef CONTROLPANEL_H
// #define CONTROLPANEL_H
//
// // ControlPanel.h
// #pragma once
// #include <QWidget>
// #include <QSlider>
// #include <QDoubleSpinBox>
// #include <QCheckBox>
// #include "Emitter.h"
//
// class ControlPanel : public QWidget {
//     Q_OBJECT
// public:
//     explicit ControlPanel(Emitter* emitter, QWidget* parent = nullptr);
//
//     private slots:
//         void onGravityChanged(double value);
//     void onViscosityChanged(double value);
//     void onParticleSizeChanged(double value);
//     void onSimulateToggled(bool checked);
//
// private:
//     Emitter* m_emitter;
//
//     QDoubleSpinBox* m_gravitySpinBox;
//     QDoubleSpinBox* m_viscositySpinBox;
//     QDoubleSpinBox* m_particleSizeSpinBox;
//     QCheckBox* m_simulateCheckBox;
// };
//
// #endif //CONTROLPANEL_H
