// //
// // Created by s5616052 on 14/05/25.
// //
//
// // ControlPanel.cpp
// #include "ControlPanel.h"
// #include <QVBoxLayout>
// #include <QLabel>
//
// ControlPanel::ControlPanel(Emitter* emitter, QWidget* parent)
//     : QWidget(parent), m_emitter(emitter) {
//
//     auto layout = new QVBoxLayout(this);
//
//     // Gravity Control
//     auto gravityLabel = new QLabel("Gravity:", this);
//     m_gravitySpinBox = new QDoubleSpinBox(this);
//     m_gravitySpinBox->setRange(-20.0, 20.0);
//     m_gravitySpinBox->setValue(-9.81);
//     connect(m_gravitySpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
//             this, &ControlPanel::onGravityChanged);
//
//     // Viscosity Control
//     auto viscosityLabel = new QLabel("Viscosity:", this);
//     m_viscositySpinBox = new QDoubleSpinBox(this);
//     m_viscositySpinBox->setRange(0.0, 100.0);
//     m_viscositySpinBox->setValue(10.0);
//     connect(m_viscositySpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
//             this, &ControlPanel::onViscosityChanged);
//
//     // Particle Size Control
//     auto sizeLabel = new QLabel("Particle Size:", this);
//     m_particleSizeSpinBox = new QDoubleSpinBox(this);
//     m_particleSizeSpinBox->setRange(0.1, 10.0);
//     m_particleSizeSpinBox->setValue(4.0);
//     connect(m_particleSizeSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
//             this, &ControlPanel::onParticleSizeChanged);
//
//     // Simulation Toggle
//     m_simulateCheckBox = new QCheckBox("Enable Simulation", this);
//     m_simulateCheckBox->setChecked(true);
//     connect(m_simulateCheckBox, &QCheckBox::toggled,
//             this, &ControlPanel::onSimulateToggled);
//
//     // Add widgets to layout
//     layout->addWidget(gravityLabel);
//     layout->addWidget(m_gravitySpinBox);
//     layout->addWidget(viscosityLabel);
//     layout->addWidget(m_viscositySpinBox);
//     layout->addWidget(sizeLabel);
//     layout->addWidget(m_particleSizeSpinBox);
//     layout->addWidget(m_simulateCheckBox);
// }
//
// // void ControlPanel::onGravityChanged(double value) {
// //     if (m_emitter && m_emitter->getPhysics()) {
// //         m_emitter->getPhysics()->setGravity(ngl::Vec3(0.0f, static_cast<float>(value), 0.0f));
// //     }
// // }
// //
// // void ControlPanel::onViscosityChanged(double value) {
// //     if (m_emitter && m_emitter->getPhysics()) {
// //         m_emitter->getPhysics()->setViscosity(static_cast<float>(value));
// //     }
// // }
// //
// // void ControlPanel::onParticleSizeChanged(double value) {
// //     if (m_emitter) {
// //         m_emitter->setParticleSize(static_cast<float>(value));
// //     }
// // }
// //
// // void ControlPanel::onSimulateToggled(bool checked) {
// //     if (m_emitter) {
// //         m_emitter->setSimulationEnabled(checked);
// //     }
// // }
